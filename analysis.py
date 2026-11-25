"""Analyze four public day-0/day-12 PBMC libraries from GSE211865."""
import os

# Set thread counts before importing Scanpy and NumPy.
os.environ.setdefault("NUMBA_NUM_THREADS", "2")

from pathlib import Path
import shutil
import tarfile
import time

import anndata
import numpy as np
import pandas as pd
import requests
import scanpy as sc

from common import DATA, RAW, RESULTS, check_hash, plot_style, save_summary

SEED = 42
CAP = 3000
MARKERS = {
    "T cells": ("CD3D", "CD3E", "TRAC", "IL7R"),
    "B cells": ("MS4A1", "CD79A", "CD74", "HLA-DRA"),
    "NK cells": ("NKG7", "GNLY", "KLRD1", "PRF1"),
    "Monocytes": ("LYZ", "LST1", "S100A8", "FCN1"),
    "Plasma cells": ("MZB1", "JCHAIN", "XBP1", "TNFRSF17"),
    "Dendritic cells": ("FCER1A", "CD1C", "CLEC10A"),
}


class ArchiveStream:
    """Read an NCBI gzip tar in byte ranges and stop after the filtered matrix."""

    BLOCK = 32 * 1024 * 1024

    def __init__(self, session, url, total):
        self.session = session
        self.url = url
        self.total = total
        self.offset = 0
        self.view = memoryview(b"")

    def read(self, size):
        output = bytearray()
        while len(output) < size:
            if not self.view:
                if self.offset >= self.total:
                    break
                stop = min(self.total, self.offset + self.BLOCK) - 1
                expected = f"bytes {self.offset}-{stop}/{self.total}"
                for attempt in range(3):
                    try:
                        response = self.session.get(
                            self.url, params={"offset": self.offset}, headers={"Range": f"bytes={self.offset}-{stop}",
                                               "Accept-Encoding": "identity"},
                            stream=True, timeout=(20, 150),
                        )
                        with response:
                            response.raise_for_status()
                            if response.status_code != 206 or response.headers.get("Content-Range") != expected:
                                raise ValueError(f"Unexpected range response for {self.url}")
                            block = response.content
                        if len(block) != stop - self.offset + 1:
                            raise ValueError("Truncated archive block")
                        break
                    except (requests.RequestException, ValueError):
                        if attempt == 2:
                            raise
                        time.sleep(2**attempt)
                self.offset = stop + 1
                self.view = memoryview(block)
            take = min(len(self.view), size - len(output))
            output.extend(self.view[:take])
            self.view = self.view[take:]
        return bytes(output)


def matrix_path(sample):
    RAW.mkdir(exist_ok=True)
    target = RAW / f"{sample.sample}.h5"
    if target.is_file():
        check_hash(target, target.name)
        return target
    archive = RAW / sample.archive
    temporary = target.with_suffix(".part")

    def extract(source):
        with tarfile.open(fileobj=source, mode="r|gz") as tar:
            for member in tar:
                if member.isfile() and Path(member.name).name == "filtered_feature_bc_matrix.h5":
                    with tar.extractfile(member) as raw, temporary.open("wb") as out:
                        shutil.copyfileobj(raw, out)
                    return
        raise ValueError(f"No filtered matrix in {sample.archive}")

    if archive.is_file():
        with archive.open("rb") as stream:
            extract(stream)
    else:
        with requests.Session() as session:
            extract(ArchiveStream(session, sample.url, sample.archive_bytes))
    check_hash(temporary, target.name)
    temporary.replace(target)
    return target


def preprocess(meta):
    rng = np.random.default_rng(SEED)
    objects, audit = [], []
    for row in meta.itertuples():
        obj = sc.read_10x_h5(matrix_path(row), gex_only=True)
        obj.var_names_make_unique()
        obj.obs_names = [f"{row.sample}_{x}" for x in obj.obs_names]
        obj.var["mt"] = obj.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(obj, qc_vars=["mt"], inplace=True, percent_top=None)
        before = obj.n_obs
        included = obj.obs.n_genes_by_counts.between(300, 6500) & obj.obs.pct_counts_mt.lt(20)
        included &= obj.obs.total_counts.le(45000)
        obj = obj[included].copy()
        after = obj.n_obs
        if after < 100:
            raise ValueError(f"Too few cells pass QC in {row.sample}")
        if after > CAP:
            obj = obj[np.sort(rng.choice(after, CAP, replace=False))].copy()
        obj.obs["sample"] = row.sample
        obj.obs["visit"] = f"Day {row.day}"
        obj.obs["library"] = row.library
        audit.append({"sample": row.sample, "donor": row.donor, "day": row.day,
                      "library": row.library, "input_cells": before, "after_qc": after,
                      "sampled_cells": obj.n_obs})
        objects.append(obj)
    return anndata.concat(objects, join="inner"), pd.DataFrame(audit)


def annotate(obj):
    present = sorted(set().union(*MARKERS.values()).intersection(obj.var_names))
    matrix = obj[:, present].X.toarray()
    z = (matrix - matrix.mean(axis=0)) / (matrix.std(axis=0) + 1e-6)
    average = pd.DataFrame(z, index=obj.obs_names, columns=present).groupby(
        obj.obs["cluster"], observed=True
    ).mean()
    labels = pd.DataFrame({name: average[[g for g in genes if g in present]].mean(axis=1)
                           for name, genes in MARKERS.items()})
    margins = np.sort(labels.to_numpy(), axis=1)
    best = labels.idxmax(axis=1).where((margins[:, -1] - margins[:, -2]) >= 0.15, "Unresolved")
    obj.obs["marker_label"] = obj.obs["cluster"].map(best).astype("category")
    return labels.assign(marker_label=best), average


def main():
    RESULTS.mkdir(exist_ok=True)
    obj, counts = preprocess(pd.read_csv(DATA / "samples.csv"))
    sc.pp.filter_genes(obj, min_cells=10)
    sc.pp.normalize_total(obj, target_sum=10000)
    sc.pp.log1p(obj)
    sc.pp.highly_variable_genes(obj, n_top_genes=1800, batch_key="sample", flavor="seurat")
    sc.tl.pca(obj, n_comps=30, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(obj, n_neighbors=15, n_pcs=25, random_state=SEED)
    sc.tl.leiden(obj, resolution=0.65, flavor="igraph", n_iterations=2,
                 directed=False, key_added="cluster", random_state=SEED)
    sc.tl.umap(obj, random_state=SEED)
    scores, genes = annotate(obj)
    counts.to_csv(RESULTS / "qc.csv", index=False)
    scores.to_csv(RESULTS / "cluster_marker_scores.csv", float_format="%.5f")
    genes.to_csv(RESULTS / "cluster_marker_genes.csv", float_format="%.5f")
    composition = pd.crosstab(obj.obs["sample"], obj.obs["marker_label"], normalize="index")
    composition = composition.reindex(counts["sample"])
    composition.to_csv(RESULTS / "composition.csv", float_format="%.5f")

    plt = plot_style()
    coords = obj.obsm["X_umap"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.75), constrained_layout=True)
    colors = dict(zip(composition.columns, plt.get_cmap("tab10").colors))
    for name, group in obj.obs.groupby("marker_label", observed=True, sort=False):
        xy = coords[obj.obs_names.get_indexer(group.index)]
        axes[0].scatter(xy[:, 0], xy[:, 1], s=3, alpha=0.58, edgecolor="none",
                        color=colors[name], label=name, rasterized=True)
    axes[0].set(xlabel="UMAP 1", ylabel="UMAP 2", title="Marker-based cluster labels")
    axes[0].legend(markerscale=3, frameon=False, fontsize=8, loc="best")
    composition.plot.bar(stacked=True, ax=axes[1], color=[colors[c] for c in composition])
    axes[1].legend().remove()
    axes[1].set(xlabel="", ylabel="Fraction of analyzed cells", ylim=(0, 1),
                title="Composition within each sample")
    axes[1].set_xticklabels(
        [f"{r.donor} · day {r.day}\n{r.library}" for r in counts.itertuples()], rotation=0, fontsize=7
    )
    fig.savefig(RESULTS / "umap_and_composition.png")
    plt.close(fig)
    save_summary({
        "donors": int(counts.donor.nunique()), "libraries": len(counts), "analyzed_cells": obj.n_obs,
        "input_cells": int(counts.input_cells.sum()), "retained_features": obj.n_vars,
        "clusters": int(obj.obs["cluster"].nunique()),
        "marker_labels": sorted(obj.obs["marker_label"].unique().tolist())
    })


if __name__ == "__main__":
    main()
