import unittest

import numpy as np
import pandas as pd

from common import DATA, RESULTS


class SingleCellAtlasTests(unittest.TestCase):
    def test_each_donor_has_two_matched_visits(self):
        samples = pd.read_csv(DATA / "samples.csv")
        self.assertEqual(samples.groupby("donor").day.apply(set).tolist(), [{0, 12}, {0, 12}])
        self.assertTrue(samples.groupby("donor").library.nunique().eq(1).all())
        self.assertEqual(samples.library.nunique(), 2)

    def test_composition_and_quality_control(self):
        samples = pd.read_csv(RESULTS / "qc.csv").set_index("sample")
        table = pd.read_csv(RESULTS / "composition.csv", index_col=0)
        self.assertEqual(samples.index.tolist(), table.index.tolist())
        self.assertTrue((samples.input_cells >= samples.after_qc).all())
        self.assertTrue((samples.after_qc >= samples.sampled_cells).all())
        np.testing.assert_allclose(table.sum(axis=1), 1.0, atol=1e-4)
