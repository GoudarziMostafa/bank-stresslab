from dataclasses import replace
from pathlib import Path
import unittest

import numpy as np

from bank_stresslab.models import fit_annual_counts, match_student, pair_probability


class AnnualCountTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1] / "data/gems_private_1994_2023.csv"
        cls.counts = np.genfromtxt(path, delimiter=",", names=True, dtype=int)

    def test_public_totals_and_distinct_pair_fit(self):
        n = self.counts["active_counterparties"]
        d = self.counts["defaults"]
        np.testing.assert_array_equal(self.counts["year"], np.arange(1994, 2024))
        self.assertEqual((d.sum(), n.sum()), (1996, 55993))
        model, pair, lag_cov, phi = fit_annual_counts(n, d)
        self.assertAlmostEqual(model.pd, np.mean(d/n))
        self.assertGreater(model.pd, d.sum()/n.sum())
        self.assertAlmostEqual(pair, np.mean(d*(d-1)/(n*(n-1))))
        self.assertAlmostEqual(pair_probability(model), pair, places=12)
        self.assertAlmostEqual(pair_probability(replace(model, rho=model.rho*phi)),
                               model.pd**2 + lag_cov, places=12)
        with self.assertRaisesRegex(ValueError, "floor"):
            match_student(model, 10)
        self.assertAlmostEqual(pair_probability(match_student(model, 100)), pair, places=12)

    def test_negative_lag_and_inadmissible_counts(self):
        later = self.counts[self.counts["year"] >= 2009]
        _, _, lag_cov, phi = fit_annual_counts(later["active_counterparties"], later["defaults"])
        self.assertLess(lag_cov, 0)
        self.assertIsNone(phi)
        for counts in (([10, 1], [2, 0]), ([10, 12], [12, 0]), ([10, 12], [1.5, 0])):
            with self.assertRaises(ValueError):
                fit_annual_counts(*counts)
