import unittest

import numpy as np
from numpy.testing import assert_allclose

from bank_stresslab.models import Model
from bank_stresslab.persistence import CountForecast
from bank_stresslab.reverse_stress import count_forecasts, decision


class ReverseStressTests(unittest.TestCase):
    def test_one_borrower_and_strict_realized_limit(self):
        grid = CountForecast(Model(.08, .25), 1, 12)
        loss = np.array([0., .45])
        rows = count_forecasts(grid, .6, [0, 1], loss, .99, .45)
        for r in rows:
            self.assertEqual(r["realized_breach_probability"], 0)
        exact = next(r for r in rows if r["information"] == "exact" and r["count"] == 1)
        alert = next(r for r in rows if r["information"] == "alert" and r["count"] == 1)
        assert_allclose([r["observation_probability"] for r in (exact, alert)], [.08, .08])
        self.assertAlmostEqual(exact["es"], alert["es"])

    def test_alerts_and_exact_disagree_only_in_integer_band(self):
        grid = CountForecast(Model(.12, .35), 20, 12)
        loss = .45*np.arange(21)/20
        rows = count_forecasts(grid, .6, np.arange(21), loss, .99, .27)
        out = decision(rows, grid.stationary, .27)
        exact = np.array([r["es"] >= .27 for r in rows if r["information"] == "exact"])
        alert = np.arange(21) >= out["alert_trigger"]
        self.assertLess(out["alert_trigger"], out["exact_trigger"])
        self.assertAlmostEqual(out["exact_frequency"], grid.stationary[exact].sum())
        self.assertAlmostEqual(out["disagreement_annual"], grid.stationary[alert & ~exact].sum())
        self.assertAlmostEqual(out["disagreement_given_alert"],
                               grid.stationary[alert & ~exact].sum()/grid.stationary[alert].sum())

    def test_limit_equality_and_no_trigger(self):
        grid = CountForecast(Model(.01, .25), 5, 12)
        rows = count_forecasts(grid, 0., np.arange(6), np.arange(6)/5, .99, 1.)
        absent = decision(rows, grid.stationary, 1.)
        self.assertIsNone(absent["exact_trigger"])
        self.assertIsNone(absent["breach_given_alert"])
        self.assertEqual(absent["alert_frequency"], 0.)
        boundary = [dict(r, es=.15) for r in rows]
        both = decision(boundary, grid.stationary, .15)
        self.assertEqual((both["exact_trigger"], both["alert_trigger"]), (0, 0))
        self.assertAlmostEqual(both["alert_frequency"], 1.)
        self.assertEqual(both["disagreement_given_alert"], 0.)
