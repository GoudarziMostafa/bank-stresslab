"""Check that the browser's data stays aligned with the published result tables."""
import csv
import json
from itertools import product
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def rows(name):
    with (ROOT / "results" / name).open() as stream:
        return list(csv.DictReader(stream))


class DashboardDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "dashboard" / "data.json").read_text())

    def test_anchor_values_match_result_tables(self):
        synthetic = next(r for r in rows("dependence_grid.csv")
                         if float(r["rho_gaussian"]) == .25 and float(r["df"]) == 10)
        gems = next(r for r in rows("gems_calibration.csv")
                    if r["first_year"] == "1994" and r["last_year"] == "2023" and r["df"] == "100")
        for case, key, reference in (("synthetic", "0000", synthetic), ("gems", "1111", gems)):
            for field in ("pd", "lgd", "rho_gaussian", "df", "pair_target", "pair_floor",
                          "gaussian_es_bp", "student_es_bp", "gap_bp"):
                self.assertAlmostEqual(self.data["anchors"][case][field], float(reference[field]))
                self.assertAlmostEqual(self.data["swaps"][key][field], float(reference[field]))
        for row in self.data["concentration"]:
            if row["case"] == "synthetic":
                reference = next(r for r in rows("reference.csv")
                                 if r["pd"] == "0.01" and r["alpha"] == "0.999"
                                 and r["model"] == row["model"]
                                 and r["allocation"] == row["allocation"])
                for field in ("var_bp", "es_bp"):
                    self.assertAlmostEqual(row[field], float(reference[field]))

    def test_swaps_use_all_inputs_and_flag_missing_matches(self):
        expected = {"".join(bits) for bits in product("01", repeat=4)}
        fields = ("pd", "lgd", "rho_gaussian", "df")
        self.assertEqual(self.data["swaps"].keys(), expected)
        for key, row in self.data["swaps"].items():
            for field, bit in zip(fields, key, strict=True):
                self.assertEqual(row[field], self.data["choices"][field][int(bit)])
            matched = row["pair_floor"] <= row["pair_target"]
            self.assertEqual(matched, row["status"] == "matched")
            self.assertEqual(matched, row["student_es_bp"] is not None)

    def test_intervals_keep_their_estimands(self):
        for reference in rows("gems_uncertainty.csv"):
            if reference["block_years"] == "3":
                interval = self.data["bootstrap"]["intervals"][reference["statistic"]]
                for field in ("point", "lower_95", "upper_95"):
                    self.assertEqual(interval[field], float(reference[field]))
        simulation = self.data["simulation"]
        self.assertEqual((simulation["model"], simulation["allocation"], simulation["measure"]),
                         ("t-G", "granular", "ES"))

    def test_gems_crossing_band_and_floor_presets(self):
        floor = self.data["floors"]["gems"]
        crossing = self.data["bootstrap"]["crossing"]
        self.assertLess(crossing["lower_95"], floor["crossing"])
        self.assertGreater(crossing["upper_95"], floor["crossing"])
        self.assertGreaterEqual(floor["curve"][-1][0], crossing["upper_95"])
        self.assertLess(crossing["excess_pair"]["lower_95"], crossing["excess_pair"]["point"])
        self.assertGreater(crossing["excess_pair"]["lower_95"], 0)
        self.assertEqual(crossing["no_finite_crossing"], 0)
        for scenario in self.data["floors"].values():
            self.assertEqual([p["df"] for p in scenario["presets"]],
                             [p["df"] for p in floor["presets"]])

    def test_large_name_tail_shares_and_var_tie(self):
        concentrated = [r for r in self.data["concentration"] if r["allocation"] == "concentrated"]
        for row in concentrated:
            self.assertGreater(row["large_name_tail_share"], 0)
            self.assertLess(row["large_name_tail_share"], 1)
        gaussian, student = [r for r in concentrated if r["case"] == "gems"]
        self.assertAlmostEqual(gaussian["var_bp"], student["var_bp"])
        self.assertNotAlmostEqual(gaussian["es_bp"], student["es_bp"])
        for row in (gaussian, student):
            self.assertGreater(row["large_name_tail_share"], .9999)
