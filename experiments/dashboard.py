"""Export the static dashboard from checked results and computed scenarios."""
import csv
import json
from dataclasses import replace
from itertools import product
from math import ceil, isclose

import numpy as np
from scipy.optimize import brentq

from analyze import BP, ROOT
from bank_stresslab.models import Model, match_student, pair_probability
from bank_stresslab.portfolio import PORTFOLIOS
from bank_stresslab.reference import reference_distribution
from bank_stresslab.risk import tail_risk
from sensitivity import annual_counts, bootstrap_pairs, checked_es, comparison, student_pair_floor


def result_rows(name):
    """Keep CSV numbers as numbers and unavailable values as nulls."""
    def value(s):
        if s == "":
            return None
        try:
            return float(s)
        except ValueError:
            return s

    with (ROOT / "results" / name).open() as stream:
        return [{key: value(s) for key, s in row.items()} for row in csv.DictReader(stream)]


def require_close(actual, expected):
    if not isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-9):
        raise ValueError(f"Dashboard result {actual} disagrees with CSV value {expected}")


def concentration(anchors):
    rows = []
    for case, row in anchors.items():
        gaussian = Model(row["pd"], row["rho_gaussian"])
        for model in (gaussian, match_student(gaussian, row["df"])):
            previous = None
            for order in (16, 32, 64):
                joint, _ = reference_distribution(model, PORTFOLIOS[0].n, order)
                values = [tuple(BP * x for x in tail_risk(
                    *replace(portfolio, lgd=row["lgd"]).distribution(joint), .999))
                    for portfolio in PORTFOLIOS]
                if previous and all(v[0] == p[0] and abs(v[1] - p[1]) < .05
                                    for v, p in zip(values, previous, strict=True)):
                    break
                previous = values
            else:
                raise ArithmeticError(f"Concentration estimate did not converge for {model}")
            for portfolio, (var, es) in zip(PORTFOLIOS, values, strict=True):
                entry = dict(case=case, model=model.name, allocation=portfolio.name,
                             var_bp=var, es_bp=es)
                if portfolio.name == "concentrated":
                    entry["large_name_tail_share"] = large_name_tail_share(
                        joint, replace(portfolio, lgd=row["lgd"]), var / BP)
                rows.append(entry)
    return rows


def large_name_tail_share(joint, portfolio, var):
    """Default share of the ES tail, taking a proportional slice of the VaR atom."""
    loss, mass = portfolio.distribution(joint)
    index = int(np.argmin(abs(loss - var)))
    default_mass = np.zeros_like(mass)
    default_mass[portfolio.large_units:] = joint[1]
    tail_mass = (1 - .999) * mass.sum()
    share_at_cutoff = (tail_mass - mass[index + 1:].sum()) / mass[index]
    return (default_mass[index + 1:].sum() + share_at_cutoff * default_mass[index]) / tail_mass


def crossing_draws(p, q):
    """Find the first feasible ν for draws with positive excess over independence."""
    positive = q > p**2
    results = np.full(p.shape, np.inf)
    pd, target = p[positive], q[positive]
    lower = np.full_like(pd, 2.0)
    upper = np.full_like(pd, 1e7)
    if np.any(student_pair_floor(pd, upper) >= target):
        raise ArithmeticError("Could not bracket a bootstrap crossing below 10 million degrees")
    for _ in range(44):
        mid = np.sqrt(lower * upper)
        infeasible = student_pair_floor(pd, mid) > target
        lower = np.where(infeasible, mid, lower)
        upper = np.where(infeasible, upper, mid)
    results[positive] = np.sqrt(lower * upper)
    return results


def crossing_summary(p, q, gems):
    crossings = crossing_draws(p, q)
    lo, hi = np.quantile(crossings, [.025, .975], method="inverted_cdf")
    excess = q - p**2
    excess_lo, excess_hi = np.quantile(excess, [.025, .975])
    return dict(lower_95=float(lo) if np.isfinite(lo) else None,
                upper_95=float(hi) if np.isfinite(hi) else None,
                no_finite_crossing=int(np.isinf(crossings).sum()),
                excess_pair=dict(point=gems["pair_target"] - gems["pd"]**2,
                                 lower_95=float(excess_lo), upper_95=float(excess_hi)))


def floor_case(pd, rho, dfs, maximum=200):
    target = pair_probability(Model(pd, rho))

    def floor(df):
        return pair_probability(Model(pd, 0, df))

    crossing = brentq(lambda df: floor(df) - target, 3, 1000)
    # The boundary and every highlighted ν are evaluated directly.
    points = sorted(set(np.geomspace(4, maximum, 96)) | set(dfs) | {crossing})
    return dict(pd=pd, rho=rho, target=target, crossing=crossing,
                curve=[(df, floor(df)) for df in points],
                presets=[dict(df=df, floor=floor(df)) for df in dfs])


def build():
    grid = result_rows("dependence_grid.csv")
    synthetic = next(r for r in grid if r["rho_gaussian"] == .25 and r["df"] == 10)
    basel = next(r for r in grid if .19 < r["rho_gaussian"] < .20 and r["df"] == 10)
    gems = next(r for r in result_rows("gems_calibration.csv")
                if r["first_year"] == 1994 and r["last_year"] == 2023 and r["df"] == 100)
    anchors = {"synthetic": synthetic, "gems": gems}
    choices = {key: [anchors[case][key] for case in anchors]
               for key in ("pd", "lgd", "rho_gaussian", "df")}
    swaps = {}
    for indices in product((0, 1), repeat=len(choices)):
        pd, lgd, rho, df = (values[i] for values, i in zip(choices.values(), indices, strict=True))
        gaussian = Model(pd, rho)
        row = comparison(gaussian, df, lgd)
        fixed = Model(pd, rho, df)
        row.update(fixed_rho_pair=pair_probability(fixed), fixed_rho_es_bp=checked_es(fixed, lgd))
        swaps["".join(str(i) for i in indices)] = row

    concentration_rows = concentration(anchors)
    reference = {(r["model"], r["allocation"]): r for r in result_rows("reference.csv")
                 if r["pd"] == synthetic["pd"] and r["alpha"] == .999}
    for row in concentration_rows:
        if row["case"] == "synthetic":
            expected = reference[row["model"], row["allocation"]]
            for field in ("var_bp", "es_bp"):
                require_close(row[field], expected[field])
    for key, expected in (("0000", synthetic), ("1111", gems)):
        for field in ("pair_target", "pair_floor", "rho_student", "gaussian_es_bp",
                      "student_es_bp", "gap_bp"):
            require_close(swaps[key][field], expected[field])
    for field in ("fixed_rho_pair", "fixed_rho_es_bp"):
        require_close(swaps["0000"][field], synthetic[field])

    p, q, floor, _, _ = bootstrap_pairs(annual_counts(), 3)
    crossing = crossing_summary(p, q, gems)
    bound = crossing["upper_95"]
    maximum = max(200, 100 * ceil(bound / 100)) if bound is not None else 500
    dfs = [10, 20, 30, 50, 60, 100, 200]
    floors = {label: floor_case(row["pd"], row["rho_gaussian"], dfs,
                                maximum=maximum if label == "gems" else 200)
              for label, row in dict(synthetic=synthetic, basel=basel, gems=gems).items()}
    if any([point["df"] for point in scenario["presets"]] != dfs for scenario in floors.values()):
        raise ValueError("Floor scenarios must share the same degrees-of-freedom presets")
    intervals = {r["statistic"]: r for r in result_rows("gems_uncertainty.csv")
                 if r["block_years"] == 3}
    for label, draws in (("pair_probability", q), ("student_floor_minus_pair", floor - q)):
        lo, hi = np.quantile(draws, [.025, .975])
        require_close(lo, intervals[label]["lower_95"])
        require_close(hi, intervals[label]["upper_95"])

    gems_model = Model(gems["pd"], gems["rho_gaussian"])
    gems_floor = floors["gems"]["crossing"]
    boundary_model = Model(gems["pd"], 0, gems_floor)
    require_close(pair_probability(boundary_model), gems["pair_target"])
    boundary_gap = checked_es(boundary_model, gems["lgd"]) - gems["gaussian_es_bp"]

    return dict(anchors=anchors, choices=choices, swaps=swaps, floors=floors,
                gems_boundary_gap_bp=boundary_gap,
                gems_range=[comparison(gems_model, df, gems["lgd"]) for df in (60, 100, 200)],
                concentration=concentration_rows,
                bootstrap=dict(intervals=intervals, crossing=crossing, correlation=np.corrcoef(q, floor)[0, 1],
                               points=list(zip(q[::20].tolist(), floor[::20].tolist(), strict=True))),
                simulation=next(r for r in result_rows("monte_carlo.csv")
                                if r["alpha"] == .999 and r["scenarios"] == 6400000
                                and r["model"] == "t-G" and r["allocation"] == "granular"
                                and r["measure"] == "ES"))


if __name__ == "__main__":
    path = ROOT / "dashboard" / "data.json"
    path.write_text(json.dumps(build(), separators=(",", ":"), allow_nan=False) + "\n")
    print(f"Wrote {path.relative_to(ROOT)}")
