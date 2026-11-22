"""Run the finite-portfolio reverse stress experiment."""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtr
from scipy.stats import binom, norm

from bank_stresslab.persistence import CountForecast
from bank_stresslab.reference import panel_rule
from bank_stresslab.reverse_stress import count_forecasts, decision
from persistence import ALPHA, MODEL, PERSISTENCE, PORT, ROOT, write_table

LIMIT = .15
LOSS = PORT.lgd * np.arange(PORT.n+1) / PORT.n


def reversed_breach_probability(phi, kind, count):
    """Integrate from next year's factor back to the report, using binomial survival."""
    z, width = panel_rule(np.linspace(-12, 12, 49), 12)
    weights = width * norm.pdf(z)
    # For stationary Gaussian AR(1), previous | next has mean phi*z and variance 1-phi².
    previous = phi*z[:, None] + np.sqrt(1-phi**2)*z[None, :]
    h = MODEL.threshold
    loading, residual = np.sqrt(MODEL.rho), np.sqrt(1-MODEL.rho)
    previous_pd = ndtr((h-loading*previous)/residual)
    likelihood = (binom.pmf(count, PORT.n, previous_pd) if kind == "exact"
                  else binom.sf(count-1, PORT.n, previous_pd))
    report_mass = weights * (likelihood @ weights)
    last_safe_count = np.flatnonzero(LOSS <= LIMIT)[-1]
    risk = binom.sf(last_safe_count, PORT.n, ndtr((h-loading*z)/residual))
    return report_mass @ risk / report_mass.sum()


def run(order):
    grid = CountForecast(MODEL, PORT.n, order)
    grouped = {phi: [] for phi in PERSISTENCE}
    for start in range(0, PORT.n+1, 32):
        counts = np.arange(start, min(start+32, PORT.n+1))
        reached = []
        for phi, rows in grouped.items():
            block = count_forecasts(grid, phi, counts, LOSS, ALPHA, LIMIT)
            rows.extend(block)
            exact_at_end = next(r for r in block if r["information"] == "exact" and r["count"] == counts[-1])
            reached.append(phi == 0 or exact_at_end["es"] >= LIMIT)
        # Count forecasts increase under positive persistence. Keep the same
        # displayed count range for each model, including the earlier 43-count case.
        if all(reached) and counts[-1] >= 43:
            break
    curves, decisions = [], []
    for phi, rows in grouped.items():
        for kind in ("exact", "alert"):
            values = [r["es"] for r in rows if r["information"] == kind]
            if min(np.diff(values)) < -1e-10:
                raise ArithmeticError("Unexpected decrease in conditional ES")
        if phi == 0:
            np.testing.assert_allclose([r["es"] for r in rows], .1048254559136992, atol=1e-12, rtol=0)
        curves.extend(rows)
        decisions.append(decision(rows, grid.stationary, LIMIT))
    return curves, decisions


def draw(rows, decisions):
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.7), sharex=True, sharey=True, layout="constrained")
    for ax, decision_row in zip(axes.flat, decisions, strict=True):
        phi = decision_row["phi"]
        ax.axhline(100*LIMIT, color="#444444", linestyle=":", label="15% ES limit")
        for kind, color, label in (("exact", "#006485", "Exact count"), ("alert", "#9a5128", "At least this count")):
            selected = [r for r in rows if r["phi"] == phi and r["information"] == kind]
            ax.plot([r["count"] for r in selected], [100*r["es"] for r in selected],
                    color=color, linestyle="-" if kind == "exact" else "--", label=label)
            trigger = decision_row[f"{kind}_trigger"]
            if trigger is not None:
                boundary = next(r for r in selected if r["count"] == trigger)
                ax.plot(trigger, 100*boundary["es"], marker="o", color=color, markersize=5)
                offset = -15 if kind == "exact" else 8
                ax.annotate(str(trigger), (trigger, 100*boundary["es"]), xytext=(0, offset),
                            textcoords="offset points", ha="center", color=color, fontsize=9)
        ax.set_title(f"Annual persistence {phi:.1f}")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.18)
    for ax in axes[-1]:
        ax.set_xlabel("Defaults this year (out of 1,000)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Next-year 99.9% ES (% of EAD)")
    axes[0, 0].legend(frameon=False, fontsize=9)
    dest = ROOT / "figures"
    dest.mkdir(exist_ok=True)
    fig.savefig(dest / "reverse_stress.png", dpi=190)
    plt.close(fig)


def main():
    coarse, _ = run(16)
    fine, decisions = run(32)
    for old, new in zip(coarse, fine, strict=True):
        assert (old["phi"], old["information"], old["count"]) == (new["phi"], new["information"], new["count"])
        assert abs(old["es"]-new["es"]) < 1e-5 and (old["es"] >= LIMIT) == (new["es"] >= LIMIT)
    with (ROOT / "results" / "persistence_forecasts.csv").open() as source:
        historical = {(float(r["phi"]),r["observation"]):float(r["es_bp"])/10000
                      for r in csv.DictReader(source) if r["order"] == "32"}
    for r in fine:
        if r["count"] == 43:
            kind = "above_cutoff" if r["information"] == "alert" else "first_qualifying_count"
            np.testing.assert_allclose(r["es"], historical[r["phi"],kind], atol=1e-12, rtol=0)
    for row in decisions:
        for kind in ("exact", "alert"):
            k = row[f"{kind}_trigger"]
            if k is not None:
                boundary = {r["count"]:r for r in fine if r["phi"] == row["phi"] and r["information"] == kind}
                assert boundary[k]["es"] >= LIMIT and (k == 0 or boundary[k-1]["es"] < LIMIT)
                check = reversed_breach_probability(row["phi"], kind, k)
                np.testing.assert_allclose(row[f"breach_given_{kind}"], check, rtol=1e-7, atol=1e-12)
    write_table("reverse_stress_curves.csv", fine)
    write_table("reverse_stress_decisions.csv", decisions)
    draw(fine, decisions)
    for row in decisions:
        print(row["phi"], row["exact_trigger"], row["alert_trigger"])


if __name__ == "__main__":
    main()
