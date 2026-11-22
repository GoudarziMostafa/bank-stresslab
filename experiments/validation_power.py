"""Measure size and power for two tests of credit-factor persistence."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom, binomtest, norm

from bank_stresslab.models import Model
from bank_stresslab.persistence import CountForecast
from bank_stresslab.validation import OrderedLikelihood, permutation_pvalues, simulate_histories

ROOT = Path(__file__).resolve().parents[1]
MODEL = Model(.01, .25)
PHI = .6
YEARS = (10, 20, 40)
SEED = 20261122
SIZE = .05
CALIBRATION = 40000
EVALUATION = 12000


def rng(stream):
    return np.random.default_rng(np.random.SeedSequence([SEED, stream]))


def write_rows(filename, rows):
    path = ROOT / "results" / filename
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def cutoffs(scores):
    ordered = np.sort(scores, axis=0)
    n = len(ordered)
    rank = int(np.ceil((1-SIZE)*(n+1)))
    lower, upper = binom.ppf([.025, .975], n, 1-SIZE).astype(int) + [0, 1]
    return ordered[[rank-1, lower-1, upper-1]]


def rate_row(true_phi, years, test, rejected):
    hits = int(rejected.sum())
    ci = binomtest(hits, len(rejected)).proportion_ci(method="wilson")
    return dict(phi=true_phi, years=years, test=test, histories=len(rejected), rejections=hits,
                rejection_rate=rejected.mean(), low95=ci.low, high95=ci.high,
                nonrejection_rate=1-rejected.mean())


def draw(rows):
    fig, ax = plt.subplots(figsize=(7.0, 4.2), layout="constrained")
    for test, label, color, marker in (("permutation", "Lag-one permutation test", "#95512a", "s"),
                                        ("likelihood", "Ordered likelihood ratio", "#006985", "o")):
        data = [r for r in rows if r["phi"] == PHI and r["test"] == test]
        power = np.array([r["rejection_rate"] for r in data])*100
        lower = power - np.array([r["low95"] for r in data])*100
        upper = np.array([r["high95"] for r in data])*100 - power
        ax.errorbar(YEARS, power, yerr=[lower, upper], label=label, color=color, marker=marker,
                    capsize=4)
    ax.set(xlabel="Years of observed default counts", ylabel="Probability of detecting persistence (%)",
           xticks=YEARS, ylim=(0, 100))
    ax.grid(axis="y", alpha=.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    directory = ROOT / "figures"
    directory.mkdir(exist_ok=True)
    fig.savefig(directory / "validation_power.png", dpi=190)
    plt.close(fig)


def main():
    grid = CountForecast(MODEL, 1000, 8)
    likelihood = OrderedLikelihood(grid, PHI)
    finer = OrderedLikelihood(CountForecast(MODEL, 1000, 16), PHI)
    calibration = simulate_histories(MODEL, 1000, 0., max(YEARS), CALIBRATION, rng(0))
    thresholds, low, high = cutoffs(likelihood.log_ratios(calibration, YEARS))
    rows, diagnostics, matched_cutoffs = [], [], []
    for index, phi in enumerate((0., PHI)):
        counts = simulate_histories(MODEL, 1000, phi, max(YEARS), EVALUATION, rng(index+1))
        scores = likelihood.log_ratios(counts, YEARS)
        # Include the most extreme count seen, in addition to a fixed held-out prefix.
        audit = np.r_[np.arange(256), np.argmax(counts.max(axis=1))]
        deviations = np.max(np.abs(finer.log_ratios(counts[audit], YEARS) - scores[audit]), axis=0)
        if np.max(deviations) > .001:
            raise ArithmeticError(f"Log likelihood ratios require finer quadrature: {deviations}")
        for j, years in enumerate(YEARS):
            pvalue = permutation_pvalues(counts[:, :years], 199, rng(10+index*3+j))
            practical = pvalue <= SIZE
            oracle = scores[:, j] > thresholds[j]
            for name, hits in (("permutation", practical), ("likelihood", oracle)):
                rows.append(rate_row(phi, years, name, hits))
            # Supplementary comparison: on the evaluation null histories set
            # a threshold matching the permutation test's realized false alarms.
            # Evaluate that fixed threshold on independent alternative histories.
            if phi == 0.:
                rank = len(scores) - practical.sum()
                matched_cutoffs.append(np.sort(scores[:, j])[rank-1])
            matched_rate = np.mean(scores[:, j] > matched_cutoffs[j])
            difference = oracle.astype(int) - practical.astype(int)
            margin = norm.ppf(.975)*difference.std(ddof=1)/np.sqrt(len(counts))
            diagnostics.append(dict(phi=phi, years=years, log_lr_cutoff=thresholds[j],
                                    cutoff_low95=low[j], cutoff_high95=high[j],
                                    rate_at_high_cutoff=np.mean(scores[:, j] > high[j]),
                                    rate_at_low_cutoff=np.mean(scores[:, j] > low[j]),
                                    matched_cutoff=matched_cutoffs[j], rate_at_matched_cutoff=matched_rate,
                                    paired_gain=difference.mean(), gain_low95=difference.mean()-margin,
                                    gain_high95=difference.mean()+margin, max_audit_log_error=deviations[j]))
    write_rows("validation_power.csv", rows)
    write_rows("validation_diagnostics.csv", diagnostics)
    draw(rows)
    for row in rows:
        print(f"phi={row['phi']}, {row['years']} years, {row['test']}: {row['rejection_rate']:.3%}")


if __name__ == "__main__":
    main()
