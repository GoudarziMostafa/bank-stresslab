"""Can a validator separate a changing grade mix from a persistent credit factor?"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest, norm

from bank_stresslab.borrower_mix import AttributionFilter, GradeBook, simulate_book
from validation_power import ROOT, cutoffs, write_rows

BOOK = GradeBook()
SEED = 20261123
CYCLE = .6
MIX = .8
VALIDATORS = ("totals_only", "grade_aware")
SCENARIOS = (("neither", 0., 0.), ("mix_only", 0., MIX),
             ("cycle_only", CYCLE, 0.), ("both", CYCLE, MIX))


def histories(phi, mix, size, stream):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, stream]))
    return simulate_book(BOOK, phi, mix, size, 20, rng)


def check_stationary_law(engine):
    pmf = engine.stationary
    k = np.arange(BOOK.n+1)
    np.testing.assert_allclose(pmf.sum(), 1, atol=1e-12, rtol=0)
    np.testing.assert_allclose(k @ pmf, 10., atol=1e-10, rtol=0)
    # Independent calculation of the count PGF without convolutions.
    lo, hi = BOOK.probabilities(engine.nodes)
    for s in (.2, .8, .98):
        by_state = [(1-lo+lo*s)**(BOOK.n-n) * (1-hi+hi*s)**n for n in BOOK.high_sizes]
        np.testing.assert_allclose(s**k @ pmf, engine.prior @ np.mean(by_state, axis=0),
                                   atol=1e-12, rtol=0)


def draw(rows):
    fig, axes = plt.subplots(1, 2, figsize=(8.7, 3.9), layout="constrained")
    for ax, selected, title in zip(axes, (("neither", "mix_only"), ("cycle_only", "both")),
                                   ("False cycle attribution", "Cycle detection"), strict=True):
        for j, (name, label, color) in enumerate(zip(VALIDATORS, ("Totals only", "Grades observed"),
                                                     ("#9a5128", "#006485"), strict=True)):
            subset = [r for r in rows if r["scenario"] in selected and r["validator"] == name
                      and r["calibration"] == "neither"]
            y = 100*np.array([r["rate"] for r in subset])
            error = [y-100*np.array([r["low95"] for r in subset]),
                     100*np.array([r["high95"] for r in subset])-y]
            x = np.arange(2)+(j-.5)*.33
            ax.bar(x, y, .31, label=label, color=color)
            ax.errorbar(x, y, yerr=error, fmt="none", color="#34393c", capsize=3)
        ax.set(xticks=[0, 1], xticklabels=["Independent mix", "Persistent mix"], title=title)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_ylabel("Histories flagged (%)")
        ax.grid(axis="y", alpha=.15)
    axes[0].axhline(5, linestyle=":", color="#444444", label="5% target")
    axes[0].set_ylim(0, max(10, axes[0].get_ylim()[1]*1.35))
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[1].set_ylim(0, 100)
    destination = ROOT / "figures"
    destination.mkdir(exist_ok=True)
    fig.savefig(destination / "borrower_mix.png", dpi=190)
    plt.close(fig)


def main():
    engine = AttributionFilter(BOOK)
    fine = AttributionFilter(BOOK, order=16)
    check_stationary_law(engine)
    thresholds = {}
    for stream, (name, phi, mix) in enumerate(SCENARIOS[:2]):
        thresholds[name] = cutoffs(engine.scores(histories(phi, mix, 12000, stream)))
        print(f"Calibrated on {name}", flush=True)
    rows, paired, audits, scores = [], [], [], {}
    for stream, (name, phi, mix) in enumerate(SCENARIOS, start=2):
        data = histories(phi, mix, 8000, stream)
        score = engine.scores(data)
        scores[name] = score
        totals = data[:, :, 1:].sum(axis=2)
        audit = np.unique(np.r_[np.arange(96), np.argmax(totals.max(axis=1)),
                                 np.argmax(data[:, :, 1].max(axis=1)), np.argmax(data[:, :, 2].max(axis=1))])
        error = np.abs(fine.scores(data[audit]) - score[audit]).max(axis=0)
        if max(error) > .002:
            raise ArithmeticError(f"Unresolved log ratios for {name}: {error}")
        audits.append(dict(scenario=name, audit_histories=len(audit), max_log_error_total=error[0],
                           max_log_error_grade=error[1], mean_annual_defaults=totals.mean(),
                           mean_high_share=np.take(BOOK.high_sizes, data[:, :, 0]).mean()/BOOK.n))
        for calibration, (cut, lower, upper) in thresholds.items():
            rejected = score > cut
            for j, validator in enumerate(VALIDATORS):
                hits = int(rejected[:, j].sum())
                ci = binomtest(hits, len(data)).proportion_ci(method="wilson")
                rows.append(dict(scenario=name, cycle=phi, mix=mix, calibration=calibration,
                                 validator=validator, histories=len(data), flagged=hits,
                                 rate=hits/len(data), low95=ci.low, high95=ci.high, log_cutoff=cut[j],
                                 rate_at_high_cutoff=(score[:, j]>upper[j]).mean(),
                                 rate_at_low_cutoff=(score[:, j]>lower[j]).mean()))
            diff = rejected[:, 0].astype(int)-rejected[:, 1].astype(int)
            margin = norm.ppf(.975)*diff.std(ddof=1)/np.sqrt(len(diff))
            paired.append(dict(scenario=name, calibration=calibration, total_minus_grade=diff.mean(),
                               low95=diff.mean()-margin, high95=diff.mean()+margin,
                               only_total=int((diff==1).sum()), only_grade=int((diff==-1).sum())))
        print(f"Evaluated {name}", flush=True)
    # Descriptive equal-size comparison. These null observations set the cutoffs;
    # only the corresponding alternative histories supply the power estimates.
    equal_size = []
    for null, alternative in (("neither", "cycle_only"), ("mix_only", "both")):
        null_scores = scores[null]
        allowed_flags = round(.05*len(null_scores))
        cut = np.sort(null_scores, axis=0)[len(null_scores)-allowed_flags-1]
        for j, validator in enumerate(VALIDATORS):
            null_rate = (scores[null][:, j]>cut[j]).mean()
            if null_rate != allowed_flags/len(null_scores):
                raise ArithmeticError("Scoring ties prevent exact size matching")
            equal_size.append(dict(null=null, alternative=alternative, validator=validator,
                                   null_rate=null_rate, detection=(scores[alternative][:,j]>cut[j]).mean()))
    write_rows("borrower_mix_rates.csv", rows)
    write_rows("borrower_mix_paired.csv", paired)
    write_rows("borrower_mix_equal_size.csv", equal_size)
    write_rows("borrower_mix_checks.csv", audits)
    draw(rows)


if __name__ == "__main__":
    main()
