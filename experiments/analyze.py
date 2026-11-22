"""Validate the reference and paired simulations and generate all public tables and figures."""
import csv
import platform
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import binom, t

from bank_stresslab.models import Model, match_student, pair_probability
from bank_stresslab.portfolio import PORTFOLIOS
from bank_stresslab.reference import reference_distribution
from bank_stresslab.risk import tail_risk

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
LEVELS = (.99, .999)
BP = 10000


def write_rows(name, rows):
    with (RESULTS / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def deterministic():
    gaussian = Model(.01, .25)
    models = (gaussian, match_student(gaussian, 10))
    risks, checks, curves, calibration = [], [], {}, []
    for pd in (.01, .02):
        for original in models:
            model = replace(original, pd=pd)
            q = pair_probability(model)
            previous = {}
            for order in (8, 16, 32):
                joint, (m1, m2, m3) = reference_distribution(model, 1000, order)
                for portfolio in PORTFOLIOS:
                    loss, mass = portfolio.distribution(joint)
                    el, variance = portfolio.moments(pd, q)
                    for alpha in LEVELS:
                        key = (portfolio.name, alpha)
                        var, es = tail_risk(loss, mass, alpha)
                        i = np.searchsorted(loss, var)
                        # Distances in probability from alpha to either side of the VaR jump.
                        below = mass[i:].sum() - (1 - alpha)
                        above = (1 - alpha) - mass[i+1:].sum()
                        checks.append(dict(pd=pd, model=model.name, allocation=portfolio.name,
                                           alpha=alpha, order=order, mass_error=mass.sum()-1,
                                           pd_error=m1-pd, pair_error=m2-q,
                                           el_error_bp=(loss @ mass - el)*BP,
                                           variance_error=((loss-el)**2 @ mass - variance),
                                           var_bp=var*BP, es_bp=es*BP,
                                           cdf_margin_left=below, cdf_margin_right=above))
                        if order == 32:
                            old_var, old_es, old_mass = previous[key]
                            sf_error = np.max(np.abs(np.cumsum((mass-old_mass)[::-1])))
                            assert abs(mass.sum()-1) < 1e-10 and np.all(mass >= 0)
                            assert max(abs(m1-pd), abs(m2-q)) < 1e-9
                            assert abs(loss @ mass - el) < 1e-10
                            assert abs((loss-el)**2 @ mass - variance) < 1e-10
                            assert old_var == var and abs(old_es-es)*BP < .1
                            assert min(below, above) > sf_error
                            risks.append(dict(pd=pd, model=model.name, allocation=portfolio.name,
                                              alpha=alpha, el_bp=el*BP, sd_bp=np.sqrt(variance)*BP,
                                              var_bp=var*BP, es_bp=es*BP, capital_bp=(var-el)*BP))
                            if pd == .01:
                                curves[model.name, portfolio.name] = (loss, mass)
                        previous[key] = var, es, mass
            calibration.append(dict(pd=pd, model=model.name, latent_rho=model.rho,
                                    df=model.df or "Gaussian", pair_pd=q,
                                    default_corr=(q-pd**2)/(pd*(1-pd)), third_given_two=m3/q))
    write_rows("reference.csv", risks)
    write_rows("calibration.csv", calibration)
    write_rows("quadrature_checks.csv", checks)
    return curves


def estimates(batch):
    values = {}
    for model, joint in zip(("Gaussian", "Student-t"), batch, strict=True):
        for portfolio in PORTFOLIOS:
            loss, counts = portfolio.distribution(joint)
            for alpha in LEVELS:
                var, es = tail_risk(loss, counts, alpha)
                values[model, portfolio.name, alpha, "ES"] = es*BP
                values[model, portfolio.name, alpha, "VaR"] = var*BP
    return add_contrasts(values)


def add_contrasts(values):
    for alpha in LEVELS:
        for portfolio in PORTFOLIOS:
            for metric in ("ES", "VaR"):
                values["t-G", portfolio.name, alpha, metric] = (
                    values["Student-t", portfolio.name, alpha, metric]
                    - values["Gaussian", portfolio.name, alpha, metric])
        values["t-G", "concentrated-minus-granular", alpha, "ES"] = (
            values["t-G", "concentrated", alpha, "ES"] - values["t-G", "granular", alpha, "ES"])
    return values


def simulation_tables(curves):
    files = sorted(RESULTS.glob("paired_*.npz"), key=lambda f: int(f.stem.split("_")[1]))
    if not files:
        raise FileNotFoundError("Run experiments/simulate.py first")
    rows = []
    refs = estimates_for_reference(curves)
    for file in files:
        with np.load(file) as data:
            counts = data["counts"]
            n = int(data["scenarios"])
        assert np.all(counts.sum(axis=(2, 3)) == n // len(counts))
        pooled = counts.sum(axis=0)
        values = estimates(pooled)
        batches = [estimates(batch) for batch in counts]
        for key, estimate in values.items():
            model, allocation, alpha, metric = key
            lo, hi, method = "", "", ""
            if metric == "ES":
                se = np.std([b[key] for b in batches], ddof=1)/np.sqrt(len(batches))
                width = t.ppf(.975, len(batches)-1)*se
                lo, hi, method = estimate-width, estimate+width, "paired batches"
            elif model != "t-G":
                port = next(p for p in PORTFOLIOS if p.name == allocation)
                loss, hist = port.distribution(pooled[int(model == "Student-t")])
                # Binomial central interval inverted into 1-based order statistics.
                ranks = binom.ppf([.025, .975], n, alpha).astype(int) + [0, 1]
                bounds = loss[np.searchsorted(hist.cumsum(), np.clip(ranks, 1, n))]*BP
                lo, hi, method = *bounds, "order statistics"
            rows.append(dict(scenarios=n, model=model, allocation=allocation, alpha=alpha,
                             measure=metric, estimate_bp=estimate, reference_bp=refs[key],
                             difference_bp=estimate-refs[key], low95_bp=lo, high95_bp=hi, interval=method))
    write_rows("monte_carlo.csv", rows)
    return rows


def estimates_for_reference(curves):
    # Reconstruct the same named contrasts as the simulation table, without using its histograms.
    values = {}
    for (model, allocation), (loss, mass) in curves.items():
        for alpha in LEVELS:
            var, es = tail_risk(loss, mass, alpha)
            values[model, allocation, alpha, "VaR"] = var*BP
            values[model, allocation, alpha, "ES"] = es*BP
    return add_contrasts(values)


def draw(curves):
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.3), sharex=True, sharey=True, layout="constrained")
    for ax, portfolio, title in zip(axes, PORTFOLIOS, ("Equal exposures", "One borrower at 10%"), strict=True):
        for model, color, style in (("Gaussian", "#005b83", "-"), ("Student-t", "#ae4417", "--")):
            loss, mass = curves[model, portfolio.name]
            survival = np.r_[np.cumsum(mass[:0:-1])[::-1], 0.]
            ax.step(loss*100, np.maximum(survival, 1e-12), where="post", label=model,
                    linewidth=1.8, color=color, linestyle=style)
        ax.set(yscale="log", ylim=(1e-5, 1), xlim=(0, 22), title=title,
               xlabel="Loss threshold (% of total exposure)")
        ax.grid(axis="y", which="major", alpha=.22)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Probability of exceeding threshold")
    axes[0].legend(frameon=False)
    directory = ROOT / "figures"
    directory.mkdir(exist_ok=True)
    fig.savefig(directory / "loss_exceedance.png", dpi=190)
    plt.close(fig)


def draw_precision(rows):
    primary = [r for r in rows if r["model"] == "t-G" and r["allocation"] == "granular"
               and r["alpha"] == .999 and r["measure"] == "ES"]
    estimates_bp = np.array([r["estimate_bp"] for r in primary])
    widths = np.array([r["high95_bp"] for r in primary]) - estimates_bp
    positions = np.arange(len(primary))
    fig, ax = plt.subplots(figsize=(6.8, 3.7), layout="constrained")
    ax.axhline(primary[0]["reference_bp"], color="#ae4417", linestyle="--",
               label="Finite-portfolio reference")
    ax.errorbar(positions, estimates_bp, yerr=widths, fmt="o", markersize=6, capsize=5,
                color="#005b83", label="Paired simulation, 95% interval")
    ax.set_xticks(positions, [f"{r['scenarios']/1e6:g} million" for r in primary])
    ax.set(xlim=(-.35, len(primary)-.65), xlabel="Cumulative scenarios",
           ylabel="Student-t minus Gaussian ES (bp)", title="Granular portfolio, 99.9% Expected Shortfall")
    ax.grid(axis="y", alpha=.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(ROOT / "figures" / "simulation_precision.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    RESULTS.mkdir(exist_ok=True)
    curves = deterministic()
    draw(curves)
    draw_precision(simulation_tables(curves))
    print("Verified quadrature and generated tables and figure.")
    print(platform.python_version(), np.__version__, scipy.__version__, matplotlib.__version__)
