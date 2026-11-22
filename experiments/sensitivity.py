"""Matched-pair sensitivity and an annual-count calibration using GEMs private lending."""
from dataclasses import replace
from functools import lru_cache
from math import expm1

import numpy as np

from analyze import BP, ROOT, write_rows
from bank_stresslab.models import Model, fit_annual_counts, match_student, pair_probability
from bank_stresslab.portfolio import PORTFOLIOS
from bank_stresslab.reference import reference_distribution
from bank_stresslab.risk import tail_risk


def corporate_rho(pd):
    weight = expm1(-50 * pd) / expm1(-50)
    return .12 * weight + .24 * (1 - weight)


@lru_cache(maxsize=None)
def checked_es(model, lgd):
    portfolio = replace(PORTFOLIOS[0], lgd=lgd)
    q = pair_probability(model)
    previous = None
    for order in (16, 32, 64):
        joint, (p1, p2, _) = reference_distribution(model, portfolio.n, order)
        loss, mass = portfolio.distribution(joint)
        var, es = tail_risk(loss, mass, .999)
        if previous is not None:
            stable = var == previous[0] and abs(es - previous[1]) * BP < .05
            error = max(abs(mass.sum() - 1), abs(p1 - model.pd), abs(p2 - q))
            if stable and error < 1e-10:
                return es * BP
        previous = var, es
    raise ArithmeticError(f"Unresolved 99.9% ES for {model}")


def comparison(gaussian, df, lgd):
    target = pair_probability(gaussian)
    floor = pair_probability(Model(gaussian.pd, 0, df))
    gaussian_es = checked_es(gaussian, lgd)
    row = dict(pd=gaussian.pd, lgd=lgd, rho_gaussian=gaussian.rho, df=df,
               pair_target=target, pair_floor=floor, gaussian_es_bp=gaussian_es,
               status="infeasible", rho_student=None, student_es_bp=None, gap_bp=None)
    if floor <= target:
        student = match_student(gaussian, df)
        if abs(pair_probability(student) - target) > 1e-12:
            raise ArithmeticError("Student pair probability misses calibration target")
        es = checked_es(student, lgd)
        row.update(status="matched", rho_student=student.rho, student_es_bp=es,
                   gap_bp=es - gaussian_es)
    return row


def sensitivity_grid():
    rows = []
    for rho in (.10, corporate_rho(.01), .25, .35, .50):
        gaussian = Model(.01, rho)
        for df in (4, 5, 8, 10, 20, 50):
            row = comparison(gaussian, df, .45)
            fixed = Model(gaussian.pd, rho, df)
            fixed_es = checked_es(fixed, .45)
            row.update(fixed_rho_pair=pair_probability(fixed), fixed_rho_es_bp=fixed_es,
                       fixed_rho_gap_bp=fixed_es - row["gaussian_es_bp"])
            rows.append(row)
    return rows


def plot_grid(rows):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D

    rhos = sorted({r["rho_gaussian"] for r in rows})
    dfs = sorted({r["df"] for r in rows})
    values = np.full((len(dfs), len(rhos)), np.nan)
    for row in rows:
        if row["status"] == "matched":
            values[dfs.index(row["df"]), rhos.index(row["rho_gaussian"])] = row["gap_bp"]

    fig, ax = plt.subplots(figsize=(8, 4.6), layout="constrained")
    palette = plt.get_cmap("Blues").reversed().with_extremes(bad="#ececea")
    ax.imshow(values, cmap=palette, vmin=np.nanmin(values), vmax=0, aspect="auto")
    for i, j in np.ndindex(values.shape):
        value = values[i, j]
        label = "×" if np.isnan(value) else f"{value:+.0f}".replace("-", "−")
        color = "white" if value < np.nanmin(values) * .55 else "#183246"
        ax.text(j, i, label, ha="center", va="center", color=color, fontsize=12)
    ax.add_patch(Rectangle((rhos.index(.25)-.49, dfs.index(10)-.48), .98, .96,
                           fill=False, edgecolor="#a44118", linewidth=2))
    ax.set_xticks(range(len(rhos)), [f"{r:.3f}" if r == corporate_rho(.01)
                                      else f"{r:.2f}" for r in rhos])
    ax.set_yticks(range(len(dfs)), dfs)
    ax.set_xlabel("Gaussian latent correlation")
    ax.set_ylabel("Student-t degrees of freedom")
    ax.set_title("Student-t minus Gaussian 99.9% expected shortfall\nMatched default pairs · basis points of total exposure",
                 loc="left", fontweight="bold", pad=12)
    baseline = Rectangle((0, 0), 1, 1, fill=False, edgecolor="#a44118", linewidth=2)
    infeasible = Line2D([], [], linestyle="none", marker="x", color="#183246")
    fig.legend([baseline, infeasible], ["Synthetic baseline (ρ = 0.25, ν = 10)", "No feasible match"],
               loc="lower center", ncol=2, frameon=False)
    fig.get_layout_engine().set(rect=(0, .09, 1, .91))
    fig.savefig(ROOT / "figures/matched-es-sensitivity.png", dpi=160, bbox_inches="tight", pad_inches=.16)
    plt.close(fig)


def student_pair_floor(pd, df):
    """Vectorized zero-correlation radial integral for bootstrap PDs and ν."""
    from scipy.special import roots_legendre, stdtrit

    p, degrees = np.broadcast_arrays(pd, df)
    nodes, weights = roots_legendre(48)
    cos_squared = np.cos(np.pi * (nodes + 3) / 8)**2
    h = stdtrit(degrees, p)
    radial = np.exp(-(degrees[..., None] / 2) * np.log1p(
        h[..., None]**2 / (degrees[..., None] * cos_squared)))
    return radial @ weights / 8


def bootstrap_pairs(counts, length):
    """Resample complete years jointly using circular blocks."""
    n, d = counts["active_counterparties"], counts["defaults"]
    rates = d / n
    pairs = rates * (d - 1) / (n - 1)
    exact = pair_probability(Model(rates.mean(), 0, 10))
    if abs(student_pair_floor(rates.mean(), 10) - exact) > 1e-12:
        raise ArithmeticError("Bootstrap floor disagrees with reference integral")
    rng = np.random.default_rng([20261120, length])
    starts = rng.integers(n.size, size=(12000, (n.size + length - 1)//length))
    idx = ((starts[:, :, None] + np.arange(length)) % n.size).reshape(12000, -1)[:, :n.size]
    p, q = rates[idx].mean(axis=1), pairs[idx].mean(axis=1)
    return p, q, student_pair_floor(p, 10), pairs.mean(), exact


def gems_uncertainty(counts):
    rows = []
    for length in (1, 3, 5):
        _, q, floor, observed, minimum = bootstrap_pairs(counts, length)
        for name, estimate, samples in (("pair_probability", observed, q),
                                        ("student_floor_minus_pair", minimum - observed, floor - q)):
            lo, hi = np.quantile(samples, [.025, .975])
            rows.append(dict(block_years=length, replicates=len(samples), statistic=name,
                             point=estimate, lower_95=lo, upper_95=hi))
    return rows


def annual_counts():
    return np.genfromtxt(ROOT / "data/gems_private_1994_2023.csv", delimiter=",",
                         names=True, dtype=int)


def empirical_rows(counts):
    rows = []
    for period in (counts, counts[:15], counts[15:]):
        n, d = period["active_counterparties"], period["defaults"]
        gaussian, pair, lag_cov, phi = fit_annual_counts(n, d)
        fit = dict(first_year=int(period["year"][0]), last_year=int(period["year"][-1]),
                   defaults=int(d.sum()), counterparty_years=int(n.sum()),
                   pooled_pd=d.sum()/n.sum(), observed_pair=pair, lag_cov=lag_cov,
                   persistence=phi)
        # Table 2 gives the full-sample average recovery on private contracts.
        for df in (10, 100):
            rows.append({**fit, **comparison(gaussian, df, 1 - .722)})
    return rows


if __name__ == "__main__":
    grid = sensitivity_grid()
    write_rows("dependence_grid.csv", grid)
    plot_grid(grid)
    counts = annual_counts()
    write_rows("gems_calibration.csv", empirical_rows(counts))
    write_rows("gems_uncertainty.csv", gems_uncertainty(counts))
