"""Forecast next-year losses after an observed high-default year."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtr
from scipy.stats import t

from bank_stresslab.models import Model
from bank_stresslab.persistence import CountForecast
from bank_stresslab.portfolio import PORTFOLIOS
from bank_stresslab.reference import reference_distribution
from bank_stresslab.risk import tail_risk

ROOT = Path(__file__).resolve().parents[1]
PORT = PORTFOLIOS[0]
MODEL = Model(.01, .25)
PERSISTENCE = (0., .3, .6, .9)
MAIN_PHI = .6
ALPHA = .999
LOSS_BP = 10000 * PORT.lgd * np.arange(PORT.n + 1) / PORT.n
SEED = 20261121
BATCHES = 32
PATHS_PER_BATCH = 400000


def write_table(name, rows):
    path = ROOT / "results" / name
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def forecast_rows(order):
    grid = CountForecast(MODEL, PORT.n, order)
    cutoff = int(np.searchsorted(grid.stationary.cumsum(), .95))
    stress = grid.likelihood(cutoff)
    likelihoods = [np.ones_like(stress), stress, grid.likelihood(cutoff+1, exact=True), 1-stress]
    weights = np.asarray(likelihoods) * grid.prior
    probabilities = weights.sum(axis=1)
    posterior = weights / probabilities[:, None]
    means = posterior @ grid.nodes
    variances = (posterior * (grid.nodes - means[:, None])**2).sum(axis=1)
    baseline_es = tail_risk(LOSS_BP, grid.stationary, ALPHA)[1]
    joint, _ = reference_distribution(MODEL, PORT.n, order)
    _, original = PORT.distribution(joint)
    rows = []
    names = ("unconditional", "above_cutoff", "first_qualifying_count", "at_or_below_cutoff")
    for phi in PERSISTENCE:
        laws = grid.predict(phi, likelihoods)
        combined = probabilities[1]*laws[1] + probabilities[3]*laws[3]
        np.testing.assert_allclose(laws[0], original, rtol=0, atol=1e-11)
        np.testing.assert_allclose(combined, original, rtol=0, atol=1e-11)
        if phi == 0:
            np.testing.assert_allclose(laws, np.broadcast_to(original, laws.shape), rtol=0, atol=1e-11)
        for name, mass, probability, mean, variance in zip(
            names, laws, probabilities, means, variances, strict=True
        ):
            var, es = tail_risk(LOSS_BP, mass, ALPHA)
            rows.append(dict(order=order, phi=phi, observation=name, cutoff=cutoff,
                             event_probability=probability, next_factor_mean=phi*mean,
                             next_factor_sd=np.sqrt(1-phi**2 + phi**2*variance),
                             el_bp=mass @ LOSS_BP, var_bp=var,
                             es_bp=es, es_revision_bp=es-baseline_es))
    return rows


def simulate(cutoff):
    """Simulate unconditioned histories and select qualifying years by their observed count."""
    hist = np.zeros((BATCHES, 3, PORT.n+1), dtype=np.int64)
    threshold = MODEL.threshold
    loading, residual = np.sqrt(MODEL.rho), np.sqrt(1-MODEL.rho)
    for batch, seed in enumerate(np.random.SeedSequence(SEED).spawn(BATCHES)):
        rng = np.random.default_rng(seed)
        current = rng.standard_normal(PATHS_PER_BATCH)
        defaults = rng.binomial(PORT.n, ndtr((threshold - loading*current)/residual))
        qualifying = defaults > cutoff
        selected = current[qualifying]
        shock = rng.standard_normal(len(selected))
        for index, phi in enumerate((0., MAIN_PHI)):
            future = phi * selected + np.sqrt(1-phi**2)*shock
            count = rng.binomial(PORT.n, ndtr((threshold - loading*future)/residual))
            hist[batch, index] = np.bincount(count, minlength=PORT.n+1)
        # Additional check of the high-persistence exact-count result. These
        # paths also occur in the stress event and use the same factor shocks.
        exact = defaults[qualifying] == cutoff+1
        phi = PERSISTENCE[-1]
        future = phi * selected[exact] + np.sqrt(1-phi**2)*shock[exact]
        count = rng.binomial(PORT.n, ndtr((threshold - loading*future)/residual))
        hist[batch, 2] = np.bincount(count, minlength=PORT.n+1)
    return hist


def monte_carlo_rows(hist, reference):
    def measures(pair):
        independent, persistent, exact = [tail_risk(LOSS_BP, counts, ALPHA)[1] for counts in pair]
        return np.array([independent, persistent, persistent-independent, exact])

    n = len(hist)
    total = hist.sum(axis=0)
    pooled = measures(total)
    # Delete whole paired batches. Even the rare exact-count check keeps nearly
    # the full pool of extreme losses in each quantile and ES calculation.
    leave_one_out = np.array([measures(total - batch) for batch in hist])
    se = np.sqrt((n-1)/n * ((leave_one_out - leave_one_out.mean(axis=0))**2).sum(axis=0))
    half_width = t.ppf(.975, n-1) * se
    refs = [r["es_bp"] for r in reference if r["observation"] == "above_cutoff" and r["phi"] in (0., MAIN_PHI)]
    refs.append(refs[1]-refs[0])
    exact_ref = next(r["es_bp"] for r in reference
                     if r["observation"] == "first_qualifying_count" and r["phi"] == PERSISTENCE[-1])
    refs.append(exact_ref)
    labels = ("phi_0", "phi_0.6", "paired_difference", "phi_0.9_exact_count")
    totals = hist.sum(axis=(0, 2))
    paths = (totals[0], totals[1], totals[0], totals[2])
    return [dict(estimate=label, seed=SEED, paths=n*PATHS_PER_BATCH, qualifying_paths=selected,
                 estimate_bp=value, low95_bp=value-width, high95_bp=value+width, reference_bp=ref)
            for label, selected, value, width, ref in zip(labels, paths, pooled, half_width, refs, strict=True)]



def draw(rows):
    fig, ax = plt.subplots(figsize=(7.1, 4.25), layout="constrained")
    first_count = rows[0]["cutoff"] + 1
    for name, label, style, color in (
        ("unconditional", "Before observing defaults", "--", "#455563"),
        ("above_cutoff", f"{first_count} or more defaults", "o-", "#995022"),
        ("first_qualifying_count", f"Exactly {first_count} defaults", "s-", "#006482"),
    ):
        subset = [row for row in rows if row["observation"] == name]
        ax.plot(PERSISTENCE, [r["es_bp"]/100 for r in subset], style, label=label, color=color)
    ax.set(xlabel="Annual persistence of the common factor", ylabel="Next-year 99.9% ES (% of exposure)",
           xticks=PERSISTENCE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=.2)
    ax.legend(frameon=False)
    directory = ROOT / "figures"
    directory.mkdir(exist_ok=True)
    fig.savefig(directory / "persistence_forecast.png", dpi=190)
    plt.close(fig)


def main():
    coarse = forecast_rows(16)
    fine = forecast_rows(32)
    for a, b in zip(coarse, fine, strict=True):
        if a["cutoff"] != b["cutoff"] or a["var_bp"] != b["var_bp"] or abs(a["es_bp"]-b["es_bp"]) >= .1:
            raise ArithmeticError("Stress forecast did not converge")
    simulation = monte_carlo_rows(simulate(fine[0]["cutoff"]), fine)
    write_table("persistence_forecasts.csv", coarse + fine)
    write_table("persistence_simulation.csv", simulation)
    draw(fine)
    for row in fine:
        if row["observation"] in ("above_cutoff", "first_qualifying_count"):
            print(f"{row['observation']}, persistence {row['phi']:.1f}, ES {row['es_bp']:.2f} bp")
    print(f"Monte Carlo paired ES difference: {simulation[2]['estimate_bp']:.2f} bp")


if __name__ == "__main__":
    main()
