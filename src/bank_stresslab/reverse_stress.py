"""Reverse an ES limit for an exact count and a binary default-count alert."""
import numpy as np
from scipy.stats import binom

from .risk import tail_risk


def count_forecasts(grid, phi, counts, loss, alpha, limit):
    """Evaluate both information sets at each integer in counts."""
    counts = np.asarray(counts)
    if (counts.ndim != 1 or not np.issubdtype(counts.dtype, np.integer)
            or np.any((counts < 0) | (counts > grid.n))):
        raise ValueError("Counts must be integers in the portfolio support")
    exact = grid.count_law[:, counts].T
    alert = binom.sf(counts[:, None]-1, grid.n, grid.pd)
    forecasts = grid.predict(phi, np.vstack((exact, alert)))
    records = []
    for kind, likelihoods, laws in zip(
        ("exact", "alert"), (exact, alert), np.split(forecasts, 2), strict=True
    ):
        probabilities = likelihoods @ grid.prior
        for count, mass, probability in zip(counts, laws, probabilities, strict=True):
            records.append(dict(phi=phi, information=kind, count=int(count),
                                observation_probability=probability,
                                es=tail_risk(loss, mass, alpha)[1],
                                realized_breach_probability=mass[loss > limit].sum()))
    return records


def decision(forecasts, count_mass, limit):
    """Use stationary count mass for event frequencies and per-alert disagreement."""
    crossings = {}
    for kind in ("exact", "alert"):
        selected = [r for r in forecasts if r["information"] == kind and r["es"] >= limit]
        crossings[kind] = min(selected, key=lambda r: r["count"]) if selected else None
    # A missing trigger occurs in none of the possible years.
    end = len(count_mass)
    e, a = [end if crossings[kind] is None else crossings[kind]["count"] for kind in ("exact", "alert")]
    if a > e:
        raise ArithmeticError("Coarsened threshold exceeds exact-count threshold")
    alert_frequency = count_mass[a:].sum()
    disagreement = count_mass[a:e].sum()
    return dict(phi=forecasts[0]["phi"], exact_trigger=e if e < end else None,
                alert_trigger=a if a < end else None, exact_frequency=count_mass[e:].sum(),
                alert_frequency=alert_frequency, disagreement_annual=disagreement,
                disagreement_given_alert=disagreement/alert_frequency if alert_frequency else None,
                breach_given_exact=crossings["exact"]["realized_breach_probability"] if e < end else None,
                breach_given_alert=crossings["alert"]["realized_breach_probability"] if a < end else None)
