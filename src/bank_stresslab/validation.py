"""Tests for persistence in ordered annual default counts."""
import numpy as np
from scipy.special import ndtr
from scipy.stats import norm


def simulate_histories(model, n, phi, years, size, rng):
    """Independent stationary histories with fixed annual portfolio composition."""
    if model.df is not None or not 0 <= phi < 1 or min(n, years, size) < 1:
        raise ValueError("Expected a Gaussian model, valid persistence and positive dimensions")
    z = rng.standard_normal(size)
    scale = np.sqrt(1-phi**2)
    threshold = model.threshold
    counts = np.empty((size, years), dtype=int)
    for year in range(years):
        if year:
            z = phi*z + scale*rng.standard_normal(size)
        pd = ndtr((threshold - np.sqrt(model.rho)*z) / np.sqrt(1-model.rho))
        counts[:, year] = rng.binomial(n, pd)
    return counts


class OrderedLikelihood:
    """Integrate a stationary AR(1) factor on an existing CountForecast grid."""

    def __init__(self, grid, phi):
        if not 0 <= phi < 1:
            raise ValueError("Persistence must be in [0, 1)")
        self.grid = grid
        scale = np.sqrt(1 - phi**2)
        distance = (grid.nodes[None, :] - phi*grid.nodes[:, None])/scale
        self.transition = norm.pdf(distance) * grid.widths / scale
        self.emission = np.ascontiguousarray(grid.count_law.T)

    def log_ratios(self, counts, horizons, chunk=256):
        """One alternative-vs-independent log likelihood ratio per history and horizon."""
        counts = np.asarray(counts)
        if (counts.ndim != 2 or not np.issubdtype(counts.dtype, np.integer)
                or np.any((counts < 0) | (counts > self.grid.n))):
            raise ValueError("Expected a matrix of valid annual counts")
        horizons = tuple(horizons)
        if (not horizons or len(set(horizons)) != len(horizons)
                or any(not isinstance(h, (int, np.integer)) for h in horizons)
                or min(horizons) < 1 or max(horizons) > counts.shape[1] or chunk < 1):
            raise ValueError("Invalid history lengths or chunk size")
        result = np.empty((len(counts), len(horizons)))
        for start in range(0, len(counts), chunk):
            paths = counts[start:start+chunk]
            predicted = self.grid.prior
            score = np.zeros(len(paths))
            for year in range(max(horizons)):
                k = paths[:, year]
                joint = predicted * self.emission[k]
                density = joint.sum(axis=1)
                if np.any(density <= 0):
                    raise ArithmeticError("Conditional count probability unresolved")
                score += np.log(density) - np.log(self.grid.stationary[k])
                if year+1 in horizons:
                    result[start:start+len(paths), horizons.index(year+1)] = score
                if year+1 < max(horizons):
                    predicted = (joint / density[:, None]) @ self.transition
        return result


def lag_product(centered):
    """Successive centered-count products, reduced over the last axis."""
    return (centered[..., :-1]*centered[..., 1:]).sum(axis=-1)


def permutation_pvalues(counts, draws, rng, chunk=64):
    """One-sided Monte Carlo permutation test with ties included."""
    counts = np.asarray(counts)
    if (counts.ndim != 2 or not np.issubdtype(counts.dtype, np.integer)
            or counts.shape[1] < 2 or draws < 1 or chunk < 1):
        raise ValueError("Expected annual integer counts, at least two years and positive draws")
    result = np.empty(len(counts))
    for start in range(0, len(counts), chunk):
        block = counts[start:start+chunk].astype(np.int64)
        # Multiply centered counts by T so tied statistics compare exactly.
        centered = block*counts.shape[1] - block.sum(axis=1, keepdims=True)
        observed = lag_product(centered)
        replicated = np.broadcast_to(centered[:, None, :], (len(block), draws, counts.shape[1]))
        reordered = rng.permuted(replicated, axis=-1)
        extreme = (lag_product(reordered) >= observed[:, None]).sum(axis=1)
        result[start:start+len(block)] = (1+extreme)/(1+draws)
    return result
