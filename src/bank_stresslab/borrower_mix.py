"""Credit-cycle likelihoods with observed or hidden beginning-of-year borrower grades."""
from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import binom, norm

from .reference import panel_rule


@dataclass(frozen=True)
class GradeBook:
    n: int = 1000
    high_sizes: tuple[int, int] = (50, 350)
    p_low: float = .005
    p_high: float = .03
    rho: float = .25

    def __post_init__(self):
        if not (0 <= self.high_sizes[0] < self.high_sizes[1] <= self.n
                and 0 < self.p_low <= self.p_high < 1 and 0 <= self.rho < 1):
            raise ValueError("Invalid grade sizes, default probabilities or asset correlation")

    def probabilities(self, z):
        offset = np.sqrt(self.rho)*z
        residual = np.sqrt(1-self.rho)
        return (ndtr((ndtri(self.p_low)-offset)/residual),
                ndtr((ndtri(self.p_high)-offset)/residual))


def simulate_book(book, phi, mix_persistence, size, years, rng):
    """Return state, low defaults and high defaults for the same simulated histories."""
    if not (0 <= phi < 1 and 0 <= mix_persistence < 1 and min(size, years) > 0):
        raise ValueError("Invalid dynamics or history dimensions")
    state = rng.integers(2, size=size)
    z = rng.standard_normal(size)
    data = np.empty((size, years, 3), dtype=int)
    for year in range(years):
        if year:
            state = state ^ (rng.random(size) < (1-mix_persistence)/2)
            z = phi*z + np.sqrt(1-phi**2)*rng.standard_normal(size)
        n_high = np.take(book.high_sizes, state)
        p_low, p_high = book.probabilities(z)
        data[:, year, 0] = state
        data[:, year, 1] = rng.binomial(book.n-n_high, p_low)
        data[:, year, 2] = rng.binomial(n_high, p_high)
    return data


def _binomial_table(size, pd):
    return binom.pmf(np.arange(size+1)[:, None], size, pd)


class AttributionFilter:
    """Marginalize the hidden cycle for both validators and the hidden mix for totals."""

    def __init__(self, book, phi=.6, mix_persistence=.8, order=8):
        if not (0 <= phi < 1 and 0 <= mix_persistence < 1 and order >= 1):
            raise ValueError("Invalid filter dynamics or order")
        self.book = book
        self.stay = (1+mix_persistence)/2
        self.nodes, widths = panel_rule(np.linspace(-10, 10, 41), order)
        self.prior = widths*norm.pdf(self.nodes)
        scale = np.sqrt(1-phi**2)
        distance = (self.nodes[None, :] - phi*self.nodes[:, None])/scale
        self.transition = norm.pdf(distance)*widths/scale
        p_low, p_high = book.probabilities(self.nodes)
        self.grade_tables = []
        totals = []
        for n_high in book.high_sizes:
            low = _binomial_table(book.n-n_high, p_low)
            high = _binomial_table(n_high, p_high)
            self.grade_tables.append((low, high))
            # Given Z, the two grade counts are independent binomials.
            totals.append(np.stack([np.convolve(l, h) for l, h in zip(low.T, high.T)], axis=1))
        self.totals = np.ascontiguousarray(np.stack(totals, axis=1))
        self.by_state = self.totals @ self.prior
        self.stationary = self.by_state.mean(axis=1)

    def grade_emissions(self, observation):
        states, low_counts, high_counts = observation.T
        emission = np.empty((len(states), len(self.nodes)))
        for state, (low, high) in enumerate(self.grade_tables):
            use = states == state
            emission[use] = low[low_counts[use]]*high[high_counts[use]]
        return emission

    def scores(self, data, chunk=128):
        """Totals-only and grade-aware log likelihood ratios, in that order."""
        data = np.asarray(data)
        if (data.ndim != 3 or data.shape[2] != 3 or data.shape[1] < 1 or chunk < 1
                or not np.issubdtype(data.dtype, np.integer) or np.any(data < 0)
                or np.any(data[:, :, 0] > 1)):
            raise ValueError("Expected a panel of mix state, low defaults and high defaults")
        sizes = np.take(self.book.high_sizes, data[:, :, 0])
        if np.any(data[:, :, 1] > self.book.n-sizes) or np.any(data[:, :, 2] > sizes):
            raise ValueError("Observed defaults exceed starting grade sizes")
        out = np.empty((len(data), 2))
        for start in range(0, len(data), chunk):
            block = data[start:start+chunk]
            out[start:start+len(block)] = self._score_block(block)
        return out

    def _score_block(self, block):
        batch, years, _ = block.shape
        # Totals-only filters relative to the baseline with both processes independent.
        iid_mix = self.prior
        persistent_mix = self.prior[None, None, :]/2
        null_mix = np.full((batch, 2), .5)
        grade_cycle = self.prior
        log_iid_mix, log_persistent_mix, log_null_mix, log_grade = np.zeros((4, batch))
        for year in range(years):
            obs = block[:, year]
            k = obs[:, 1]+obs[:, 2]
            marginal = self.stationary[k]
            emission = self.totals[k] / marginal[:, None, None]
            iid_joint = iid_mix * emission.mean(axis=1)
            both_joint = persistent_mix * emission
            null_joint = null_mix * (self.by_state[k]/marginal[:, None])
            graded = self.grade_emissions(obs)
            grade_joint = grade_cycle * graded
            iid_mass = iid_joint.sum(axis=1)
            both_mass = both_joint.sum(axis=(1, 2))
            null_mass = null_joint.sum(axis=1)
            grade_mass = grade_joint.sum(axis=1)
            log_iid_mix += np.log(iid_mass)
            log_persistent_mix += np.log(both_mass)
            log_null_mix += np.log(null_mass)
            log_grade += np.log(grade_mass) - np.log(graded @ self.prior)
            if year+1 < years:
                iid_mix = (iid_joint/iid_mass[:, None]) @ self.transition
                grade_cycle = (grade_joint/grade_mass[:, None]) @ self.transition
                post = both_joint/both_mass[:, None, None]
                shifted = (post.reshape(2*batch, -1) @ self.transition).reshape(post.shape)
                persistent_mix = self.stay*shifted + (1-self.stay)*shifted[:, ::-1]
                null_post = null_joint/null_mass[:, None]
                null_mix = self.stay*null_post + (1-self.stay)*null_post[:, ::-1]
        log_total = np.logaddexp(log_iid_mix, log_persistent_mix) - np.logaddexp(0, log_null_mix)
        scores = np.column_stack((log_total, log_grade))
        if not np.isfinite(scores).all():
            raise ArithmeticError("Count likelihood unresolved on the factor grid")
        return scores
