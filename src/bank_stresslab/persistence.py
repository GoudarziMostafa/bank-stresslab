"""One-step forecasts from a stationary Gaussian credit factor and a count observation."""
import numpy as np
from scipy.special import ndtr
from scipy.stats import binom, norm

from .reference import panel_rule


class CountForecast:
    """Quadrature for the annual default count in a fixed-composition portfolio."""

    def __init__(self, model, borrowers, order):
        if model.df is not None or borrowers < 1 or order < 1:
            raise ValueError("A Gaussian model, positive portfolio size and quadrature order are required")
        self.n = borrowers
        self.nodes, self.widths = panel_rule(np.linspace(-12, 12, 49), order)
        self.prior = self.widths * norm.pdf(self.nodes)
        self.pd = ndtr((model.threshold - np.sqrt(model.rho)*self.nodes) / np.sqrt(1-model.rho))
        self.count_law = binom.pmf(np.arange(borrowers+1), borrowers, self.pd[:, None])
        self.stationary = self.prior @ self.count_law

    def likelihood(self, cutoff, exact=False):
        """Probability of the given count, or a count strictly above the cutoff."""
        if not 0 <= cutoff <= self.n:
            raise ValueError("Count must be in the portfolio's support")
        if exact:
            return binom.pmf(cutoff, self.n, self.pd)
        return binom.sf(cutoff, self.n, self.pd)

    def predict(self, persistence, likelihoods):
        """Each row describes P(observation | current factor) on the quadrature grid."""
        if not 0 <= persistence < 1:
            raise ValueError("Persistence must be in [0, 1)")
        likelihoods = np.atleast_2d(likelihoods)
        if likelihoods.shape[1] != len(self.nodes) or not np.all((0 <= likelihoods) & (likelihoods <= 1)):
            raise ValueError("Each likelihood must give probabilities on the current factor grid")
        weighted = likelihoods * self.prior
        evidence = weighted.sum(axis=1)
        if np.any(evidence == 0):
            raise ValueError("The observation has zero numerical probability")
        posterior = weighted / evidence[:, None]
        # Integrate the transition density from today's posterior to tomorrow's factor.
        scale = np.sqrt(1 - persistence**2)
        transition = norm.pdf((self.nodes[None, :] - persistence*self.nodes[:, None])/scale)
        factor_mass = (posterior @ transition) * self.widths / scale
        predicted = factor_mass @ self.count_law
        if not np.allclose(predicted.sum(axis=1), 1, atol=1e-10, rtol=0):
            raise ArithmeticError("Predictive mass is unresolved; increase quadrature order")
        return predicted
