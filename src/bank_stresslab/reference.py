"""Deterministic integration of the exact finite-portfolio conditional law."""
import numpy as np
from scipy.special import ndtr, roots_hermitenorm, roots_legendre
from scipy.stats import binom, chi, norm


def panel_rule(edges, order):
    nodes, weights = roots_legendre(order)
    half = np.diff(edges)[:, None] / 2
    mid = edges[:-1, None] + half
    return (mid + half * nodes).ravel(), (half * weights).ravel()


def mixing_rule(model, order):
    """Integrate the density of X = Phi^{-1}(conditional PD) on a composite grid."""
    if model.df is None and model.rho == 0:
        return np.array([model.pd]), np.ones(1)
    h = model.threshold
    residual = np.sqrt(1 - model.rho)
    spread = np.sqrt(model.rho) / residual
    if model.df is None:
        # Integrate the underlying standard normal directly so small rho is well resolved.
        z, weights = panel_rule(np.linspace(-12, 12, 49), order)
        return ndtr(h / residual + spread * z), weights * norm.pdf(z)
    # S = sqrt(W/df). For fixed S the probit X is normal with this
    # standard deviation and mean h*S/residual. A 1D mixture-density grid
    # avoids recomputing a thousand binomial masses for every (S,Z) pair.
    df = model.df
    upper = chi.isf(1e-16, df) / np.sqrt(df)
    scales, scale_weights = panel_rule(np.linspace(0, upper, 13), order)
    scale_weights *= np.sqrt(df) * chi.pdf(scales * np.sqrt(df), df)
    means = h * scales / residual
    if model.rho == 0:
        return ndtr(means), scale_weights
    # For weak loadings, the normal kernels are narrower than the spacing of
    # their centers. Integrate conditional probabilities directly over (S,Z).
    if spread < np.max(np.abs(np.diff(means))):
        z, normal_weights = roots_hermitenorm(order)
        u = ndtr(means[:, None] + spread*z)
        weights = scale_weights[:, None] * (normal_weights/np.sqrt(2*np.pi))
        return u.ravel(), weights.ravel()
    xmin, xmax = means.min() - 10 * spread, 10 * spread
    panels = int(np.ceil((xmax - xmin) / .4))
    x, weights = panel_rule(np.linspace(xmin, xmax, panels + 1), order)
    density = norm.pdf((x[:, None] - means) / spread) @ scale_weights / spread
    return ndtr(x), weights * density


def reference_distribution(model, n, order):
    u, weights = mixing_rule(model, order)
    k = np.arange(n)
    joint = np.zeros((2, n))
    for start in range(0, len(u), 128):
        section = slice(start, start + 128)
        probability = u[section]
        w = weights[section]
        conditional = binom.pmf(k, n - 1, probability[:, None])
        joint[0] += (w * (1 - probability)) @ conditional
        joint[1] += (w * probability) @ conditional
    moments = [weights @ (u**power) for power in (1, 2, 3)]
    return joint, moments
