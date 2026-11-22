"""Common-PD Gaussian and Student threshold models."""
from dataclasses import dataclass
from math import acos, cos, exp, pi

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm, t


@dataclass(frozen=True)
class Model:
    pd: float
    rho: float
    df: float | None = None

    def __post_init__(self):
        if not (0 < self.pd < 0.5 and 0 <= self.rho < 1):
            raise ValueError("This lower-tail experiment requires 0 < PD < 0.5 and 0 <= rho < 1")
        if self.df is not None and not (np.isfinite(self.df) and self.df > 2):
            raise ValueError("Student degrees of freedom must be finite and exceed two")

    @property
    def threshold(self):
        return norm.ppf(self.pd) if self.df is None else t.ppf(self.pd, self.df)

    @property
    def name(self):
        return "Gaussian" if self.df is None else "Student-t"


def pair_probability(model):
    """Bivariate lower-tail probability from the elliptical radial survival law."""
    h = model.threshold
    df = model.df

    def radial_survival(angle):
        radius_squared = (h / cos(angle)) ** 2
        if df is None:
            return exp(-radius_squared / 2)
        return (1 + radius_squared / df) ** (-df / 2)

    return quad(radial_survival, acos(model.rho) / 2, pi / 2, epsabs=1e-14)[0] / pi


def match_student(gaussian, df):
    if gaussian.df is not None:
        raise ValueError("A Gaussian benchmark is required")
    target = pair_probability(gaussian)

    def residual(rho):
        return pair_probability(Model(gaussian.pd, rho, df)) - target

    if residual(0) > 0:
        raise ValueError("Student pair-default floor exceeds the target")
    rho = brentq(residual, 0, np.nextafter(1., 0.), xtol=1e-14)
    return Model(gaussian.pd, rho, df)


def fit_annual_counts(active, defaults):
    """Fit an equal-year PD, distinct-pair probability and nonnegative Gaussian AR(1) factor."""
    n, d = np.asarray(active, dtype=float), np.asarray(defaults, dtype=float)
    if n.ndim != 1 or n.shape != d.shape or n.size < 2:
        raise ValueError("Expected consecutive annual risk sets and default counts")
    if not np.all((n >= 2) & (d >= 0) & (d <= n) & (n == np.floor(n)) & (d == np.floor(d))
                  & np.isfinite(n) & np.isfinite(d)):
        raise ValueError("Invalid annual counts")
    rate = d / n
    pd = rate.mean()
    pair = np.mean(rate * ((d - 1) / (n - 1)))

    def gaussian_pair(rho):
        if rho == 0:
            return pd**2
        return pair_probability(Model(pd, rho))

    if not pd**2 <= pair < pd:
        raise ValueError("Observed pair probability is outside the Gaussian factor model range")
    rho = brentq(lambda r: gaussian_pair(r) - pair, 0., np.nextafter(1., 0.))
    # With independent borrower shocks across years, latent lag correlation is rho * phi.
    centered = rate - pd
    lag_cov = np.mean(centered[1:] * centered[:-1])
    phi = None
    if rho > 0 and 0 <= lag_cov < pair - pd**2:
        phi = brentq(lambda v: gaussian_pair(rho * v) - pd**2 - lag_cov, 0., 1.)
    return Model(pd, rho), pair, lag_cov, phi
