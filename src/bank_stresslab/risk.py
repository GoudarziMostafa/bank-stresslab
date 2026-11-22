"""Upper-tail risk for finitely supported losses in ascending order."""
import numpy as np


def tail_risk(loss, mass, alpha):
    if not 0 < alpha < 1:
        raise ValueError("Confidence level must lie between zero and one")
    loss = np.asarray(loss)
    mass = np.asarray(mass)
    if loss.ndim != 1 or loss.shape != mass.shape or not loss.size:
        raise ValueError("Loss and mass must be nonempty vectors of the same length")
    if not np.all(np.isfinite(loss)) or np.any(loss[1:] < loss[:-1]):
        raise ValueError("Losses must be finite and in ascending order")
    total = mass.sum()
    if not np.isfinite(total) or total <= 0 or np.any(mass < 0):
        raise ValueError("Mass must be nonnegative with a positive finite total")
    if 1-alpha > np.sqrt(np.finfo(float).eps):
        # Compare the CDF directly at ordinary levels, including atoms at alpha.
        i = np.searchsorted(np.cumsum(mass, dtype=float)/total, alpha)
    else:
        survival = np.r_[np.cumsum(mass[:0:-1], dtype=float)[::-1], 0.]
        target = (1-alpha)*total
        # Only allow rounding at the scale of the remaining tail mass.
        i = np.searchsorted(-survival, -np.nextafter(target, np.inf))
    var = loss[i]
    excess = np.dot(loss[i+1:] - var, mass[i+1:]) / total
    return var, var + excess / (1 - alpha)
