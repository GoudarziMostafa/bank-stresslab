"""Paired scenario simulation using borrower-level independent normals."""
import numpy as np


def simulate_block(gaussian, student, n, size, seed, chunk=512):
    rng = np.random.default_rng(seed)
    counts = np.zeros((2, 2, n), dtype=np.int64)
    h_g, h_t = gaussian.threshold, student.threshold
    for offset in range(0, size, chunk):
        m = min(chunk, size - offset)
        z = rng.standard_normal((m, 1))
        scale = np.sqrt(rng.chisquare(student.df, size=(m, 1)) / student.df)
        errors = rng.standard_normal((m, n))
        for model, boundary, output in zip(
            (gaussian, student), (h_g, h_t * scale), counts, strict=True
        ):
            # Threshold the same borrower draws in both models, without evaluating a CDF.
            cutoff = (boundary - np.sqrt(model.rho) * z) / np.sqrt(1 - model.rho)
            defaults = errors <= cutoff
            indices = defaults[:, 0] * n + defaults[:, 1:].sum(axis=1)
            output += np.bincount(indices, minlength=2*n).reshape(2, n)
    return counts
