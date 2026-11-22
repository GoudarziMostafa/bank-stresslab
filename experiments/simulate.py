"""Produce independent paired batch histograms, with cumulative precision checkpoints."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import t

from bank_stresslab.models import Model, match_student
from bank_stresslab.portfolio import PORTFOLIOS
from bank_stresslab.risk import tail_risk
from bank_stresslab.simulation import simulate_block

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261120
BATCHES = 32
BLOCK = 50000


def run_batch(task):
    batch, block = task
    gaussian = Model(.01, .25)
    student = match_student(gaussian, 10)
    return simulate_block(gaussian, student, 1000, BLOCK, [SEED, batch, block])


def half_width(histograms):
    portfolio = PORTFOLIOS[0]
    gaps = []
    for batch in histograms:
        es = [tail_risk(*portfolio.distribution(joint), .999)[1] for joint in batch]
        gaps.append((es[1] - es[0]) * 10000)
    return t.ppf(.975, BATCHES-1) * np.std(gaps, ddof=1) / np.sqrt(BATCHES)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    path = ROOT / "results"
    path.mkdir(exist_ok=True)
    histograms = np.zeros((BATCHES, 2, 2, 1000), dtype=np.int64)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for block in range(4):
            tasks = [(batch, block) for batch in range(BATCHES)]
            for batch, result in enumerate(pool.map(run_batch, tasks)):
                histograms[batch] += result
            if block + 1 in (1, 2, 4):
                scenarios = (block + 1) * BLOCK * BATCHES
                np.savez_compressed(path / f"paired_{scenarios}.npz", counts=histograms,
                                    seed=SEED, scenarios=scenarios)
                width = half_width(histograms)
                print(f"{scenarios:,} scenarios, primary half-width {width:.2f} bp", flush=True)
                if width <= 2:
                    break


if __name__ == "__main__":
    main()
