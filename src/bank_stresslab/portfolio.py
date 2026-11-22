"""Integer exposure units preserve exact ties in the loss distribution."""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Portfolio:
    name: str
    n: int
    large_units: int
    lgd: float = 0.45

    def __post_init__(self):
        if self.n < 1 or self.large_units < 1 or not 0 < self.lgd <= 1:
            raise ValueError("Positive borrower count, exposure units and LGD are required")

    @property
    def total_units(self):
        return self.n - 1 + self.large_units

    def moments(self, pd, pair):
        hhi = (self.large_units**2 + self.n - 1) / self.total_units**2
        return self.lgd * pd, self.lgd**2 * ((pair - pd**2) + (pd - pair) * hhi)

    def distribution(self, joint):
        """Rows give joint mass or counts for designated survivor/default and k others."""
        if joint.shape != (2, self.n):
            raise ValueError("Expected two rows and one column per possible other-default count")
        mass = np.zeros(self.total_units + 1, dtype=joint.dtype)
        mass[:self.n] += joint[0]
        mass[self.large_units:] += joint[1]
        loss = self.lgd * np.arange(len(mass)) / self.total_units
        return loss, mass


PORTFOLIOS = (Portfolio("granular", 1000, 1), Portfolio("concentrated", 1000, 111))
