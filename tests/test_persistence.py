import unittest

import numpy as np
from numpy.testing import assert_allclose

from bank_stresslab.models import Model, pair_probability
from bank_stresslab.persistence import CountForecast


class PersistenceTests(unittest.TestCase):
    def test_single_obligor_forecast_against_bivariate_normal(self):
        # Correlation of the latent borrower draws across years equals rho * phi.
        model = Model(.04, .25)
        grid = CountForecast(model, 1, 16)
        for phi in (0., .6, .9):
            predicted = grid.predict(phi, grid.likelihood(1, exact=True))[0]
            joint = pair_probability(Model(model.pd, model.rho*phi))
            self.assertAlmostEqual(predicted[1], joint/model.pd, places=11)

    def test_exact_counts_pool_to_event_and_stationary_law(self):
        grid = CountForecast(Model(.1, .25), 8, 12)
        exact = grid.count_law.T
        predicted = grid.predict(.65, exact)
        assert_allclose(grid.stationary @ predicted, grid.stationary, atol=2e-13)
        cutoff = 2
        tail_weights = grid.stationary[cutoff+1:]
        pooled = (tail_weights @ predicted[cutoff+1:])/tail_weights.sum()
        event = grid.predict(.65, grid.likelihood(cutoff))[0]
        assert_allclose(event, pooled, atol=2e-13)

    def test_independent_years_and_uninformative_observations(self):
        grid = CountForecast(Model(.01, .25), 100, 16)
        observed = grid.likelihood(4)
        assert_allclose(grid.predict(0., observed)[0], grid.stationary, atol=2e-13)
        assert_allclose(grid.predict(.9, np.ones_like(observed))[0], grid.stationary, atol=2e-13)

    def test_absent_credit_factor_and_impossible_observation(self):
        grid = CountForecast(Model(.03, 0.), 10, 12)
        assert_allclose(grid.predict(.9, grid.likelihood(1))[0], grid.stationary, atol=2e-13)
        with self.assertRaisesRegex(ValueError, "zero"):
            grid.predict(.6, grid.likelihood(10))


if __name__ == "__main__":
    unittest.main()
