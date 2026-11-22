from itertools import permutations, product
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from bank_stresslab.models import Model, pair_probability
from bank_stresslab.persistence import CountForecast
from bank_stresslab.validation import OrderedLikelihood, permutation_pvalues, simulate_histories


class OrderedCountTests(unittest.TestCase):
    def test_two_year_single_borrower_closed_form(self):
        model = Model(.08, .25)
        phi = .6
        grid = CountForecast(model, 1, 12)
        seq = np.array(list(product((0, 1), repeat=2)))
        q = pair_probability(Model(model.pd, model.rho*phi))
        expected = [1-2*model.pd+q, model.pd-q, model.pd-q, q]
        scores = OrderedLikelihood(grid, phi).log_ratios(seq, (1, 2))
        assert_allclose(scores[:, 0], 0, atol=2e-13)
        assert_allclose(np.exp(scores[:, 1])*grid.stationary[seq].prod(axis=1), expected, atol=2e-13)

    def test_full_joint_law_and_reversibility(self):
        grid = CountForecast(Model(.06, .25), 3, 8)
        obs = np.array(list(product(range(4), repeat=3)))
        alternative = OrderedLikelihood(grid, .6)
        scores = alternative.log_ratios(obs, (3,), chunk=17)
        independent = grid.stationary[obs].prod(axis=1)
        self.assertAlmostEqual(np.dot(np.exp(scores[:, 0]), independent), 1., places=12)
        assert_allclose(scores, alternative.log_ratios(obs[:, ::-1], (3,)), atol=2e-13)
        assert_allclose(OrderedLikelihood(grid, 0.).log_ratios(obs, (2, 3)), 0, atol=2e-13)

    def test_simulated_histories_reproduce_and_prefix_likelihood_agrees(self):
        model = Model(.01, .25)
        paths = simulate_histories(model, 1000, .6, 10, 30, np.random.default_rng(128))
        again = simulate_histories(model, 1000, .6, 10, 30, np.random.default_rng(128))
        assert_array_equal(paths, again)
        likelihood = OrderedLikelihood(CountForecast(model, 1000, 8), .6)
        multiple = likelihood.log_ratios(paths, (10, 3))
        assert_allclose(multiple[:, 1], likelihood.log_ratios(paths[:, :3], (3,))[:, 0])
        with self.assertRaises(ValueError):
            likelihood.log_ratios(paths, (3, 3))
        with self.assertRaises(ValueError):
            likelihood.log_ratios(np.array([[-1, 2]]), (2,))


class PermutationTests(unittest.TestCase):
    def test_full_enumeration_including_ties(self):
        paths = np.array([[3, 3, 1, 0], [5, 5, 5, 5]])

        class AllOtherOrderings:
            def permuted(self, rows, axis):
                return np.stack([np.array(list(permutations(history[0]))[1:]) for history in rows])

        observed = permutation_pvalues(paths, math.factorial(4)-1, AllOtherOrderings())
        expected = []
        for path in paths:
            def statistic(seq):
                centered = np.array(seq)-path.mean()
                return np.dot(centered[:-1], centered[1:])
            expected.append(np.mean([statistic(seq) >= statistic(path) for seq in permutations(path)]))
        assert_allclose(observed, expected)
        self.assertEqual(observed[1], 1)


if __name__ == "__main__":
    unittest.main()
