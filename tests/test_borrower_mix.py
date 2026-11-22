from itertools import product
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.stats import binom

from bank_stresslab.borrower_mix import AttributionFilter, GradeBook, simulate_book


class BorrowerMixTests(unittest.TestCase):
    def setUp(self):
        self.book = GradeBook(n=7, high_sizes=(1, 5), p_low=.04, p_high=.2, rho=.25)
        self.engine = AttributionFilter(self.book, order=4)

    def test_annual_distribution_and_history_reversal(self):
        scenarios = [simulate_book(self.book, phi, mix, 50, 6, np.random.default_rng(22))
                     for phi, mix in product((0., .6), (0., .8))]
        for data in scenarios:
            assert_array_equal(data[:, 0], scenarios[0][:, 0])
            assert_allclose(self.engine.scores(data[:, :1]), 0, atol=1e-12)
            assert_allclose(self.engine.scores(data), self.engine.scores(data[:, ::-1]), atol=1e-10)

    def test_aggregate_emissions_by_explicit_grade_count_sum(self):
        book, eng = self.book, self.engine
        p_low, p_high = book.probabilities(eng.nodes)
        for state, high in enumerate(book.high_sizes):
            for k in range(book.n+1):
                mass = sum(binom.pmf(h, high, p_high)*binom.pmf(k-h, book.n-high, p_low)
                           for h in range(high+1))
                assert_allclose(eng.totals[k, state], mass, atol=1e-13)

    def test_two_year_likelihood_by_enumerating_mix_and_factor_pairs(self):
        eng = self.engine
        obs = np.array([[[0, 2, 1], [1, 0, 1]]])
        total = [3, 1]
        prior_pair = eng.prior[:, None] * eng.prior[None, :]
        cycle_pair = eng.prior[:, None] * eng.transition
        model_likelihoods = []
        for factors in (prior_pair, cycle_pair):
            marginal = 0.
            for stay in (.5, eng.stay):
                for a, b in product((0, 1), repeat=2):
                    emissions = np.outer(eng.totals[total[0], a], eng.totals[total[1], b])
                    # Half for each candidate mix law and half for each initial mix state.
                    marginal += .25*(stay if a==b else 1-stay)*np.sum(factors*emissions)
            graded = np.outer(eng.grade_emissions(obs[:, 0])[0], eng.grade_emissions(obs[:, 1])[0])
            model_likelihoods.append([marginal, np.sum(factors*graded)])
        expected = np.log(np.array(model_likelihoods[1])/model_likelihoods[0])
        assert_allclose(eng.scores(obs)[0], expected, atol=2e-12)

    def test_homogeneous_grades_and_unobserved_grade_labels(self):
        book = GradeBook(n=7, high_sizes=(1, 5), p_low=.1, p_high=.1)
        engine = AttributionFilter(book, order=4)
        data = simulate_book(book, .6, .8, 40, 10, np.random.default_rng(12))
        total, grade = engine.scores(data).T
        assert_allclose(total, grade, atol=1e-12)
        two_versions = np.array([[[0, 2, 1], [0, 2, 0]], [[1, 0, 3], [1, 1, 1]]])
        scores = self.engine.scores(two_versions)
        self.assertAlmostEqual(scores[0, 0], scores[1, 0])
        self.assertNotAlmostEqual(scores[0, 1], scores[1, 1])
        with self.assertRaises(ValueError):
            self.engine.scores(np.array([[[0, 0, 2]]]))
