import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.stats import binom

from bank_stresslab.models import Model, match_student, pair_probability
from bank_stresslab.portfolio import PORTFOLIOS, Portfolio
from bank_stresslab.reference import mixing_rule, reference_distribution
from bank_stresslab.risk import tail_risk
from bank_stresslab.simulation import simulate_block


class CalibrationTests(unittest.TestCase):
    def test_matched_pair_matches_independent_conditioning_benchmark(self):
        gaussian = Model(.01, .25)
        student = match_student(gaussian, 10)
        self.assertAlmostEqual(pair_probability(gaussian), .0004375151261025, places=12)
        self.assertAlmostEqual(student.rho, .03949298796509, places=9)
        self.assertAlmostEqual(pair_probability(student), pair_probability(gaussian), places=13)
        for model in (gaussian, student):
            u, w = mixing_rule(model, 24)
            self.assertAlmostEqual(w @ u, .01, places=11)
            self.assertAlmostEqual(w @ (u**2), pair_probability(model), places=11)

    def test_near_zero_student_loading_preserves_mass_and_calibration(self):
        for rho in (1e-8, .0001):
            model = Model(.01, rho, 10)
            u, w = mixing_rule(model, 32)
            assert_allclose([w.sum(), w @ u, w @ (u*u)],
                            [1, model.pd, pair_probability(model)], rtol=0, atol=1e-11)

    def test_shared_student_scale_and_feasibility(self):
        model = Model(.01, 0, 5)
        self.assertGreater(pair_probability(model), model.pd**2)
        u, w = mixing_rule(model, 16)
        self.assertAlmostEqual(w @ (u*u), pair_probability(model), places=11)
        with self.assertRaisesRegex(ValueError, "floor"):
            match_student(Model(.01, .25), 5)


class DistributionTests(unittest.TestCase):
    def test_independent_gaussian_exact_binomial(self):
        n = 17
        joint, _ = reference_distribution(Model(.1, 0), n, 12)
        _, mass = Portfolio("equal", n, 1).distribution(joint)
        assert_allclose(mass, binom.pmf(np.arange(n+1), n, .1), atol=2e-15)

    def test_exchangeable_split_against_direct_binomial(self):
        model = match_student(Model(.01, .25), 10)
        n = 41
        u, w = mixing_rule(model, 24)
        direct = w @ binom.pmf(np.arange(n+1), n, u[:, None])
        joint, _ = reference_distribution(model, n, 24)
        _, mass = Portfolio("equal", n, 1).distribution(joint)
        assert_allclose(mass, direct, atol=3e-14)

    def test_one_borrower_and_atom_at_var(self):
        joint, _ = reference_distribution(Model(.01, .25), 1, 16)
        loss, mass = Portfolio("single", 1, 1).distribution(joint)
        assert_allclose(tail_risk(loss, mass, .98), [0, .225], atol=1e-12)
        assert_allclose(tail_risk(loss, mass, .999), [.45, .45], atol=1e-12)
        assert_allclose(tail_risk(loss, np.array([99, 1]), .99), [0, .45])

    def test_tiny_tail_atom_and_invalid_mass(self):
        loss = np.array([0., 1.])
        assert_allclose(tail_risk(loss, [1-2**-51, 2**-51], 1-2**-52), [1, 1])
        assert_allclose(tail_risk(loss, [0., 1.], 1e-16), [1, 1])
        for weights in ([93, 7], [93., 7.], [.93, .07]):
            self.assertEqual(tail_risk(loss, weights, .93)[0], 0)
        with self.assertRaises(ValueError):
            tail_risk(loss, [0., 0.], .99)

    def test_tied_loss_units_and_normalization(self):
        p = Portfolio("overlap", 4, 2)
        _, counts = p.distribution(np.array([[4, 3, 2, 1], [1, 2, 3, 4]]))
        assert_array_equal(counts, [4, 3, 3, 3, 3, 4])
        for port in PORTFOLIOS:
            el, variance = port.moments(.01, .00043751512615615)
            self.assertAlmostEqual(el, .0045)
            target = 83.835085873 if port.large_units == 1 else 94.48857156
            self.assertAlmostEqual(np.sqrt(variance)*10000, target, places=6)


class SimulationTests(unittest.TestCase):
    def test_reproducibility_and_marginal_laws(self):
        gaussian = Model(.08, .25)
        student = match_student(gaussian, 10)
        args = (gaussian, student, 12, 30000, [7, 2])
        counts = simulate_block(*args)
        assert_array_equal(counts, simulate_block(*args))
        for model, observed in zip((gaussian, student), counts, strict=True):
            expected, _ = reference_distribution(model, 12, 24)
            self.assertEqual(observed.sum(), 30000)
            deviations = abs(observed - 30000*expected)
            self.assertTrue(np.all(deviations < 6*np.sqrt(30000*expected*(1-expected)) + 2))


if __name__ == "__main__":
    unittest.main()
