import math
import unittest

import numpy as np

from collider_extra_dimension import (
    SimplifiedLikelihood,
    add_radius_and_gap,
    expand_qualifiers,
    first_number,
    normal_cdf,
)


class ParsingTests(unittest.TestCase):
    def test_expand_qualifiers(self):
        table = {
            "qualifiers": {
                "year": [
                    {"group": 0, "colspan": 2, "value": "2017"},
                    {"group": 2, "colspan": 1, "value": "2018"},
                ]
            }
        }
        self.assertEqual(expand_qualifiers(table)[1]["year"], "2017")
        self.assertEqual(expand_qualifiers(table)[2]["year"], "2018")

    def test_first_number(self):
        self.assertEqual(first_number("10.7 TeV"), 10.7)

    def test_normal_cdf(self):
        self.assertAlmostEqual(normal_cdf(0.0), 0.5)
        self.assertAlmostEqual(normal_cdf(1.9599639845), 0.975, places=6)

    def test_add_radius_and_gap_are_reciprocal(self):
        radius_m, gap_ev = add_radius_and_gap(2, 10.0)
        self.assertGreater(radius_m, 0.0)
        self.assertAlmostEqual(radius_m * gap_ev, 1.973269804e-7, places=16)


class LikelihoodTests(unittest.TestCase):
    def test_background_asimov_has_zero_signal_mle(self):
        likelihood = SimplifiedLikelihood(
            background=np.array([100.0]),
            observed=np.array([100.0]),
            covariance=np.array([[25.0]]),
        )
        mu_hat, _ = likelihood.best_fit_mu(np.array([15.0]))
        self.assertLess(mu_hat, 1e-3)

    def test_stronger_signal_has_smaller_expected_cls(self):
        likelihood = SimplifiedLikelihood(
            background=np.array([100.0]),
            observed=np.array([100.0]),
            covariance=np.array([[25.0]]),
        )
        _, weak = likelihood.cls(np.array([5.0]))
        _, strong = likelihood.cls(np.array([25.0]))
        self.assertLess(strong["expected_cls"], weak["expected_cls"])
        self.assertTrue(math.isfinite(strong["q_asimov"]))


if __name__ == "__main__":
    unittest.main()
