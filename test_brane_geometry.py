import unittest

from brane_geometry import (
    CONSTITUENT_QUARK_MASS_EV,
    HBAR_C_EV_M,
    conditional_geometry_limit,
    coupling_5d_noncompact_per_m,
    geometric_coupling_from_mixing,
    mixing_energy_limit_from_probability,
    mixing_from_geometric_coupling,
    separation_lower_bound_5d_noncompact,
    swapping_probability,
)


class BraneGeometryTests(unittest.TestCase):
    def test_probability_inversion_round_trip(self):
        p = 3.1e-11
        epsilon = mixing_energy_limit_from_probability(p, 2e3, 1e5)
        self.assertAlmostEqual(swapping_probability(epsilon, 2e3, 1e5), p)

    def test_stereo_braneworld_benchmark(self):
        epsilon = mixing_energy_limit_from_probability(3.1e-11, 2e3)
        self.assertAlmostEqual(epsilon, 7.874007874e-3, places=11)

    def test_geometric_mixing_round_trip(self):
        coupling = 1e-4
        vector_potential = 2e9
        epsilon = mixing_from_geometric_coupling(coupling, vector_potential)
        self.assertAlmostEqual(
            geometric_coupling_from_mixing(epsilon, vector_potential), coupling
        )

    def test_5d_exponential(self):
        scale = 1e25
        compton_length = HBAR_C_EV_M / CONSTITUENT_QUARK_MASS_EV
        g0 = coupling_5d_noncompact_per_m(0.0, scale)
        g1 = coupling_5d_noncompact_per_m(compton_length, scale)
        self.assertAlmostEqual(g1 / g0, 1.0 / 2.718281828459045)

    def test_separation_inverse(self):
        scale = 1e25
        distance = 2e-15
        coupling = coupling_5d_noncompact_per_m(distance, scale)
        inferred = separation_lower_bound_5d_noncompact(coupling, scale)
        self.assertIsNotNone(inferred)
        self.assertAlmostEqual(inferred, distance, places=27)

    def test_no_distance_reach_when_limit_above_maximum(self):
        result = conditional_geometry_limit(
            probability_limit=3.1e-11,
            effective_detuning_ev=2e3,
            vector_potential_difference_tesla_m=2e9,
            brane_scale_ev=1.22089e28,
        )
        self.assertFalse(result.separation_is_constrained)
        self.assertIsNone(result.conditional_separation_lower_bound_m)


if __name__ == "__main__":
    unittest.main()
