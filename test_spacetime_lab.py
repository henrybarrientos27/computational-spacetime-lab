import math
import unittest

from spacetime_lab import (
    C,
    G0,
    alcubierre_bubble,
    alcubierre_variational_bound,
    casimir_energy_density,
    compact_radius_for_kk_gap,
    cruise_mission,
    lorentz_gamma,
    kaluza_klein_mode,
    minimum_alcubierre_for_cavity,
    morris_thorne_wormhole,
    proper_acceleration_trip,
    synthetic_dimension_transfer,
)


class SpacetimeLabTests(unittest.TestCase):
    def test_lorentz_gamma(self):
        self.assertAlmostEqual(lorentz_gamma(0.6), 1.25, places=12)

    def test_cruise_clock_ratio_and_energy(self):
        mission = cruise_mission(1000.0, 0.6, 4.0)
        self.assertAlmostEqual(
            mission.earth_cruise_years / mission.traveler_cruise_years,
            mission.gamma,
        )
        self.assertAlmostEqual(mission.kinetic_energy_j, 250.0 * C**2)

    def test_one_g_alpha_centauri_trip(self):
        trip = proper_acceleration_trip(4.2465, G0)
        self.assertTrue(5.8 < trip.earth_years < 6.1)
        self.assertTrue(3.4 < trip.traveler_years < 3.7)
        self.assertTrue(0.94 < trip.peak_beta < 0.97)

    def test_warp_energy_is_negative_and_scales_as_beta_squared(self):
        slow = alcubierre_bubble(100.0, 10.0, 0.5, 4000)
        fast = alcubierre_bubble(100.0, 10.0, 1.0, 4000)
        self.assertLess(fast.negative_energy_j, 0.0)
        self.assertAlmostEqual(
            fast.negative_energy_j / slow.negative_energy_j, 4.0, places=10
        )
        self.assertGreater(fast.flat_interior_radius_m, 0.0)

    def test_warp_cavity_optimization_respects_flat_radius(self):
        result = minimum_alcubierre_for_cavity(5.0, beta=1.0, samples=20)
        self.assertGreaterEqual(result.bubble_radius_m, 5.0)
        self.assertLess(result.negative_energy_j, 0.0)

    def test_warp_variational_bound(self):
        infinite = alcubierre_variational_bound(20.0, 1.0)
        finite = alcubierre_variational_bound(20.0, 1.0, 100.0)
        self.assertAlmostEqual(infinite.minimum_radial_integral_m, 20.0)
        self.assertAlmostEqual(finite.minimum_radial_integral_m, 25.0)
        self.assertGreater(
            abs(finite.negative_energy_upper_bound_j),
            abs(infinite.negative_energy_upper_bound_j),
        )
        tanh = alcubierre_bubble(40.0, 8.7, 1.0, 4000)
        bound = alcubierre_variational_bound(tanh.flat_interior_radius_m, 1.0)
        self.assertGreaterEqual(
            abs(tanh.negative_energy_j), abs(bound.negative_energy_upper_bound_j)
        )

    def test_wormhole_nec_and_tidal_limit(self):
        result = morris_thorne_wormhole(100.0, 1e-8)
        self.assertLess(result.throat_nec_j_m3, 0.0)
        at_limit = morris_thorne_wormhole(
            100.0, result.maximum_beta_for_one_g_tides
        )
        self.assertAlmostEqual(
            at_limit.lateral_tidal_acceleration_m_s2, G0, places=8
        )
        larger = morris_thorne_wormhole(200.0, 1e-8)
        self.assertAlmostEqual(
            larger.two_sided_proper_volume_nec_integral_j
            / result.two_sided_proper_volume_nec_integral_j,
            2.0,
        )

    def test_casimir_scaling(self):
        micron = casimir_energy_density(1e-6)
        two_micron = casimir_energy_density(2e-6)
        self.assertLess(micron, 0.0)
        self.assertAlmostEqual(micron / two_micron, 16.0, places=12)

    def test_kaluza_klein_gap_and_inverse(self):
        mode = kaluza_klein_mode(30e-6, 1)
        self.assertAlmostEqual(mode.kk_gap_energy_ev, 0.006577566, places=9)
        radius = compact_radius_for_kk_gap(1.5e12)
        self.assertAlmostEqual(radius, 1.3155e-19, delta=1e-23)

    def test_synthetic_dimension_perfect_transfer(self):
        beginning = synthetic_dimension_transfer(9, 0.0)
        middle = synthetic_dimension_transfer(9, math.pi / 4)
        end = synthetic_dimension_transfer(9, math.pi / 2)
        self.assertAlmostEqual(beginning.sender_probability, 1.0)
        self.assertAlmostEqual(sum(middle.probabilities), 1.0)
        self.assertAlmostEqual(middle.mean_synthetic_site, 4.0)
        self.assertAlmostEqual(end.receiver_probability, 1.0)


if __name__ == "__main__":
    unittest.main()
