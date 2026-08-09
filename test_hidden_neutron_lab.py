import math
import unittest

from hidden_neutron_lab import (
    HBAR_EV_S,
    asimov_discovery_time,
    beamline_result,
    coherent_probability,
    field_detuning_ev,
    medium_swapping_probability,
    optimize_beamline,
    pi_pulse_time,
)


class HiddenNeutronLabTests(unittest.TestCase):
    def test_resonant_pi_pulse_gives_full_swap(self):
        tau = 2.0
        self.assertAlmostEqual(
            coherent_probability(pi_pulse_time(tau), tau), 1.0, places=12
        )

    def test_detuning_suppresses_swap(self):
        tau = 1.0
        time = 0.5
        resonant = coherent_probability(time, tau)
        detuned = coherent_probability(time, tau, field_detuning_ev(1e-6))
        self.assertLess(detuned, resonant)

    def test_medium_probability_high_detuning_limit(self):
        epsilon = 1e-6
        detuning = 2.0
        actual = medium_swapping_probability(epsilon, detuning, 0.0, 0.0)
        expected = 2.0 * epsilon**2 / detuning**2
        self.assertAlmostEqual(actual, expected, places=20)

    def test_two_zone_probability_is_square(self):
        result = beamline_result(
            oscillation_time_s=10.0,
            interaction_time_s=1.0,
            neutron_velocity_m_s=5.0,
            neutron_flux_hz=1e5,
            detection_efficiency=0.5,
            background_rate_hz=0.0,
            field_range_tesla=1e-6,
        )
        self.assertAlmostEqual(
            result.through_wall_probability_worst_case,
            result.single_swap_probability_worst_case**2,
        )

    def test_slower_mixing_requires_longer_search(self):
        fast = optimize_beamline(10.0)
        slow = optimize_beamline(100.0)
        self.assertGreater(slow.total_scan_time_s, fast.total_scan_time_s)
        self.assertTrue(math.isfinite(slow.total_scan_time_s))

    def test_epsilon_tau_conversion(self):
        self.assertAlmostEqual(HBAR_EV_S / 1.0, 6.582119569e-16)

    def test_asimov_time_is_stable_for_tiny_signal(self):
        signal = 1e-20
        background = 1e-4
        expected_small_signal_time = 25.0 * background / signal**2
        actual = asimov_discovery_time(signal, background)
        self.assertAlmostEqual(actual / expected_small_signal_time, 1.0)


if __name__ == "__main__":
    unittest.main()
