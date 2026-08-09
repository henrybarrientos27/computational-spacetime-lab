import math
import unittest

from beamline_realism import (
    ensemble_through_wall_probability,
    flux_targets,
    gravity_limited_length_m,
    optimize_realistic_beamline,
    required_signal_rate_for_exposure,
)
from hidden_neutron_lab import (
    NEUTRON_LIFETIME_S,
    asimov_discovery_time,
    coherent_probability,
)


class BeamlineRealismTests(unittest.TestCase):
    def test_zero_spread_matches_monochromatic_result(self):
        length = 1.0
        velocity = 5.0
        tau = 2.0
        interaction_time = length / velocity
        expected = coherent_probability(interaction_time, tau) ** 2
        expected *= math.exp(-2.0 * interaction_time / NEUTRON_LIFETIME_S)
        actual = ensemble_through_wall_probability(
            length, velocity, 0.0, tau, velocity_samples=1, field_samples=1
        )
        self.assertAlmostEqual(actual, expected)

    def test_gravity_limit_hits_aperture(self):
        velocity = 5.0
        spread = 0.05
        aperture = 0.10
        length = gravity_limited_length_m(velocity, spread, aperture)
        slow_velocity = velocity * (1.0 - 3.0 * spread)
        sag = 0.5 * 9.80665 * (length / slow_velocity) ** 2
        self.assertAlmostEqual(sag, aperture)

    def test_faster_beam_has_longer_gravity_limit(self):
        slow = gravity_limited_length_m(5.0, 0.05, 0.10)
        fast = gravity_limited_length_m(500.0, 0.05, 0.10)
        self.assertAlmostEqual(fast / slow, 100.0)

    def test_field_nonuniformity_reduces_pi_pulse_peak(self):
        tau = 1.0
        velocity = 5.0
        length = math.pi * tau * velocity / 2.0
        perfect = ensemble_through_wall_probability(length, velocity, 0.0, tau)
        spread = ensemble_through_wall_probability(
            length, velocity, 0.0, tau, field_sigma_tesla=1e-6
        )
        self.assertLess(spread, perfect)

    def test_optimizer_respects_gravity(self):
        result = optimize_realistic_beamline(10.0, mean_velocity_m_s=5.0)
        self.assertLessEqual(
            result.converter_length_m,
            result.gravity_length_limit_m * (1.0 + 1e-12),
        )
        self.assertLessEqual(result.three_sigma_gravity_sag_m, 0.10 * (1 + 1e-12))

    def test_required_rate_meets_target_exposure(self):
        target = 86400.0
        rate = required_signal_rate_for_exposure(target, 1e-4)
        self.assertLessEqual(asimov_discovery_time(rate, 1e-4), target)
        self.assertGreater(asimov_discovery_time(rate * 0.999, 1e-4), target)

    def test_blind_scan_needs_more_flux_than_known_resonance(self):
        result = optimize_realistic_beamline(100.0, mean_velocity_m_s=5.0)
        targets = flux_targets(result, 365.25)
        self.assertGreater(
            targets["blind_scan_required_incident_flux_hz"],
            targets["known_resonance_required_incident_flux_hz"],
        )


if __name__ == "__main__":
    unittest.main()
