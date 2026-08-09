#!/usr/bin/env python3
"""Add beam spread, field nonuniformity, and gravity to a portal beamline.

The apparatus is still conditional on nonzero neutron/hidden-neutron mixing.
This module asks a narrower engineering question: if the two-state Hamiltonian
is real, can a horizontal two-converter beamline preserve enough coherent
amplitude to test it?
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import NormalDist

from hidden_neutron_lab import (
    NEUTRON_LIFETIME_S,
    asimov_discovery_time,
    coherent_probability,
    field_detuning_ev,
)


STANDARD_GRAVITY_M_S2 = 9.80665


def normal_quantile_samples(mean: float, sigma: float, count: int) -> list[float]:
    if count < 1 or sigma < 0:
        raise ValueError("sample count must be positive and sigma nonnegative")
    if sigma == 0:
        return [mean]
    normal = NormalDist()
    samples = [
        mean + sigma * normal.inv_cdf((index + 0.5) / count)
        for index in range(count)
    ]
    positive = [sample for sample in samples if sample > 0]
    if not positive:
        raise ValueError("velocity distribution has no positive samples")
    return positive


def gravity_limited_length_m(
    mean_velocity_m_s: float,
    fractional_velocity_sigma: float,
    vertical_half_aperture_m: float,
    sigma_safety: float = 3.0,
    gravity_m_s2: float = STANDARD_GRAVITY_M_S2,
) -> float:
    """Maximum horizontal zone length before a slow neutron exits aperture."""
    if mean_velocity_m_s <= 0 or not 0 <= fractional_velocity_sigma < 1:
        raise ValueError("velocity must be positive and fractional spread in [0,1)")
    if vertical_half_aperture_m <= 0 or sigma_safety < 0 or gravity_m_s2 <= 0:
        raise ValueError("aperture and gravity must be positive")
    slow_velocity = mean_velocity_m_s * (
        1.0 - sigma_safety * fractional_velocity_sigma
    )
    if slow_velocity <= 0:
        raise ValueError("sigma safety range includes zero velocity")
    return slow_velocity * math.sqrt(
        2.0 * vertical_half_aperture_m / gravity_m_s2
    )


def ensemble_through_wall_probability(
    converter_length_m: float,
    mean_velocity_m_s: float,
    fractional_velocity_sigma: float,
    oscillation_time_s: float,
    field_offset_tesla: float = 0.0,
    field_sigma_tesla: float = 0.0,
    velocity_samples: int = 17,
    field_samples: int = 11,
) -> float:
    """Average P(n->n')P(n'->n), including decay across both zones.

    The two zones have the same statistical field quality but independent local
    field errors. Velocity is common to both zones, so the average is
    E_v[(E_B P(v,B))^2], not (E_{v,B} P)^2.
    """
    if converter_length_m <= 0 or mean_velocity_m_s <= 0:
        raise ValueError("length and velocity must be positive")
    if not 0 <= fractional_velocity_sigma < 1 or field_sigma_tesla < 0:
        raise ValueError("spreads must be nonnegative")

    velocities = normal_quantile_samples(
        mean_velocity_m_s,
        mean_velocity_m_s * fractional_velocity_sigma,
        velocity_samples,
    )
    field_errors = normal_quantile_samples(0.0, field_sigma_tesla, field_samples)
    total = 0.0
    for velocity in velocities:
        interaction_time = converter_length_m / velocity
        mean_single = sum(
            coherent_probability(
                interaction_time,
                oscillation_time_s,
                field_detuning_ev(field_offset_tesla + field_error),
            )
            for field_error in field_errors
        ) / len(field_errors)
        survival = math.exp(-2.0 * interaction_time / NEUTRON_LIFETIME_S)
        total += mean_single**2 * survival
    return total / len(velocities)


def ensemble_field_fwhm_tesla(
    converter_length_m: float,
    mean_velocity_m_s: float,
    fractional_velocity_sigma: float,
    oscillation_time_s: float,
    field_sigma_tesla: float,
) -> float:
    peak = ensemble_through_wall_probability(
        converter_length_m,
        mean_velocity_m_s,
        fractional_velocity_sigma,
        oscillation_time_s,
        field_sigma_tesla=field_sigma_tesla,
    )
    if peak <= 1e-30:
        raise ValueError("ensemble conversion peak is too small")
    target = peak / 2.0
    interaction_time = converter_length_m / mean_velocity_m_s
    high = max(1e-12, field_sigma_tesla, 1e-8 / max(interaction_time, 1e-9))
    for _ in range(80):
        value = ensemble_through_wall_probability(
            converter_length_m,
            mean_velocity_m_s,
            fractional_velocity_sigma,
            oscillation_time_s,
            field_offset_tesla=high,
            field_sigma_tesla=field_sigma_tesla,
        )
        if value <= target:
            break
        high *= 2.0
    else:
        raise RuntimeError("failed to bracket ensemble FWHM")

    low = 0.0
    for _ in range(55):
        middle = 0.5 * (low + high)
        value = ensemble_through_wall_probability(
            converter_length_m,
            mean_velocity_m_s,
            fractional_velocity_sigma,
            oscillation_time_s,
            field_offset_tesla=middle,
            field_sigma_tesla=field_sigma_tesla,
        )
        if value > target:
            low = middle
        else:
            high = middle
    return low + high


@dataclass(frozen=True)
class RealisticBeamlineResult:
    oscillation_time_s: float
    mean_velocity_m_s: float
    fractional_velocity_sigma: float
    field_sigma_microtesla: float
    vertical_half_aperture_m: float
    gravity_length_limit_m: float
    converter_length_m: float
    three_sigma_gravity_sag_m: float
    mean_interaction_time_s: float
    ensemble_peak_through_wall_probability: float
    ensemble_worst_sampled_probability: float
    field_fwhm_microtesla: float
    field_step_microtesla: float
    scan_settings: int
    signal_rate_hz: float
    dwell_time_s: float
    total_scan_time_s: float
    known_resonance_time_s: float


def evaluate_beamline(
    oscillation_time_s: float,
    converter_length_m: float,
    mean_velocity_m_s: float,
    fractional_velocity_sigma: float,
    field_sigma_microtesla: float,
    vertical_half_aperture_m: float,
    neutron_flux_hz: float,
    detection_efficiency: float,
    background_rate_hz: float,
    field_range_microtesla: float,
    oversampling: float = 2.0,
    field_settling_time_s: float = 2.0,
) -> RealisticBeamlineResult:
    gravity_limit = gravity_limited_length_m(
        mean_velocity_m_s,
        fractional_velocity_sigma,
        vertical_half_aperture_m,
    )
    if converter_length_m > gravity_limit * (1.0 + 1e-12):
        raise ValueError("converter exceeds the three-sigma gravity aperture")
    if not 0 < detection_efficiency <= 1 or oversampling < 1:
        raise ValueError("invalid detector efficiency or oversampling")

    field_sigma_tesla = field_sigma_microtesla * 1e-6
    fwhm_tesla = ensemble_field_fwhm_tesla(
        converter_length_m,
        mean_velocity_m_s,
        fractional_velocity_sigma,
        oscillation_time_s,
        field_sigma_tesla,
    )
    step_tesla = fwhm_tesla / oversampling
    settings = max(1, math.ceil(field_range_microtesla * 1e-6 / step_tesla))
    peak_probability = ensemble_through_wall_probability(
        converter_length_m,
        mean_velocity_m_s,
        fractional_velocity_sigma,
        oscillation_time_s,
        field_sigma_tesla=field_sigma_tesla,
    )
    worst_probability = ensemble_through_wall_probability(
        converter_length_m,
        mean_velocity_m_s,
        fractional_velocity_sigma,
        oscillation_time_s,
        field_offset_tesla=step_tesla / 2.0,
        field_sigma_tesla=field_sigma_tesla,
    )
    peak_rate = neutron_flux_hz * detection_efficiency * peak_probability
    signal_rate = neutron_flux_hz * detection_efficiency * worst_probability
    known_time = asimov_discovery_time(peak_rate, background_rate_hz)
    dwell = asimov_discovery_time(signal_rate, background_rate_hz)
    slow_velocity = mean_velocity_m_s * (1.0 - 3.0 * fractional_velocity_sigma)
    sag = 0.5 * STANDARD_GRAVITY_M_S2 * (
        converter_length_m / slow_velocity
    ) ** 2
    return RealisticBeamlineResult(
        oscillation_time_s=oscillation_time_s,
        mean_velocity_m_s=mean_velocity_m_s,
        fractional_velocity_sigma=fractional_velocity_sigma,
        field_sigma_microtesla=field_sigma_microtesla,
        vertical_half_aperture_m=vertical_half_aperture_m,
        gravity_length_limit_m=gravity_limit,
        converter_length_m=converter_length_m,
        three_sigma_gravity_sag_m=sag,
        mean_interaction_time_s=converter_length_m / mean_velocity_m_s,
        ensemble_peak_through_wall_probability=peak_probability,
        ensemble_worst_sampled_probability=worst_probability,
        field_fwhm_microtesla=fwhm_tesla * 1e6,
        field_step_microtesla=step_tesla * 1e6,
        scan_settings=settings,
        signal_rate_hz=signal_rate,
        dwell_time_s=dwell,
        total_scan_time_s=settings * (dwell + field_settling_time_s),
        known_resonance_time_s=known_time,
    )


def optimize_realistic_beamline(
    oscillation_time_s: float,
    mean_velocity_m_s: float = 5.0,
    fractional_velocity_sigma: float = 0.05,
    field_sigma_microtesla: float = 0.01,
    vertical_half_aperture_m: float = 0.10,
    engineering_length_limit_m: float = 50.0,
    neutron_flux_hz: float = 5e5,
    detection_efficiency: float = 0.30,
    background_rate_hz: float = 1e-4,
    field_min_microtesla: float = 50.0,
    field_max_microtesla: float = 1100.0,
) -> RealisticBeamlineResult:
    gravity_limit = gravity_limited_length_m(
        mean_velocity_m_s,
        fractional_velocity_sigma,
        vertical_half_aperture_m,
    )
    max_length = min(engineering_length_limit_m, gravity_limit)
    min_length = min(0.01, max_length / 100.0)
    best: RealisticBeamlineResult | None = None
    for index in range(90):
        fraction = index / 89.0
        length = math.exp(
            math.log(min_length) + fraction * math.log(max_length / min_length)
        )
        try:
            candidate = evaluate_beamline(
                oscillation_time_s,
                length,
                mean_velocity_m_s,
                fractional_velocity_sigma,
                field_sigma_microtesla,
                vertical_half_aperture_m,
                neutron_flux_hz,
                detection_efficiency,
                background_rate_hz,
                field_max_microtesla - field_min_microtesla,
            )
        except ValueError:
            continue
        if best is None or candidate.total_scan_time_s < best.total_scan_time_s:
            best = candidate
    if best is None:
        raise RuntimeError("no valid realistic beamline point")
    return best


def required_signal_rate_for_exposure(
    exposure_time_s: float,
    background_rate_hz: float,
) -> float:
    """Smallest signal rate satisfying the lab's discovery rule in a time."""
    if exposure_time_s <= 0 or background_rate_hz < 0:
        raise ValueError("exposure must be positive and background nonnegative")
    low = 0.0
    high = max(1.0 / exposure_time_s, background_rate_hz, 1e-30)
    for _ in range(100):
        if asimov_discovery_time(high, background_rate_hz) <= exposure_time_s:
            break
        high *= 10.0
    else:
        raise RuntimeError("failed to bracket required signal rate")
    for _ in range(100):
        middle = 0.5 * (low + high)
        if asimov_discovery_time(middle, background_rate_hz) > exposure_time_s:
            low = middle
        else:
            high = middle
    return high


def flux_targets(
    result: RealisticBeamlineResult,
    target_days: float,
    detection_efficiency: float = 0.30,
    background_rate_hz: float = 1e-4,
    field_settling_time_s: float = 2.0,
) -> dict[str, float]:
    """Incident rates needed for targeted and blind tests in a fixed time."""
    if target_days <= 0 or not 0 < detection_efficiency <= 1:
        raise ValueError("target days and efficiency must be positive")
    total_time = target_days * 86400.0
    known_signal = required_signal_rate_for_exposure(
        total_time, background_rate_hz
    )
    per_setting_time = total_time / result.scan_settings - field_settling_time_s
    if per_setting_time <= 0:
        blind_flux = math.inf
    else:
        blind_signal = required_signal_rate_for_exposure(
            per_setting_time, background_rate_hz
        )
        blind_flux = blind_signal / (
            detection_efficiency * result.ensemble_worst_sampled_probability
        )
    known_flux = known_signal / (
        detection_efficiency * result.ensemble_peak_through_wall_probability
    )
    return {
        "target_days": target_days,
        "known_resonance_required_incident_flux_hz": known_flux,
        "blind_scan_required_incident_flux_hz": blind_flux,
        "known_resonance_required_beam_rejection": (
            background_rate_hz / known_flux
        ),
        "blind_scan_required_beam_rejection": (
            0.0 if not math.isfinite(blind_flux) else background_rate_hz / blind_flux
        ),
    }


def write_sweep(path: Path) -> None:
    rows = []
    for velocity in (5.0, 50.0, 500.0, 1000.0):
        for tau in (1.0, 10.0, 100.0):
            rows.append(asdict(optimize_realistic_beamline(tau, velocity)))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    optimize = subparsers.add_parser("optimize")
    optimize.add_argument("--tau", type=float, required=True)
    optimize.add_argument("--velocity", type=float, default=5.0)
    optimize.add_argument("--velocity-spread", type=float, default=0.05)
    optimize.add_argument("--field-sigma-microtesla", type=float, default=0.01)
    optimize.add_argument("--aperture", type=float, default=0.10)
    optimize.add_argument("--length-limit", type=float, default=50.0)
    optimize.add_argument("--flux", type=float, default=5e5)

    sweep = subparsers.add_parser("sweep")
    sweep.add_argument(
        "--output", type=Path, default=Path("realistic-beamline-scan.csv")
    )

    target = subparsers.add_parser("target", help="calculate facility-rate targets")
    target.add_argument("--tau", type=float, required=True)
    target.add_argument("--velocity", type=float, default=5.0)
    target.add_argument("--target-days", type=float, default=365.25)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "optimize":
        result = optimize_realistic_beamline(
            oscillation_time_s=args.tau,
            mean_velocity_m_s=args.velocity,
            fractional_velocity_sigma=args.velocity_spread,
            field_sigma_microtesla=args.field_sigma_microtesla,
            vertical_half_aperture_m=args.aperture,
            engineering_length_limit_m=args.length_limit,
            neutron_flux_hz=args.flux,
        )
        print(json.dumps(asdict(result), indent=2))
    elif args.command == "sweep":
        write_sweep(args.output)
        print(json.dumps({"output": str(args.output)}, indent=2))
    elif args.command == "target":
        result = optimize_realistic_beamline(args.tau, args.velocity)
        output = asdict(result)
        output.update(flux_targets(result, args.target_days))
        print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
