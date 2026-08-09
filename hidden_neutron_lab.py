#!/usr/bin/env python3
"""Model reversible neutron/hidden-neutron transfer and a wall experiment.

The model is the optimistic two-state Hamiltonian used in hidden-neutron and
braneworld searches. It can design a falsifiable neutron-shining-through-a-wall
experiment; it does not assume that a hidden state has already been observed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


HBAR_EV_S = 6.582119569e-16
NEUTRON_MAGNETIC_MOMENT_EV_T = 6.030774e-8
NEUTRON_LIFETIME_S = 878.4


def coherent_probability(
    interaction_time_s: float,
    oscillation_time_s: float,
    detuning_ev: float = 0.0,
) -> float:
    """Probability for n -> n' under a constant two-state Hamiltonian.

    The off-diagonal mixing energy is epsilon = hbar/tau. At resonance the
    result reduces to sin^2(t/tau).
    """
    if interaction_time_s < 0 or oscillation_time_s <= 0:
        raise ValueError("times must be nonnegative and tau must be positive")
    epsilon_ev = HBAR_EV_S / oscillation_time_s
    splitting_ev = math.hypot(detuning_ev, 2.0 * epsilon_ev)
    amplitude = 4.0 * epsilon_ev**2 / splitting_ev**2
    phase = splitting_ev * interaction_time_s / (2.0 * HBAR_EV_S)
    return max(0.0, min(1.0, amplitude * math.sin(phase) ** 2))


def medium_swapping_probability(
    epsilon_ev: float,
    delta_energy_ev: float,
    fermi_potential_ev: float,
    collision_rate_hz: float,
) -> float:
    """Collision-averaged swap probability per projection from STEREO."""
    if epsilon_ev < 0 or collision_rate_hz < 0:
        raise ValueError("epsilon and collision rate must be nonnegative")
    denominator = (
        (delta_energy_ev + fermi_potential_ev) ** 2
        + 4.0 * epsilon_ev**2
        + (HBAR_EV_S * collision_rate_hz / 2.0) ** 2
    )
    return 2.0 * epsilon_ev**2 / denominator


def pi_pulse_time(oscillation_time_s: float) -> float:
    return math.pi * oscillation_time_s / 2.0


def resonance_field_tesla(delta_mass_ev: float) -> float:
    return delta_mass_ev / NEUTRON_MAGNETIC_MOMENT_EV_T


def field_detuning_ev(field_offset_tesla: float) -> float:
    return NEUTRON_MAGNETIC_MOMENT_EV_T * field_offset_tesla


def field_fwhm_tesla(interaction_time_s: float, oscillation_time_s: float) -> float:
    """FWHM of the central conversion peak as magnetic-field detuning."""
    peak = coherent_probability(interaction_time_s, oscillation_time_s, 0.0)
    if peak <= 1e-30:
        raise ValueError("interaction time lies at a conversion node")
    target = peak / 2.0
    low = 0.0
    high = max(
        HBAR_EV_S / max(interaction_time_s, 1e-30),
        HBAR_EV_S / oscillation_time_s,
    ) / NEUTRON_MAGNETIC_MOMENT_EV_T
    for _ in range(100):
        probability = coherent_probability(
            interaction_time_s,
            oscillation_time_s,
            field_detuning_ev(high),
        )
        if probability <= target:
            break
        high *= 2.0
    else:
        raise RuntimeError("failed to bracket resonance half maximum")

    for _ in range(100):
        middle = 0.5 * (low + high)
        probability = coherent_probability(
            interaction_time_s,
            oscillation_time_s,
            field_detuning_ev(middle),
        )
        if probability > target:
            low = middle
        else:
            high = middle
    return 2.0 * 0.5 * (low + high)


def asimov_discovery_time(
    signal_rate_hz: float,
    background_rate_hz: float,
    significance: float = 5.0,
    minimum_signal_events: float = 10.0,
) -> float:
    """Exposure needed for an Asimov counting significance and event floor."""
    if signal_rate_hz <= 0:
        return math.inf
    event_floor_time = minimum_signal_events / signal_rate_hz
    if background_rate_hz <= 0:
        return event_floor_time
    ratio = signal_rate_hz / background_rate_hz
    if ratio < 1e-4:
        # Stable expansion of 2*b*((1+x)log(1+x)-x). The direct expression
        # catastrophically cancels when a portal rate is far below background.
        information_rate = background_rate_hz * (
            ratio**2
            - ratio**3 / 3.0
            + ratio**4 / 6.0
            - ratio**5 / 10.0
        )
    else:
        information_rate = 2.0 * background_rate_hz * (
            (1.0 + ratio) * math.log1p(ratio) - ratio
        )
    if information_rate <= 0:
        return math.inf
    significance_time = significance**2 / information_rate
    return max(event_floor_time, significance_time)


@dataclass(frozen=True)
class BeamlineResult:
    oscillation_time_s: float
    interaction_time_s: float
    converter_length_m: float
    field_fwhm_microtesla: float
    field_step_microtesla: float
    scan_settings: int
    single_swap_probability_worst_case: float
    through_wall_probability_worst_case: float
    neutron_survival: float
    signal_rate_hz: float
    dwell_time_s: float
    field_settling_time_s: float
    total_scan_time_s: float
    known_resonance_time_s: float


def beamline_result(
    oscillation_time_s: float,
    interaction_time_s: float,
    neutron_velocity_m_s: float,
    neutron_flux_hz: float,
    detection_efficiency: float,
    background_rate_hz: float,
    field_range_tesla: float,
    oversampling: float = 2.0,
    field_settling_time_s: float = 2.0,
    significance: float = 5.0,
    minimum_signal_events: float = 10.0,
) -> BeamlineResult:
    if not 0 < detection_efficiency <= 1:
        raise ValueError("detection efficiency must lie in (0, 1]")
    if oversampling < 1:
        raise ValueError("oversampling must be at least one")

    fwhm = field_fwhm_tesla(interaction_time_s, oscillation_time_s)
    step = fwhm / oversampling
    settings = max(1, math.ceil(field_range_tesla / step))
    worst_offset = step / 2.0
    single_swap = coherent_probability(
        interaction_time_s,
        oscillation_time_s,
        field_detuning_ev(worst_offset),
    )
    through_wall = single_swap**2
    survival = math.exp(-2.0 * interaction_time_s / NEUTRON_LIFETIME_S)
    signal_rate = neutron_flux_hz * detection_efficiency * through_wall * survival
    dwell = asimov_discovery_time(
        signal_rate,
        background_rate_hz,
        significance=significance,
        minimum_signal_events=minimum_signal_events,
    )

    peak_swap = coherent_probability(
        interaction_time_s, oscillation_time_s, 0.0
    )
    peak_rate = (
        neutron_flux_hz
        * detection_efficiency
        * peak_swap**2
        * survival
    )
    known_resonance_time = asimov_discovery_time(
        peak_rate,
        background_rate_hz,
        significance=significance,
        minimum_signal_events=minimum_signal_events,
    )
    return BeamlineResult(
        oscillation_time_s=oscillation_time_s,
        interaction_time_s=interaction_time_s,
        converter_length_m=interaction_time_s * neutron_velocity_m_s,
        field_fwhm_microtesla=fwhm * 1e6,
        field_step_microtesla=step * 1e6,
        scan_settings=settings,
        single_swap_probability_worst_case=single_swap,
        through_wall_probability_worst_case=through_wall,
        neutron_survival=survival,
        signal_rate_hz=signal_rate,
        dwell_time_s=dwell,
        field_settling_time_s=field_settling_time_s,
        total_scan_time_s=settings * (dwell + field_settling_time_s),
        known_resonance_time_s=known_resonance_time,
    )


def optimize_beamline(
    oscillation_time_s: float,
    max_converter_length_m: float = 10.0,
    neutron_velocity_m_s: float = 5.0,
    neutron_flux_hz: float = 5e5,
    detection_efficiency: float = 0.30,
    background_rate_hz: float = 1e-4,
    field_min_microtesla: float = 50.0,
    field_max_microtesla: float = 1100.0,
    oversampling: float = 2.0,
    field_settling_time_s: float = 2.0,
) -> BeamlineResult:
    max_time = max_converter_length_m / neutron_velocity_m_s
    min_time = min(1e-3, max_time / 1000.0)
    best: BeamlineResult | None = None
    samples = 500
    for index in range(samples):
        fraction = index / (samples - 1)
        interaction_time = math.exp(
            math.log(min_time)
            + fraction * math.log(max_time / min_time)
        )
        # Avoid exact conversion nodes, where the resonance peak vanishes.
        if coherent_probability(interaction_time, oscillation_time_s) <= 1e-20:
            continue
        candidate = beamline_result(
            oscillation_time_s=oscillation_time_s,
            interaction_time_s=interaction_time,
            neutron_velocity_m_s=neutron_velocity_m_s,
            neutron_flux_hz=neutron_flux_hz,
            detection_efficiency=detection_efficiency,
            background_rate_hz=background_rate_hz,
            field_range_tesla=(field_max_microtesla - field_min_microtesla) * 1e-6,
            oversampling=oversampling,
            field_settling_time_s=field_settling_time_s,
        )
        if best is None or candidate.total_scan_time_s < best.total_scan_time_s:
            best = candidate
    if best is None:
        raise RuntimeError("beamline optimizer found no valid point")
    return best


def sweep(output: Path) -> list[BeamlineResult]:
    results = [
        optimize_beamline(oscillation_time_s=tau)
        for tau in (1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 414.0)
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(results[0])))
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    coherent = subparsers.add_parser("coherent")
    coherent.add_argument("--tau", type=float, required=True)
    coherent.add_argument("--time", type=float, required=True)
    coherent.add_argument("--field-offset-microtesla", type=float, default=0.0)

    medium = subparsers.add_parser("medium")
    medium.add_argument("--epsilon-ev", type=float, required=True)
    medium.add_argument("--delta-ev", type=float, required=True)
    medium.add_argument("--fermi-ev", type=float, default=0.0)
    medium.add_argument("--collision-rate", type=float, required=True)

    optimize = subparsers.add_parser("optimize")
    optimize.add_argument("--tau", type=float, required=True)
    optimize.add_argument("--max-length", type=float, default=10.0)
    optimize.add_argument("--velocity", type=float, default=5.0)
    optimize.add_argument("--flux", type=float, default=5e5)
    optimize.add_argument("--efficiency", type=float, default=0.30)
    optimize.add_argument("--background", type=float, default=1e-4)
    optimize.add_argument("--field-min", type=float, default=50.0)
    optimize.add_argument("--field-max", type=float, default=1100.0)
    optimize.add_argument("--settling-time", type=float, default=2.0)

    scan = subparsers.add_parser("sweep")
    scan.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "hidden-neutron-beamline-scan.csv",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "coherent":
        probability = coherent_probability(
            args.time,
            args.tau,
            field_detuning_ev(args.field_offset_microtesla * 1e-6),
        )
        result = {
            "probability": probability,
            "pi_pulse_time_s": pi_pulse_time(args.tau),
            "field_fwhm_microtesla": field_fwhm_tesla(args.time, args.tau) * 1e6,
        }
    elif args.command == "medium":
        result = {
            "swapping_probability_per_projection": medium_swapping_probability(
                args.epsilon_ev,
                args.delta_ev,
                args.fermi_ev,
                args.collision_rate,
            )
        }
    elif args.command == "optimize":
        result = asdict(
            optimize_beamline(
                oscillation_time_s=args.tau,
                max_converter_length_m=args.max_length,
                neutron_velocity_m_s=args.velocity,
                neutron_flux_hz=args.flux,
                detection_efficiency=args.efficiency,
                background_rate_hz=args.background,
                field_min_microtesla=args.field_min,
                field_max_microtesla=args.field_max,
                field_settling_time_s=args.settling_time,
            )
        )
    else:
        result = {"results": [asdict(item) for item in sweep(args.output)]}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
