#!/usr/bin/env python3
"""Numerical checks for relativistic travel and speculative spacetime metrics.

This program deliberately separates three levels of confidence:

* Relativistic travel uses experimentally verified special relativity.
* The Alcubierre calculation evaluates the stress-energy implied by a chosen
  metric; it is not a construction method for a warp drive.
* The Morris-Thorne calculation evaluates a pedagogical zero-redshift
  wormhole; it assumes the topology already exists and does not establish
  stability or manufacturability.

All calculations use SI units internally and only Python's standard library.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


C = 299_792_458.0
G = 6.67430e-11
G0 = 9.80665
HBAR = 1.054_571_817e-34
PROTON_MASS = 1.672_621_923_69e-27
ELECTRON_VOLT = 1.602_176_634e-19
JULIAN_YEAR = 365.25 * 86_400.0
LIGHT_YEAR = C * JULIAN_YEAR
SOLAR_MASS = 1.98847e30


def _require_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and greater than zero")


def _require_beta(beta: float, *, subluminal: bool = True) -> None:
    if not math.isfinite(beta) or beta < 0:
        raise ValueError("beta must be finite and non-negative")
    if subluminal and beta >= 1:
        raise ValueError("this calculation requires 0 <= beta < 1")


def lorentz_gamma(beta: float) -> float:
    """Return the Lorentz factor for a subluminal speed beta = v/c."""
    _require_beta(beta)
    return 1.0 / math.sqrt(1.0 - beta * beta)


@dataclass(frozen=True)
class CruiseMission:
    beta: float
    gamma: float
    mass_kg: float
    distance_ly: float
    earth_cruise_years: float
    traveler_cruise_years: float
    kinetic_energy_j: float
    ideal_antimatter_for_acceleration_kg: float
    proton_impact_energy_ev: float
    dust_impact_energy_j: float


def cruise_mission(
    mass_kg: float,
    beta: float,
    distance_ly: float,
    dust_mass_kg: float = 1e-9,
) -> CruiseMission:
    """Constant-speed cruise, excluding acceleration, braking, and inefficiency."""
    _require_positive("mass_kg", mass_kg)
    _require_positive("distance_ly", distance_ly)
    _require_positive("dust_mass_kg", dust_mass_kg)
    _require_beta(beta)
    if beta == 0:
        raise ValueError("beta must be greater than zero for a cruise mission")

    gamma = lorentz_gamma(beta)
    earth_years = distance_ly / beta
    traveler_years = earth_years / gamma
    kinetic_energy = (gamma - 1.0) * mass_kg * C**2
    proton_energy = (gamma - 1.0) * PROTON_MASS * C**2 / ELECTRON_VOLT
    dust_energy = (gamma - 1.0) * dust_mass_kg * C**2

    # Perfect annihilation releases 2*m*c^2 from equal antimatter and matter.
    antimatter_mass = kinetic_energy / (2.0 * C**2)
    return CruiseMission(
        beta=beta,
        gamma=gamma,
        mass_kg=mass_kg,
        distance_ly=distance_ly,
        earth_cruise_years=earth_years,
        traveler_cruise_years=traveler_years,
        kinetic_energy_j=kinetic_energy,
        ideal_antimatter_for_acceleration_kg=antimatter_mass,
        proton_impact_energy_ev=proton_energy,
        dust_impact_energy_j=dust_energy,
    )


@dataclass(frozen=True)
class ProperAccelerationTrip:
    distance_ly: float
    proper_acceleration_m_s2: float
    earth_years: float
    traveler_years: float
    peak_beta: float
    peak_gamma: float


def proper_acceleration_trip(
    distance_ly: float, proper_acceleration_m_s2: float = G0
) -> ProperAccelerationTrip:
    """Symmetric accelerate-halfway/decelerate-halfway relativistic trip."""
    _require_positive("distance_ly", distance_ly)
    _require_positive("proper_acceleration_m_s2", proper_acceleration_m_s2)

    distance_m = distance_ly * LIGHT_YEAR
    half_distance = distance_m / 2.0
    rapidity = math.acosh(
        1.0 + proper_acceleration_m_s2 * half_distance / C**2
    )
    half_ship_s = C * rapidity / proper_acceleration_m_s2
    half_earth_s = C * math.sinh(rapidity) / proper_acceleration_m_s2
    return ProperAccelerationTrip(
        distance_ly=distance_ly,
        proper_acceleration_m_s2=proper_acceleration_m_s2,
        earth_years=2.0 * half_earth_s / JULIAN_YEAR,
        traveler_years=2.0 * half_ship_s / JULIAN_YEAR,
        peak_beta=math.tanh(rapidity),
        peak_gamma=math.cosh(rapidity),
    )


def _sech_squared(value: float) -> float:
    if abs(value) > 40.0:
        return 0.0
    cosh_value = math.cosh(value)
    return 1.0 / (cosh_value * cosh_value)


def alcubierre_shape(radius_m: float, bubble_radius_m: float, wall_scale_m: float) -> float:
    """Alcubierre's smooth top-hat shape f(r), with sigma = 1/wall_scale."""
    _require_positive("bubble_radius_m", bubble_radius_m)
    _require_positive("wall_scale_m", wall_scale_m)
    if radius_m < 0:
        raise ValueError("radius_m cannot be negative")
    sigma = 1.0 / wall_scale_m
    denominator = 2.0 * math.tanh(sigma * bubble_radius_m)
    return (
        math.tanh(sigma * (radius_m + bubble_radius_m))
        - math.tanh(sigma * (radius_m - bubble_radius_m))
    ) / denominator


def alcubierre_shape_derivative(
    radius_m: float, bubble_radius_m: float, wall_scale_m: float
) -> float:
    _require_positive("bubble_radius_m", bubble_radius_m)
    _require_positive("wall_scale_m", wall_scale_m)
    if radius_m < 0:
        raise ValueError("radius_m cannot be negative")
    sigma = 1.0 / wall_scale_m
    denominator = 2.0 * math.tanh(sigma * bubble_radius_m)
    return sigma * (
        _sech_squared(sigma * (radius_m + bubble_radius_m))
        - _sech_squared(sigma * (radius_m - bubble_radius_m))
    ) / denominator


def _simpson_integral(function, lower: float, upper: float, intervals: int) -> float:
    if intervals < 2:
        raise ValueError("intervals must be at least two")
    if intervals % 2:
        intervals += 1
    width = (upper - lower) / intervals
    total = function(lower) + function(upper)
    for index in range(1, intervals):
        coefficient = 4.0 if index % 2 else 2.0
        total += coefficient * function(lower + index * width)
    return total * width / 3.0


@dataclass(frozen=True)
class AlcubierreBubble:
    beta: float
    bubble_radius_m: float
    wall_scale_m: float
    negative_energy_j: float
    mass_equivalent_kg: float
    mass_equivalent_solar: float
    peak_negative_energy_density_j_m3: float
    radial_integral_m: float
    flat_interior_radius_m: float


def alcubierre_flat_interior_radius(
    bubble_radius_m: float,
    wall_scale_m: float,
    threshold: float = 0.99,
) -> float:
    """Radius over which the Alcubierre shape remains at least `threshold`."""
    _require_positive("bubble_radius_m", bubble_radius_m)
    _require_positive("wall_scale_m", wall_scale_m)
    if not 0.0 < threshold < 1.0:
        raise ValueError("threshold must lie strictly between zero and one")

    lower = 0.0
    upper = bubble_radius_m + 40.0 * wall_scale_m
    for _ in range(100):
        middle = (lower + upper) / 2.0
        if alcubierre_shape(middle, bubble_radius_m, wall_scale_m) >= threshold:
            lower = middle
        else:
            upper = middle
    return lower


def alcubierre_bubble(
    bubble_radius_m: float,
    wall_scale_m: float,
    beta: float,
    integration_intervals: int = 20_000,
) -> AlcubierreBubble:
    """Evaluate negative Eulerian energy on one Alcubierre spatial slice.

    This integrates the standard Eulerian density for the original metric:

        rho = -(c^4/G) beta^2 sin(theta)^2 [f'(r)]^2 / (32*pi)

    Angular integration gives E = -(c^4/G) beta^2 I/12, where
    I = integral r^2 [f'(r)]^2 dr. This is not a formation-energy estimate.
    """
    _require_positive("bubble_radius_m", bubble_radius_m)
    _require_positive("wall_scale_m", wall_scale_m)
    _require_beta(beta, subluminal=False)

    upper = bubble_radius_m + 20.0 * wall_scale_m

    def integrand(radius: float) -> float:
        derivative = alcubierre_shape_derivative(
            radius, bubble_radius_m, wall_scale_m
        )
        return radius * radius * derivative * derivative

    radial_integral = _simpson_integral(
        integrand, 0.0, upper, integration_intervals
    )
    energy = -(C**4 / G) * beta * beta * radial_integral / 12.0

    # Search the wall for the equatorial peak of the negative density.
    peak_derivative_squared = 0.0
    samples = 4000
    for index in range(samples + 1):
        radius = upper * index / samples
        derivative = alcubierre_shape_derivative(
            radius, bubble_radius_m, wall_scale_m
        )
        peak_derivative_squared = max(peak_derivative_squared, derivative**2)
    peak_density = -(C**4 / G) * beta * beta * peak_derivative_squared / (
        32.0 * math.pi
    )
    mass_equivalent = abs(energy) / C**2
    return AlcubierreBubble(
        beta=beta,
        bubble_radius_m=bubble_radius_m,
        wall_scale_m=wall_scale_m,
        negative_energy_j=energy,
        mass_equivalent_kg=mass_equivalent,
        mass_equivalent_solar=mass_equivalent / SOLAR_MASS,
        peak_negative_energy_density_j_m3=peak_density,
        radial_integral_m=radial_integral,
        flat_interior_radius_m=alcubierre_flat_interior_radius(
            bubble_radius_m, wall_scale_m
        ),
    )


@dataclass(frozen=True)
class AlcubierreOptimization:
    requested_flat_radius_m: float
    flatness_threshold: float
    beta: float
    bubble_radius_m: float
    wall_scale_m: float
    characteristic_outer_radius_m: float
    negative_energy_j: float
    mass_equivalent_kg: float
    mass_equivalent_solar: float
    peak_negative_energy_density_j_m3: float


@dataclass(frozen=True)
class AlcubierreVariationalBound:
    cavity_radius_m: float
    outer_radius_m: float
    beta: float
    minimum_radial_integral_m: float
    negative_energy_upper_bound_j: float
    minimum_mass_equivalent_kg: float
    minimum_mass_equivalent_solar: float


def alcubierre_variational_bound(
    cavity_radius_m: float,
    beta: float = 1.0,
    outer_radius_m: float = math.inf,
) -> AlcubierreVariationalBound:
    """Least |negative slice energy| for any spherical Alcubierre shape.

    For f=1 through radius a and f=0 at radius b, minimizing
    integral(r^2 f'(r)^2 dr) gives r^2 f'=constant and
    I_min=ab/(b-a). In the limit b->infinity, I_min=a. Smooth profiles can
    approach the bound. It does not remove the negative sign or solve formation.
    """
    _require_positive("cavity_radius_m", cavity_radius_m)
    _require_beta(beta, subluminal=False)
    if not math.isinf(outer_radius_m):
        _require_positive("outer_radius_m", outer_radius_m)
        if outer_radius_m <= cavity_radius_m:
            raise ValueError("outer_radius_m must exceed cavity_radius_m")
        radial_integral = (
            cavity_radius_m
            * outer_radius_m
            / (outer_radius_m - cavity_radius_m)
        )
    else:
        radial_integral = cavity_radius_m

    energy = -(C**4 / G) * beta**2 * radial_integral / 12.0
    mass = abs(energy) / C**2
    return AlcubierreVariationalBound(
        cavity_radius_m=cavity_radius_m,
        outer_radius_m=outer_radius_m,
        beta=beta,
        minimum_radial_integral_m=radial_integral,
        negative_energy_upper_bound_j=energy,
        minimum_mass_equivalent_kg=mass,
        minimum_mass_equivalent_solar=mass / SOLAR_MASS,
    )


def minimum_alcubierre_for_cavity(
    requested_flat_radius_m: float,
    beta: float = 1.0,
    threshold: float = 0.99,
    samples: int = 120,
) -> AlcubierreOptimization:
    """Grid-search the original shape for the least |negative energy|.

    The search requires the nominal bubble radius to be at least the requested
    cabin radius, and requires f(r) >= threshold throughout that cabin. Wall
    scales from 0.01 to 10 times the cabin radius are explored. This optimizes
    only one spatial-slice integral, not formation, stability, or control.
    """
    _require_positive("requested_flat_radius_m", requested_flat_radius_m)
    _require_beta(beta, subluminal=False)
    if not 0.0 < threshold < 1.0:
        raise ValueError("threshold must lie strictly between zero and one")
    if samples < 10:
        raise ValueError("samples must be at least ten")

    best: AlcubierreBubble | None = None
    for index in range(samples):
        exponent = -2.0 + 3.0 * index / (samples - 1)
        wall = requested_flat_radius_m * 10.0**exponent
        lower_radius = requested_flat_radius_m
        if (
            alcubierre_flat_interior_radius(lower_radius, wall, threshold)
            >= requested_flat_radius_m
        ):
            radius = lower_radius
        else:
            upper_radius = requested_flat_radius_m + 50.0 * wall
            while (
                alcubierre_flat_interior_radius(upper_radius, wall, threshold)
                < requested_flat_radius_m
            ):
                upper_radius *= 2.0
            for _ in range(80):
                middle = (lower_radius + upper_radius) / 2.0
                if (
                    alcubierre_flat_interior_radius(middle, wall, threshold)
                    >= requested_flat_radius_m
                ):
                    upper_radius = middle
                else:
                    lower_radius = middle
            radius = upper_radius
        candidate = alcubierre_bubble(radius, wall, beta, 4000)
        if best is None or abs(candidate.negative_energy_j) < abs(best.negative_energy_j):
            best = candidate

    assert best is not None
    return AlcubierreOptimization(
        requested_flat_radius_m=requested_flat_radius_m,
        flatness_threshold=threshold,
        beta=beta,
        bubble_radius_m=best.bubble_radius_m,
        wall_scale_m=best.wall_scale_m,
        characteristic_outer_radius_m=best.bubble_radius_m
        + 6.0 * best.wall_scale_m,
        negative_energy_j=best.negative_energy_j,
        mass_equivalent_kg=best.mass_equivalent_kg,
        mass_equivalent_solar=best.mass_equivalent_solar,
        peak_negative_energy_density_j_m3=best.peak_negative_energy_density_j_m3,
    )


@dataclass(frozen=True)
class MorrisThorneWormhole:
    throat_radius_m: float
    traversal_beta: float
    body_length_m: float
    throat_energy_density_j_m3: float
    throat_radial_pressure_pa: float
    throat_nec_j_m3: float
    lateral_tidal_acceleration_m_s2: float
    maximum_beta_for_one_g_tides: float
    two_sided_proper_volume_nec_integral_j: float


def morris_thorne_wormhole(
    throat_radius_m: float,
    traversal_beta: float,
    body_length_m: float = 2.0,
) -> MorrisThorneWormhole:
    """Evaluate b(r)=r0^2/r and Phi=0 at the throat.

    The geometry has no redshift horizon, but it requires negative energy and
    this function does not address creation or dynamical stability.
    """
    _require_positive("throat_radius_m", throat_radius_m)
    _require_positive("body_length_m", body_length_m)
    _require_beta(traversal_beta)
    gamma = lorentz_gamma(traversal_beta)

    density = -(C**4 / G) / (8.0 * math.pi * throat_radius_m**2)
    radial_pressure = density
    nec = 2.0 * density
    lateral_tidal = (
        C**2
        * gamma**2
        * traversal_beta**2
        * body_length_m
        / throat_radius_m**2
    )

    # Solve gamma^2 beta^2 = q, where q = g*r0^2/(c^2*L).
    q = G0 * throat_radius_m**2 / (C**2 * body_length_m)
    maximum_beta = math.sqrt(q / (1.0 + q))
    # Integral of (rho+p_r) over the proper spatial volume of both sides.
    integrated_nec = -math.pi * (C**4 / G) * throat_radius_m
    return MorrisThorneWormhole(
        throat_radius_m=throat_radius_m,
        traversal_beta=traversal_beta,
        body_length_m=body_length_m,
        throat_energy_density_j_m3=density,
        throat_radial_pressure_pa=radial_pressure,
        throat_nec_j_m3=nec,
        lateral_tidal_acceleration_m_s2=lateral_tidal,
        maximum_beta_for_one_g_tides=maximum_beta,
        two_sided_proper_volume_nec_integral_j=integrated_nec,
    )


def casimir_energy_density(plate_separation_m: float) -> float:
    """Ideal parallel-plate electromagnetic Casimir energy density."""
    _require_positive("plate_separation_m", plate_separation_m)
    return -(math.pi**2 * HBAR * C) / (720.0 * plate_separation_m**4)


@dataclass(frozen=True)
class KaluzaKleinMode:
    compact_radius_m: float
    circumference_m: float
    mode_number: int
    zero_mode_mass_energy_ev: float
    extra_dimension_momentum_kg_m_s: float
    kk_gap_energy_ev: float
    total_mode_mass_energy_ev: float
    excitation_above_zero_mode_ev: float
    equivalent_total_mass_kg: float
    gap_frequency_hz: float


def kaluza_klein_mode(
    compact_radius_m: float,
    mode_number: int = 1,
    zero_mode_mass_energy_ev: float = 0.0,
) -> KaluzaKleinMode:
    """Energy of a field mode carrying momentum around a compact S1 dimension.

    The quantized extra-dimensional momentum is p_y=n*hbar/R and the apparent
    four-dimensional mass obeys E_n^2=E_0^2+(n*hbar*c/R)^2. This calculation
    assumes the chosen field is allowed to propagate in the bulk. In common
    brane models, ordinary Standard Model matter is not.
    """
    _require_positive("compact_radius_m", compact_radius_m)
    if not isinstance(mode_number, int) or mode_number < 1:
        raise ValueError("mode_number must be a positive integer")
    if not math.isfinite(zero_mode_mass_energy_ev) or zero_mode_mass_energy_ev < 0:
        raise ValueError("zero_mode_mass_energy_ev must be finite and non-negative")

    momentum = mode_number * HBAR / compact_radius_m
    gap_ev = momentum * C / ELECTRON_VOLT
    total_ev = math.hypot(zero_mode_mass_energy_ev, gap_ev)
    excitation_ev = total_ev - zero_mode_mass_energy_ev
    total_mass_kg = total_ev * ELECTRON_VOLT / C**2
    gap_frequency_hz = gap_ev * ELECTRON_VOLT / (2.0 * math.pi * HBAR)
    return KaluzaKleinMode(
        compact_radius_m=compact_radius_m,
        circumference_m=2.0 * math.pi * compact_radius_m,
        mode_number=mode_number,
        zero_mode_mass_energy_ev=zero_mode_mass_energy_ev,
        extra_dimension_momentum_kg_m_s=momentum,
        kk_gap_energy_ev=gap_ev,
        total_mode_mass_energy_ev=total_ev,
        excitation_above_zero_mode_ev=excitation_ev,
        equivalent_total_mass_kg=total_mass_kg,
        gap_frequency_hz=gap_frequency_hz,
    )


def compact_radius_for_kk_gap(gap_energy_ev: float, mode_number: int = 1) -> float:
    """Return R=n*hbar*c/E for a requested Kaluza-Klein gap."""
    _require_positive("gap_energy_ev", gap_energy_ev)
    if not isinstance(mode_number, int) or mode_number < 1:
        raise ValueError("mode_number must be a positive integer")
    return mode_number * HBAR * C / (gap_energy_ev * ELECTRON_VOLT)


@dataclass(frozen=True)
class SyntheticDimensionTransfer:
    sites: int
    normalized_time: float
    probabilities: tuple[float, ...]
    engineered_couplings: tuple[float, ...]
    sender_probability: float
    receiver_probability: float
    intermediate_probability: float
    mean_synthetic_site: float
    entropy_bits: float


def synthetic_dimension_transfer(
    sites: int = 9, normalized_time: float = math.pi / 4.0
) -> SyntheticDimensionTransfer:
    """Perfect-state-transfer chain representing a finite synthetic dimension.

    Adjacent synthetic sites have engineered dimensionless couplings
    sqrt((j+1)(N-1-j)). Starting at site zero, the probability distribution is
    binomial and reaches site N-1 exactly at normalized time pi/2. The sites can
    represent photon frequencies, atomic spin states, or other internal modes;
    they are not evidence of an additional physical direction in space.
    """
    if not isinstance(sites, int) or sites < 2:
        raise ValueError("sites must be an integer of at least two")
    if not math.isfinite(normalized_time) or not 0.0 <= normalized_time <= math.pi / 2.0:
        raise ValueError("normalized_time must lie between zero and pi/2")

    trials = sites - 1
    sin_squared = math.sin(normalized_time) ** 2
    cos_squared = math.cos(normalized_time) ** 2
    probabilities = tuple(
        math.comb(trials, site)
        * sin_squared**site
        * cos_squared ** (trials - site)
        for site in range(sites)
    )
    normalization = sum(probabilities)
    probabilities = tuple(value / normalization for value in probabilities)
    couplings = tuple(
        math.sqrt((site + 1) * (sites - 1 - site))
        for site in range(sites - 1)
    )
    mean_site = sum(index * value for index, value in enumerate(probabilities))
    entropy = -sum(value * math.log2(value) for value in probabilities if value > 0)
    return SyntheticDimensionTransfer(
        sites=sites,
        normalized_time=normalized_time,
        probabilities=probabilities,
        engineered_couplings=couplings,
        sender_probability=probabilities[0],
        receiver_probability=probabilities[-1],
        intermediate_probability=sum(probabilities[1:-1]),
        mean_synthetic_site=mean_site,
        entropy_bits=entropy,
    )


def _scientific_json(value):
    if isinstance(value, float):
        return float(f"{value:.12g}")
    if isinstance(value, dict):
        return {key: _scientific_json(item) for key, item in value.items()}
    return value


def _print_result(result) -> None:
    print(json.dumps(_scientific_json(asdict(result)), indent=2, sort_keys=True))


def write_sweep(path: Path) -> None:
    """Write a compact cross-model parameter sweep for further analysis."""
    rows: list[dict[str, float | str]] = []
    for beta in (0.1, 0.5, 0.9, 0.99, 0.999, 0.9999):
        mission = cruise_mission(1_000.0, beta, 4.2465)
        rows.append(
            {
                "model": "relativistic_cruise",
                "parameter_1": beta,
                "parameter_2": 1_000.0,
                "primary_value": mission.traveler_cruise_years,
                "secondary_value": mission.kinetic_energy_j,
            }
        )
    for wall_scale in (0.1, 1.0, 10.0, 100.0):
        bubble = alcubierre_bubble(100.0, wall_scale, 1.0, 6000)
        rows.append(
            {
                "model": "alcubierre_slice",
                "parameter_1": wall_scale,
                "parameter_2": 100.0,
                "primary_value": bubble.negative_energy_j,
                "secondary_value": bubble.peak_negative_energy_density_j_m3,
            }
        )
    for throat_radius in (1.0, 10.0, 100.0, 1_000.0, 10_000.0):
        wormhole = morris_thorne_wormhole(throat_radius, 1e-6)
        rows.append(
            {
                "model": "morris_thorne_throat",
                "parameter_1": throat_radius,
                "parameter_2": 1e-6,
                "primary_value": wormhole.throat_nec_j_m3,
                "secondary_value": wormhole.lateral_tidal_acceleration_m_s2,
            }
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def write_extra_dimension_sweep(path: Path) -> None:
    """Write representative compact and synthetic-dimension calculations."""
    rows: list[dict[str, float | int | str]] = []
    for radius in (30e-6, 1e-6, 1e-9, compact_radius_for_kk_gap(1.5e12), 1e-22):
        mode = kaluza_klein_mode(radius)
        rows.append(
            {
                "model": "compact_S1_KK",
                "parameter_1": radius,
                "parameter_2": 1,
                "primary_value": mode.kk_gap_energy_ev,
                "secondary_value": mode.equivalent_total_mass_kg,
            }
        )
    for progress in (0.0, 0.25, 0.5, 0.75, 1.0):
        transfer = synthetic_dimension_transfer(9, progress * math.pi / 2.0)
        rows.append(
            {
                "model": "synthetic_state_transfer",
                "parameter_1": progress,
                "parameter_2": 9,
                "primary_value": transfer.receiver_probability,
                "secondary_value": transfer.intermediate_probability,
            }
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    mission = subparsers.add_parser("mission", help="constant-speed cruise")
    mission.add_argument("--mass", type=float, default=1_000.0, help="ship mass, kg")
    mission.add_argument("--beta", type=float, default=0.99, help="speed divided by c")
    mission.add_argument(
        "--distance", type=float, default=4.2465, help="distance, light-years"
    )
    mission.add_argument(
        "--dust-mass", type=float, default=1e-9, help="impacting grain mass, kg"
    )

    accelerated = subparsers.add_parser(
        "accelerated", help="accelerate halfway, then decelerate"
    )
    accelerated.add_argument("--distance", type=float, default=4.2465)
    accelerated.add_argument("--acceleration", type=float, default=G0)

    warp = subparsers.add_parser("warp", help="Alcubierre spatial-slice energy")
    warp.add_argument("--radius", type=float, default=100.0, help="bubble radius, m")
    warp.add_argument("--wall", type=float, default=10.0, help="wall scale 1/sigma, m")
    warp.add_argument("--beta", type=float, default=1.0, help="bubble coordinate speed/c")
    warp.add_argument("--intervals", type=int, default=20_000)

    optimize_warp = subparsers.add_parser(
        "optimize-warp", help="minimize slice energy for a flat cabin radius"
    )
    optimize_warp.add_argument(
        "--cavity-radius", type=float, default=20.0, help="99%-flat cabin radius, m"
    )
    optimize_warp.add_argument("--beta", type=float, default=1.0)
    optimize_warp.add_argument("--threshold", type=float, default=0.99)
    optimize_warp.add_argument("--samples", type=int, default=120)

    warp_bound = subparsers.add_parser(
        "warp-bound", help="variational energy bound for any spherical wall"
    )
    warp_bound.add_argument("--cavity-radius", type=float, default=20.0)
    warp_bound.add_argument("--outer-radius", type=float, default=math.inf)
    warp_bound.add_argument("--beta", type=float, default=1.0)

    wormhole = subparsers.add_parser(
        "wormhole", help="zero-redshift Morris-Thorne throat"
    )
    wormhole.add_argument("--radius", type=float, default=100.0, help="throat radius, m")
    wormhole.add_argument("--beta", type=float, default=1e-6, help="traversal speed/c")
    wormhole.add_argument("--body-length", type=float, default=2.0, help="body length, m")

    casimir = subparsers.add_parser("casimir", help="ideal plate energy density")
    casimir.add_argument("--separation", type=float, default=1e-6, help="plate gap, m")

    kk = subparsers.add_parser("kk", help="compact extra-dimension mode threshold")
    kk.add_argument("--radius", type=float, default=30e-6, help="compact radius, m")
    kk.add_argument("--mode", type=int, default=1, help="positive KK mode number")
    kk.add_argument(
        "--zero-mass-ev",
        type=float,
        default=0.0,
        help="zero-mode rest energy, eV",
    )

    synthetic = subparsers.add_parser(
        "synthetic", help="state transfer along a synthetic dimension"
    )
    synthetic.add_argument("--sites", type=int, default=9)
    synthetic.add_argument(
        "--tau", type=float, default=math.pi / 4.0, help="normalized time, 0 to pi/2"
    )

    sweep = subparsers.add_parser("sweep", help="write representative CSV sweep")
    sweep.add_argument("--output", type=Path, default=Path("spacetime-sweep.csv"))
    extra_sweep = subparsers.add_parser(
        "extra-sweep", help="write compact and synthetic-dimension CSV sweep"
    )
    extra_sweep.add_argument(
        "--output", type=Path, default=Path("extra-dimension-sweep.csv")
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "mission":
        _print_result(cruise_mission(args.mass, args.beta, args.distance, args.dust_mass))
    elif args.command == "accelerated":
        _print_result(proper_acceleration_trip(args.distance, args.acceleration))
    elif args.command == "warp":
        _print_result(alcubierre_bubble(args.radius, args.wall, args.beta, args.intervals))
    elif args.command == "optimize-warp":
        _print_result(
            minimum_alcubierre_for_cavity(
                args.cavity_radius, args.beta, args.threshold, args.samples
            )
        )
    elif args.command == "warp-bound":
        _print_result(
            alcubierre_variational_bound(
                args.cavity_radius, args.beta, args.outer_radius
            )
        )
    elif args.command == "wormhole":
        _print_result(morris_thorne_wormhole(args.radius, args.beta, args.body_length))
    elif args.command == "casimir":
        print(json.dumps({"energy_density_j_m3": casimir_energy_density(args.separation)}, indent=2))
    elif args.command == "kk":
        _print_result(kaluza_klein_mode(args.radius, args.mode, args.zero_mass_ev))
    elif args.command == "synthetic":
        _print_result(synthetic_dimension_transfer(args.sites, args.tau))
    elif args.command == "sweep":
        write_sweep(args.output)
        print(args.output.resolve())
    elif args.command == "extra-sweep":
        write_extra_dimension_sweep(args.output)
        print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
