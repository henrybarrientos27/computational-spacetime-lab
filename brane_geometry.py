#!/usr/bin/env python3
"""Conditional inference from neutron swapping to a two-brane geometry.

This module deliberately separates three layers that are often conflated:

1. an observed or bounded swapping probability p;
2. a model-dependent mixing energy epsilon;
3. a still more model-dependent interbrane coupling g and separation d.

No experiment has observed a nonzero value at layer 1.  The geometric map at
layer 3 is implemented only for the non-compact 5D DGP-brane ansatz

    g = (m_q**2 / M_B) exp(-m_q d)

in natural units.  It is not a generic prediction of extra dimensions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


HBAR_EV_S = 6.582119569e-16
HBAR_C_EV_M = 1.973269804e-7
NEUTRON_MAGNETIC_MOMENT_EV_T = 6.030774e-8
CONSTITUENT_QUARK_MASS_EV = 340e6
PLANCK_ENERGY_EV = 1.22089e28


def mixing_energy_limit_from_probability(
    probability: float,
    effective_detuning_ev: float,
    collision_rate_hz: float = 0.0,
) -> float:
    """Invert the collision-averaged STEREO probability for epsilon.

    ``effective_detuning_ev`` is Delta E + V_F.  The probability model is

        p = 2 epsilon^2 /
            (effective_detuning^2 + 4 epsilon^2 + (hbar Gamma/2)^2).

    A probability upper limit therefore gives an epsilon upper limit only
    after the detuning and collision rate have been specified.
    """
    if not 0.0 < probability < 0.5:
        raise ValueError("probability must lie between zero and one half")
    if collision_rate_hz < 0:
        raise ValueError("collision rate must be nonnegative")
    projection_width_ev = HBAR_EV_S * collision_rate_hz / 2.0
    scale_squared = effective_detuning_ev**2 + projection_width_ev**2
    return math.sqrt(probability * scale_squared / (2.0 - 4.0 * probability))


def swapping_probability(
    mixing_energy_ev: float,
    effective_detuning_ev: float,
    collision_rate_hz: float = 0.0,
) -> float:
    """Collision-averaged probability per projection."""
    if mixing_energy_ev < 0 or collision_rate_hz < 0:
        raise ValueError("mixing energy and collision rate must be nonnegative")
    projection_width_ev = HBAR_EV_S * collision_rate_hz / 2.0
    denominator = (
        effective_detuning_ev**2
        + 4.0 * mixing_energy_ev**2
        + projection_width_ev**2
    )
    if denominator == 0:
        return 0.0
    return 2.0 * mixing_energy_ev**2 / denominator


def geometric_coupling_from_mixing(
    mixing_energy_ev: float,
    vector_potential_difference_tesla_m: float,
) -> float:
    """Infer g in m^-1 under epsilon = |mu_n| g |A_+ - A_-|."""
    if mixing_energy_ev < 0:
        raise ValueError("mixing energy must be nonnegative")
    if vector_potential_difference_tesla_m <= 0:
        raise ValueError("vector-potential difference must be positive")
    return mixing_energy_ev / (
        NEUTRON_MAGNETIC_MOMENT_EV_T
        * vector_potential_difference_tesla_m
    )


def mixing_from_geometric_coupling(
    coupling_per_m: float,
    vector_potential_difference_tesla_m: float,
) -> float:
    if coupling_per_m < 0 or vector_potential_difference_tesla_m <= 0:
        raise ValueError("coupling must be nonnegative and vector potential positive")
    return (
        NEUTRON_MAGNETIC_MOMENT_EV_T
        * coupling_per_m
        * vector_potential_difference_tesla_m
    )


def coupling_5d_noncompact_per_m(
    separation_m: float,
    brane_scale_ev: float,
    quark_mass_ev: float = CONSTITUENT_QUARK_MASS_EV,
) -> float:
    """DGP-brane coupling g=(m_q^2/M_B) exp(-m_q d), in m^-1."""
    if separation_m < 0 or brane_scale_ev <= 0 or quark_mass_ev <= 0:
        raise ValueError("separation must be nonnegative and energy scales positive")
    exponent = -quark_mass_ev * separation_m / HBAR_C_EV_M
    zero_separation_coupling_per_m = (
        quark_mass_ev**2 / brane_scale_ev / HBAR_C_EV_M
    )
    if exponent < -745.0:
        return 0.0
    return zero_separation_coupling_per_m * math.exp(exponent)


def separation_lower_bound_5d_noncompact(
    coupling_upper_limit_per_m: float,
    brane_scale_ev: float,
    quark_mass_ev: float = CONSTITUENT_QUARK_MASS_EV,
) -> float | None:
    """Conditional lower bound on d; None means the g limit has no reach."""
    if coupling_upper_limit_per_m <= 0:
        raise ValueError("coupling limit must be positive")
    zero_coupling = coupling_5d_noncompact_per_m(
        0.0, brane_scale_ev, quark_mass_ev
    )
    if coupling_upper_limit_per_m >= zero_coupling:
        return None
    return (
        HBAR_C_EV_M
        / quark_mass_ev
        * math.log(zero_coupling / coupling_upper_limit_per_m)
    )


@dataclass(frozen=True)
class ConditionalGeometry:
    swapping_probability_limit: float
    assumed_effective_detuning_ev: float
    assumed_collision_rate_hz: float
    mixing_energy_limit_ev: float
    assumed_vector_potential_difference_tesla_m: float
    geometric_coupling_limit_per_m: float
    assumed_brane_scale_ev: float
    zero_separation_coupling_per_m: float
    separation_is_constrained: bool
    conditional_separation_lower_bound_m: float | None
    conditional_separation_lower_bound_fm: float | None


def conditional_geometry_limit(
    probability_limit: float,
    effective_detuning_ev: float,
    vector_potential_difference_tesla_m: float,
    brane_scale_ev: float,
    collision_rate_hz: float = 0.0,
    quark_mass_ev: float = CONSTITUENT_QUARK_MASS_EV,
) -> ConditionalGeometry:
    epsilon = mixing_energy_limit_from_probability(
        probability_limit,
        effective_detuning_ev,
        collision_rate_hz,
    )
    coupling = geometric_coupling_from_mixing(
        epsilon, vector_potential_difference_tesla_m
    )
    zero_coupling = coupling_5d_noncompact_per_m(
        0.0, brane_scale_ev, quark_mass_ev
    )
    distance = separation_lower_bound_5d_noncompact(
        coupling, brane_scale_ev, quark_mass_ev
    )
    return ConditionalGeometry(
        swapping_probability_limit=probability_limit,
        assumed_effective_detuning_ev=effective_detuning_ev,
        assumed_collision_rate_hz=collision_rate_hz,
        mixing_energy_limit_ev=epsilon,
        assumed_vector_potential_difference_tesla_m=(
            vector_potential_difference_tesla_m
        ),
        geometric_coupling_limit_per_m=coupling,
        assumed_brane_scale_ev=brane_scale_ev,
        zero_separation_coupling_per_m=zero_coupling,
        separation_is_constrained=distance is not None,
        conditional_separation_lower_bound_m=distance,
        conditional_separation_lower_bound_fm=(
            None if distance is None else distance / 1e-15
        ),
    )


def write_degeneracy_sweep(path: Path, probability_limit: float) -> None:
    """Write equally data-compatible geometric interpretations."""
    detunings = (1e-8, 1.0, 2e3)
    vector_potentials = (2e9, 2e12)
    brane_scales = (1e12, 1e25, PLANCK_ENERGY_EV)
    rows = []
    for detuning in detunings:
        for vector_potential in vector_potentials:
            for brane_scale in brane_scales:
                rows.append(
                    asdict(
                        conditional_geometry_limit(
                            probability_limit,
                            detuning,
                            vector_potential,
                            brane_scale,
                        )
                    )
                )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    infer = subparsers.add_parser("infer", help="make one conditional inference")
    infer.add_argument("--probability", type=float, default=3.1e-11)
    infer.add_argument("--detuning-ev", type=float, default=2e3)
    infer.add_argument("--collision-rate", type=float, default=0.0)
    infer.add_argument("--vector-potential-tm", type=float, default=2e9)
    infer.add_argument("--brane-scale-ev", type=float, default=PLANCK_ENERGY_EV)

    forward = subparsers.add_parser("forward", help="predict p for a 5D point")
    forward.add_argument("--separation-m", type=float, required=True)
    forward.add_argument("--detuning-ev", type=float, default=2e3)
    forward.add_argument("--collision-rate", type=float, default=0.0)
    forward.add_argument("--vector-potential-tm", type=float, default=2e9)
    forward.add_argument("--brane-scale-ev", type=float, default=PLANCK_ENERGY_EV)

    sweep = subparsers.add_parser("sweep", help="write an assumption sweep")
    sweep.add_argument("--probability", type=float, default=3.1e-11)
    sweep.add_argument(
        "--output", type=Path, default=Path("brane-geometry-degeneracy.csv")
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "infer":
        result = conditional_geometry_limit(
            args.probability,
            args.detuning_ev,
            args.vector_potential_tm,
            args.brane_scale_ev,
            args.collision_rate,
        )
        print(json.dumps(asdict(result), indent=2))
    elif args.command == "forward":
        coupling = coupling_5d_noncompact_per_m(
            args.separation_m, args.brane_scale_ev
        )
        epsilon = mixing_from_geometric_coupling(
            coupling, args.vector_potential_tm
        )
        probability = swapping_probability(
            epsilon, args.detuning_ev, args.collision_rate
        )
        print(
            json.dumps(
                {
                    "coupling_per_m": coupling,
                    "mixing_energy_ev": epsilon,
                    "swapping_probability": probability,
                },
                indent=2,
            )
        )
    elif args.command == "sweep":
        write_degeneracy_sweep(args.output, args.probability)
        print(json.dumps({"output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
