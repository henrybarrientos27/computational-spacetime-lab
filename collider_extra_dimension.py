#!/usr/bin/env python3
"""Reinterpret the CMS EXO-20-004 monojet search for ADD extra dimensions.

This program downloads the official HEPData tables, constructs the CMS
simplified likelihood (Poisson counts plus correlated Gaussian background
nuisances), reproduces approximate limits on the ADD fundamental scale, and
projects two transparent luminosity-scaling scenarios.

It is an independent monojet-only reinterpretation, not the official CMS fit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np


RECORD_ID = 106115
VERSION = 2
BASE_LUMINOSITY_FB = 137.0
LATE_LUMINOSITY_FB = 101.0
REDUCED_PLANCK_GEV = 2.4e18
HBAR_C_GEV_M = 1.973269804e-16
DEFAULT_DATA_DIR = Path(__file__).resolve().parent / "data" / "cms-exo-20-004"

# These are the stable internal table IDs exposed by the HEPData record.
TABLES = {
    "add-monojet.json": 1272790,
    "monojet-covariance.json": 1272812,
    "monojet-yields.json": 1272814,
    "add-md-limits.json": 1272859,
}

LIMIT_BRACKETS_TEV = {
    2: (8.0, 15.0),
    3: (6.0, 12.0),
    4: (5.0, 10.0),
    5: (4.5, 9.0),
    6: (4.0, 8.0),
    7: (4.0, 7.0),
}


def normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def add_radius_and_gap(dimensions: int, md_tev: float) -> tuple[float, float]:
    """Return ADD radius (m) and first KK gap (eV) in the PDG M_D convention.

    Uses M_P^2 = R^d M_D^(d+2), with reduced M_P = 2.4e18 GeV. This is a
    model-specific toroidal equal-radius interpretation, not a generic bound on
    every possible extra-dimensional geometry.
    """
    md_gev = md_tev * 1e3
    radius_gev_inverse = (
        REDUCED_PLANCK_GEV**2 / md_gev ** (dimensions + 2)
    ) ** (1.0 / dimensions)
    radius_m = radius_gev_inverse * HBAR_C_GEV_M
    gap_ev = 1e9 / radius_gev_inverse
    return radius_m, gap_ev


def fetch_tables(data_dir: Path) -> list[Path]:
    """Download the four official HEPData tables used by this analysis."""
    data_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for filename, table_id in TABLES.items():
        url = (
            f"https://www.hepdata.net/record/data/"
            f"{RECORD_ID}/{table_id}/{VERSION}/0"
        )
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "spacetime-lab/1.0 (HEPData reinterpretation)"},
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
        decoded = json.loads(payload)
        if "values" not in decoded:
            raise ValueError(f"HEPData response for {filename} contains no values")
        path = data_dir / filename
        path.write_text(json.dumps(decoded, indent=2) + "\n", encoding="utf-8")
        written.append(path)
    return written


def ensure_tables(data_dir: Path) -> None:
    missing = [name for name in TABLES if not (data_dir / name).exists()]
    if missing:
        fetch_tables(data_dir)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def expand_qualifiers(table: dict) -> dict[int, dict[str, str]]:
    """Expand HEPData colspan qualifiers to a mapping for every y group."""
    result: dict[int, dict[str, str]] = {}
    for qualifier_name, entries in table.get("qualifiers", {}).items():
        for entry in entries:
            start = int(entry["group"])
            stop = start + int(entry["colspan"])
            for group in range(start, stop):
                result.setdefault(group, {})[qualifier_name] = str(entry["value"])
    return result


def first_number(text: str) -> float:
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", text)
    if not match:
        raise ValueError(f"No number found in {text!r}")
    return float(match.group())


@dataclass(frozen=True)
class SignalVector:
    values: np.ndarray
    mask: np.ndarray
    years: tuple[int, ...]
    extrapolated: bool


class SignalLibrary:
    """ADD signal templates with binwise log interpolation in M_D."""

    def __init__(self, table: dict):
        attributes = expand_qualifiers(table)
        self.bin_count = len(table["values"])
        self.templates: dict[tuple[int, int, float], np.ndarray] = {}
        for group, attr in attributes.items():
            year = int(attr["Data-taking period"])
            dimensions = int(attr["d"])
            md_tev = first_number(attr["$M_{D}$"])
            values = np.array(
                [float(row["y"][group]["value"]) for row in table["values"]],
                dtype=float,
            )
            self.templates[(year, dimensions, md_tev)] = values

    def available_masses(self, year: int, dimensions: int) -> list[float]:
        return sorted(
            md
            for y, d, md in self.templates
            if y == year and d == dimensions
        )

    def interpolate(
        self, year: int, dimensions: int, md_tev: float
    ) -> tuple[np.ndarray | None, bool]:
        points = sorted(
            (md, values)
            for (y, d, md), values in self.templates.items()
            if y == year and d == dimensions
        )
        if not points:
            return None, False

        masses = np.array([point[0] for point in points], dtype=float)
        values = np.array([point[1] for point in points], dtype=float)
        exact = np.flatnonzero(np.isclose(masses, md_tev, rtol=0, atol=1e-12))
        if exact.size:
            return values[int(exact[0])].copy(), False

        # Outside the simulated grid, preserve the nearest published shape and
        # use the leading ADD normalization s ~ M_D^-(d+2). Extrapolating two
        # finite-Monte-Carlo bins in log space can otherwise amplify sampling
        # noise into an unphysical rising cross section.
        if md_tev < masses[0]:
            scale = (masses[0] / md_tev) ** (dimensions + 2)
            return values[0] * scale, True
        if md_tev > masses[-1]:
            scale = (masses[-1] / md_tev) ** (dimensions + 2)
            return values[-1] * scale, True

        insertion = int(np.searchsorted(masses, md_tev))
        low = max(0, min(len(masses) - 2, insertion - 1))
        high = low + 1
        fraction = math.log(md_tev / masses[low]) / math.log(
            masses[high] / masses[low]
        )
        floor = 1e-12
        interpolated = np.exp(
            (1.0 - fraction) * np.log(np.maximum(values[low], floor))
            + fraction * np.log(np.maximum(values[high], floor))
        )
        return interpolated, False

    def vector(self, dimensions: int, md_tev: float) -> SignalVector:
        blocks: list[np.ndarray] = []
        mask: list[bool] = []
        years: list[int] = []
        extrapolated = False
        for year in (2016, 2017, 2018):
            block, outside = self.interpolate(year, dimensions, md_tev)
            if block is None:
                blocks.append(np.zeros(self.bin_count, dtype=float))
                mask.extend([False] * self.bin_count)
            else:
                blocks.append(block)
                mask.extend([True] * self.bin_count)
                years.append(year)
                extrapolated = extrapolated or outside
        return SignalVector(
            values=np.concatenate(blocks),
            mask=np.array(mask, dtype=bool),
            years=tuple(years),
            extrapolated=extrapolated,
        )


def load_yields(table: dict) -> tuple[list[str], np.ndarray, np.ndarray]:
    names = [str(row["x"][0]["value"]) for row in table["values"]]
    background = np.array(
        [float(row["y"][0]["value"]) for row in table["values"]], dtype=float
    )
    observed = np.array(
        [float(row["y"][1]["value"]) for row in table["values"]], dtype=float
    )
    return names, background, observed


def load_covariance(table: dict, names: list[str]) -> np.ndarray:
    index = {name: position for position, name in enumerate(names)}
    covariance = np.zeros((len(names), len(names)), dtype=float)
    for row in table["values"]:
        first = str(row["x"][0]["value"])
        second = str(row["x"][1]["value"])
        covariance[index[first], index[second]] = float(row["y"][0]["value"])
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues = np.linalg.eigvalsh(covariance)
    if eigenvalues[0] <= 0:
        raise ValueError("Published covariance matrix is not positive definite")
    return covariance


def load_official_limits(table: dict) -> dict[int, dict[str, float]]:
    labels = {
        int(entry["group"]): str(entry["value"])
        for entry in table["qualifiers"]["Quantile"]
    }
    result: dict[int, dict[str, float]] = {}
    for row in table["values"]:
        dimensions = int(row["x"][0]["value"])
        result[dimensions] = {
            labels[group]: float(value["value"])
            for group, value in enumerate(row["y"])
        }
    return result


class SimplifiedLikelihood:
    """Poisson likelihood with a multivariate Gaussian background nuisance."""

    def __init__(
        self,
        background: np.ndarray,
        observed: np.ndarray,
        covariance: np.ndarray,
    ):
        self.background = np.asarray(background, dtype=float)
        self.observed = np.asarray(observed, dtype=float)
        self.covariance = np.asarray(covariance, dtype=float)
        self.inverse_covariance = np.linalg.inv(self.covariance)

    def profile_nll(
        self,
        mu: float,
        signal: np.ndarray,
        data: np.ndarray | None = None,
    ) -> float:
        """Profile the additive correlated background nuisance with Newton steps."""
        signal = np.asarray(signal, dtype=float)
        counts = self.observed if data is None else np.asarray(data, dtype=float)
        delta = np.zeros_like(self.background)

        def objective(candidate: np.ndarray) -> float:
            if not np.all(np.isfinite(candidate)) or np.max(np.abs(candidate)) > 1e20:
                return math.inf
            expectation = self.background + mu * signal + candidate
            if np.min(expectation) <= 0:
                return math.inf
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                value = float(
                    np.sum(expectation - counts * np.log(expectation))
                    + 0.5 * candidate @ self.inverse_covariance @ candidate
                )
            return value if math.isfinite(value) else math.inf

        current = objective(delta)
        for _ in range(100):
            expectation = self.background + mu * signal + delta
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                gradient = (
                    1.0
                    - counts / expectation
                    + self.inverse_covariance @ delta
                )
                hessian = self.inverse_covariance + np.diag(
                    counts / np.square(expectation)
                )
            if not np.all(np.isfinite(gradient)) or not np.all(np.isfinite(hessian)):
                break
            step = np.linalg.solve(hessian, gradient)
            if np.max(np.abs(step) / (1.0 + np.abs(delta))) < 1e-10:
                break

            rate = 1.0
            accepted = False
            while rate > 1e-12:
                candidate = delta - rate * step
                value = objective(candidate)
                if value < current - 1e-10:
                    delta = candidate
                    current = value
                    accepted = True
                    break
                rate *= 0.5
            if not accepted:
                break
        return current

    def best_fit_mu(self, signal: np.ndarray) -> tuple[float, float]:
        objective = lambda mu: self.profile_nll(mu, signal)
        upper = 2.0
        mu_hat, minimum = golden_minimize(objective, 0.0, upper)
        while mu_hat > 0.98 * upper and upper < 32.0:
            upper *= 2.0
            mu_hat, minimum = golden_minimize(objective, 0.0, upper)
        return mu_hat, minimum

    def cls(self, signal: np.ndarray) -> tuple[float, dict[str, float]]:
        """Return the asymptotic modified-frequentist CLs for signal strength 1."""
        mu_hat, best_nll = self.best_fit_mu(signal)
        fixed_nll = self.profile_nll(1.0, signal)
        q_observed = 0.0 if mu_hat > 1.0 else max(0.0, 2.0 * (fixed_nll - best_nll))

        asimov = self.background
        nll_a0 = self.profile_nll(0.0, signal, data=asimov)
        nll_a1 = self.profile_nll(1.0, signal, data=asimov)
        q_asimov = max(0.0, 2.0 * (nll_a1 - nll_a0))

        root_q = math.sqrt(q_observed)
        root_qa = math.sqrt(q_asimov)
        numerator = 1.0 - normal_cdf(root_q)
        denominator = max(normal_cdf(root_qa - root_q), 1e-15)
        cls_value = min(1.0, numerator / denominator)
        expected_cls = min(1.0, 2.0 * (1.0 - normal_cdf(root_qa)))
        q_zero = max(
            0.0,
            2.0 * (self.profile_nll(0.0, signal) - best_nll),
        )
        return cls_value, {
            "expected_cls": expected_cls,
            "mu_hat": mu_hat,
            "q_observed": q_observed,
            "q_asimov": q_asimov,
            "local_z": math.sqrt(q_zero),
        }


def golden_minimize(
    function: Callable[[float], float],
    low: float,
    high: float,
    tolerance: float = 1e-5,
) -> tuple[float, float]:
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = high - ratio * (high - low)
    right = low + ratio * (high - low)
    f_left = function(left)
    f_right = function(right)
    for _ in range(100):
        if high - low <= tolerance * (1.0 + abs(low) + abs(high)):
            break
        if f_left < f_right:
            high, right, f_right = right, left, f_left
            left = high - ratio * (high - low)
            f_left = function(left)
        else:
            low, left, f_left = left, right, f_right
            right = low + ratio * (high - low)
            f_right = function(right)
    optimum = 0.5 * (low + high)
    return optimum, function(optimum)


@dataclass
class AnalysisInputs:
    signals: SignalLibrary
    background: np.ndarray
    observed: np.ndarray
    covariance: np.ndarray
    official_limits: dict[int, dict[str, float]]


def load_inputs(data_dir: Path) -> AnalysisInputs:
    ensure_tables(data_dir)
    signal_table = read_json(data_dir / "add-monojet.json")
    yield_table = read_json(data_dir / "monojet-yields.json")
    covariance_table = read_json(data_dir / "monojet-covariance.json")
    limit_table = read_json(data_dir / "add-md-limits.json")
    names, background, observed = load_yields(yield_table)
    return AnalysisInputs(
        signals=SignalLibrary(signal_table),
        background=background,
        observed=observed,
        covariance=load_covariance(covariance_table, names),
        official_limits=load_official_limits(limit_table),
    )


def likelihood_for_signal(
    inputs: AnalysisInputs,
    vector: SignalVector,
    luminosity_scale: float = 1.0,
    covariance_power: float = 1.0,
    asimov_observed: bool = False,
) -> tuple[SimplifiedLikelihood, np.ndarray]:
    mask = vector.mask
    background = inputs.background[mask] * luminosity_scale
    observed = (
        background.copy()
        if asimov_observed
        else inputs.observed[mask] * luminosity_scale
    )
    covariance = (
        inputs.covariance[np.ix_(mask, mask)]
        * luminosity_scale**covariance_power
    )
    signal = vector.values[mask] * luminosity_scale
    return SimplifiedLikelihood(background, observed, covariance), signal


def cls_at_mass(
    inputs: AnalysisInputs,
    dimensions: int,
    md_tev: float,
    expected: bool,
    luminosity_scale: float = 1.0,
    covariance_power: float = 1.0,
) -> tuple[float, dict[str, float], SignalVector]:
    vector = inputs.signals.vector(dimensions, md_tev)
    likelihood, signal = likelihood_for_signal(
        inputs,
        vector,
        luminosity_scale=luminosity_scale,
        covariance_power=covariance_power,
        asimov_observed=expected,
    )
    value, diagnostics = likelihood.cls(signal)
    if expected:
        value = diagnostics["expected_cls"]
    return value, diagnostics, vector


def solve_limit(
    inputs: AnalysisInputs,
    dimensions: int,
    expected: bool,
    luminosity_scale: float = 1.0,
    covariance_power: float = 1.0,
) -> tuple[float, dict[str, float], SignalVector]:
    low, high = LIMIT_BRACKETS_TEV[dimensions]

    def evaluate(md_tev: float):
        return cls_at_mass(
            inputs,
            dimensions,
            md_tev,
            expected=expected,
            luminosity_scale=luminosity_scale,
            covariance_power=covariance_power,
        )

    low_cls, _, _ = evaluate(low)
    high_cls, _, _ = evaluate(high)
    while low_cls > 0.05 and low > 1.0:
        low *= 0.8
        low_cls, _, _ = evaluate(low)
    while high_cls < 0.05 and high < 100.0:
        high *= 1.25
        high_cls, _, _ = evaluate(high)
    if not (low_cls <= 0.05 <= high_cls):
        raise RuntimeError(
            f"Could not bracket d={dimensions} limit: "
            f"CLs({low:.3g})={low_cls:.3g}, CLs({high:.3g})={high_cls:.3g}"
        )

    for _ in range(32):
        middle = 0.5 * (low + high)
        middle_cls, _, _ = evaluate(middle)
        if middle_cls < 0.05:
            low = middle
        else:
            high = middle
    limit = 0.5 * (low + high)
    _, diagnostics, vector = evaluate(limit)
    return limit, diagnostics, vector


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_analysis(data_dir: Path, output: Path, projection_output: Path) -> dict:
    inputs = load_inputs(data_dir)
    rows = []
    projections = []
    for dimensions in range(2, 8):
        observed_limit, observed_diag, vector = solve_limit(
            inputs, dimensions, expected=False
        )
        expected_limit, expected_diag, _ = solve_limit(
            inputs, dimensions, expected=True
        )
        official = inputs.official_limits[dimensions]
        radius_m, kk_gap_ev = add_radius_and_gap(dimensions, observed_limit)
        rows.append(
            {
                "dimensions": dimensions,
                "bins_used": int(np.count_nonzero(vector.mask)),
                "years_used": "+".join(str(year) for year in vector.years),
                "computed_observed_md_tev": f"{observed_limit:.5f}",
                "official_observed_md_tev": f"{official['Observed']:.5f}",
                "observed_difference_percent": f"{100.0 * (observed_limit / official['Observed'] - 1.0):.3f}",
                "computed_expected_md_tev": f"{expected_limit:.5f}",
                "official_expected_md_tev": f"{official['Median Expected']:.5f}",
                "expected_difference_percent": f"{100.0 * (expected_limit / official['Median Expected'] - 1.0):.3f}",
                "best_fit_signal_strength_at_limit": f"{observed_diag['mu_hat']:.4f}",
                "local_z_at_limit": f"{observed_diag['local_z']:.4f}",
                "add_radius_upper_m": f"{radius_m:.8e}",
                "kk_gap_lower_ev": f"{kk_gap_ev:.8e}",
                "interpolation_extrapolated": str(vector.extrapolated).lower(),
            }
        )

        base_luminosity = (
            BASE_LUMINOSITY_FB if 2016 in vector.years else LATE_LUMINOSITY_FB
        )
        scale = 3000.0 / base_luminosity
        for scenario, covariance_power in (
            ("statistics-like covariance", 1.0),
            ("fixed fractional background uncertainty", 2.0),
        ):
            projected_limit, _, projected_vector = solve_limit(
                inputs,
                dimensions,
                expected=True,
                luminosity_scale=scale,
                covariance_power=covariance_power,
            )
            projections.append(
                {
                    "dimensions": dimensions,
                    "scenario": scenario,
                    "target_luminosity_fb": "3000",
                    "projected_expected_md_tev": f"{projected_limit:.5f}",
                    "template_extrapolated": str(projected_vector.extrapolated).lower(),
                }
            )

    write_csv(
        output,
        [
            "dimensions",
            "bins_used",
            "years_used",
            "computed_observed_md_tev",
            "official_observed_md_tev",
            "observed_difference_percent",
            "computed_expected_md_tev",
            "official_expected_md_tev",
            "expected_difference_percent",
            "best_fit_signal_strength_at_limit",
            "local_z_at_limit",
            "add_radius_upper_m",
            "kk_gap_lower_ev",
            "interpolation_extrapolated",
        ],
        rows,
    )
    write_csv(
        projection_output,
        [
            "dimensions",
            "scenario",
            "target_luminosity_fb",
            "projected_expected_md_tev",
            "template_extrapolated",
        ],
        projections,
    )
    return {"reproduction": rows, "projection": projections}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("fetch", help="download the official HEPData tables")

    analyze = subparsers.add_parser(
        "analyze", help="reproduce current limits and project 3 ab^-1 sensitivity"
    )
    analyze.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "collider-limit-reproduction.csv",
    )
    analyze.add_argument(
        "--projection-output",
        type=Path,
        default=Path(__file__).resolve().parent / "collider-hllhc-projection.csv",
    )

    benchmark = subparsers.add_parser(
        "benchmark", help="evaluate one ADD model point"
    )
    benchmark.add_argument("--dimensions", type=int, choices=range(2, 8), required=True)
    benchmark.add_argument("--md", type=float, required=True, help="M_D in TeV")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "fetch":
        paths = fetch_tables(args.data_dir)
        print(json.dumps({"downloaded": [str(path) for path in paths]}, indent=2))
        return

    inputs = load_inputs(args.data_dir)
    if args.command == "benchmark":
        cls_value, diagnostics, vector = cls_at_mass(
            inputs, args.dimensions, args.md, expected=False
        )
        print(
            json.dumps(
                {
                    "dimensions": args.dimensions,
                    "md_tev": args.md,
                    "cls": cls_value,
                    "excluded_at_95_percent": cls_value < 0.05,
                    "bins_used": int(np.count_nonzero(vector.mask)),
                    "years_used": vector.years,
                    "template_extrapolated": vector.extrapolated,
                    **diagnostics,
                },
                indent=2,
            )
        )
        return

    result = run_analysis(
        args.data_dir,
        output=args.output,
        projection_output=args.projection_output,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
