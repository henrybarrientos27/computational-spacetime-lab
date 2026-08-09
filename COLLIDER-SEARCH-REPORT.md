# Collider bounds on extra spatial dimensions

## Bottom line

We cannot currently make ordinary matter leave three-dimensional space. There
is no experimentally established extra spatial dimension, no known portal that
couples a payload to one, and no demonstrated mechanism that detaches Standard
Model matter from our three-dimensional brane and returns it intact.

The strongest physically real version of the idea is nevertheless testable:
produce a particle that is allowed to propagate in a higher-dimensional bulk
and look for its momentum leaving the detector. In the ADD model that particle
is a Kaluza–Klein graviton. This is not a metaphor—the missing-momentum channel
is the closest existing experiment to sending a physical excitation into an
extra spatial dimension.

I built a local, reproducible analysis of that experiment. It finds no evidence
of passage. It does reproduce the published exclusion boundary well enough to
define the next search frontier.

## Two corrections that change the engineering problem

1. Antimatter is not negative mass or negative gravitational energy. The
   ALPHA-g experiment observed antihydrogen behaving consistently with
   attraction toward Earth and ruled out repulsive gravity of magnitude 1g.
   Antimatter can release enormous positive energy when it annihilates, but it
   is not a known way to curve spacetime with the opposite sign. See the
   [ALPHA-g paper](https://www.nature.com/articles/s41586-023-06527-1).
2. Curving four-dimensional spacetime and moving along another spatial
   coordinate are different operations. General relativity supplies curvature;
   an extra-dimensional theory must additionally supply the coordinate, its
   geometry, a field that can enter it, and a coupling that moves a controlled
   state off our brane.

In the ADD scenario, Standard Model fields are localized on a 3-brane and
gravity propagates in the bulk. Therefore an ADD graviton can be a bulk
messenger, but an atom, probe, or person cannot follow it under the model's own
assumptions. The current [Particle Data Group review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-extra-dimensions.pdf)
states this setup explicitly and describes missing transverse momentum from
escaping gravitons as its collider signature.

## What the computer analysis actually does

`collider_extra_dimension.py` retrieves four official HEPData tables from the
CMS EXO-20-004 search:

- 66 observed monojet bins from 2016, 2017, and 2018;
- the corresponding Standard Model background prediction;
- the full 66 × 66 correlated background covariance matrix;
- simulated ADD graviton signal spectra at multiple values of the number of
  extra dimensions `d` and fundamental scale `M_D`;
- CMS's published observed and expected limits for comparison.

For every candidate point it profiles the likelihood

```text
L(mu, delta) = product_i Poisson(n_i | mu*s_i + b_i + delta_i)
               × Normal(delta | 0, covariance)
```

and computes an asymptotic modified-frequentist CLs value. Signal templates are
interpolated bin by bin inside the simulated grid. Beyond that grid, the nearest
published shape is scaled with the leading ADD relation
`signal ∝ M_D^-(d+2)` so Monte Carlo noise cannot create a rising cross section.

This is a monojet-only simplified-likelihood reinterpretation. The official CMS
result combines more categories and nuisance information, so agreement—not
identity—is the appropriate validation target. Source data are at the
[HEPData record](https://www.hepdata.net/record/ins1894408), and the official
analysis summary is on the [CMS EXO-20-004 page](https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-20-004/).

## Present result

The independent calculation reproduces the official 95% confidence limits to
within 1.4% for `d=2–6`. For `d=7`, only 2017–2018 ADD signal templates are
public, so the calculation uses 44 rather than 66 bins and lands 2.4% below the
official observed limit.

| d | Computed observed M_D lower limit | Official | Computed median expected | Local fit Z* | ADD radius upper estimate | First KK gap lower estimate |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 10.638 TeV | 10.692 TeV | 12.358 TeV | 1.94 | 4.19 µm | 0.047 eV |
| 3 | 8.025 TeV | 8.035 TeV | 9.091 TeV | 1.81 | 0.110 nm | 1.79 keV |
| 4 | 6.728 TeV | 6.748 TeV | 7.417 TeV | 1.80 | 0.554 pm | 356 keV |
| 5 | 6.010 TeV | 6.005 TeV | 6.564 TeV | 1.88 | 0.0227 pm | 8.68 MeV |
| 6 | 5.478 TeV | 5.501 TeV | 5.970 TeV | 1.82 | 2.74 fm | 72.1 MeV |
| 7 | 5.037 TeV | 5.159 TeV | 5.397 TeV | 1.68 | 0.612 fm | 322 MeV |

\*The local `Z` is the one-template background-only likelihood ratio at the
reported boundary. It is not a global discovery significance. The official
analysis reports no significant excess. Values below about 2σ, before a
look-elsewhere correction, are ordinary fluctuations—not evidence of a bulk.

The radius and gap columns use the PDG ADD convention
`M_P² = R^d M_D^(d+2)` with reduced `M_P = 2.4×10^18 GeV`. They apply only to a
flat equal-radius ADD compactification. They are not generic geometric bounds
on every higher-dimensional theory.

## What 3 ab⁻¹ could test

I projected the same analysis to 3000 fb⁻¹ using two deliberately separated
bookends:

- **Statistics-like covariance:** the background covariance grows linearly
  with luminosity, so relative uncertainty continues to improve.
- **Fixed fractional uncertainty:** the covariance grows as luminosity squared,
  so present fractional background systematics do not improve.

| d | Expected M_D reach, statistics-like | Expected M_D reach, fixed fractional uncertainty |
|---:|---:|---:|
| 2 | 18.37 TeV | 14.08 TeV |
| 3 | 12.40 TeV | 10.10 TeV |
| 4 | 9.65 TeV | 8.11 TeV |
| 5 | 8.20 TeV | 7.11 TeV |
| 6 | 7.30 TeV | 6.41 TeV |
| 7 | 6.67 TeV | 5.76 TeV |

These are sensitivity scenarios, not official HL-LHC forecasts. They do not
model the energy increase to 14 TeV, trigger or detector upgrades, a new event
selection, or a decomposition of present covariance into statistical and
systematic sources. Several points extrapolate beyond the published signal
grid. Their purpose is to reveal the governing bottleneck: better control of
correlated backgrounds is worth almost as much as the added collision count.

## The passage problem, decomposed

A real transport device needs every link below. Detecting an extra dimension
would satisfy only the first.

| Gate | Required observation | Current status | Stop condition |
|---|---|---|---|
| 0. Existence | A reproducible spectrum requiring an extra spatial coordinate | Not observed | No global discovery in independent channels |
| 1. Bulk access | A known field demonstrably leaves and re-enters the brane | ADD predicts graviton escape only; no return control | Missing energy alone gives no destination or recovery |
| 2. Payload coupling | A tunable interaction transfers a non-gravitational state off-brane | No known coupling for ordinary matter | Standard Model localization forbids payload motion in minimal ADD |
| 3. Coherence | Identity or quantum state survives an out-and-back transfer | Not tested because Gate 2 is absent | Decoherence, decay, or irreversible leakage |
| 4. Conservation and targeting | Energy, momentum, charge, and endpoint are measured on both sides | No protocol | Any unexplained nonconservation or unaddressable endpoint |
| 5. Scaling | The process works from one excitation to atoms to macroscopic matter | No mechanism | Energy, tidal force, radiation, or instability diverges |

The key conceptual trap is treating missing momentum as transportation. A
graviton that irreversibly escapes a detector is closer to throwing away an
untracked bit than moving a package through a shortcut.

## The shortest honest route forward

1. **Keep Gate 0 falsifiable.** Re-run this likelihood whenever new public
   monojet data or covariance models appear. Add photon-plus-missing-momentum,
   dilepton/diphoton virtual-graviton, and short-range gravity constraints to a
   single joint parameter scan.
2. **Demand cross-channel geometry.** A candidate must yield one consistent
   `d`, `M_D`, KK spectrum, and rate in at least two independent production
   channels. A generic invisible particle can mimic one monojet excess.
3. **Require discovery evidence before portal engineering.** Use a global
   significance of at least 5σ, independent replication, detector-systematic
   closure, and a predictive spectrum in held-out data.
4. **If and only if Gates 0–1 pass, search for reversible coupling.** The first
   device target should be coherent transfer of a single controllable quantum
   state out and back, with complete state tomography—not a macroscopic object.
5. **Do not substitute a synthetic dimension for a physical one.** Synthetic
   dimensions are valuable quantum simulators and can test control protocols,
   but their “sites” are internal modes in ordinary space.

## Reproduce it

```bash
cd computational-spacetime-lab
python3 -m unittest -v
python3 collider_extra_dimension.py fetch
python3 collider_extra_dimension.py analyze
python3 collider_extra_dimension.py benchmark --dimensions 2 --md 10.7
```

Generated results are in `collider-limit-reproduction.csv` and
`collider-hllhc-projection.csv`. The cached source tables remain under
`data/cms-exo-20-004/` so every reported number can be traced back to an
official bin and covariance entry.
