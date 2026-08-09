# Computational Spacetime Lab

This is a small, auditable numerical laboratory for testing proposed “4D
travel” configurations against known physics. It does not claim to create a
warp drive or wormhole. Its purpose is to identify the first failed physical
requirement and quantify how large the failure is.

See [FINDINGS.md](FINDINGS.md) for the first parameter search and its strongest
results.

The lab contains these calculations:

- Experimentally established relativistic time dilation and kinetic energy.
- A symmetric constant-proper-acceleration interstellar trip.
- The negative Eulerian energy implied by the original Alcubierre metric.
- Throat stress-energy and traveler tidal force for a simple Morris-Thorne
  wormhole with `b(r)=r0²/r` and zero redshift function.
- Ideal parallel-plate Casimir energy density for scale comparison.
- Kaluza-Klein excitation thresholds for a compact fourth spatial dimension.
- Perfect state transfer through a finite synthetic dimension.
- A data-driven monojet search for ADD gravitons using the official CMS
  EXO-20-004 HEPData yields, covariance matrix, and signal templates.
- Coherent neutron/hidden-neutron oscillations and an optimized two-converter
  neutron-shining-through-a-wall experiment.
- A conditional inference gate from a swapping probability to an interbrane
  coupling and distance, with explicit degeneracy across geometric assumptions.
- A beamline engineering audit with velocity spread, magnetic-field
  nonuniformity, neutron decay, and gravitational aperture.

## Run it

Python 3.10 or newer is recommended. Install the one numerical dependency
before running the full test suite:

```bash
cd computational-spacetime-lab
python3 -m pip install -r requirements.txt
python3 -m unittest -v
python3 spacetime_lab.py mission --mass 1000 --beta 0.99 --distance 4.2465
python3 spacetime_lab.py accelerated --distance 4.2465 --acceleration 9.80665
python3 spacetime_lab.py warp --radius 100 --wall 10 --beta 1
python3 spacetime_lab.py optimize-warp --cavity-radius 20 --beta 1
python3 spacetime_lab.py warp-bound --cavity-radius 20 --outer-radius 100 --beta 1
python3 spacetime_lab.py wormhole --radius 100 --beta 0.000001
python3 spacetime_lab.py casimir --separation 0.000001
python3 spacetime_lab.py sweep --output spacetime-sweep.csv
python3 spacetime_lab.py kk --radius 1.3155e-19 --mode 1
python3 spacetime_lab.py synthetic --sites 9 --tau 1.57079632679
python3 spacetime_lab.py extra-sweep --output extra-dimension-sweep.csv
```

The collider analysis uses NumPy:

```bash
python3 collider_extra_dimension.py fetch
python3 collider_extra_dimension.py analyze
python3 collider_extra_dimension.py benchmark --dimensions 2 --md 10.7
python3 hidden_neutron_lab.py coherent --tau 10 --time 1
python3 hidden_neutron_lab.py optimize --tau 100
python3 hidden_neutron_lab.py sweep
python3 brane_geometry.py infer
python3 brane_geometry.py sweep
python3 beamline_realism.py optimize --tau 10 --velocity 5
python3 beamline_realism.py target --tau 100 --target-days 365.25
python3 beamline_realism.py sweep
```

See [COLLIDER-SEARCH-REPORT.md](COLLIDER-SEARCH-REPORT.md) for the meaning and
limits of this reinterpretation.

See [HIDDEN-NEUTRON-TRANSFER.md](HIDDEN-NEUTRON-TRANSFER.md) for the closest
known reversible matter-transfer hypothesis and the experiment needed to test
it.

See [BRANE-GEOMETRY-GATE.md](BRANE-GEOMETRY-GATE.md) for the test that prevents
a generic invisible neutron state from being mislabeled as travel through a
fourth spatial dimension.

See [REALISTIC-PORTAL-BEAMLINE.md](REALISTIC-PORTAL-BEAMLINE.md) for the
gravity, beam-spread, and field-uniformity audit of the reversible apparatus.

Outputs are JSON so they can be piped into other analysis tools. The sweep is a
CSV containing representative points from all three travel models.

The Kaluza-Klein command assumes the selected field is allowed to enter a
compact bulk. This is not true for ordinary matter in common gravity-only brane
models. The synthetic command models a buildable transfer between internal
states or photonic modes; it does not move matter out of physical 3D space.

## What the numbers mean

### Relativistic mission

The cruise model intentionally excludes acceleration, braking, propulsion
efficiency, shielding mass, and the relativistic rocket equation. Consequently
its energy is an absolute lower bound, not a vehicle estimate. “Ideal
antimatter” means perfect conversion of equal masses of antimatter and matter;
real production, storage, and propulsion would require much more.

The constant-proper-acceleration model accelerates for half the coordinate
distance and reverses acceleration for the other half. It reports elapsed time
on Earth and aboard the ship plus the midpoint speed.

### Alcubierre bubble

The program evaluates the standard negative energy density measured by
Eulerian observers on a spatial slice of Alcubierre’s original metric. It uses
the smooth top-hat shape from the paper and numerically integrates the negative
energy over the slice. This is not the energy necessary to form, steer, or
destroy the bubble. A metric is a target spacetime geometry, not a mechanism.

`wall` is the characteristic scale `1/sigma` of the hyperbolic-tangent wall;
it is not a universally defined physical wall thickness. Values with bubble
radius much smaller than the wall scale are mathematical explorations rather
than useful vehicle geometries.

`optimize-warp` searches wall scales and bubble radii for the smallest negative
slice energy that still leaves the requested cabin radius at least 99% flat.
It is useful for ruling out the idea that a simple parameter adjustment solves
the energy problem; it still does not model bubble formation or stability.

`warp-bound` goes one step further: it uses a variational minimum to calculate
the least possible negative spatial-slice energy for *any* spherical shape
function in the original Alcubierre metric, given a flat cabin radius and an
optional outer support radius. This makes the result independent of the chosen
hyperbolic-tangent wall profile.

### Morris-Thorne throat

The wormhole calculation assumes the throat and its nontrivial topology already
exist. The chosen model avoids a redshift horizon, but has negative energy at
the throat and violates the null energy condition. It does not test dynamical
stability, quantum backreaction, mouth creation, causality, or radiation.
The reported two-sided NEC integral includes pressure as well as energy density;
it is an exoticity measure, not a reservoir of extractable energy.

### Casimir comparison

The Casimir formula assumes ideal, perfectly conducting parallel plates. At
very small gaps, finite conductivity, surface roughness, temperature, and
material structure matter. The value is included only to compare an accepted
negative-energy phenomenon with the stress-energy scales implied by the toy
metrics.

## Scientific stopping rule

A candidate shortcut should not be called physically possible until one model
simultaneously provides:

1. A realizable stress-energy tensor.
2. A causal formation history from ordinary initial conditions.
3. Stability under classical and quantum perturbations.
4. An external and internal control procedure.
5. Human-safe tidal forces and radiation.
6. A topology and endpoint consistent with conservation laws.

The current warp and wormhole examples fail the first items before ordinary
vehicle engineering begins. Relativistic future travel does not, but its energy
and shielding requirements remain extreme.

## Verification and status

The repository contains 36 deterministic tests across the relativistic,
warp-metric, wormhole, collider, hidden-neutron, brane-geometry, and beamline
models. The calculations are research and educational models, not a published
new physical law or evidence that a portal exists.

AI-assisted development contributed substantially. See [AI_USAGE.md](AI_USAGE.md)
for the authorship boundary.
