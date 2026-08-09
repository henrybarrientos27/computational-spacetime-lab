# Computational Findings

These results were generated on July 15, 2026 with `spacetime_lab.py`. They are
tests of published geometries, not evidence that those geometries can be built.

## Result 1: forward time travel remains the viable route

A one-gravity trip to Alpha Centauri, accelerating for half of 4.2465
light-years and decelerating for the other half, gives:

- Earth elapsed time: **5.873 years**
- Traveler elapsed time: **3.542 years**
- Midpoint speed: **0.94965 c**
- Midpoint Lorentz factor: **3.1918**

This is a real spacetime effect and does not require exotic matter. The model
does not yet include propulsion efficiency, reaction mass, heat rejection,
radiation, or dust shielding.

At a constant `0.99 c`, a 1,000 kg probe crossing 4.27 light-years has a kinetic
energy lower bound of `5.47e20 J`. A one-microgram grain impacts with `5.47e8 J`,
and an interstellar proton carries `5.71 GeV` in the probe frame. Those impact
numbers make shielding a first-class problem rather than an afterthought.

## Result 2: optimizing the original warp shape does not rescue it

The optimizer required a 20 m-radius cabin in which Alcubierre's shape function
remains at least 99% of its central value. Searching bubble radii and wall
scales produced this minimum for a coordinate speed of `1 c`:

- Nominal bubble radius: **40.04 m**
- Wall scale `1/sigma`: **8.76 m**
- Characteristic outer radius `R + 6/sigma`: **92.60 m**
- Integrated negative energy on one spatial slice: **-6.25e44 J**
- Mass equivalent: **6.95e27 kg**, or **0.00350 solar masses**
- Peak negative energy density: **-3.92e39 J/m³**

The mass equivalent is about 3.7 Jupiter masses, all with the wrong sign. This
is already an optimized lower value for the chosen metric slice; it excludes
formation energy, stability, steering, and shutdown.

The numerical scaling exposes why parameter tweaking fails. For geometrically
similar optimized bubbles, the integrated negative energy scales approximately
as

```text
|E_negative| proportional to beta² × usable-cavity-radius.
```

A 1 m flat cavity still requires `3.12e43 J`, with a mass equivalent of
`3.48e26 kg`. For the 20 m cavity to have only a 1,000 kg mass equivalent, its
coordinate speed would need to be about `3.79e-13 c`, or `0.000114 m/s`. That is
not a useful warp drive, and superluminal values make the negative energy grow.

There is also a profile-independent lower bound inside the original spherical
Alcubierre class. Minimizing `integral r²(f')² dr` for a flat cabin radius `a`
and an outer wall radius `b` gives

```text
I_min = ab/(b-a),       E_min = -(c⁴/G) beta² I_min / 12.
```

Even allowing the wall to extend without limit gives `I_min -> a`, not zero.
For a 20 m cabin at `1 c`, absolutely any spherical wall profile therefore
requires at least **-2.02e44 J**, with a mass-equivalent magnitude of
**2.24e27 kg**—about 1.2 Jupiter masses. Confining the wall inside 100 m raises
the bound to **-2.52e44 J**. The optimized smooth tanh profile is roughly three
times above the unlimited-support mathematical bound.

## Result 3: making a wormhole larger trades one failure for another

For the zero-redshift Morris-Thorne example `b(r)=r0²/r`, a 100 m throat crossed
at `1e-6 c` gives:

- Throat `rho + p_r`: **-9.63e38 J/m³**
- Lateral tidal acceleration across 2 m: **17.98 m/s²**, or **1.83 g**
- Maximum traversal speed for one-g lateral tides: **221 m/s**
- Two-sided proper-volume integral of `rho + p_r`: **-3.80e46 J**

Increasing the throat radius helps local passage conditions: energy density and
tidal acceleration fall as `1/r0²`. But the integrated null-energy violation
grows linearly with `r0`. A larger throat therefore makes local conditions
gentler while worsening the global exoticity measure.

For comparison, ideal conducting plates separated by 1 micrometer have a
Casimir energy density of only `-4.33e-4 J/m³`. A wormhole throat would need a
radius of about **15.8 million light-years** merely to reduce its local NEC
violation to that magnitude. Even comparing with the ideal formula at a 1 nm
gap gives a required throat radius of about **15.8 light-years**, while the
integrated violation increases with that radius. Real plates add material and
finite-conductivity limitations, so this comparison is generous.

## What the computer ruled out

Within the original Alcubierre metric and this simple Morris-Thorne geometry,
no macroscopic parameter adjustment reaches laboratory negative-energy scales.
The obstacle is not a missing amount of ordinary positive energy or antimatter;
it is the sign and distribution of the required stress-energy, followed by
formation and stability.

## Best next computations

1. Add a relativistic shield model using measured interstellar gas and dust
   distributions, including secondary particle showers and heat deposition.
2. Add propulsion models—beamed sail, fusion, antimatter-catalyzed fusion, and
   magnetic braking—and optimize payload fraction and travel time.
3. Generalize the metric engine so a numerical optimizer searches stress-energy
   tensors while enforcing energy conditions, flat cabin size, tidal limits,
   asymptotic flatness, and causal initial data simultaneously.
4. Add linear perturbation evolution. A candidate should be rejected if small
   field or matter perturbations grow uncontrollably.
5. Compare quantum-energy-inequality bounds directly with every negative-energy
   pulse used by a proposed metric.

## Reference starting points

- [Alcubierre's original warp metric](https://arxiv.org/abs/gr-qc/0009013)
- [Generic warp drives and energy-condition violations](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.105.064038)
- [Morris, Thorne, and Yurtsever on wormholes and causality](https://authors.library.caltech.edu/records/m644f-tbz27)
- [Quantum-energy-inequality restrictions on wormholes](https://arxiv.org/abs/2405.05963)
- [ALPHA measurement of antihydrogen gravity](https://www.nature.com/articles/s41586-023-06527-1)

The current actionable direction is therefore a relativistic time-skip vehicle,
while warp and wormhole calculations serve as falsification tests for proposed
new physics.
