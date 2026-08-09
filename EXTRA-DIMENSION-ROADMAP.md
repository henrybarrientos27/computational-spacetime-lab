# Realistic Fourth-Spatial-Dimension Roadmap

## What “sending something through” would mean

If there is one compact fourth spatial direction `w` with radius `R`, spacetime
would be five-dimensional: three familiar spatial directions, `w`, and time.
Motion around the compact direction is quantized:

```text
p_w = n hbar/R
E_n² = E_0² + (n hbar c/R)².
```

To observers who cannot resolve `w`, that momentum looks like an additional
particle mass. The excitations are called Kaluza-Klein modes. Consequently, the
realistic action is not opening a geometric door. It is creating a particle
state with nonzero `p_w`—but only if that particle's field is allowed to exist
away from our three-dimensional brane.

The 2025 Particle Data Group review describes this mode tower and the crucial
model split. In ADD-type models, gravity propagates through the higher-
dimensional bulk while Standard Model fields are localized on a 3-brane. In
universal-extra-dimension models, Standard Model fields can have KK modes, but
collider constraints put the compactification energy above roughly 1.4–1.5
TeV in representative minimal models. [PDG extra-dimensions review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-extra-dimensions.pdf)

## What the computer says

For the first KK mode,

```text
E_gap = hbar c/R.
```

- At the model-specific short-range-gravity limit `R = 30 micrometers`, the gap
  is only **0.00658 eV**. This does not make an optical doorway: in the relevant
  ADD model, photons and atoms are still brane-confined, and the bulk messenger
  is gravitational.
- A representative universal-extra-dimension bound `1/R > 1.5 TeV` corresponds
  to `R < 1.32e-19 m`. The first mode then requires at least **1.5 TeV** and is a
  particle excitation produced in a collider, not a transported macroscopic
  object.
- A spacecraft cannot be assigned a KK mode unless every field binding its
  atoms is a bulk field. No experimentally supported theory provides that.

## The experiment ladder

### 1. Computer experiment—available now

Simulate transfer through a synthetic coordinate. The included model uses an
engineered chain of `N` internal states with adjacent couplings

```text
J_j proportional to sqrt((j+1)(N-1-j)).
```

A state initialized at synthetic site zero reaches site `N-1` with 100%
idealized probability at normalized time `pi/2`. This is an exact, testable
state-transfer protocol.

### 2. Tabletop synthetic-dimension experiment—physically buildable

Encode the synthetic coordinate in photon frequencies, waveguide modes, or
atomic internal states. Prepare a photon or atom at `w=0`, apply the engineered
mode couplings, and read out population at `w=N-1`. Photonic waveguide arrays
and ultracold atoms have already reproduced responses associated with four-
dimensional quantum Hall systems. [Photonic experiment](https://www.nature.com/articles/nature25011), [ultracold-atom experiment](https://www.nature.com/articles/nature25000)

This genuinely sends a quantum state through a synthetic dimension, but the
carrier remains inside ordinary space. A 2025 potassium-atom experiment also
observed a four-dimensional Anderson transition using three synthetic
dimensions, demonstrating how far this emulation can be pushed. [4D synthetic-dimension phase transition](https://www.nature.com/articles/s41467-025-57396-3)

### 3. Search for an actual bulk messenger

The most defensible candidate is a graviton or another presently unknown field
that is not brane-confined. Two relevant signatures are:

- A deviation from the inverse-square gravitational potential below the compact
  radius. For a two-torus ADD example, the PDG quotes `R < 30 micrometers` at
  95% confidence for the referenced short-range test.
- Collision events with missing transverse momentum consistent with emission
  of a KK graviton into the bulk, or a reproducible tower of resonances with
  masses proportional to `n/R`.

Neither experiment is feasible with a normal home computer alone. The computer
can analyze public data and calculate signals, but producing the relevant
collisions or measuring gravity at tens of micrometers requires specialized
equipment.

## Evidence standard

A disappearing signal is not enough. A genuine extra dimension should produce
several mutually consistent observations:

1. Multiple KK levels with the same `1/R` spacing.
2. Energy and momentum conservation including the inferred bulk momentum.
3. The same `R` inferred from collider and short-range-gravity measurements.
4. Dependence on particle species matching which fields the model places in the
   bulk.
5. Replication with detector backgrounds and ordinary invisible particles
   excluded.

## Bottom line

We can realistically send information through a **synthetic** fourth coordinate
today. Sending a real particle through a physical fourth spatial direction
requires first discovering a bulk-capable field. In the leading gravity-only
model, the candidate payload is gravitational radiation—not atoms or objects.
The next useful computer project is therefore a missing-momentum/KK-tower
analysis, not a macroscopic portal design.
