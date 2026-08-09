# Can neutron transfer prove travel through a fourth spatial dimension?

## Gate result

Not from a swapping probability alone.

A neutron disappearance/reappearance signal would be an extraordinary and
potentially reversible portal, but the observable does not uniquely specify a
second brane, an interbrane distance, or even an extra spatial dimension. The
inference has three separate layers:

```text
measured p
  + assumed energy detuning and collision width
    -> mixing energy epsilon
      + assumed interbrane vector-potential difference
        -> geometric coupling g
          + assumed brane scale, bulk dimension, metric, and topology
            -> physical separation d
```

Every arrow after the first introduces parameters that the swapping count does
not measure. `brane_geometry.py` now implements this chain explicitly and
returns no distance when the chosen experiment/model combination has no
geometric reach.

## First layer: what an experiment can measure

The STEREO two-state treatment gives

```text
p = 2 epsilon² /
    ((Delta E + V_F)² + 4 epsilon² + (hbar Gamma/2)²).
```

Consequently, `p` does not determine `epsilon` unless the effective detuning
`Delta E + V_F` and collision rate `Gamma` are known. STEREO found
`p < 3.1×10^-11` at 95% confidence. Under the specific two-brane benchmark
`Delta E = 2 keV`, its limit becomes `epsilon < 7.9 meV`. That numerical
agreement is reproduced by the program.

The important distinction is that neither `p` nor `epsilon` yet says that a
particle occupied a different position along an extra coordinate. The same
two-level Hamiltonian describes mirror neutrons, sterile neutrons, and other
hidden sectors in ordinary effective field theory.

## Second layer: converting mixing into an interbrane coupling

One published two-brane model uses

```text
epsilon = |mu_n| B_perp
B_perp = g |A_+ - A_-|,
```

so the inferred geometric coupling is

```text
g = epsilon / (|mu_n| |A_+ - A_-|).
```

The vector-potential difference between the visible and hidden branes is not
measured. The model literature discusses astrophysical magnitudes from roughly
`10^9` to `10^12 T m`, a factor of one thousand. The same measured `epsilon`
therefore permits a factor of one thousand in `g` before a bulk geometry is
even selected.

## Third layer: converting coupling into distance

For one particular non-compact 5D DGP-brane ansatz, the published map is

```text
g = (m_q² / M_B) exp(-m_q d)                 (hbar = c = 1),
```

where the model fixes the constituent-quark mass to `m_q = 340 MeV`, `M_B` is
an effective brane energy/thickness scale, and `d` is the physical separation
in the bulk. The inverse used by the lab is

```text
d = (hbar c / m_q) ln[(m_q² / (M_B hbar c)) / g].
```

This expression is not universal. The same paper derives different functions
for a compact `S1/Z2` bulk, a non-compact 6D bulk, and a two-dimensional torus.
Some compact brane locations even make `g` exactly zero regardless of nearby
distance. Its conclusion is therefore that a null disappearance experiment
cannot rule out hidden branes universally.

## Reproduced Planck-scale benchmark

Using all of the following assumptions at once:

- `p < 3.1×10^-11`;
- effective detuning `Delta E = 2 keV`;
- no collision broadening;
- `|A_+ - A_-| = 2×10^9 T m`;
- non-compact 5D DGP form;
- `M_B = 1.22089×10^28 eV` (Planck energy);
- constituent-quark mass `340 MeV`;

the program obtains

```text
epsilon_limit = 7.874×10^-3 eV
g_limit       = 6.528×10^-5 m^-1
g(d=0)        = 4.798×10^-5 m^-1
```

The experimental upper limit is above the largest coupling this specific
Planck-scale model predicts, so **there is no distance constraint**. This is a
useful negative result: the present null measurement has not yet reached even
the zero-separation edge of that benchmark.

Changing only the unknown vector-potential difference to `2×10^12 T m` makes
the same count limit appear to imply `d > 3.83 fm`. That is not new evidence;
it is the consequence of an unmeasured assumption. In addition, the cited
phenomenological model treats neutron exchange beyond roughly `0.5 fm` as
effectively precluded, so a computed bound beyond that range should be read as
“this parameter choice is excluded,” not as a ruler locating a hidden brane.

The full 18-point assumption table is in
`brane-geometry-degeneracy.csv`. It shows orders-of-magnitude changes in
`epsilon`, `g`, and inferred `d` while keeping the exact same probability
limit.

## What would establish geometry rather than a generic portal

A credible fourth-spatial-dimension claim would require several independent
signatures that share one geometric parameter set:

1. Simultaneous upstream disappearance and downstream regeneration with the
   complete two-state resonance line shape, fixing both nonzero mixing and
   detuning.
2. Reproducible phase evolution with controlled interaction time, establishing
   coherent transport rather than an unknown neutron background.
3. A second particle species or bound system with coupling ratios predicted by
   the same localization model; an arbitrary sterile-neutron interaction would
   not be enough.
4. An independent geometric spectrum—such as multiple Kaluza-Klein masses with
   integer spacing—giving the same compactification scale.
5. Consistency with short-range gravity and collider missing-momentum limits,
   using one bulk dimension, metric, and brane scale rather than refitting each
   experiment separately.

Only after those checks could `d` be treated as a physical displacement. Before
them, “the neutron went through the fourth dimension” is one interpretation of
a hidden-state transition, not an experimental conclusion.

## Consequence for transport

The neutron route remains the strongest realistic first portal because its
ideal Hamiltonian permits a full reversible swap. It does not yet transport an
atom. A complete object would need coherent partners and matched transfer for
protons, electrons, binding energy, charge, spin, and every environmental
degree of freedom without decoherence. No model or experiment currently
provides that map.

Antimatter does not remove this obstacle. Antimatter has positive inertial
energy and is not a known source of the negative stress-energy required by
traversable-wormhole or warp geometries. Annihilation can supply energy, but it
does not choose a direction “upward” out of three-dimensional space.

## Run the geometry gate

```bash
cd /Users/henrybarrientos/Documents/Codex/2026-07-15/i-w/outputs/spacetime-lab
python3 brane_geometry.py infer
python3 brane_geometry.py forward --separation-m 1e-16
python3 brane_geometry.py sweep
python3 -m unittest -v
```

Primary sources:

- [STEREO hidden-neutron search](https://arxiv.org/abs/2111.01519)
- [2016 neutron-passing-through-walls experiment and two-brane interpretation](https://arxiv.org/abs/1604.07861)
- [Bulk-geometry-dependent neutron/hidden-neutron coupling](https://arxiv.org/abs/2009.12149)
