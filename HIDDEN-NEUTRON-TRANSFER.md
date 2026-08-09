# Reversible hidden-neutron transfer: the closest real portal experiment

## Result

There is a mathematically valid, experimentally testable mechanism that is
closer to “send something through another dimension and get it back” than
missing momentum at a collider: neutron–hidden-neutron oscillation.

In a two-brane interpretation, an ordinary neutron state `n` is localized on
our brane and a hidden state `n'` is localized on another brane. A weak mixing
term can convert one into the other. A two-converter apparatus could therefore
perform this sequence:

```text
visible neutron
    -> controlled conversion region A
hidden neutron
    -> ordinary neutron absorber / wall
hidden neutron
    -> controlled conversion region B
visible neutron in a shielded detector
```

That is reversible single-particle transfer in the model. It is not known to
occur in nature. A hidden neutron could also be an internal sterile or mirror
state rather than a particle moving along a literal spatial coordinate, so a
positive experiment would establish a portal first and extra-dimensional
geometry only after additional tests.

## The controlled transfer equation

The optimistic spin-independent two-state Hamiltonian is

```text
H = [[E_n, epsilon],
     [epsilon, E_n']]
```

For mixing energy `epsilon`, detuning `Delta = E_n - E_n'`, and coherent time
`t`, the conversion probability is

```text
P(n -> n') = 4 epsilon^2 / (Delta^2 + 4 epsilon^2)
              * sin^2(sqrt(Delta^2 + 4 epsilon^2) t / (2 hbar)).
```

Writing `epsilon = hbar/tau`, exact resonance gives
`P = sin²(t/tau)`. A complete swap is mathematically possible at
`t = pi*tau/2`. Two equal resonant regions then give a through-wall probability
`P²`, which can also reach one in the ideal model.

Matter changes the problem. The collision-averaged probability per quantum
projection used by the STEREO analysis is

```text
p = 2 epsilon² /
    ((Delta E + V_F)² + 4 epsilon² + hbar² Gamma²/4),
```

where `V_F` is the material Fermi potential and `Gamma` is the collision rate.
Collisions can create a large incoherent source, but they also reset coherent
amplitude. This is why a reactor is excellent for disappearance/regeneration
statistics while a controlled full swap needs a long collision-free flight.

## What experiments currently say

- The STEREO reactor experiment found no hidden-neutron signal and constrained
  the per-projection swapping probability to
  `p < 3.1×10^-11` at 95% confidence. Its analysis says a background-free
  version of the same exposure could reach `2×10^-12`; the present result was
  limited primarily by reactor-induced and cosmic backgrounds. See the
  [STEREO paper](https://arxiv.org/abs/2111.01519).
- A 2023 ultracold-neutron beam experiment scanned magnetic fields from 50 to
  1100 µT and obtained the conservative limit `tau > 1 s` for mass splittings
  from 2 to 69 peV. It found no significant signal. The experiment used a 5 m
  solenoid, approximately 500,000 UCN counts/s, an average 32.2 ms free-flight
  time, and about 26 wall collisions per detected neutron. See the
  [PRL experiment and open data record](https://repo.scoap3.org/records/81314/).
- The newest dedicated PSI result, accepted in June 2026, scanned 5–109 µT and
  found no anomalous losses. It excludes 99.98% of the solid-angle parameter
  space associated with previously claimed mirror-neutron anomalies. See the
  [2026 PSI paper](https://arxiv.org/abs/2602.23487).

No neutron has been demonstrated to leave our sector or brane and return.

## The computer-designed two-zone scan

`hidden_neutron_lab.py` calculates coherent transfer, matter-induced swapping,
the magnetic resonance width, a two-zone through-wall rate, and the exposure
needed for a counting discovery. It optimizes converter interaction time while
accounting for scan bandwidth and field-settling overhead.

The baseline is deliberately explicit:

- UCN velocity: 5 m/s;
- incident rate: 500,000/s;
- downstream efficiency: 30%;
- downstream background: `10^-4/s`;
- magnetic search interval: 50–1100 µT;
- two samples across each resonance FWHM;
- 2 s field-settling overhead per setting;
- two equal conversion zones with at most 10 m equivalent coherent path each;
- discovery requirement: 5σ and at least ten regenerated neutrons.

| Mixing time tau | Optimized equivalent path per zone | FWHM | Settings | Full scan time | Time if resonance were known |
|---:|---:|---:|---:|---:|---:|
| 1 s | 0.55 m | 0.557 µT | 3,771 | 0.12 day | 0.48 s |
| 3 s | 1.63 m | 0.186 µT | 11,290 | 0.35 day | 0.48 s |
| 10 s | 5.44 m | 0.0558 µT | 37,607 | 1.16 days | 0.48 s |
| 30 s | 10 m | 0.0304 µT | 69,146 | 5.38 days | 3.4 s |
| 100 s | 10 m | 0.0304 µT | 69,136 | 467 days | 7.0 min |
| 300 s | 10 m | 0.0304 µT | 69,135 | 191 years | 14.0 h |
| 414 s | 10 m | 0.0304 µT | 69,135 | 1,884 years | 5.44 days |

These times are optimistic lower bounds, not a construction estimate. A 5 m/s
neutron falls about 19.6 m during a two-second horizontal free flight. A real
10 m *coherent* path therefore needs a material-free magnetic/gravitational
trajectory or fountain geometry. Material guide collisions cannot simply be
ignored: the experiments treat them as state projections that reset the
oscillation.

The `tau=414 s` row is an illustrative zero-field mixing bound, not a uniform
limit over the entire 50–1100 µT interval. The table exposes the scaling rather
than claiming every listed point remains allowed.

## What the optimization discovered

The weak coupling is only half the problem. The unknown detuning is equally
important.

If the resonance field were known, one setting could accumulate a test signal
in minutes or days even for weak mixing. If it is unknown, longer coherent
interaction increases conversion but narrows the resonance. The apparatus must
then scan tens of thousands of field values. Once the converter length reaches
its practical maximum, the two-zone probability falls approximately as
`tau^-4`; search time becomes enormous.

This produces a precise present blocker:

> No nonzero neutron–hidden-neutron mixing has been observed, and the unknown
> energy detuning prevents long-coherence resonant transfer from being targeted
> without an increasingly expensive high-resolution scan.

That statement can be changed only by experimental data, not additional
algebra or computer time.

## Strongest next experiment

The next apparatus should not merely count unexplained downstream neutrons. It
should demand the complete resonance fingerprint:

1. A tunable upstream conversion region and independently tunable downstream
   regeneration region.
2. An absorber thick enough to reduce the ordinary neutron beam below the
   downstream background while leaving the hypothetical sterile state
   unaffected.
3. Simultaneous upstream disappearance and downstream regeneration at the same
   resonance energy.
4. Field reversal, source-off, absorber-removed, and deliberately detuned
   control runs.
5. Signal scaling with both conversion times as predicted by the two-state
   Hamiltonian.
6. Repetition with a different absorber and detector technology.
7. A material-free flight strategy—magnetogravitational trap, vertical UCN
   fountain, or a faster neutron beam with a longer magnetic region—to prevent
   wall projections from erasing coherence.

The most valuable immediate improvement to the reactor method is less exotic:
reduce reactor/cosmic backgrounds. The published STEREO estimate says that
doing so could improve its probability sensitivity from `3.1×10^-11` to about
`2×10^-12` without increasing the existing exposure.

## What a positive result would and would not prove

A valid discovery would establish reversible coupling between a neutron and an
invisible propagating state. To establish a *spatial* extra dimension, the next
tests would still need to recover geometric information—such as a KK spectrum,
distance-dependent propagation law, or a common compactification scale also
seen in collider or short-range-gravity data.

It would be the first genuine gate toward transport, not yet a transporter.
Atoms contain protons and electrons as well as neutrons; their charges and
binding fields have no demonstrated hidden-state partners. Scaling from a free
neutron to an intact atom therefore remains a separate unsolved coupling
problem.

## Run it

```bash
cd /Users/henrybarrientos/Documents/Codex/2026-07-15/i-w/outputs/spacetime-lab
python3 hidden_neutron_lab.py coherent --tau 10 --time 1
python3 hidden_neutron_lab.py medium --epsilon-ev 1e-6 --delta-ev 2 --collision-rate 1e5
python3 hidden_neutron_lab.py optimize --tau 100
python3 hidden_neutron_lab.py sweep
python3 -m unittest -v
```

The complete baseline scan is in `hidden-neutron-beamline-scan.csv`.
