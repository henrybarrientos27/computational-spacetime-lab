# Engineering audit of a reversible neutron portal beamline

## Result

The ideal two-state transfer mechanism survives velocity spread and finite
field uniformity, but a straight horizontal ultracold-neutron beamline loses
the long coherent interaction time that made weak mixing look accessible.
Gravity is the dominant new constraint.

`beamline_realism.py` adds:

- a Gaussian neutron-velocity distribution;
- independent field nonuniformity in the two conversion zones;
- neutron beta-decay survival;
- a three-sigma gravitational aperture constraint;
- ensemble resonance width and scan sampling;
- stable discovery-time statistics even when signal is far below background.

This is still a conditional apparatus design. It does not assume that
neutron/hidden-neutron mixing exists.

## Baseline

The comparison holds these assumptions fixed:

- 5% RMS velocity spread;
- 0.01 microtesla RMS field variation in each conversion zone;
- 0.10 m vertical half-aperture;
- at most 50 m per conversion zone;
- 500,000 incident neutrons/s;
- 30% downstream detection efficiency;
- `10^-4` background counts/s;
- 50–1100 microtesla resonance search interval;
- two samples per ensemble FWHM and 2 s settling per field setting;
- discovery requires 5 sigma and at least ten regenerated neutrons.

The 10 nT field-uniformity assumption is already demanding. Fixed flux across
all velocities is an apples-to-apples sensitivity comparison, not a claim that
real facilities provide identical beams.

## Straight horizontal result

For the slowest three-sigma neutron, the zone length must obey

```text
L <= v_slow sqrt(2 a / g),
```

where `a` is the vertical half-aperture. At a 5 m/s mean speed and 5% spread,
the limit is only `0.607 m`, corresponding to `0.121 s` at the mean velocity.

| Mean speed | Mixing time tau | Optimized zone | Peak wall probability | Settings | Full scan | Known resonance |
|---:|---:|---:|---:|---:|---:|---:|
| 5 m/s | 1 s | 0.52 m | `1.18×10^-4` | 5,259 | 3.9 h | 0.56 s |
| 5 m/s | 10 s | 0.607 m | `2.22×10^-8` | 6,181 | 259 d | 50 min |
| 5 m/s | 100 s | 0.607 m | `2.22×10^-12` | 6,181 | 6.4 million y | 718 y |
| 500 m/s | 1 s | 50 m | `1.02×10^-4` | 5,052 | 3.9 h | 0.66 s |
| 500 m/s | 10 s | 50 m | `1.02×10^-8` | 5,050 | 458 d | 1.81 h |
| 500 m/s | 100 s | 50 m | `1.02×10^-12` | 5,050 | 24.5 million y | 3,370 y |
| 1000 m/s | 10 s | 50 m | `6.39×10^-10` | 2,477 | 37.5 y | 4.04 d |

The table reveals a time-bandwidth trap:

- slow neutrons provide interaction time but fall out of a narrow straight
  apparatus;
- fast neutrons stay in the apparatus but cross it too quickly;
- longer interaction narrows the magnetic resonance and increases scan points;
- after the zone reaches its length/aperture limit, two-zone regeneration falls
  approximately as `tau^-4`.

The ensemble model makes the earlier ideal table more honest. For `tau=10 s`,
the ideal 10 m UCN path predicted a 1.16-day scan; a straight buildable aperture
changes that to 259 days. For `tau=100 s`, the change is from about 467 days to
millions of years.

## Can gravity be engineered around?

Gravity is not a no-go theorem, but every workaround has a measurable cost.

- A freely falling/parabolic vacuum and magnet path can preserve coherence.
  A 5 m/s neutron observed for 2 s drops about 19.6 m, so the earlier 10 m
  horizontal-span concept becomes a roughly building-scale vertical apparatus,
  not a tabletop beamline.
- A vertical UCN fountain reduces floor space, but conversion and regeneration
  must occur on opposite sides of an absorber while the hidden state's own
  gravitational acceleration is unknown. That unknown is part of the physics
  being tested.
- Material neutron guides restore the path but wall collisions project the
  state and reset coherent oscillation.
- Magnetic support requires field gradients. For a neutron, balancing gravity
  is of order tesla per metre; that is incompatible with a resonance whose
  useful width can be tens of nanotesla unless a substantially more elaborate
  spin/detuning compensation scheme is demonstrated.
- A space or free-fall experiment removes sag but greatly increases experiment
  complexity before any nonzero portal coupling has been found.

The rational sequence is therefore staged: first search for nonzero mixing in
an existing high-flux, collision-based or short-coherence experiment; only
after a resonance is located should a long ballistic two-zone regenerator be
built. Knowing the resonance deletes thousands of scan settings and changes
the engineering decision completely.

## Strongest presently testable point

Under the explicit baseline, `tau` near 1 s remains accessible in hours and
could produce the full disappearance/regeneration fingerprint. Existing UCN
experiments already report no significant signal and a conservative
`tau > 1 s` limit over part of the scanned detuning range. The beamline is
therefore a way to test remaining parameter space and systematics, not a
prediction that a portal lies at `tau = 1 s`.

For weaker mixing, the strongest computer-derived design rule is:

> Locate the resonance with a broad, high-flux discovery stage before spending
> coherence time on a narrow, reversible two-zone transfer stage.

Trying to do discovery and near-complete coherent transfer with the same
apparatus is what creates the prohibitive scan.

## Rate targets turn elapsed time into an engineering specification

The long baseline times are not fundamental if incident rate can be increased
without increasing downstream background. Inverting the same 5-sigma counting
model for a one-year campaign gives:

| Mixing time tau | Targeted test if resonance is known | Blind 50–1100 microtesla scan |
|---:|---:|---:|
| 10 s | `1.36×10^3` neutrons/s | `3.54×10^5` neutrons/s |
| 30 s | `1.10×10^5` neutrons/s | `2.87×10^7` neutrons/s |
| 100 s | `1.36×10^7` neutrons/s | `3.54×10^9` neutrons/s |

These rates apply to the fixed 5 m/s, gravity-limited geometry and the stated
efficiency/background assumptions. The `tau=100 s` blind target replaces a
6.4-million-year low-rate scan with one calendar year, but it creates a severe
shielding requirement: keeping direct ordinary-neutron leakage below
`10^-4/s` at `3.54×10^9/s` means a raw beam rejection below about
`2.8×10^-14`, before accounting for capture gammas and secondary particles.
This rejection, not computational uncertainty, becomes the apparatus gate.

## Run it

```bash
cd /Users/henrybarrientos/Documents/Codex/2026-07-15/i-w/outputs/spacetime-lab
python3 beamline_realism.py optimize --tau 10 --velocity 5
python3 beamline_realism.py target --tau 100 --target-days 365.25
python3 beamline_realism.py sweep
python3 -m unittest -v
```

The full result table is `realistic-beamline-scan.csv`.
