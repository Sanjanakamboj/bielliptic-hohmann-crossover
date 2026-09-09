# Bi-elliptic / Hohmann Crossover

**When does a three-impulse bi-elliptic transfer cost less total delta-v than a
two-impulse Hohmann transfer, for raising a spacecraft from one circular orbit to a
much higher circular one?**

**Status: Milestone 1 complete** — problem definition, analytical derivation, hand
calculations and verification plan. No solver, sweep, optimizer or plotting code is
implemented yet; that is Milestone 2 onward.

The full derivation, every hand calculation, the verification matrix and the
regression targets live in **[DESIGN.md](DESIGN.md)**.

---

## Objective

Determine the radius ratio at which bi-elliptic transfers become cheaper than
Hohmann transfers, express the result in a form independent of the central body and
the starting orbit, and turn it into an engineering recommendation rather than a
bare minimization.

Portfolio deliverable: a crossover radius-ratio plot plus an engineering
recommendation.

---

## Definitions

Radius is measured from the centre of the Earth; altitude from its surface. This
project is strict about the difference — **every ratio below is a ratio of radii.**

| Symbol | Definition | Constraint |
|---|---|---|
| `r1` | initial circular orbit **radius** | — |
| `r2` | final circular orbit **radius** | `r2 > r1` |
| `rb` | bi-elliptic intermediate apoapsis **radius** | `rb > r2` |
| `R` | `r2 / r1` — the final radius ratio | `R > 1` |
| `B` | `rb / r1` — the intermediate apoapsis ratio | `B > R` |
| `v1` | `sqrt(mu / r1)` — initial circular speed | — |
| `dv_bar` | `dv / v1` — normalized delta-v | — |

## Dimensionless formulation

Because `mu` and `r1` enter every speed only through `v1 = sqrt(mu/r1)`, normalized
delta-v depends on `R` and `B` alone. **The crossover is a pure radius-ratio result,
independent of the central body and of the starting orbit.**

Hohmann (two impulses, both prograde):

```
dv_bar_H(R) = sqrt(2R/(1+R)) - 1  +  1/sqrt(R) - sqrt(2/(R*(1+R)))
```

Bi-elliptic (three impulses; burns 1 and 2 prograde, **burn 3 retrograde**):

```
dv_bar_B(R,B) = sqrt(2B/(1+B)) - 1                        (burn 1, at r1)
              + sqrt(2R/(B*(R+B))) - sqrt(2/(B*(1+B)))    (burn 2, at rb)
              + sqrt(2B/(R*(R+B))) - 1/sqrt(R)            (burn 3, at r2)
```

Two exact limits anchor the whole analysis:

```
dv_bar_B(R, R)   = dv_bar_H(R)                      (B -> R+ IS the Hohmann transfer)
dv_bar_B_inf(R)  = (sqrt(2) - 1) * (1 + 1/sqrt(R))  (B -> infinity, derived in closed form)
```

## Hohmann vs. bi-elliptic — the physical picture

The Hohmann transfer raises apoapsis to exactly the target radius, then circularizes
there with a **prograde** burn. It is the minimum-delta-v *two-impulse* transfer.

The bi-elliptic transfer overshoots: it raises apoapsis far beyond the target, raises
periapsis to the target radius while moving slowly at that distant apoapsis (where
velocity changes are cheap), then arrives at the target radius moving *faster* than
circular and brakes with a **retrograde** burn.

The trade is that burn 2 is cheap precisely because the spacecraft is far away and
slow — but getting far away and slow costs extra on burn 1 and costs a great deal of
time. For small `R` the extra cost wins; for large `R` the cheap plane-of-attack at
apoapsis wins.

## There are TWO thresholds, not one

This is the point most treatments blur. The two thresholds answer **different
questions** and are computed from **opposite ends** of the `B` domain.

| | Threshold #1 | Threshold #2 |
|---|---|---|
| Value (recomputed here) | `R1* = 11.938765472645870716` | `R2* = 15.581718738763179213` |
| Condition | `dv_bar_H(R) = dv_bar_B_inf(R)` | `∂dv_bar_B/∂B = 0` at `B = R+` |
| Endpoint examined | `B -> infinity` | `B -> R+` |
| Exact form | `u³ − (1+2√2)u² + u + 1 = 0`, `u = √R` | `R³ − 15R² − 9R − 1 = 0` |
| Question | Can **any** bi-elliptic beat Hohmann? | Do **all** bi-elliptics beat Hohmann? |
| Quantifier | **∃** `B` | **∀** `B` |

**`R1*` is where the bi-elliptic family first contains a winner. `R2*` is where the
bi-elliptic family contains nothing but winners.**

Both values were derived and solved independently here — three or four separate
formulations each, including exact algebraic characterizations — and only then
compared against the classical textbook values (~11.94 and ~15.58) as a sanity
check. **Literature numbers are never used as production constants.**

### The three regions

| Region | Range | Behaviour |
|---|---|---|
| **A** | `1 < R < 11.9388` | Hohmann beats **every** admissible bi-elliptic transfer. |
| **B** | `11.9388 < R < 15.5817` | Bi-elliptic wins **only** for `B > B_crit(R)`. `B_crit -> ∞` as `R -> R1*`; `B_crit -> R` as `R -> R2*`. |
| **C** | `R > 15.5817` | Bi-elliptic wins for **every** `B > R`. |

A structural result established in M1: `dv_bar_B(R, ·)` has **no interior local
minimum** in `B`, so the best achievable bi-elliptic delta-v is always
`min(dv_bar_H, dv_bar_B_inf)` — an endpoint value, never a finite interior optimum.
See [DESIGN.md §6](DESIGN.md) for the proof sketch and the numerical scan, including
a correction to a plausible-but-empty alternative definition of Region B.

## Representative Earth example

`r1 = R_Earth + 300 km = 6678.1363 km`, `v1 = 7.7257606 km/s`. Final **altitude** is
a conversion from the universal **radius-ratio** result: `h2 = R*r1 - R_Earth`.

| `R` | final altitude `h2` [km] | `dv_H` [km/s] | best bi-elliptic (`B -> ∞`) | max saving | region |
|---|---|---|---|---|---|
| 5 | 27013 | 3.7084 | 4.6313 | −922.8 m/s | A |
| 10 | 60403 | 4.0930 | 4.2121 | −119.1 m/s | A |
| 12 | 73759 | 4.1269 | 4.1239 | **+3.0 m/s** | B |
| 16 | 100472 | 4.1429 | 4.0001 | +142.7 m/s | C |
| 20 | 127185 | 4.1312 | 3.9157 | +215.5 m/s | C |
| 50 | 327529 | 3.9687 | 3.6527 | +316.0 m/s | C |

Threshold radii for this case: `R1*` → `r2 = 79728.7 km` (altitude 73350.6 km);
`R2*` → `r2 = 104056.8 km` (altitude 97678.7 km).

**GEO from a 300 km parking orbit is `R = 6.31`** — deep in Region A. For the most
commonly flown high-orbit transfer, Hohmann simply wins.

## Delta-v is not the whole story

`B -> infinity` is a **mathematical asymptote with infinite transfer time**, never a
flight recommendation. Transfer duration scales as `B^(3/2)`:

| `R` | Hohmann time | `B_crit` | time at break-even | ratio |
|---|---|---|---|---|
| 12 | 0.52 d | 815.8 | **524 d** | **1006×** |
| 16 | 0.78 d | 16.0 | 2.79 d | 3.58× |
| 20 | 1.07 d | 20.0 | 3.88 d | 3.63× |
| 50 | 4.05 d | 50.0 | 15.16 d | 3.75× |

At `R = 12` the entire theoretical prize is 3 m/s, and merely breaking even costs a
1006× increase in transfer time. The final recommendation must therefore be an
**engineering trade**, not "choose the minimum delta-v".

## Roadmap

| Milestone | Content | Status |
|---|---|---|
| **M1** | Problem definition, analytical derivation, hand calculations, verification plan | **complete** |
| M2 | Production transfer solver (Hohmann + bi-elliptic), normalized and dimensional, with the DESIGN.md §9 verification suite | pending |
| M3 | Parameter sweeps in `R` and `B`; `B_crit(R)`; independent root-finding for both thresholds | pending |
| M4 | Crossover radius-ratio plot and supporting figures | pending |
| M5 | Transfer-time trade analysis and combined delta-v / time figures | pending |
| M6 | Engineering recommendation and portfolio write-up | pending |

## Scope and limitations

Coplanar circular initial and final orbits; ideal two-body central gravity;
impulsive burns; no finite-burn losses; no plane change; no J2; no third-body
perturbations; no atmospheric drag; no lunar or sphere-of-influence effects for
enormous apoapses; no launch constraints; no radiation model; no collision or
debris-environment model; no mission-specific operational optimization.
`B -> infinity` is mathematical only.

**No flight-operations or mission-optimality claims are made.** Full list in
[DESIGN.md §12](DESIGN.md).

## Repository layout

```
README.md                            this file
DESIGN.md                            M1 derivation, hand calculations, verification plan
pyproject.toml                       src-layout package, pytest configured for src/
.gitignore
src/bielliptic_crossover/__init__.py package scaffold (version only)
tests/test_placeholder.py            M1 scaffold tests
scripts/                             (empty) M2+ analysis drivers
figures/                             (empty) M4+ generated plots
results/                             (empty) M3+ numerical outputs
```

## Development

```bash
pip install -e ".[dev]"
pytest -W error
```

Version `0.1.0`.
