# Bi-elliptic / Hohmann Crossover

**When does a three-impulse bi-elliptic transfer cost less total delta-v than a
two-impulse Hohmann transfer, for raising a spacecraft from one circular orbit to a
much higher circular one?**

**Status: Milestone 3 complete** — the verified solver is now a full trade study:
a global crossover map, the break-even apoapsis locus, the delta-v/time trade, and
a scoped engineering recommendation, backed by **1026 tests** passing under
`pytest -W error`.

![Crossover map](figures/m3_crossover_map.png)

The derivation, verification matrix and implementation notes live in
**[DESIGN.md](DESIGN.md)**; reproducible numbers in
**[results/](results/)**.

---

## Objective

Determine the radius ratio at which bi-elliptic transfers become cheaper than
Hohmann transfers, express the result in a form independent of the central body and
the starting orbit, and turn it into an engineering recommendation rather than a
bare minimization.

Portfolio deliverable: a crossover radius-ratio plot plus an engineering
recommendation.

```python
>>> from bielliptic_crossover import threshold_R1, threshold_R2, break_even_B
>>> round(threshold_R1(), 10), round(threshold_R2(), 10)
(11.9387654726, 15.5817187388)
>>> break_even_B(12.0).B_crit
815.8202504753092
>>> from bielliptic_crossover.trade import recommendation_for
>>> recommendation_for(6.3137).headline          # 300 km LEO to GEO
'Use the Hohmann transfer.'
```

## The two transfers

A **Hohmann** transfer raises apoapsis to exactly the target radius, then
circularizes there with a prograde burn. It is the minimum-delta-v *two-impulse*
transfer between coplanar circular orbits.

A **bi-elliptic** transfer deliberately overshoots: it raises apoapsis far beyond
the target, raises periapsis to the target radius while moving slowly at that
distant apoapsis (where velocity changes are cheap), then arrives at the target
moving *faster* than circular and brakes with a retrograde burn.

The trade is that the middle burn is cheap precisely because the spacecraft is far
away and slow — but getting far away and slow costs extra on the first burn, and
costs a great deal of time.

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
| Value (recomputed here) | `R1* = 11.9387654726458923` | `R2* = 15.5817187387631790` |
| Independent second path | `11.9387654726458692` (cubic) | `15.5817187387631879` (cubic) |
| Condition | `dv_bar_H(R) = dv_bar_B_inf(R)` | `∂dv_bar_B/∂B = 0` at `B = R+` |
| Endpoint examined | `B -> infinity` | `B -> R+` |
| Exact form | `u³ − (1+2√2)u² + u + 1 = 0`, `u = √R` | `R³ − 15R² − 9R − 1 = 0` |
| Question | Can **any** bi-elliptic beat Hohmann? | Do **all** bi-elliptics beat Hohmann? |
| Quantifier | **∃** `B` | **∀** `B` |

**`R1*` is where the bi-elliptic family first contains a winner. `R2*` is where the
bi-elliptic family contains nothing but winners.**

Both values are **recomputed at run time by two independent solver paths each** —
root-finding on the transfer equations, and the exact polynomial characterization —
which agree to `2.3e-14` and `8.9e-15` respectively. They were compared against the
classical textbook values (~11.94 and ~15.58) only after independent derivation.
**No literature number appears as a production constant**; a test asserts that
`constants.py` contains no crossover value at all.

### The three regions

| Region | Range | Behaviour |
|---|---|---|
| **A** | `1 < R < 11.9388` | Hohmann beats **every** admissible bi-elliptic transfer. |
| **B** | `11.9388 < R < 15.5817` | Bi-elliptic wins **only** for `B > B_crit(R)`. `B_crit -> ∞` as `R -> R1*`; `B_crit -> R` as `R -> R2*`. |
| **C** | `R > 15.5817` | Bi-elliptic wins for **every** `B > R`. |

Three statements this project does **not** make: that bi-elliptic is always better
above 11.94; that 15.58 is where bi-elliptic first wins; that `B -> ∞` is optimal.

A structural result established in M1 and **re-verified numerically in M2**:
`dv_bar_B(R, ·)` has **no interior local minimum** in `B`, so the best achievable
bi-elliptic delta-v is always `min(dv_bar_H, dv_bar_B_inf)` — an endpoint value,
never a finite interior optimum. Dense sweeps over `R` from 1.5 to 200 found
**zero** interior minima; every interior turning point present is a *maximum*.
An optimizer reporting a finite interior optimum is a bug, not a discovery.

### Break-even apoapsis ratio in Region B

In Region B the bi-elliptic transfer wins only once the intermediate apoapsis is
large enough. `break_even_B(R)` solves for that threshold, explicitly excluding the
trivial root at `B = R` (where equality holds by construction):

| `R` | `B_crit` | `rb_crit / r2` | `rb_crit` altitude (Earth ref) | `t_B/t_H` at break-even |
|---|---|---|---|---|
| 12 | 815.8202504753 | 68.0 | 5.44e+06 km | **1006** |
| 12.5 | 90.7509442091 | 7.26 | 6.00e+05 km | 38.9 |
| 13 | 48.9048433284 | 3.76 | 3.20e+05 km | 16.0 |
| 14 | 26.1046112824 | 1.86 | 1.68e+05 km | 8.0 |
| 15 | 18.1902815122 | 1.21 | 1.15e+05 km | 4.3 |
| 15.5 | 15.8968710115 | 1.03 | 9.98e+04 km | 3.6 |

Measured asymptotic laws across Region B:
`B_crit/R ∝ (R − R1*)^−1` and `t_B/t_H ∝ (R − R1*)^−3/2`.

Region A returns "no bi-elliptic win"; Region C returns no finite root either,
because every `B > R` already wins and the boundary `B = R` is open.

## Representative Earth examples

`r1 = R_Earth + 300 km = 6678.1363 km`, `v1 = 7.7257606 km/s`. Final **altitude** is
a conversion from the universal **radius-ratio** result: `h2 = R*r1 - R_Earth`.
These large target radii are dimensional illustrations of a normalized result, not
proposed mission designs.

| case | `R` | region | `h2` [km] | `dv_H` [km/s] | `B→∞` infimum | best possible saving | `B=5R` saving |
|---|---|---|---|---|---|---|---|
| **GEO** | 6.3137 | A | 35786 | 3.8926 | 4.4737 | **0 m/s** | **−486 m/s** |
| | 10 | A | 60403 | 4.0930 | 4.2121 | 0 m/s | −125 m/s |
| | 12 | B | 73759 | 4.1269 | 4.1239 | +3.0 m/s | **−27 m/s** |
| | 13 | B | 80438 | 4.1355 | 4.0877 | +47.9 m/s | +8.5 m/s |
| | 15 | B | 93794 | 4.1427 | 4.0264 | +116.3 m/s | +63.9 m/s |
| | 16 | C | 100472 | 4.1429 | 4.0001 | +142.7 m/s | +85.4 m/s |
| | 20 | C | 127185 | 4.1312 | 3.9157 | +215.5 m/s | +145.3 m/s |
| | 50 | C | 327529 | 3.9687 | 3.6527 | +316.0 m/s | +233.9 m/s |

Two rows carry the message. **GEO** sits firmly in Region A — the best possible
saving is exactly zero, and a representative bi-elliptic *costs 486 m/s more*.
**`R = 12`** is in Region B, yet `B = 5R` still costs 27 m/s more than Hohmann,
because `5R = 60` is far below `B_crit = 816`.

## Why delta-v savings alone are insufficient

`B -> infinity` is a **mathematical asymptote with infinite transfer time**, never a
flight recommendation. Transfer duration scales as `B^(3/2)`:

| `R` | Hohmann time | `B_crit` | time at break-even | ratio |
|---|---|---|---|---|
| 12 | 0.52 d | 815.8 | **524 d** | **1006×** |
| 16 | 0.78 d | 16.0 | 2.79 d | 3.58× |
| 20 | 1.07 d | 20.0 | 3.88 d | 3.63× |
| 50 | 4.05 d | 50.0 | 15.16 d | 3.75× |

![Delta-v vs. time trade](figures/m3_dv_time_trade.png)

At `R = 12` (Region B) the numbers are stark, and come straight from production code:

| strategy | `B` | saving | `t_B/t_H` | `r_b` / lunar distance |
|---|---|---|---|---|
| `B/R = 2` | 24 | **−39.0 m/s** | 7.3 | 0.42 |
| `B_crit` | 815.8 | 0 m/s | **1006** | 14.2 |
| `2 B_crit` | 1631.6 | +1.5 m/s | **2829** | 28.3 |
| `B → ∞` | — | +3.04 m/s (infimum) | ∞ | ∞ |

Even at twice the break-even apoapsis, the saving is half the theoretical prize for
a 2829× time penalty at 28 lunar distances. A 3 m/s saving sits inside the noise of
launch dispersion, navigation and finite-burn losses.

At `R = 16` (Region C) the benefit appears far more readily — `B/R = 2` already
saves 31.8 m/s — but still costs 6.9× the transfer time, and every design saving
more than ~85 m/s is already beyond the Moon.

Just above `R2*`, a `B = 1.001R` transfer saves **0.0005 m/s** while taking 3.58×
as long: *"every bi-elliptic beats Hohmann"* is a statement about sign, not
magnitude.

## Model-validity warning

Whenever the intermediate apoapsis `r_b` approaches or exceeds the **lunar distance
(384 400 km)**, an Earth-only two-body model is no longer physically appropriate;
beyond the Earth Hill radius (~1.5e6 km) the spacecraft is not meaningfully
Earth-bound at all. Every result table and both trade figures carry this flag.

Solving `r_b,crit` = lunar distance gives **`R = 12.8349`**: below that — the lower
quarter of Region B — **the break-even apoapsis lies beyond the Moon**. Numbers
quoted there are mathematical extrapolations of an idealized model, not trajectory
designs.

## Recommendation

- **`R < 11.9388`** — Hohmann is delta-v superior to every admissible bi-elliptic
  transfer. Its rival's longer transfer time buys nothing.
- **`11.9388 < R < 15.5817`** — bi-elliptic beats Hohmann *only* if `r_b` exceeds
  `B_crit(R)·r_1`. Near the lower threshold the required apoapsis and duration are
  absurd for a vanishing saving.
- **`R > 15.5817`** — every `B > R` is delta-v better, but the practical size of the
  saving must still be weighed against the transfer-time penalty.
- **`B → ∞` is a mathematical lower bound**, never an engineering recommendation and
  never a realizable transfer.
- **For a 300 km LEO to GEO transfer (`R = 6.3137`), Hohmann remains firmly
  preferred** within this idealized coplanar impulsive two-body model.

![Break-even penalty](figures/m3_break_even_time_penalty.png)

### A subtlety the code makes explicit

Delta-v and transfer time degenerate **differently** at `B = R`. The delta-v
reduces exactly to Hohmann, but the time does not:
`t_B(R,R) = t_H(R) + π·R^{3/2}·t_star`, i.e. the Hohmann time *plus half the period
of the final circular orbit*, because ellipse 2 becomes that circular orbit and the
path still coasts half a revolution of it. See
`bielliptic_time_excess_at_B_equals_R`.

## Roadmap

| Milestone | Content | Status |
|---|---|---|
| M1 | Problem definition, analytical derivation, hand calculations, verification plan | complete |
| M2 | Production transfer equations, dimensional cross-checks, crossover and break-even solvers, verification suite | complete |
| **M3** | Crossover map, break-even locus, delta-v/time trade, portfolio figures, scoped engineering recommendation | **complete** |
| M4 | Final repository packaging and polish (licence, CI, release) | pending |

The crossover plot, transfer-time trade and engineering recommendation originally
sketched for M4–M6 were delivered in M3; M4 is now final packaging and polish only.

## Scope and limitations

Coplanar circular initial and final orbits; ideal two-body central gravity;
impulsive burns; no finite-burn losses; no plane change; no J2; no third-body
perturbations; no atmospheric drag; no lunar or sphere-of-influence effects for
enormous apoapses; no launch constraints; no radiation model; no collision or
debris-environment model; no mission-specific operational optimization.
`B -> infinity` is mathematical only.

**No flight-operations or mission-optimality claims are made.** Full list in
[DESIGN.md §12](DESIGN.md).

## Verification and result artifacts

**1026 tests** pass under `pytest -W error` with no warnings. Every result file is
regenerated from production code and is byte-identical across runs (checked by
tests). Thresholds are recomputed by two independent solver paths each and are
never hardcoded — a test asserts `constants.py` contains no crossover value at all.

M3 (trade study):
- **[results/m3_summary.json](results/m3_summary.json)** — thresholds, Earth
  reference, `B_crit` examples, GEO, R=12/R=16 trades, recommendation
- **[results/m3_radius_trade.csv](results/m3_radius_trade.csv)** — envelope over the
  1090-point authoritative grid
- **[results/m3_break_even_curve.csv](results/m3_break_even_curve.csv)** — 219 solved
  `B_crit(R)` roots
- **[results/m3_finite_b_strategies.csv](results/m3_finite_b_strategies.csv)**,
  **[results/m3_earth_examples.csv](results/m3_earth_examples.csv)**

M2 (equation verification): **[results/m2_verification_report.txt](results/m2_verification_report.txt)**,
plus the diagnostic figures
[m2_fig1](figures/m2_fig1_excess_vs_B.png) and [m2_fig2](figures/m2_fig2_time_trade_R12.png).

An independent cross-check recomputes selected points with a scratch vis-viva
implementation that does not import the trade layer; worst relative residual
`5.4e-13`.

## Repository layout

```
README.md                             this file
DESIGN.md                             M1 derivation + M2 implementation notes
pyproject.toml                        src-layout package, pytest configured for src/
src/bielliptic_crossover/
    __init__.py                       public API, version
    constants.py                      physical constants (NO crossover constants)
    _stable.py                        cancellation-free sqrt(1+x)-1
    hohmann.py                        two-impulse transfer
    bielliptic.py                     three-impulse transfer + infinite-B limit
    timing.py                         transfer times, B=R degeneracy
    dimensional.py                    independent direct vis-viva path
    crossover.py                      thresholds, break-even B, classifier, structure scan
    trade.py                          M3 trade study: map, locus, metrics, recommendation
tests/                                1026 tests across 11 files
scripts/m2_verification_report.py     regenerates the M2 numerical report
scripts/m2_diagnostic_figures.py      regenerates the M2 diagnostic figures
scripts/m3_trade_analysis.py          regenerates the M3 CSV/JSON artifacts
scripts/m3_figures.py                 regenerates the three M3 portfolio figures
figures/                              M2 diagnostic + M3 portfolio figures
results/                              M2 report + M3 trade artifacts
```

## Development

```bash
pip install -e ".[dev,figures]"
pytest -W error
python scripts/m2_verification_report.py
python scripts/m2_diagnostic_figures.py
python scripts/m3_trade_analysis.py
python scripts/m3_figures.py
```

Version `0.3.0`. Runtime dependencies: `numpy`, `scipy` (`matplotlib` for figures only).
