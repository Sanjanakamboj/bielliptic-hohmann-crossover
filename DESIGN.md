# DESIGN — status

**Technical project complete (v1.0.0).**

| Milestone | Role | Authority |
|---|---|---|
| M1 | Analytical derivation, hand calculations, verification plan | derivation history |
| M2 | Production transfer equations, solvers, verification suite | verification history |
| **M3** | Crossover map, break-even locus, Δv/time trade, recommendation | **authoritative for all final results** |
| M4 | Audit, packaging, reproducibility, CI, release | audit record only — no physics changes |

**M3 is authoritative** for the crossover results, the break-even locus and the
engineering recommendation. M1 and M2 are preserved unchanged as derivation and
verification history; where an earlier section's prose disagrees with M3 output,
M3 wins. M4 changed **no physics** — it corrected documentation transcription
errors, added packaging, and recorded the audit (see the M4 section at the end).

The single source of truth for headline numbers is
`results/m4_final_summary.md`, generated from production code.

---

# DESIGN — Milestone 1

**Bi-elliptic / Hohmann crossover: problem definition, analytical derivation, hand
calculations, and verification plan.**

Status: **M1 complete.** This milestone is design and derivation only. No transfer
solver, sweep, optimizer, root-finding module or plotting code is implemented.
Every number below was derived and checked here with a private scratch calculation
that is **not** part of this repository; Milestone 2 will implement these results as
tested production code and must reproduce the values recorded in §10.

---

## 1. Problem definition and normalization

### 1.1 Model

An ideal, Earth-centred, two-body, coplanar, impulsive-burn transfer between two
circular orbits. A spacecraft starts on a circular orbit of radius `r1` and must
reach a circular orbit of larger radius `r2`. We compare two manoeuvre families:

- the two-impulse **Hohmann** transfer, and
- the three-impulse **bi-elliptic** transfer via an intermediate apoapsis radius `rb`.

### 1.2 Radius vs. altitude — terminology

This document is strict about the distinction, because the entire result is a
*radius* ratio and mixing in altitude silently changes the answer:

- **radius** `r` is measured from the centre of the Earth;
- **altitude** `h` is measured from the Earth's mean equatorial radius, `h = r - R_Earth`.

All of `r1`, `r2`, `rb`, `R` and `B` are **radii or ratios of radii**, never altitudes.
Altitudes appear only in §7, purely to give a dimensional reader some intuition, and
are always labelled as a conversion.

### 1.3 Nondimensional variables

Final-to-initial **radius** ratio:

```
R = r2 / r1,        R > 1
```

Intermediate-apoapsis **radius** ratio:

```
B = rb / r1,        B > R
```

The constraint `B > R` is what makes the manoeuvre a genuine bi-elliptic transfer:
the intermediate apoapsis must lie beyond the target radius, so that the second
ellipse descends to `r2` at its periapsis. `B = R` is the degenerate limit and,
as shown in §4.1, reduces *exactly* to the Hohmann transfer.

Initial circular speed and normalized delta-v:

```
v1 = sqrt(mu / r1)
dv_bar = dv / v1
```

Because `mu` and `r1` enter every speed only through the common factor `sqrt(mu/r1)`,
`dv_bar` depends on `R` and `B` alone. **The crossover is therefore a pure radius-ratio
result: it is independent of `mu`, of `r1`, and of which central body is used.**
This is verified numerically in §9 item C.

### 1.4 Dimensional reference case (intuition only)

```
mu_Earth = 398600.4418 km^3/s^2
R_Earth  =   6378.1363 km
h1       =    300      km          (initial ALTITUDE)
r1       = R_Earth + h1 = 6678.1363 km      (initial RADIUS)
v1       = sqrt(mu/r1)  = 7.72576063698292 km/s
```

This case is used only to convert normalized results into m/s and days. It never
enters the derivation of the thresholds.

---

## 2. Hohmann transfer — derivation from vis-viva

The vis-viva equation for a Keplerian orbit of semi-major axis `a`:

```
v^2 = mu * (2/r - 1/a)
```

The Hohmann transfer ellipse has periapsis `r1` and apoapsis `r2`, hence

```
a_H = (r1 + r2) / 2
```

### 2.1 First impulse (at r1, periapsis of the transfer ellipse)

```
v_p^2 = mu * (2/r1 - 1/a_H)
      = mu * (2/r1 - 2/(r1 + r2))
      = (mu/r1) * (2 - 2*r1/(r1 + r2))
      = (mu/r1) * 2*r2/(r1 + r2)
```

Dividing by `v1^2 = mu/r1` and substituting `r2 = R*r1`:

```
v_p_bar = sqrt(2R / (1 + R))
```

The spacecraft is on a circular orbit at `v_bar = 1` and must speed up to `v_p_bar`,
so the burn is **prograde** with magnitude

```
dv_bar_H1 = sqrt(2R/(1+R)) - 1                                             (H1)
```

### 2.2 Second impulse (at r2, apoapsis of the transfer ellipse)

```
v_a^2 = mu * (2/r2 - 2/(r1 + r2))
      = mu * 2*[(r1 + r2) - r2] / (r2*(r1 + r2))
      = 2*mu*r1 / (r2*(r1 + r2))
```

Dividing by `mu/r1` and substituting `r2 = R*r1`, `r1 + r2 = r1*(1+R)`:

```
v_a_bar^2 = 2*r1^2 / (r2*(r1+r2)) = 2 / (R*(1+R))
v_a_bar   = sqrt(2/(R*(1+R)))
```

The required circular speed at `r2` is `v_bar = 1/sqrt(R)`. Since `R > 1` implies
`2/(1+R) < 1`, we have `v_a_bar < 1/sqrt(R)`: the spacecraft arrives **slower** than
circular, so this burn is also **prograde**:

```
dv_bar_H2 = 1/sqrt(R) - sqrt(2/(R*(1+R)))                                  (H2)
```

### 2.3 Total

```
dv_bar_H(R) = sqrt(2R/(1+R)) - 1 + 1/sqrt(R) - sqrt(2/(R*(1+R)))           (H)
```

A more compact equivalent form, obtained by collecting the two `sqrt(2/(1+R))` terms,
is used later in §5.1:

```
dv_bar_H(R) = sqrt(2/(1+R)) * (R-1)/sqrt(R) + 1/sqrt(R) - 1                (H')
```

### 2.4 Required checks

- **`R -> 1` gives `dv_bar_H -> 0`.** At `R = 1`, (H) evaluates term by term to
  `sqrt(2/2) - 1 + 1 - sqrt(2/2) = 1 - 1 + 1 - 1 = 0`. Numerically,
  `dv_bar_H(1) = 0` exactly and `dv_bar_H(1 + 1e-9) = 5.0e-10`, i.e. the leading
  behaviour is `dv_bar_H ≈ (R-1)/2` — smooth and vanishing, as required.
- **Both burns positive for `R > 1`.** (H1) is positive because `2R/(1+R) > 1`
  whenever `R > 1`. (H2) is positive because `2/(1+R) < 1` whenever `R > 1`.
- **Dimensional / nondimensional agreement.** Computing the two impulses directly
  from dimensional vis-viva with the Earth reference case and dividing by `v1`
  reproduces (H) to `3.4e-41` at 40-digit precision across
  `R ∈ {1.5, 2, 5, 10, 11.9, 12, 15, 15.6, 16, 20, 50, 200}`.

---

## 3. Bi-elliptic transfer — derivation from vis-viva

Three impulses, two transfer ellipses:

```
ellipse 1:  periapsis r1, apoapsis rb    ->  a_1 = (r1 + rb)/2
ellipse 2:  periapsis r2, apoapsis rb    ->  a_2 = (r2 + rb)/2
```

### 3.1 Burn 1 — at `r1`, circular -> ellipse 1

Identical in form to §2.1 with `r2` replaced by `rb`, hence with `R` replaced by `B`:

```
v_bar(r1, ellipse 1) = sqrt(2B/(1+B))
dv_bar_B1 = sqrt(2B/(1+B)) - 1                                             (B1)
```

**Prograde.** Positive for all `B > 1`.

### 3.2 Burn 2 — at `rb`, ellipse 1 -> ellipse 2

Speed at the apoapsis `rb` of ellipse 1 (same algebra as §2.2 with `r2 -> rb`):

```
v_bar(rb, ellipse 1) = sqrt(2/(B*(1+B)))
```

Speed at the apoapsis `rb` of ellipse 2:

```
v^2 = mu*(2/rb - 2/(r2 + rb)) = 2*mu*r2 / (rb*(r2 + rb))
```

Dividing by `mu/r1` with `r2 = R*r1`, `rb = B*r1`:

```
v_bar(rb, ellipse 2)^2 = 2*r2*r1 / (rb*(r2+rb)) = 2R / (B*(R+B))
v_bar(rb, ellipse 2)   = sqrt(2R/(B*(R+B)))
```

Both are apoapsis speeds at the same radius `rb`; ellipse 2 has the larger
semi-major axis (`r2 > r1`), so it is the faster orbit. The burn is **prograde**:

```
dv_bar_B2 = sqrt(2R/(B*(R+B))) - sqrt(2/(B*(1+B)))                         (B2)
```

**Positivity.** `dv_bar_B2 > 0` iff `2R/(B(R+B)) > 2/(B(1+B))`, i.e.
`R(1+B) > (R+B)`, i.e. `R + RB > R + B`, i.e. `RB > B`, i.e. **`R > 1`**.
So burn 2 is positive for every admissible `B` whenever `R > 1`.

### 3.3 Burn 3 — at `r2`, ellipse 2 -> circular

`r2` is the **periapsis** of ellipse 2:

```
v^2 = mu*(2/r2 - 2/(r2 + rb)) = 2*mu*rb / (r2*(r2 + rb))
```

Dividing by `mu/r1`:

```
v_bar(r2, ellipse 2)^2 = 2*rb*r1 / (r2*(r2+rb)) = 2B / (R*(R+B))
v_bar(r2, ellipse 2)   = sqrt(2B/(R*(R+B)))
```

The required circular speed at `r2` is `1/sqrt(R)`. Now

```
2B/(R+B) > 1   <=>   2B > R + B   <=>   B > R
```

which holds by construction. So the spacecraft arrives at `r2` **faster** than
circular and must **decelerate**. This is the key qualitative difference from the
Hohmann transfer, whose second burn is prograde:

```
dv_bar_B3 = sqrt(2B/(R*(R+B))) - 1/sqrt(R)                                 (B3)
```

**Burn 3 is retrograde.** `dv_bar_B3` as written is the (positive) *magnitude* of
that retrograde impulse, positive for every `B > R`.

### 3.4 Total

```
dv_bar_B(R,B) = sqrt(2B/(1+B)) - 1
              + sqrt(2R/(B*(R+B))) - sqrt(2/(B*(1+B)))
              + sqrt(2B/(R*(R+B))) - 1/sqrt(R)                             (B)
```

### 3.5 Verification of §3

The expressions given as the target for this milestone were **not** assumed. Each
burn was rederived above from vis-viva, and then checked independently:

- **Physical speeds before/after every burn** were computed from dimensional
  vis-viva and compared against the normalized closed forms.
- **Burn directions** were determined by comparing the two speeds at each burn
  radius, not by assumption: burn 1 prograde, burn 2 prograde, burn 3 **retrograde**.
- **Positivity for `B > R > 1`** is proved algebraically above for each of the three
  burns, and confirmed numerically over a wide sample grid.
- **Dimensional agreement.** `dv_bar_B(R,B)` matches
  `(|dv1| + |dv2| + |dv3|)/v1` computed from dimensional vis-viva to `3.4e-41`
  at 40-digit precision, over `R ∈ {1.5 … 200}` and `B/R ∈ {1.0001, 1.5, 2, 5, 50, 1000}`.

The target expressions were found to be **correct**; no correction was required.
One sign-convention slip was made and caught during this milestone's own scratch
checking — burn 3's *magnitude* is `v_transfer - v_circular > 0` while its
*direction* is retrograde. Both facts are recorded above so the M2 implementation
cannot repeat the confusion.

---

## 4. Limiting behaviour

### 4.1 `B -> R+` reduces exactly to Hohmann

Setting `B = R` in (B):

```
dv_bar_B1 -> sqrt(2R/(1+R)) - 1                      = dv_bar_H1
dv_bar_B2 -> sqrt(2R/(R*2R)) - sqrt(2/(R*(1+R)))
           = 1/sqrt(R) - sqrt(2/(R*(1+R)))           = dv_bar_H2
dv_bar_B3 -> sqrt(2R/(R*2R)) - 1/sqrt(R)
           = 1/sqrt(R) - 1/sqrt(R)                   = 0
```

so

```
dv_bar_B(R, R) = dv_bar_H(R)                                               (LIM-1)
```

Physically: ellipse 2 degenerates to the circular target orbit, burn 3 vanishes,
and the manoeuvre *is* the Hohmann transfer. Verified numerically to 20 significant
figures at `R = 2, 10, 50`. **(LIM-1) is the structural fact that makes the second
threshold well posed** — see §5.2.

### 4.2 `B -> infinity` — derived analytically, not sampled

Take the limit of (B) term by term. Each term is elementary, so no numerical
"large B" substitute is needed:

```
lim sqrt(2B/(1+B)) - 1                  = sqrt(2) - 1
 B->inf

lim sqrt(2R/(B*(R+B))) - sqrt(2/(B*(1+B))) = 0 - 0 = 0      (both decay as 1/B)
 B->inf

lim sqrt(2B/(R*(R+B))) - 1/sqrt(R)      = sqrt(2/R) - 1/sqrt(R)
 B->inf                                  = (sqrt(2) - 1)/sqrt(R)
```

Therefore

```
dv_bar_B_inf(R) = (sqrt(2) - 1) * (1 + 1/sqrt(R))                          (BINF)
```

This has a clean interpretation: the first burn becomes a **parabolic escape**
impulse from `r1` (`sqrt(2) - 1` is exactly the escape-minus-circular increment),
the middle burn becomes free, and the third burn is the mirror-image parabolic
capture increment at `r2`, scaled by `1/sqrt(R)`.

Verified: `dv_bar_B(R, 1e30)` agrees with (BINF) to 20 significant figures at
`R = 2, 12, 20`.

### 4.3 Large-`B` approach direction

Expanding (B) for large `B`:

```
dv_bar_B(R,B) = dv_bar_B_inf(R) + (sqrt(R) - 3)/(sqrt(2)*B) + O(1/B^2)     (ASY)
```

Confirmed numerically: at `B = 1e8` the predicted and actual residuals agree to
8 significant figures for `R = 4, 16, 25`, and the residual vanishes at `R = 9`.

So `R = 9` is where the asymptote is approached from below (`R < 9`) rather than
from above (`R > 9`). **`R = 9` is not a crossover threshold** — it never involves
comparing against the Hohmann transfer. It is a *shape* number, and it is recorded
here only because it explains why `dv_bar_B(R, ·)` stops being monotone (§6).

---

## 5. The two crossover thresholds

There are **two distinct radius-ratio thresholds**, and they answer **two different
questions**. Confusing them is the single most common error in this problem.

### 5.1 Threshold #1 — `R1*`: when the infinite-`B` limit first beats Hohmann

**Mathematical condition:**

```
dv_bar_H(R) = dv_bar_B_inf(R)                                              (T1)
```

Substituting (H') and (BINF) and simplifying (multiply through by `sqrt(R)`, then
divide by `sqrt(2)`), (T1) reduces to the compact form

```
(R - 1) / sqrt(R + 1) = sqrt(R) + 1 - sqrt(2)                              (T1')
```

Squaring (T1') and writing `u = sqrt(R)` yields an exact algebraic characterization:

```
u^3 - (1 + 2*sqrt(2))*u^2 + u + 1 = 0,     R1* = u^2  (largest real root)   (T1'')
```

**Independently computed value** (four independent routes: Newton on (T1); Newton on
(T1'); bisection on the raw *dimensional* difference; and the cubic (T1'')):

```
R1* = 11.938765472645870715530055180410
```

The literature value was **not** used as an input; see §9 item I.

**What it means physically.** For `R > R1*`, a bi-elliptic transfer with a
sufficiently distant intermediate apoapsis costs less total delta-v than the
Hohmann transfer. Below `R1*`, **no** bi-elliptic transfer of any `B` can beat
Hohmann (proved in §6). `R1*` is therefore the ratio at which the bi-elliptic
family *first becomes competitive at all*.

**`R1*` is an asymptote, not a manoeuvre.** Achieving `dv_bar_B_inf` requires
`B -> infinity`, i.e. **unbounded transfer time** (§8). It is a mathematical
delta-v bound, not a flight recommendation. Just above `R1*` the delta-v prize is
vanishingly small while the required `B` — and hence the transfer duration — is
enormous: at `R = 12`, only 3.0 m/s is available even at infinite `B`, and the
smallest `B` that breaks even is `B ≈ 816`, a 524-day transfer against 0.52 days
for Hohmann (§8.3).

### 5.2 Threshold #2 — `R2*`: when *every* bi-elliptic beats Hohmann

**Mathematical condition.** By (LIM-1), `dv_bar_B(R, B) = dv_bar_H(R)` at `B = R`.
So whether a *marginally* bi-elliptic transfer (an apoapsis barely above the target
radius) helps or hurts is decided by the **sign of the slope at that endpoint**:

```
d/dB [ dv_bar_B(R,B) ] |_(B = R+)  =  0                                    (T2)
```

Differentiating (B) term by term and evaluating at `B = R` (full working reproduced
in the scratch derivation; each of the four terms contributes a closed-form
expression in `R`), the whole slope collapses to

```
R^(3/2) * d(dv_bar_B)/dB |_(B=R+)  =  (1 + 3R) / (sqrt(2) * (1+R)^(3/2)) - 1/2
```

Setting the bracket to zero gives `(1+R)^(3/2) = sqrt(2)*(1+3R)`, and squaring gives
`(1+R)^3 = 2*(1+3R)^2`, i.e. an **exact cubic**:

```
R^3 - 15*R^2 - 9*R - 1 = 0                                                 (T2')
```

**Independently computed value** (three independent routes: the cubic (T2'); the
pre-squaring form; and a purely numerical derivative of `dv_bar_B` in `B` at `B=R+`,
which uses none of the hand algebra):

```
R2* = 15.581718738763179213242444820820
```

The cubic has exactly one positive real root (the other two are
`-0.4338` and `-0.1480`), so `R2*` is unambiguous.

**What it means physically.** For `R > R2*`, the slope at `B = R+` is negative:
raising the intermediate apoapsis *even infinitesimally* above the target radius
already reduces total delta-v. Every admissible bi-elliptic transfer beats Hohmann.
For `R < R2*` the slope is positive, so a *modest* bi-elliptic transfer is worse
than Hohmann even where a very distant one would win.

### 5.3 Why the two thresholds are different — stated so they cannot be confused

| | Threshold #1 | Threshold #2 |
|---|---|---|
| Symbol | `R1*` | `R2*` |
| Value | `11.938765472645870716` | `15.581718738763179213` |
| Condition | `dv_bar_H(R) = dv_bar_B_inf(R)` | `∂dv_bar_B/∂B = 0` at `B = R+` |
| Endpoint examined | `B -> infinity` | `B -> R+` |
| Exact form | `u^3-(1+2√2)u^2+u+1=0`, `u=√R` | `R^3-15R^2-9R-1=0` |
| Question answered | Can **any** bi-elliptic beat Hohmann? | Do **all** bi-elliptics beat Hohmann? |
| Quantifier | **∃** `B` | **∀** `B` |
| Crossing at | the far endpoint of the `B` domain | the near endpoint of the `B` domain |

In one sentence: **`R1*` is where the bi-elliptic family first contains a winner;
`R2*` is where the bi-elliptic family contains nothing but winners.** They differ
because they test opposite ends of the `B` domain, and `R1* < R2*` because the far
endpoint becomes favourable before the near endpoint does.

---

## 6. Region structure — and a correction to the assumed region definitions

### 6.1 What the delta-v curve actually looks like in `B`

For fixed `R`, `dv_bar_B(R, ·)` on `B ∈ (R, ∞)` starts at `dv_bar_H(R)` (by LIM-1)
and ends at `dv_bar_B_inf(R)` (by BINF). Its shape is governed by the two endpoint
slopes: the `B = R+` slope changes sign at `R2*` (§5.2), and the large-`B` slope
changes sign at `R = 9` (§4.3). This gives three shapes:

| `R` range | slope at `B=R+` | slope at large `B` | shape |
|---|---|---|---|
| `1 < R < 9` | `+` | `+` | monotonically **increasing** |
| `9 < R < R2*` | `+` | `-` | **up, then down** (interior maximum) |
| `R > R2*` | `-` | `-` | monotonically **decreasing** |

**A dense brute-force scan found zero interior local minima** — 3000 logarithmically
spaced `B` samples per `R`, over `R ∈ [1.5, 370]`, `B/R` up to `1e14`. The curve is
increasing, decreasing, or up-then-down; it never dips and recovers.

**Consequence (important for M2):**

```
inf over B of dv_bar_B(R,B) = min( dv_bar_H(R), dv_bar_B_inf(R) )          (INF)
```

The infimum is **always at an endpoint and never attained at a finite interior `B`**.
There is no finite-`B` optimum to find. A Milestone 2 "finite-`B` optimizer" will
therefore correctly report that the objective is monotone toward an endpoint over
the relevant range; §9 item D specifies the brute-force cross-check that must
confirm this rather than an interior stationary point.

### 6.2 The three regions

| Region | Range | Behaviour |
|---|---|---|
| **A** | `1 < R < R1* = 11.9388` | Hohmann is strictly better than **every** admissible bi-elliptic transfer (all `B > R`). |
| **B** | `R1* < R < R2*` (`11.9388 … 15.5817`) | Bi-elliptic beats Hohmann **only for `B > B_crit(R)`**, a finite threshold. `B_crit -> ∞` as `R -> R1*+` and `B_crit -> R` as `R -> R2*-`. |
| **C** | `R > R2* = 15.5817` | Bi-elliptic beats Hohmann for **every** `B > R`. `B_crit = R`. |

`B_crit(R)`, the smallest `B` that breaks even against Hohmann, is the practically
meaningful quantity in Region B:

| `R` | `B_crit` | `rb_crit / r2` | `rb_crit` altitude (Earth ref) |
|---|---|---|---|
| 11.94 | 40370.6 | 3381 | 2.70e8 km |
| 12 | 815.820 | 68.0 | 5.44e6 km |
| 13 | 48.9048 | 3.76 | 3.20e5 km |
| 14 | 26.1046 | 1.86 | 1.68e5 km |
| 15 | 18.1903 | 1.21 | 1.15e5 km |
| `R2*` = 15.5817 | 15.5817 | 1.000 | 9.77e4 km |
| ≥ `R2*` | `= R` | 1.000 | — |

### 6.3 Correction: the assumed "Region B" definition is empty

This milestone was specified with a candidate Region B reading *"some finite-`B`
bi-elliptic transfers can beat Hohmann, but the `B -> infinity` limit does not yet
beat Hohmann."* **That region does not exist.** It is empty, and the reason follows
directly from (INF):

> Suppose some finite `B0 > R` satisfies `dv_bar_B(R,B0) < dv_bar_H(R)`. Since
> `dv_bar_B(R,R) = dv_bar_H(R)` and the curve has no interior local minimum, `B0`
> must lie on a decreasing branch that continues decreasing all the way to the
> `B -> ∞` endpoint. Hence `dv_bar_B_inf(R) < dv_bar_B(R,B0) < dv_bar_H(R)`.
> So a finite-`B` win **forces** an infinite-`B` win. The converse of the assumed
> Region B is impossible.

Checked numerically as well: over every `R` on a grid where
`dv_bar_B_inf(R) >= dv_bar_H(R)`, the minimum of `dv_bar_B` over 2001 logarithmically
spaced `B` values up to `B/R = 1e15` never fell below `dv_bar_H(R)`. **Zero violations.**

The region structure in §6.2 is the corrected version and is what M2 will implement.
The intended three-way distinction survives intact — it is just anchored differently:
Region B is *"only sufficiently distant bi-elliptics win"*, not *"finite wins while
infinite loses"*.

---

## 7. Representative hand calculations

All values independently recomputed at 30-digit precision.

### 7.1 Normalized total delta-v (units of `v1`)

| `R` | `dv_bar_H` | `dv_bar_B` (`B=2R`) | `dv_bar_B` (`B=5R`) | `dv_bar_B` (`B=1e6`, sanity) | `dv_bar_B_inf` (exact) | better | region |
|---|---|---|---|---|---|---|---|
| 2 | 0.28445705 | 0.46632139 | 0.60189809 | 0.70710566 | 0.70710678 | Hohmann | A |
| 5 | 0.48000915 | 0.54094283 | 0.57688720 | 0.59945496 | 0.59945550 | Hohmann | A |
| 10 | 0.52978752 | 0.54261935 | 0.54594466 | 0.54519951 | 0.54519939 | Hohmann | A |
| 12 | 0.53417987 | 0.53923048 | 0.53773598 | 0.53378705 | 0.53378672 | bi-elliptic | B |
| 15 | 0.53621819 | 0.53385750 | 0.52794811 | 0.52116366 | 0.52116304 | bi-elliptic | B |
| 16 | 0.53623939 | 0.53211454 | 0.52518686 | 0.51776766 | 0.51776695 | bi-elliptic | C |
| 20 | 0.53473136 | 0.52563061 | 0.51592650 | 0.50683557 | 0.50683453 | bi-elliptic | C |
| 50 | 0.51369584 | 0.49665076 | 0.48341568 | 0.47279508 | 0.47279221 | bi-elliptic | C |

The `B = 1e6` column is the required large-`B` numerical sanity check: it tracks
the exact `B -> ∞` column to ~6 decimal places, approaching it **from below** for
`R < 9` (see `R = 2, 5`) and **from above** for `R > 9` (see `R = 10, 12, 15, …`),
exactly as (ASY) predicts.

**`B = 2R` and `B = 5R` are arbitrary samples, not optima.** By §6.1 there is no
finite-`B` optimum. Note `R = 10`, where `dv_bar_B(B=5R) > dv_bar_B(B=1e6)`: this is
the up-then-down shape, with the interior *maximum* near `B ≈ 66.5`. Reading the
`B=2R`/`B=5R` columns as if they bracketed an optimum would be wrong.

**`R = 12` is the clearest illustration of the two-threshold structure.** It sits in
Region B: `dv_bar_H = 0.53418` beats both `B=2R` and `B=5R`, yet `dv_bar_B_inf =
0.53379` beats Hohmann. The bi-elliptic wins there only at very large `B`.

### 7.2 Earth reference case, dimensional

`r1 = 6678.1363 km` (300 km altitude), `v1 = 7.7257606 km/s`.
Final **altitude** is a conversion from the universal **radius ratio** result:
`h2 = R*r1 - R_Earth`.

| `R` | `h2` [km] (converted) | `dv_H` [km/s] | `dv_B(2R)` | `dv_B(5R)` | `dv_B_inf` | saving [m/s] | saving [%] |
|---|---|---|---|---|---|---|---|
| 2 | 6978.1 | 2.197647 | 3.602687 | 4.650121 | 5.462938 | −3265.3 | −148.6 |
| 5 | 27012.5 | 3.708436 | 4.179195 | 4.456892 | 4.631250 | −922.8 | −24.9 |
| 10 | 60403.2 | 4.093012 | 4.192147 | 4.217838 | 4.212080 | −119.1 | −2.9 |
| 12 | 73759.5 | 4.126946 | 4.165966 | 4.154419 | 4.123908 | **+3.0** | +0.07 |
| 15 | 93793.9 | 4.142693 | 4.124455 | 4.078801 | 4.026381 | +116.3 | +2.8 |
| 16 | 100472.0 | 4.142857 | 4.110990 | 4.057468 | 4.000144 | +142.7 | +3.4 |
| 20 | 127184.6 | 4.131206 | 4.060896 | 3.985925 | 3.915682 | +215.5 | +5.2 |
| 50 | 327528.7 | 3.968691 | 3.837005 | 3.734754 | 3.652679 | +316.0 | +8.0 |

"saving" is `dv_H - dv_B_inf`, i.e. the **maximum theoretically available** prize,
attainable only at infinite transfer time. Negative means Hohmann wins.

Threshold radii in the Earth reference case:

```
R1* = 11.9388  ->  r2 =  79728.70 km   (altitude 73350.57 km)
R2* = 15.5817  ->  r2 = 104056.84 km   (altitude 97678.71 km)
```

**Practical note:** GEO (`r2 = 42164 km`) from a 300 km parking orbit is
`R = 6.31` — deep inside Region A. For the single most common high-orbit transfer
flown, Hohmann is unambiguously better and the bi-elliptic question does not arise.

---

## 8. Transfer time

### 8.1 Hohmann

The transfer is exactly half of the transfer ellipse's period:

```
a_H = (r1 + r2)/2
t_H = pi * sqrt(a_H^3 / mu)                                                (TH)
```

### 8.2 Bi-elliptic

Two half-ellipses, each traversed once:

```
a_1 = (r1 + rb)/2
a_2 = (r2 + rb)/2
t_B = pi*sqrt(a_1^3/mu) + pi*sqrt(a_2^3/mu)                                (TB)
```

Normalized by `t_ref = sqrt(r1^3/mu)`, both are functions of `R` and `B` alone:
`t_bar_H = pi*((1+R)/2)^(3/2)` and `t_bar_B = pi*[((1+B)/2)^(3/2) + ((R+B)/2)^(3/2)]`.
Since `t_bar_B` grows as `B^(3/2)`, **`B -> ∞` implies `t_B -> ∞`**: the delta-v
asymptote (BINF) is reached only in infinite time.

### 8.3 The trade, quantified (Earth reference, days)

| `R` | `t_H` | `t_B(B=2R)` | `t_B(B=5R)` | `B_crit` | `t_B(B_crit)` | `t_B/t_H` at break-even |
|---|---|---|---|---|---|---|
| 2 | 0.0577 | 0.288 | 0.867 | never | — | — |
| 5 | 0.163 | 1.051 | 3.299 | never | — | — |
| 10 | 0.405 | 2.895 | 9.212 | never | — | — |
| 12 | 0.521 | 3.789 | 12.08 | 815.8 | **524.1** | **1006×** |
| 15 | 0.711 | 5.272 | 16.85 | 18.19 | 3.06 | 4.30× |
| 16 | 0.779 | 5.802 | 18.55 | 16.0 | 2.79 | 3.58× |
| 20 | 1.069 | 8.082 | 25.89 | 20.0 | 3.88 | 3.63× |
| 50 | 4.047 | 31.69 | 101.9 | 50.0 | 15.16 | 3.75× |

The four conclusions this milestone commits to:

1. Bi-elliptic can reduce delta-v, but only for sufficiently large `R` (`R > R1*`).
2. The price is dramatically longer transfer duration — at minimum ~3.6× the
   Hohmann time even to *break even* in Region C.
3. Increasing `B` keeps reducing delta-v in Regions B and C while making duration
   progressively impractical (`t ~ B^(3/2)`).
4. `B -> ∞` has infinite duration, so the delta-v asymptote is never attainable.

**The `R = 12` row is the decisive engineering datum:** just inside Region B, the
*entire theoretical prize* is 3.0 m/s, and merely breaking even costs a 1006×
increase in transfer time — 524 days instead of 12.5 hours. The final
recommendation must therefore be an **engineering trade**, not "choose minimum
delta-v".

---

## 9. Verification plan for later milestones

This is the acceptance matrix M2+ must satisfy. Every check is against an
*independent* computation, never against a number transcribed from this document
alone (except the regression targets in §10, which are themselves the output of
multiple independent routes).

**A. Hohmann limiting behaviour**
- `dv_bar_H(1) == 0` exactly; `dv_bar_H(1+eps) -> 0` with leading order `(R-1)/2`.
- Closed form (H) vs. an independent two-impulse vis-viva computation, over a wide
  `R` grid, to tight tolerance.
- Both impulses individually positive for all `R > 1`.

**B. Bi-elliptic limiting behaviour**
- `B -> R+`: `dv_bar_B(R,B) -> dv_bar_H(R)` and `dv_bar_B3 -> 0` (LIM-1).
- `B -> ∞`: `dv_bar_B(R,B) -> (sqrt(2)-1)(1+1/sqrt(R))` (BINF), with the residual
  matching `(sqrt(R)-3)/(sqrt(2)B)` (ASY).
- Direct dimensional vis-viva vs. the normalized equations, all three burns.
- Burn *directions* asserted, not just magnitudes: burns 1 and 2 prograde, burn 3
  retrograde for all `B > R > 1`.

**C. Scaling / invariance**
- `dv_dimensional ∝ sqrt(mu/r1)`: scale `mu` and `r1` and confirm delta-v scales
  exactly as `v1`.
- Crossover `R` values are invariant to `mu` and `r1`. *Already demonstrated:* `R1*`
  recomputed with Earth, Mars, and an artificially scaled body (`mu × 1e6`,
  `r1 × 1e6`) reproduces `11.93876547264587071553006` identically to 25 digits.

**D. Finite-`B` optimization**
- Any M2 optimizer must agree with an independent dense brute-force `B` sweep.
- The sweep must confirm §6.1: **no interior local minimum**, and
  `inf_B dv_bar_B = min(dv_bar_H, dv_bar_B_inf)` (INF). An optimizer reporting an
  interior optimum is a bug, not a discovery.
- The interior *maximum* for `9 < R < R2*` must be located consistently by
  optimizer and brute force (e.g. `B_max ≈ 66.47` at `R = 10`, `26.47` at `R = 12`).

**E. Crossover roots**
- `R1*` reproduced by at least two independent formulations: the direct condition
  (T1), the reduced form (T1'), the algebraic cubic (T1''), and a bisection on the
  raw dimensional difference. All must agree to the precision in §10.
- `R2*` reproduced by the cubic (T2'), the pre-squaring form, and a purely numerical
  `∂/∂B` at `B = R+` that uses none of the hand algebra.
- `B_crit(R)` in Region B verified by bracketing rather than by a fitted formula,
  and checked at the endpoints: `B_crit -> ∞` as `R -> R1*+`, `B_crit -> R` as
  `R -> R2*-`.

**F. Continuity**
- `dv_bar_H(R)` continuous on `R > 1`; `dv_bar_B(R,·)` continuous on `B > R`;
  `dv_bar_B(R,B)` jointly continuous. Check for NaN/overflow at large `B` and for
  catastrophic cancellation in (B2), where two nearly equal `O(1/B)` square roots
  are subtracted — this is the numerically dangerous term and needs an explicit
  high-precision or algebraically rearranged cross-check.

**G. Time of flight**
- Half-period formulas (TH), (TB) vs. an independent Kepler propagation
  (solve Kepler's equation from periapsis to apoapsis and compare elapsed time).
- Normalized and dimensional time formulations must agree.

**H. Dimensional / nondimensional round trip**
- Compute dimensionally, normalize by `v1`, and recover the normalized values;
  and the reverse. Must close to machine precision.

**I. Literature sanity check — last, never first**
- Only *after* independent derivation, compare `R1*` and `R2*` against classical
  textbook values (~11.94 and ~15.58). **Confirmed consistent.** These literature
  numbers must never appear as production constants, hardcoded thresholds, test
  fixtures, or optimizer seeds. Production code computes the roots.

---

## 10. Regression targets for Milestone 2

Recorded to more precision than any M2 test will need, so tolerances can be chosen
freely. Each was produced by at least three independent routes (§5).

```
R1* = 11.938765472645870715530055180410063729     (dv_bar_H = dv_bar_B_inf)
R2* = 15.581718738763179213242444820819518455     (d/dB dv_bar_B = 0 at B = R+)

exact characterizations:
  R1* = u^2, u the largest real root of  u^3 - (1 + 2*sqrt(2))*u^2 + u + 1 = 0
  R2* = the unique positive root of      R^3 - 15*R^2 - 9*R - 1 = 0

dv_bar_B_inf(R) = (sqrt(2) - 1) * (1 + 1/sqrt(R))          [closed form, exact]
dv_bar_B(R,R)   = dv_bar_H(R)                              [exact identity]
large-B residual: dv_bar_B - dv_bar_B_inf ~ (sqrt(R)-3)/(sqrt(2)*B)

Earth reference: r1 = 6678.1363 km, v1 = 7.72576063698292 km/s
  R1* -> r2 =  79728.7030801 km  (altitude  73350.5667801 km)
  R2* -> r2 = 104056.841526  km  (altitude  97678.7052257 km)
```

---

## 11. Engineering recommendation criteria (to be applied in later milestones)

The final deliverable is a recommendation, not a minimization. It will weigh:

1. **Total delta-v** — the raw objective, but never the sole criterion.
2. **Absolute delta-v saving in m/s** — a 3 m/s saving (Region B at `R = 12`) is
   inside the noise of launch dispersion, navigation and finite-burn losses.
3. **Percentage delta-v saving** — and its propellant-mass consequence through the
   rocket equation, which is what actually matters to a mission.
4. **Transfer duration** — absolute, and as a multiple of the Hohmann time.
5. **Intermediate apoapsis radius and altitude** — is `rb` even physically sensible?
6. **Practicality of very high apoapsis** — at `R = 12`, `B_crit ≈ 816` puts
   apoapsis at ~5.4 million km, roughly 14× the lunar distance. That is not an
   Earth orbit in any operational sense.
7. **Perturbation exposure** — at such radii the two-body model is not merely
   inaccurate, it is the wrong model (lunar/solar third-body effects dominate).
8. **Radiation and environment** — repeated slow passages through the trapped
   radiation belts, and prolonged exposure outside the magnetosphere.
9. **Mission timeline** — launch windows, commissioning schedules, spacecraft
   consumables and staffing cost over a multi-month or multi-year cruise.
10. **Whether the saving is operationally material** — a saving smaller than the
    delta-v budget margin buys nothing.

**Explicit commitment:** this project will not claim that a mathematically minimum
transfer is automatically mission-optimal. The `R = 12` case is the standing
counterexample — mathematically bi-elliptic wins, operationally it is absurd.

---

## 12. Scope and limitations

This analysis is deliberately idealized. It assumes:

- **coplanar** circular initial and final orbits;
- **ideal two-body** central gravity (point-mass `mu`);
- **impulsive** burns (instantaneous velocity changes);
- **no finite-burn losses** (no gravity losses, no steering losses, no throttle profile);
- **no plane change** (no inclination or RAAN change, and no combined-manoeuvre optimization);
- **no J2** or any other non-spherical gravity harmonic;
- **no third-body perturbations** — notably no lunar or solar gravity, which for the
  very large `rb` values in Regions B and C would dominate the dynamics entirely;
- **no atmospheric drag**;
- **no lunar / sphere-of-influence effects** for enormous apoapses, even where `rb`
  exceeds the Earth's Hill radius (~1.5e6 km) and the two-body model is simply invalid;
- **no launch constraints**, launch windows, or injection accuracy;
- **no radiation model**;
- **no collision, conjunction or debris-environment model**;
- **no mission-specific operational optimization** (consumables, staffing, ground
  contact, schedule risk).
- **`B -> ∞` is a mathematical limit only.** It has infinite transfer time and is
  never a flight recommendation.

**No flight-operations or mission-optimality claims are made anywhere in this
project.** The deliverable is a normalized delta-v crossover analysis plus an
engineering discussion of when the result could matter.

---

## 13. Milestone roadmap

| Milestone | Content | Status |
|---|---|---|
| **M1** | Problem definition, analytical derivation, hand calculations, verification plan | **complete** |
| M2 | Production transfer solver (Hohmann + bi-elliptic), normalized and dimensional, with the §9 verification suite | pending |
| M3 | Parameter sweeps in `R` and `B`; `B_crit(R)`; independent root-finding for both thresholds | pending |
| M4 | Crossover radius-ratio plot and supporting figures | pending |
| M5 | Transfer-time trade analysis and combined delta-v/time figures | pending |
| M6 | Engineering recommendation and portfolio write-up | pending |

---
---

# DESIGN — Milestone 2

**Production transfer equations, dimensional cross-checks, break-even solver, and
numerical verification.**

Status: **M2 complete.** Everything in the M1 sections above is unchanged; no M1
headline value was altered. M2 turns that derivation into tested production code
and re-derives every number from it. The final portfolio crossover figure and the
engineering recommendation remain out of scope (M4 and M6).

## M2.1 Production API

`src`-layout package `bielliptic_crossover`, version `0.2.0`.

| Module | Contents |
|---|---|
| `constants.py` | `MU_EARTH`, `R_EARTH`, `H1_REFERENCE`, `R1_REFERENCE`, `V1_REFERENCE`, `MU_MARS`, `R_MARS`, `altitude_from_radius`, `radius_from_altitude`. **Contains no crossover constants** — a test asserts the strings `11.93` / `15.58` never appear in it. |
| `_stable.py` | `sqrt1pm1(x) = sqrt(1+x) - 1` evaluated without cancellation. |
| `hohmann.py` | `hohmann_burns_normalized(R)`, `hohmann_total_normalized(R)`, `hohmann_total_dv(mu, r1, r2)`. |
| `bielliptic.py` | `bielliptic_burns_normalized(R, B)`, `bielliptic_total_normalized(R, B)`, `bielliptic_infinite_limit_normalized(R)`, `bielliptic_total_derivative_wrt_B(R, B)`, `burn_directions()`. |
| `timing.py` | `time_scale`, `hohmann_transfer_time[_normalized]`, `bielliptic_transfer_time[_normalized]`, `bielliptic_time_excess_at_B_equals_R`. |
| `dimensional.py` | Independent direct-vis-viva path: `vis_viva_speed`, `circular_speed`, `hohmann_burns_direct`, `hohmann_total_dv_direct`, `bielliptic_burns_direct`, `bielliptic_total_dv`. Exposes signed velocity changes and every apsis speed so burn direction can be checked physically. |
| `crossover.py` | `threshold_R1`, `threshold_R1_from_polynomial`, `threshold_R2`, `threshold_R2_from_cubic`, `break_even_B`, `BreakEvenResult`, `classify_radius_ratio`, `scan_B_structure`, `BStructure`, `bielliptic_infimum_normalized`, `B_RATIO_RESOLVABLE_MAX`, regime constants. |

Domain rules are enforced, not assumed: `R >= 1`, `B >= R`, positive `mu` and `r1`,
`r2 >= r1`, `rb >= r2`. NaN and infinity are rejected with explicit messages —
`bielliptic_total_normalized(R, inf)` directs the caller to the analytical limit
instead of silently returning a large-`B` approximation.

**Burn magnitudes are always returned positive; direction is documented, never
encoded as a sign.** `burn_directions()` returns
`("prograde", "prograde", "retrograde")`.

## M2.2 Cancellation-free formulation

Every burn has the shape `sqrt(1+x) - 1` for an `x` that vanishes at a physically
interesting boundary. The production code uses the exact algebraic rewrite
`sqrt(1+x) - 1 = x / (sqrt(1+x) + 1)`:

| Burn | `x` |
|---|---|
| Hohmann 1 | `(R-1)/(1+R)` |
| Hohmann 2 | `(1-R)/(1+R)`, scaled by `-1/sqrt(R)` |
| Bi-elliptic 1 | `(B-1)/(1+B)` |
| Bi-elliptic 2 | `B(R-1)/(R+B)`, scaled by `sqrt(2/(B(1+B)))` |
| Bi-elliptic 3 | `(B-R)/(R+B)`, scaled by `1/sqrt(R)` |

These are identities, not approximations. Two consequences are load-bearing:

- **`dv3` is exactly `0.0` at `B = R`** (since `sqrt1pm1(0) == 0`), so burn 3
  vanishes exactly rather than to a tolerance.
- **`dv1` at `B = R` is bitwise identical to Hohmann's first impulse**, because
  the two reduce to the same expression.

The `B = R` total identity `dv_bar_B(R,R) = dv_bar_H(R)` is exact mathematically
and agrees to **at most 1 ulp** in floating point over `R ∈ [1.5, 500]` — burn 2
is reached by two different but algebraically equal expressions, so the last bit
may differ. Tests assert `<= 2 ulp` rather than claiming bitwise equality.

## M2.3 Thresholds recomputed

No crossover constant is hardcoded. Each threshold is solved twice by paths that
share no algebra.

| | `R1*` | `R2*` |
|---|---|---|
| Path 1 (transfer equations) | `brentq` on `dv_bar_H(R) - dv_bar_B_inf(R)` → `11.9387654726458923` | `brentq` on `d(dv_bar_B)/dB` at `B=R` → `15.5817187387631790` |
| Path 2 (exact polynomial) | `u³-(1+2√2)u²+u+1=0`, `u=√R` → `11.9387654726458727` | `R³-15R²-9R-1=0` → `15.5817187387631790` |
| Difference between paths | `1.954e-14` | `0.000e+00` |

*(Path-2 values updated in M4: the polynomials are now solved by bracketed
root-finding rather than by `numpy.roots`, which was not reproducible across
BLAS builds. See M4.5b. The `R2*` paths now agree exactly.)*
| M1 accepted value | `11.938765472645870716` | `15.581718738763179213` |

Both agree with M1 to the limit of double precision. Residuals at the recomputed
roots: `dv_bar_H(R1*) - dv_bar_B_inf(R1*) = 0.000e+00`, and
`d(dv_bar_B)/dB|(R2*,R2*) = -8.674e-19`.

Separation `R2* - R1* = 3.6429532661172868`.

The analytical `d/dB` is cross-checked against **complex-step differentiation**,
which shares none of the hand algebra and suffers no subtractive cancellation;
they agree to `1e-12` relative across the tested grid.

## M2.4 Break-even `B_crit(R)`

`break_even_B(R)` returns a `BreakEvenResult` carrying the regime, the nontrivial
`B_crit`, and the infimum of the winning set.

**The trivial root at `B = R` is explicitly excluded.** `dv_bar_B(R,R) = dv_bar_H(R)`
holds by construction, so a naive root-finder would return it. The solver probes
at `B = R(1 + 1e-6)`, expands geometrically until the excess changes sign, and
solves in `log B` (the root spans orders of magnitude). It asserts `B_crit > R`
before returning.

| Region | Behaviour | `B_crit` | `winning_B_infimum` |
|---|---|---|---|
| A (`R < R1*`) | Hohmann beats every `B` | `None` | `None` |
| B (`R1* < R < R2*`) | only `B > B_crit` wins | finite root | `B_crit` |
| C (`R > R2*`) | every `B > R` wins | `None` | `R` (open boundary, not attained) |

Region C returns no finite root because none exists: the winning set is the open
interval `(R, ∞)`, whose infimum `B = R` is not attained. Manufacturing a root
there would be fiction. If `R` is so close to `R1*` that `B_crit` exceeds the
search ceiling (`1e30`), the solver **raises** rather than misreporting Region A.

Production values (residual `dv_bar_B(R,B_crit) - dv_bar_H(R)` is exactly `0.0`):

| `R` | `B_crit` | `rb_crit/r2` | `rb_crit` altitude (Earth ref) |
|---|---|---|---|
| 12 | `815.8202504753092` | 67.985021 | 5.4418e+06 km |
| 13 | `48.90484332838889` | 3.761911 | 3.2022e+05 km |
| 14 | `26.10461128235042` | 1.864615 | 1.6795e+05 km |
| 15 | `18.19028151222189` | 1.212685 | 1.1510e+05 km |

M1 quoted `B_crit(12) ≈ 815.8`. The production value `815.8202504753092` confirms
it; the M1 rounded figure was never assumed to be exact.

## M2.5 Numerical stability audit

Findings, all reproduced as tests in `tests/test_stability.py`.

**Finding 1 — a genuine defect in the naive burn-3 formula.** Evaluating
`sqrt(2B/(R(R+B))) - 1/sqrt(R)` directly loses precision as `B -> R+`: relative
error `5.8e-12` at `B/R-1 = 1e-4`, `9.2e-4` at `1e-12`, and `8.3e-2` at `1e-14`.
At exactly `B = R` it returns a **strictly negative** value — a physically
impossible negative burn magnitude — for **58 of 398** sampled `R` in `[1.5, 200]`
(e.g. `-1.110e-16` at `R = 1.5`). The production form returns exactly `0.0` in all
398 cases. This matters because the break-even solver must distinguish the trivial
`B = R` root from a real one, and sign noise there is precisely the failure mode.

**Finding 2 — `R -> 1+` is limited by input representation, not by algebra.**
Against a 50-digit reference evaluated at the *same stored double*, the production
form holds ~`1e-16` relative error at every `R-1` from `1e-2` to `1e-14`, while the
naive form degrades to `3.7e-14` and `6.7e-11`. Separately, storing `R = 1 + 1e-12`
as a double reproduces `R - 1` to only ~4 significant digits (relative error
`8.9e-5`), and since `dv ≈ (R-1)/2` the result inherits exactly that. **No
reformulation can recover information the input no longer carries**; this is a
parameterization limit, and it is documented rather than hidden.

**Finding 3 — the large-`B` tail becomes unresolvable in double precision.** The
residual `dv_bar_B - dv_bar_B_inf ≈ (sqrt(R)-3)/(sqrt(2)B)` eventually falls below
a few ulp of the value itself, after which the curve is flat to machine precision.
`B_RATIO_RESOLVABLE_MAX = 1e10` records the generic safe limit (at `R = 20` the
residual is still ~47,000 ulp there). **`R = 9` is the worst case**: the leading
coefficient `sqrt(R) - 3` vanishes identically, so the tail decays as `1/B²`
(measured log-ratio `2.000` per decade) and reaches the noise floor by `B/R ≈ 1e8`.

Consequently `scan_B_structure` takes a `noise_ulps` parameter (default 8) and
requires a turning point to clear the rounding noise floor before counting it.
With the floor disabled and a sweep to `B/R = 1e14`, `R = 9` reports **278 spurious
"minima" and 222 "maxima" that are pure rounding noise**; with the floor on, it
correctly reports zero. This was diagnosed during M2 and is the reason the
detector is noise-aware — it is a measurement limit, **not** a contradiction of the
M1 structural result.

**Finding 4 — extreme scales.** Across `mu ∈ [1e2, 1e20]` and `r1 ∈ [1e1, 1e12]`
the normalized and direct-dimensional paths agree to `9.8e-16` relative.

## M2.6 The `B = R` time degeneracy

**Delta-v and transfer time degenerate differently at `B = R`, and the code says so.**

`dv_bar_B(R,R) = dv_bar_H(R)` exactly. But ellipse 2 degenerates into the *circular
target orbit*, and the bi-elliptic path still coasts a half revolution of it before
reaching the nominal periapsis point. Hence

```
t_bar_B(R, R) = t_bar_H(R) + pi * R^(3/2)
```

— the Hohmann time **plus half the period of the final circular orbit**, verified
exactly at `R = 2, 12, 16, 50`. The extra term is loiter the Hohmann transfer never
flies. `t_bar_B` is continuous in `B` down to `B = R`; it simply does not converge
to `t_bar_H` there.

This is exposed as `bielliptic_time_excess_at_B_equals_R(R)` so the subtlety is a
testable quantity rather than a footnote, and a test asserts the inequality
`t_bar_B(R,R) > t_bar_H(R)` specifically to stop it being "fixed" into a false
equality later.

## M2.7 Structural verification: no finite interior minimum

Re-derived numerically, not assumed from M1. Across `R ∈ {2, 5, 8, 9, 10, 12, 14,
15, 16, 20, 50, 100}` plus a sweep of `R` from 1.5 to 200 in steps of 2.5, with
4001 logarithmically spaced `B` samples per `R`:

- **interior local minima found: 0** in every case;
- at most one interior stationary point, always a **maximum**, present only for
  `9 < R < R2*` (e.g. `B_max ≈ 66.45` at `R = 10`, `26.40` at `R = 12`);
- an independent brute-force argmin over 3001 samples lands at a domain **endpoint**
  in every case;
- the shape matches the M1 endpoint-slope table exactly: increasing for `R < 9`,
  up-then-down for `9 < R < R2*`, decreasing for `R > R2*`; the combination
  (negative near-slope, positive far-slope) never occurs.

Hence `inf_B dv_bar_B(R,B) = min(dv_bar_H(R), dv_bar_B_inf(R))`, exposed as
`bielliptic_infimum_normalized(R)` together with which endpoint attains it.

A regression test also re-confirms M1 Section 6.3 directly: wherever
`dv_bar_B_inf(R) >= dv_bar_H(R)`, no finite `B` beats Hohmann either — **a finite-`B`
win forces an infinite-`B` win**, so the "finite wins while infinite loses" region
is empty.

## M2.8 Verification results against the M1 plan

| M1 §9 item | M2 status |
|---|---|
| A Hohmann limits | `dv_bar_H(1) == 0.0` exactly; `(R-1)/2` leading order; both burns positive; closed form vs. direct vis-viva to `1.1e-16` |
| B Bi-elliptic limits | `B->R+` identity (`<= 1 ulp`, `dv3` exactly 0); `B->inf` limit; `(sqrt(R)-3)/(sqrt(2)B)` residual to `1e-4` relative; burn directions checked from actual speeds |
| C Scaling / invariance | `dv ∝ sqrt(mu/r1)` verified by explicit scaling law; `R1*` recomputed from dimensional transfers for 5 bodies agrees to `1e-11` |
| D Finite-`B` optimization | dense brute force confirms no interior minimum; infimum is an endpoint |
| E Crossover roots | two independent formulations per threshold; `B_crit` bracketed, endpoints checked (`B_crit -> ∞` as `R -> R1*+`, `-> R` as `R -> R2*-`) |
| F Continuity | refinement test: halving the grid step halves the largest jump; no NaN/overflow across 12 decades of `B`; the cancellation-prone term B2 covered by the stability tests |
| G Time of flight | half-period formulas vs. independent `2*pi*sqrt(a^3/mu)/2` construction; normalized/dimensional agreement to `1e-14` |
| H Dimensional round trip | normalize and recover, both directions, worst relative residual `4.7e-16` |
| I Literature check | performed only after independent derivation; consistent with classical `~11.94` and `~15.58`; **no literature number appears as a production constant** |

## M2.9 Test suite

**882 tests, all passing under `pytest -W error`** (no warnings).

| File | Focus |
|---|---|
| `tests/test_package.py` | scaffold, version, public API present, M3 API absent, no crossover constants in `constants.py` |
| `tests/test_hohmann.py` | verification A, plus M1 table regression |
| `tests/test_bielliptic.py` | verification B and C, complex-step derivative check, M1 table regression |
| `tests/test_dimensional.py` | verification G and H, scale invariance across five bodies |
| `tests/test_timing.py` | verification H, `B = R` degeneracy, `B^{3/2}` growth, M1 time-table regression |
| `tests/test_crossover.py` | verification D, E, F: thresholds, classifier, break-even, continuity |
| `tests/test_structure.py` | verification I: no interior minimum, shape table, resolution limit |
| `tests/test_stability.py` | the M2.5 conditioning findings |
| `tests/test_report.py` | report determinism and freedom from absolute paths |

The M1 file `tests/test_placeholder.py` was renamed (`git mv`) to
`tests/test_package.py` and its "no M2 API exists" assertions replaced by their M2
counterparts. That file's M1 content remains unchanged in git history.

## M2.10 Artifacts

| Path | Content |
|---|---|
| `scripts/m2_verification_report.py` | regenerates the report from production code |
| `results/m2_verification_report.txt` | 10-section numerical report; deterministic (byte-identical across runs, checked by test); relative paths only |
| `scripts/m2_diagnostic_figures.py` | regenerates both figures |
| `figures/m2_fig1_excess_vs_B.png` | `dv_bar_B - dv_bar_H` vs. `B/R`, log `x`, **symlog `y`** so no curve is clipped; Region A/B/C structure and the maximum-only turning points are directly visible; break-even crossings marked |
| `figures/m2_fig2_time_trade_R12.png` | delta-v saving and transfer time for `R = 12`, full unclipped range |

Both figures show units and normalization, use no zero suppression, and represent
`B -> infinity` only as a labelled asymptote — never as a fake finite point.

## M2.11 Headline trade, from production code

For `R = 12` (Region B), Earth reference:

| Quantity | Value |
|---|---|
| entire theoretical prize (`B -> ∞`) | **3.0374 m/s** |
| worst case at moderate `B` (`B/R ≈ 2.2`) | **−39.4 m/s** (bi-elliptic *worse*) |
| break-even `B_crit` | 815.820250 |
| break-even apoapsis radius | 5.4482e+06 km (**14.17× lunar distance**) |
| Hohmann transfer time | 0.520859 d |
| time merely to break even | **524.09 d** |
| time penalty factor | **1006.20×** |

At break-even the saving is zero by definition; the full 3 m/s requires infinite
time. This is why the deliverable must be an engineering trade rather than a
minimization — but the recommendation itself is **M6, not M2**.

## M2.12 Scope guard

M2 did **not** implement, and must not be read as implementing:

- the final crossover radius-ratio portfolio plot (M4)
- the final mission/engineering recommendation (M6)
- sensitivity or large parameter sweeps (M3)
- perturbations of any kind (no J2, no third body, no drag)
- finite-burn losses
- plane-change coupling
- radiation or environment models
- a mission-specific optimizer or operational decision tool

M2 is equation and solver verification only. The M1 scope and limitations in
Section 12 continue to apply in full, unchanged.

## M2.13 Changes to M1 material

**None.** No M1 headline value was altered; every M1 number in Sections 1–13 above
was independently reproduced by the production code as a pre-flight check and again
as regression tests. The only M1-era file modified is `tests/test_placeholder.py`,
renamed and updated as described in M2.9, because its explicit purpose was to
assert that the M2 API did not yet exist.

One clarification, not a correction: M1 Section 6.1 described the shape scan; M2
adds that the scan needs a rounding-noise floor to remain meaningful past
`B/R ~ 1e10`, and that `R = 9` is the special case where the tail decays
quadratically. The structural conclusion is unchanged.

---
---

# DESIGN — Milestone 3

**Global crossover map, break-even locus, delta-v/time trade, and the engineering
recommendation.**

Status: **M3 complete.** The M1 and M2 sections above are unchanged and no
accepted headline value was altered. M3 composes the verified M2 primitives into
the principal trade study; it contains no transfer equations of its own.

## M3.1 Production module

`src/bielliptic_crossover/trade.py`. Every delta-v and every transfer time is
obtained by calling M2 production functions — a test asserts that
`practical_trade_metrics` reproduces direct M2 calls exactly, so the trade layer
can only delegate, never re-derive.

| API | Purpose |
|---|---|
| `TradePoint`, `RadiusRatioTrade`, `EndpointInfimum`, `BreakEvenLocusPoint`, `Recommendation` | result records |
| `minimum_endpoint_transfer(R)` | mathematical best-possible delta-v and which endpoint attains it |
| `practical_trade_metrics(R, B)` | full delta-v / time / geometry metrics for one finite design |
| `break_even_trade(R)` | metrics exactly at `B_crit`, or `None` outside Region B |
| `evaluate_trade(R, B_values)` | the whole finite-`B` family at one `R` |
| `authoritative_R_grid()`, `break_even_locus_grid()` | the documented sweep grids |
| `break_even_locus()` | `B_crit(R)` across Region B, solved at every point |
| `sweep_radius_ratios()` | endpoint-infimum envelope |
| `model_validity(rb_km)` | two-body model-validity flag and note |
| `recommendation_for(R)`, `recommendation_summary()` | scoped guidance, built from recomputed thresholds |

Three ideas are kept rigorously separate throughout: **mathematical delta-v
dominance**, **finite-`B` practical saving**, and **transfer-time / model-validity
penalty**. The `B -> infinity` branch is called the *mathematical infimum*
everywhere; `EndpointInfimum.attainable` is `False` whenever it is the best
endpoint, because it is not a realizable transfer.

## M3.2 Authoritative sweep domain

`authoritative_R_grid()` returns **1090 points** on `R ∈ [1.01, 100]`:

- 600 logarithmically spaced points across the full domain;
- 241 linearly spaced points in a window of half-width 0.5 about **each**
  threshold (so both are locally resolved to ~4e-3 in `R`);
- the thresholds themselves plus offsets at `±1e-9`, `+1e-6`, `+1e-3` relative;
- deduplicated and sorted.

The domain resolves both thresholds, contains ordinary transfers (GEO sits at
`R = 6.3137`), and reaches `R = 100` where the asymptotic trend is clear.

**The grid never defines a threshold.** `R1*` and `R2*` always come from the M2
root solvers; a test re-runs the sweep at three very different resolutions and
confirms the thresholds are bit-identical.

Finite apoapsis family: `B/R ∈ {1.25, 1.5, 2, 3, 5, 10, 20, 50, 100}`, plus
`B_crit`, `1.1 B_crit` and `2 B_crit` where Region B defines them. These are
**representative bounded designs, never optima** — M1/M2 established that no
finite interior optimum exists.

## M3.3 Break-even locus `B_crit(R)`

`break_even_locus()` solves the root at each of **219** points spanning Region B,
parameterised as `R = R1* + (R2* - R1*)·s` with `s` logarithmically spaced so both
ends are resolved. Nothing is interpolated: a test asserts every locus point
satisfies `dv_bar_B(R, B_crit) - dv_bar_H(R) = 0` to better than `1e-15`, and
cross-checks a sample against an independent bisection.

| `R` | `B_crit` | `B_crit/R` | `r_b,crit` altitude | `t_B/t_H` at break-even |
|---|---|---|---|---|
| 12 | `815.8202504753092` | 67.985 | 5.442e+06 km | 1006.20 |
| 12.5 | `90.75094420913885` | 7.260 | 6.001e+05 km | 38.87 |
| 13 | `48.90484332838889` | 3.762 | 3.202e+05 km | 16.03 |
| 14 | `26.10461128235042` | 1.865 | 1.680e+05 km | 8.02 |
| 15 | `18.19028151222189` | 1.213 | 1.151e+05 km | 4.30 |
| 15.5 | `15.896871011509875` | 1.026 | 9.978e+04 km | 3.63 |

Verified asymptotics, both required and both confirmed numerically:

- **`B_crit -> infinity` as `R -> R1*+`.** At `R - R1* = 3.6e-5`, `B_crit = 1.37e6`.
  The measured law is `B_crit/R ∝ (R - R1*)^-1` (fitted log-log slope `-1.00005`).
- **`B_crit -> R` as `R -> R2*-`.** At `R2* - R = 3.6e-6`, `B_crit/R = 1.0000011`.
  The measured law is `B_crit/R - 1 ≈ 0.306·(R2* - R)` (fitted slope `1.0004`).
- `B_crit` and `B_crit/R` both **decrease strictly** across Region B (asserted over
  all 219 points).

Infinity is never faked: the locus simply runs off the top of the figure axes,
and Region C is drawn with no `B_crit` at all rather than a placeholder.

A consequence of the `-3/2` time law (`t ∝ B^{3/2}`): the break-even **time**
penalty diverges as `t_B/t_H ∝ (R - R1*)^{-3/2}` (fitted slope `-1.49983`).

## M3.4 Endpoint-optimal envelope

Because no finite interior minimum exists,

```
dv_best(R) = min( dv_bar_H(R), dv_bar_B_inf(R) )
saving_bar(R) = dv_bar_H(R) - dv_best(R)
```

- **Region A: `saving_bar` is exactly `0.0`** (not merely small) — asserted
  bit-exactly, since `dv_best` returns the Hohmann value itself.
- Regions B and C: `saving_bar > 0`, but `attainable = False` — it requires
  `B -> infinity` and therefore unbounded transfer time.

Peak achievable saving over the domain is ≈ **316 m/s near `R ≈ 50`**, falling
again toward `R = 100` (≈ 289 m/s), because both curves converge to `sqrt(2) - 1`.

## M3.5 Earth dimensional examples

`r1 = 6678.1363 km`, `v1 = 7.72576063698292 km/s`. Final altitude is a conversion
from the universal radius-ratio result: `h2 = R·r1 - R_Earth`.

| case | `R` | region | `h2` [km] | `dv_H` [km/s] | `B→∞` infimum [km/s] | infimum saving | `B=5R` saving | `t_H` [d] | `t(5R)` [d] |
|---|---|---|---|---|---|---|---|---|---|
| GEO | 6.3137 | A | 35785.9 | 3.8926 | 4.4737 | **0.00 m/s** | **−486.4 m/s** | 0.220 | 4.66 |
| — | 10 | A | 60403.2 | 4.0930 | 4.2121 | 0.00 m/s | −124.8 m/s | 0.405 | 9.21 |
| — | 12 | B | 73759.5 | 4.1269 | 4.1239 | +3.04 m/s | **−27.5 m/s** | 0.521 | 12.08 |
| — | 13 | B | 80437.6 | 4.1355 | 4.0877 | +47.87 m/s | +8.54 m/s | 0.582 | 13.61 |
| — | 15 | B | 93793.9 | 4.1427 | 4.0264 | +116.31 m/s | +63.89 m/s | 0.711 | 16.85 |
| — | 16 | C | 100472.0 | 4.1429 | 4.0001 | +142.71 m/s | +85.39 m/s | 0.779 | 18.55 |
| — | 20 | C | 127184.6 | 4.1312 | 3.9157 | +215.52 m/s | +145.28 m/s | 1.069 | 25.89 |
| — | 50 | C | 327528.7 | 3.9687 | 3.6527 | +316.01 m/s | +233.94 m/s | 4.047 | 101.93 |

Two rows carry the whole engineering message. **GEO** is firmly Region A: the
infimum saving is exactly zero and the representative bi-elliptic *costs 486 m/s
more*. **`R = 12`** is in Region B yet `B = 5R` still costs **27.5 m/s more** than
Hohmann, because `5R = 60` is far below `B_crit = 815.8`.

These enormous target radii are dimensional illustrations of a normalized result,
not proposed Earth mission designs.

## M3.6 Threshold neighbourhoods

**Around `R1*` = 11.938765472645892.** No finite win exists below it (`break_even_B`
returns "no bi-elliptic win" at `R = 11.5, 11.9`). Immediately above, `B_crit`
explodes and the prize is negligible:

| `R` | infimum saving | `B_crit` | `r_b,crit` | `t_B/t_H` at break-even |
|---|---|---|---|---|
| `R1*(1+1e-6)` | 0.0006 m/s | 4.17e+06 | 2.79e+10 km | ~1e10 |
| `R1*(1+1e-4)` | 0.060 m/s | 4.17e+04 | 2.79e+08 km | ~1e7 |
| `R1*(1+1e-2)` | 5.89 m/s | 4.19e+02 | 2.80e+06 km | ~3e3 |
| 12.0 | 3.04 m/s | 815.82 | 5.45e+06 km | 1006.20 |
| 13.0 | 47.87 m/s | 48.90 | 3.26e+05 km | 16.03 |

**Around `R2*` = 15.581718738763179.** `B_crit -> R`, and above it every `B > R`
wins — but the win can be worth almost nothing:

| `R` | region | `B_crit/R` | saving at `B = 1.01R` | `t_B/t_H` at `B = 1.01R` |
|---|---|---|---|---|
| 15.0 | B | 1.2127 | **−0.154 m/s** | 3.601 |
| 15.5 | B | 1.0256 | **−0.013 m/s** | 3.609 |
| `R2*` | boundary | → 1 | +0.009 m/s | 3.610 |
| 16.0 | C | — | +0.117 m/s | 3.616 |
| 17.0 | C | — | +0.352 m/s | 3.630 |

Just above `R2*`, a `B = 1.001R` transfer saves **+0.0005 m/s** while taking
**3.58×** as long. "Every bi-elliptic beats Hohmann" is a statement about sign,
not about magnitude.

## M3.7 Model-validity warnings

`model_validity(rb_km)` classifies each design:

| flag | condition |
|---|---|
| `two_body_reasonable` | `rb < 0.1 ×` lunar distance |
| `lunar_third_body_caution` | `rb ≥ 0.1 ×` lunar distance |
| `two_body_invalid` | `rb ≥` lunar distance (384 400 km), or `≥` Earth Hill radius (~1.5e6 km) |

The flag appears in every result table and on both trade figures.

**A headline finding.** Solving `r_b,crit = lunar distance` gives

```
R = 12.834879147      (B_crit = 57.5610)
```

**Below `R ≈ 12.835` — the lower 25 % of Region B — the break-even apoapsis lies
beyond the Moon.** At `R = 12` it is `5.45e6 km`, **14.17 lunar distances**, and
`2 B_crit` reaches 28.3 lunar distances, well past the Earth Hill radius. In that
range the Earth-only two-body model is not merely inaccurate, it is the wrong
model, and the quoted numbers are mathematical extrapolations of an idealized
model rather than trajectory designs. This project makes no claim otherwise.

## M3.8 Practical trade at the two headline radius ratios

**`R = 12` (Region B).** Entire mathematical prize `+3.04 m/s`, and only at
`B -> infinity`:

*(Values below corrected during the M4 audit — see M4.2. The production code,
CSVs and figures were always correct; three cells in this hand-typed table were
not.)*

| strategy | `B` | saving | `t_B/t_H` | `r_b` / lunar | model |
|---|---|---|---|---|---|
| `B/R=1.25` | 15.00 | −22.30 m/s | 4.36 | 0.26 | caution |
| `B/R=2` | 24.00 | **−39.02 m/s** | 7.28 | 0.42 | caution |
| `B/R=10` | 120.00 | −14.99 m/s | 60.75 | 2.08 | invalid |
| `B_crit` | 815.82 | **0.00 m/s** | **1006.20** | 14.17 | invalid |
| `1.1 B_crit` | 897.40 | +0.27 m/s | 1159.59 | 15.59 | invalid |
| `2 B_crit` | 1631.64 | +1.50 m/s | **2829.07** | 28.35 | invalid |

Even at twice the break-even apoapsis the saving is **half** the theoretical
prize, for a **2829×** time penalty at 28 lunar distances.

**`R = 16` (Region C).** Every `B > R` wins, and the benefit appears far more
readily — but still trades time, and the useful designs are already past the
lunar distance:

*(Values below corrected during the M4 audit — see M4.2.)*

| strategy | `B` | saving | `t_B/t_H` | `r_b` / lunar | model |
|---|---|---|---|---|---|
| `B/R=1.25` | 20.00 | +5.97 m/s | 4.45 | 0.35 | caution |
| `B/R=2` | 32.00 | +31.87 m/s | 7.45 | 0.56 | caution |
| `B/R=5` | 80.00 | +85.39 m/s | 23.82 | 1.39 | invalid |
| `B/R=100` | 1600.00 | +139.33 m/s | 1840.74 | 27.80 | invalid |
| `B -> ∞` | — | +142.71 m/s (infimum) | ∞ | ∞ | not a transfer |

## M3.9 Engineering recommendation

Produced by `recommendation_for(R)` and `recommendation_summary()` from the
recomputed thresholds. Three concepts are kept distinct: **mathematical delta-v
dominance**, **finite-transfer practical saving**, **time and model-validity cost**.

- **`R < 11.9388` (Region A).** Hohmann is delta-v superior to **every** admissible
  bi-elliptic transfer. No intermediate apoapsis helps; its longer transfer time
  buys nothing.
- **`11.9388 < R < 15.5817` (Region B).** A bi-elliptic transfer beats Hohmann
  **only if `r_b` exceeds `B_crit(R)·r_1`**. Near the lower threshold the required
  apoapsis and transfer time are absurd for a vanishing saving — at `R = 12`,
  3 m/s in exchange for a 1006× time penalty and an apoapsis 14 lunar distances
  out. Below `R ≈ 12.835` the break-even apoapsis is beyond the Moon, so the model
  does not apply there at all.
- **`R > 15.5817` (Region C).** Every `B > R` is delta-v better than Hohmann, but
  the **practical size** of the saving must still be weighed against transfer time.
  Just above `R2*` a marginally bi-elliptic transfer saves under 0.01 m/s for
  3.6× the duration.
- **`B -> infinity` is a mathematical lower bound, never an engineering
  recommendation** and never a realizable transfer.
- **For ordinary Earth transfers such as 300 km LEO to GEO (`R = 6.3137`),
  Hohmann remains firmly preferred** within this idealized coplanar impulsive
  two-body model.

Three claims this project explicitly does **not** make:

1. that the bi-elliptic transfer simply becomes the better choice everywhere
   above the lower threshold — in Region B it is worse for every apoapsis below
   `B_crit(R)`;
2. that the *upper* threshold is where the bi-elliptic family first becomes
   competitive — that is the *lower* threshold's role; the upper one is where
   **every** admissible apoapsis wins;
3. that `B -> infinity` is an optimum — it is a mathematical infimum requiring
   unbounded transfer time and is never a realizable transfer.

`tests/test_m4_docs.py` asserts that neither this document nor the README ever
states any of them, and that the generated recommendation text avoids them too.

## M3.10 Numerical verification and boundary findings

**Independent cross-check.** Selected points were recomputed with a scratch
mpmath vis-viva implementation that does not import `trade.py`: `R=12` at `B_crit`
and `2 B_crit`, `R=16` at `B=1.1R`, `R=20` and `R=50` at `B=2R`, and the GEO
classification. Delta-v, time, `r_b`, sign and region all agree; **worst relative
residual `5.43e-13`**. (A sign comparison exactly at `B = B_crit` is undefined,
since the saving is zero there by construction.)

**Two documented precision boundaries**, both resolution limits rather than
errors, both now covered by regression tests:

1. **Near `R2*`.** `break_even_B` steps off the trivial `B = R` root with a fixed
   relative probe of `1e-6`. Since `B_crit/R - 1 ≈ 0.306·(R2* - R)` vanishes
   linearly, it falls below that probe once `R2* - R < 3.27e-6`, and the solver
   then reports Region C while the threshold-based classifier still reports
   Region B. Inside that band the break-even apoapsis is under one part in a
   million above the target radius, so nothing physically meaningful is lost.
2. **At `R1*` exactly.** The double-rounded `threshold_R1()` sits `2.16e-14` above
   the true root, so at that exact `R` the true difference
   `dv_H - dv_B_inf = 1.39e-16` is positive (verified at 50 digits) and a finite
   break-even genuinely exists — about `4e15`, i.e. `6.9e13` lunar distances. The
   classifier's tolerance band assigns Region A, which is the sensible engineering
   reading. Both are defensible; the value is a meaningless extrapolation.

Neither boundary required any change to M1 or M2 code. **No error was found in
the accepted M1/M2 results.**

## M3.11 Convergence and determinism

- Thresholds are **bit-identical** across sweep resolutions of 80, 600 and 1500
  global points.
- All five result artifacts regenerate **byte-identically** (verified by repeated
  runs and by tests that rebuild them and compare against the committed files).
- CSV floats are written with `repr`, so every value round-trips exactly; a test
  asserts `repr(float(value)) == value` for the sampled rows.
- JSON uses `sort_keys=True` and contains no timestamps or machine paths.

## M3.12 Artifacts and figures

| Path | Content |
|---|---|
| `results/m3_radius_trade.csv` | endpoint-infimum envelope, 1090 rows |
| `results/m3_finite_b_strategies.csv` | finite-`B` metrics, 162 rows over 16 radius ratios |
| `results/m3_break_even_curve.csv` | `B_crit(R)` locus, 219 solved roots |
| `results/m3_earth_examples.csv` | 8 Earth dimensional cases including GEO |
| `results/m3_summary.json` | thresholds, Earth reference, `B_crit` examples, GEO, R12/R16 trades, recommendation |
| `figures/m3_crossover_map.png` | **primary figure**: delta-v curves with inset zoom on the crossing, best-possible saving, and the `B_crit/R` locus, all over shaded Regions A/B/C |
| `figures/m3_dv_time_trade.png` | delta-v saving vs. transfer time at `R = 12` and `R = 16`, with finite-`B` markers, `B_crit`, and the lunar-distance validity boundary |
| `figures/m3_break_even_time_penalty.png` | Region-B divergence of `B_crit/R`, `t_B/t_H` and `r_b,crit`/lunar distance, with the fitted power laws |

Every figure was visually inspected. Three defects were found and fixed during
review: the Region B/C labels were occluded by the legend in the primary figure
(replaced with a dedicated region-key legend); the infimum curve was clipped at
the top of panel 1 (y-limit extended to include its `0.826` value at `R = 1.01`);
and the `R = 16` panel of the trade figure lacked the lunar-distance warning
(added to both panels, which matters more at `R = 16` where every `B/R ≥ 5`
design is already beyond the Moon). `B -> infinity` appears only as a labelled
asymptote or an arrow direction, never as a finite plotted point.

## M3.13 Test suite

**1026 tests, all passing under `pytest -W error`.** New files:
`tests/test_trade.py` (118) and `tests/test_m3_artifacts.py` (16), covering the
required items A–O: sweep-resolution
independence, locus points being solved roots rather than interpolants,
monotonicity, both asymptotic laws, envelope equals endpoint minimum, no finite
interior optimum introduced, delegation to the M2 API, exact dimensional scaling,
GEO in Region A, time ratio always above one, `B^{3/2}` growth, the `R = 12`
penalty reproducing M2, recommendation text generated from production code, and
artifact determinism.

## M3.14 Scope guard

M3 did **not** implement, and must not be read as implementing: J2 or any gravity
harmonic; lunar or solar third-body propagation; Lambert targeting; finite burns;
plane changes; low thrust; radiation or environment models; launch or operational
design; a mission-specific optimizer; or CI/release packaging. This remains an
idealized coplanar two-body impulsive trade study. The M1 Section 12 scope and
limitations continue to apply in full.

Remaining for later milestones: final portfolio packaging and polish (M4+).

## M3.15 Changes to M1/M2 material

**No numerical change whatsoever.** Every M1/M2 source file that computes
anything -- `constants.py`, `_stable.py`, `hohmann.py`, `bielliptic.py`,
`timing.py`, `dimensional.py`, `crossover.py` -- is **byte-identical** to the M2
commit (`git diff` against `410681d` over those files is empty).

`src/bielliptic_crossover/trade.py` is purely additive and is imported explicitly
as `bielliptic_crossover.trade`, deliberately leaving the M2 public surface in
`__init__.py` exactly as verified.

The only edits to M1/M2-era files are non-numerical: the version string
(`0.2.0 -> 0.3.0` in `pyproject.toml` and `__init__.py`, following the M1->M2
precedent), the corresponding assertion in `tests/test_package.py`, and that
file's scaffold list extended with the new M3 paths. Its obsolete
`test_no_production_modules_yet` check -- written at M1 when the package was
empty -- was removed, since the M2 modules it forbade have existed since the
previous milestone and the surrounding tests already assert the current layout.

---
---

# DESIGN — Milestone 4

**Final audit, portfolio polish, reproducibility, CI and release.**

Status: **complete — released as v1.0.0.** M4 changed **no physics**. It audited
the repository, corrected documentation transcription errors, fixed one packaging
bug it introduced, and added licensing, CI and reproducibility guarantees.

## M4.1 What was audited

| Area | Method | Result |
|---|---|---|
| Accepted results | independent scratch mpmath vis-viva implementation importing no project physics | all reproduced, worst relative residual `3.21e-13` |
| Physics recheck A–O | separate scratch implementation, 15 categories | all passed, worst `2.82e-14` (excluding the finite-`B` asymptotic check at `4.06e-7`) |
| Static analysis | `pyflakes`, `ruff` | clean after fixes; lint config now pinned |
| Documentation numbers | every prose table cell compared against production output | **5 erroneous cells found and corrected** (M4.2) |
| Figures | inspected at full resolution and at ~880 px README embed width | one genuine readability issue fixed (M4.3) |
| Dependencies | imports cross-checked against `pyproject.toml` | correctly partitioned |
| Fresh environment | new venv, editable install, full regeneration | **one packaging bug found and fixed** (M4.4) |
| Artifacts | absolute paths, determinism, ordering, duplication | clean; all numerical artifacts byte-identical across environments |
| Claims/wording | forbidden-phrase scan, now enforced by tests | clean |

## M4.2 Genuine error found: documentation drift

**The only substantive defect the audit found.**

The prose tables in M3.8 were hand-typed rather than generated. Comparing every
cell against production output revealed **five wrong cells**:

| location | documented | actual |
|---|---|---|
| `R = 12`, `B/R = 10` | −15.02 m/s, `t_B/t_H` = 60.0 | **−14.99 m/s, 60.75** |
| `R = 16`, `B/R = 1.25` | +5.87 m/s, 3.86 | **+5.97 m/s, 4.45** |
| `R = 16`, `B/R = 2` | +31.79 m/s, 6.86 | **+31.87 m/s, 7.45** |
| `R = 16`, `B/R = 5` | 23.8 | **23.82** |
| `R = 16`, `B/R = 100` | +139.85 m/s, 2005 | **+139.33 m/s, 1840.74** |

**Cause.** Transcription, not computation. The production code, the CSV/JSON
artifacts and all figures were correct throughout — verified independently. Only
the hand-typed markdown was wrong, and the `R = 16` row values appear to have come
from an intermediate exploratory run rather than the final production output.

**Impact.** Documentation only. No threshold, no `B_crit`, no accepted headline
value and no figure changed. The qualitative conclusions are unaffected: `R = 16`
still shows a benefit appearing far more readily than at `R = 12`, and the time
penalty at `B/R = 2` is 7.45× rather than 6.86× — slightly *worse* for the
bi-elliptic case than previously written, so no conclusion was overstated in the
bi-elliptic's favour by the error.

**Fix.**
1. `scripts/m4_final_summary.py` now **generates** the headline tables from
   production code into `results/m4_final_summary.md`, which the README embeds
   verbatim.
2. The M3.8 tables were corrected in place and annotated.
3. `tests/test_m4_docs.py` asserts every documented strategy row still matches
   production, that the README embeds the generated table verbatim, and that the
   documented thresholds equal `repr(threshold_R1())` / `repr(threshold_R2())`.

The lesson is recorded rather than buried: **numbers that appear in prose should
be generated, not typed.**

## M4.3 Figure audit

Every figure was inspected at full resolution and at ~880 px, the width GitHub
renders README images at.

Final hierarchy:

| Tier | Figure | Role |
|---|---|---|
| **Primary 1** | `figures/m3_crossover_map.png` | the project in one glance: Δv curves with inset zoom on the crossing, best-possible saving, `B_crit/R` locus, over shaded Regions A/B/C |
| **Primary 2** | `figures/m3_dv_time_trade.png` | the practical story: Δv saving vs. transfer time at `R = 12` and `R = 16` |
| Supporting 3 | `figures/m3_break_even_time_penalty.png` | Region-B divergence with fitted power laws |
| Supporting 4 | `figures/m2_fig1_excess_vs_B.png` | structural verification: no interior minimum |
| Diagnostic 5 | `figures/m2_fig2_time_trade_R12.png` | M2-era Δv/time diagnostic, superseded by Primary 2 |

**One figure was changed**, `m3_dv_time_trade.png`: at embed width its smallest
annotations were only marginally legible. Font sizes were raised (8.5 → 10–11 pt)
and the canvas slightly narrowed. Enlarging the type then caused the
lunar-distance note to collide with the legend, so that note was moved *into* the
legend, where it cannot collide at any scale. **The underlying numerical arrays
are unchanged** — a SHA-256 digest of the plotted time/saving arrays is identical
before and after. No other figure was regenerated.

## M4.4 Packaging bug found and fixed

Adding `license-files = ["LICENSE"]` to `pyproject.toml` switched setuptools into
PEP 639 mode, under which `project.license` must be an SPDX **string**, not the
legacy `{ text = "MIT" }` table. The editable install then failed outright:

```
configuration error: `project.license` must be string
```

Caught immediately by the fresh-environment check, which is precisely its purpose.
Fixed by moving to `license = "MIT"` and requiring `setuptools>=77`. The pre-M4
form was verified to still install (deprecated but accepted), so **this bug was
introduced during M4 and fixed within it** — it never affected M1–M3.

## M4.5 Hygiene fixes

All non-numerical, all proven to cause zero drift (result artifacts re-hashed
before and after):

- removed two genuinely unused imports in `scripts/m3_figures.py`;
- `typing.Sequence` → `collections.abc.Sequence` in `trade.py` (annotation only);
- `zip(..., strict=True)` on equal-length iterations, and `itertools.pairwise`
  for the successive-pair idioms where `strict=True` would have been *wrong*;
- `chmod +x` on the four runnable scripts, which carried shebangs but were not
  executable;
- pinned `[tool.ruff]` configuration so lint results are reproducible rather than
  dependent on the reviewer's ruff version.

Deliberately **not** changed: `__all__` ordering (grouped by module, which is more
readable than alphabetical for this package) and the multi-line message strings
ruff flags as `ISC004` (they are intentional, not missing commas).

**Every M1/M2 numerical source file except `crossover.py` — `constants.py`,
`_stable.py`, `hohmann.py`, `bielliptic.py`, `timing.py`, `dimensional.py` — is
byte-identical to the M2 commit `410681d`.** `trade.py` changed only by the import
move above. `crossover.py` changed only in its two polynomial root-finders, for
the CI-proven reproducibility defect documented in M4.5b; every transfer equation
in it is untouched.

## M4.5b Second genuine error: polynomial roots were not platform-reproducible

**Found by CI, after the M4 commit was first pushed.** Both GitHub Actions jobs
failed on `test_committed_report_matches_regeneration`: the committed M2
verification report did not match the report regenerated on Linux.

The cause was `numpy.roots`, used by `threshold_R1_from_polynomial` and
`threshold_R2_from_cubic`. It solves for polynomial roots as the eigenvalues of a
companion matrix, delegating to LAPACK — and LAPACK's last bits differ between
BLAS builds. The polynomial threshold therefore came out as

```
macOS / Accelerate : R1* = 11.9387654726458692
Linux / OpenBLAS   : R1* = 11.9387654726458745
```

Both are correct to ~1e-15, so **no physical conclusion was ever affected** —
every threshold quoted throughout the project comes from the `brentq`
transfer-equation path, which was always deterministic. But the *generated
artifact* embedded the platform-dependent digits, which silently broke the
byte-reproducibility claim.

**Fix.** Both polynomials are now solved by bracketed root-finding on the
polynomial itself, over brackets justified from the known root structure
(`u ∈ [1, 10]`; `R ∈ [0, 100]`). This uses only IEEE-754 arithmetic on our own
function, so it is reproducible on any platform, and it removed a dependency on
the linear-algebra stack from a scalar root-find. The two solver paths remain
genuinely independent — polynomial versus transfer equation — and no crossover
constant is hardcoded: a bracket is a search interval, exactly as in
`threshold_R1()`.

It is also **more accurate**. Against a 40-digit reference:

| | `numpy.roots` | bracketed | change |
|---|---|---|---|
| `R1*` residual | `1.53e-15` (macOS), `3.78e-15` (Linux) | `2.02e-15` | comparable, now deterministic |
| `R2*` residual | `8.70e-15` | **`1.80e-16`** | **48× better** |

As a result the two independent `R2*` paths now agree **exactly** (difference
`0.000e+00`, previously `8.882e-15`). The M2.3 table above was updated
accordingly; this is the only M1/M2 numerical source change in the entire
project, and it is a reproducibility fix, not a physics change.

**Lesson:** cross-platform CI earned its place immediately — this defect was
invisible on the development machine and would have shipped without it.

## M4.6 Reproducibility

A fresh virtual environment was created, the package installed with
`pip install -e ".[dev,figures]"`, the suite run, and every artifact regenerated.

Versions used: Python 3.14.5, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.1,
pytest 9.1.1 (the committed figures were rendered with matplotlib 3.10.9).

| Output | Result |
|---|---|
| Test suite | 1079 passed under `pytest -W error` |
| `results/*.csv`, `*.json`, `*.md`, `*.txt` | **8 of 8 byte-identical** to the committed files |
| `figures/*.png` | bytes differ under matplotlib 3.11.1 vs 3.10.9 |

The distinction matters and is not glossed over: **numerical reproducibility is
exact and enforced in CI**; **PNG byte-reproducibility is not claimed across
matplotlib versions**, only within one. The arrays behind the figures are
deterministic; the rasterisation is not portable. Committed figures were left as
rendered by matplotlib 3.10.9 rather than regenerated to chase byte equality.

CI enforces the numerical half directly: after regenerating the artifacts it runs
`git diff --exit-code -- results/`, so any drift fails the build.

## M4.7 Final authoritative hierarchy

1. **`results/m4_final_summary.md`** — generated headline table, single source of
   truth for numbers quoted in prose.
2. **M3 artifacts** (`results/m3_*.csv`, `results/m3_summary.json`) — the full
   trade study.
3. **M2 report** (`results/m2_verification_report.txt`) — equation verification.
4. **DESIGN.md M3** — authoritative narrative for final results.
5. **DESIGN.md M1/M2** — derivation and verification history, preserved unchanged.

## M4.8 Release state

- Version **1.0.0**, consistent between `pyproject.toml` and `__init__.py`
  (asserted by a test).
- **MIT LICENSE** added, author `Sanjana` taken from the git configuration.
- **CI** (`.github/workflows/tests.yml`) on push and pull request, Python 3.11 and
  3.12, running `pytest -W error` and the artifact-determinism check.
- Runtime dependencies `numpy`, `scipy`; `matplotlib` under the `figures` extra;
  `pytest`, `ruff`, `pyflakes` under `dev`.
- **1079 tests** passing under `pytest -W error` with no warnings.

## M4.9 Remaining weaknesses

Stated plainly rather than hidden:

1. **The model is idealized and that is the dominant limitation.** Every
   conclusion holds only for coplanar, two-body, impulsive transfers. The
   regime where the bi-elliptic transfer wins materially is exactly the regime
   where an Earth-only two-body model stops being appropriate — which is a
   finding, but it also means the practical recommendation rests on a model
   invalid at the apoapses in question.
2. **No finite-burn, plane-change or perturbation modelling.** Real
   high-apoapsis transfers would be dominated by third-body effects and finite-burn
   losses, both of which plausibly shift the crossover.
3. **Threshold behaviour near `R1*`/`R2*` is resolution-limited in double
   precision** (M3.10). The bands are ~2e-14 and ~3.3e-6 wide respectively and are
   documented and tested, but a reviewer wanting exact behaviour at the thresholds
   would need extended precision.
4. **Figure rendering is not byte-portable** across matplotlib versions, so CI
   verifies numerical artifacts only, not images.
5. **Single-author verification.** The independent checks are genuinely separate
   implementations, but written by the same author against the same understanding;
   a shared conceptual error would not be caught.
6. **No experimental or mission validation.** Comparison is against classical
   textbook values only, and only after independent derivation.

## M4.10 Scope guard

M4 added **no orbital-transfer physics**. No J2, no third-body propagation, no
Lambert targeting, no finite burns, no plane changes, no low thrust, no radiation
model, no mission optimizer. The M1 Section 12 limitations continue to apply in
full and unchanged.
