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
