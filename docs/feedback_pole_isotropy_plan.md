# Feedback pole isotropy — analysis & plan

**Status: PLAN — not implemented.** The D1 tangent-warp variant was measured
and **rejected 2026-07-23**; of the D3 helpers only `pole_wrap` is in the tree.
Read "Measured outcome" before proposing any variant.

Fix two MeshFeedback pole artifacts — trails pinching toward the poles, and a
~4–5 px "lensing" cap under low-frequency noise — so the warp is isotropic over
the whole sphere, without giving back the 2026-07-20 GREEN device profile
(worst frame 54.45 ms, 8.05 ms margin; report no longer tracked).

## How the current path works

Per flush (`Filter::Pixel::Feedback` in `core/render/filter.h`):

1. `populate_warp_field` (~line 1960): for each coarse point (ds=4 grid),
   reconstruct the sphere point `v = pixel_to_vector(x, y)`, warp it on the
   sphere (`noise_warp` / `melt_warp` — both are tangent-plane slides, fully
   isotropic in 3D), then project the warped vector back to equirect and store
   **equirect pixel offsets** `(projected_x − x, projected_y − y)` as int16
   (`WARP_SCALE = 128`, ±255 px range).
2. `composite_previous_frame` (~line 2038): bilinearly upsample the offsets
   (with a 3-branch seam-unwrap per coarse cell anchored on `d00`), add to the
   pixel coordinate, and `sample_bilinear_prev` the previous frame there.
   Rows outside `[0, H)` contribute **black**.

The 3D warp itself is isotropic. Every artifact comes from the round trip
through equirect offset space and from equirect-blind sampling.

## Root causes

### R1 — degenerate pole rows in the warp grid (the lensing cap)

`pixel_to_vector(x, 0)` is the north pole `(0,1,0)` for **every** column
(`sin_phi[0] = 0`), so all coarse points of row 0 warp to the *same* distorted
vector, which projects to one fixed `(projected_x, projected_y)`. The stored
row-0 offsets therefore encode "every column samples from one longitude". All
pixel rows 0..3 interpolate against that degenerate row → the cap collapses
toward a point: exactly the 4–5 px lensing (ds=4 → one coarse cell). Host/sim
south pole (row H−1, `H_OFFSET = 0`) is the mirror case. On device
(`H_OFFSET = 3`) the bottom row stops short of the pole, so the south cap is a
milder version of R2/R3, not a true degeneracy.

### R2 — offsets stored in a quantity that blows up as 1/sin(φ) (the pinch)

The equirect x-offset of an isotropic angular displacement scales as
`1/sin(φ)`. The stored field is exact *at* coarse points, but:

- Bilinear interpolation in y linearly blends samples of a hyperbolically
  growing, rapidly rotating quantity — systematic overshoot between coarse rows
  near the poles, re-applied every frame by the feedback loop, accumulating
  into visible azimuthal drag/pinch.
- Near the poles a modest tangent slide exceeds the (vanishing) sphere distance
  between adjacent columns, so whole near-pole bands legitimately converge to
  common longitudes in offset space; interpolating and quantizing that
  representation (int16 clamp at ±255 px is only ~5 coarse cells of azimuth
  wrap at row 1) distorts it further.
- With higher-frequency noise the offsets decorrelate faster between coarse
  rows, making the y-interpolation error (and thus the pinch) more visible.

### R3 — sampling falls off the sphere at the poles

`sample_bilinear_prev` treats rows `< 0` or `≥ H` as black. Spherically, a
sample past the north pole should continue on the far side: `(x, −y)` →
`(x + W/2 mod W, +y)`. Black taps drain trail energy at the cap every frame and
forbid any content from flowing across the pole. Related: the row-0 offsets can
only ever have `y_offset ≥ 0` (projected φ ≥ 0), so the top band can only pull
content *up* toward the pole — an asymmetric bias feeding the pinch.

### Why previous attempts failed (do not revisit)

- Smootherstep on the coarse grid (reverted 57c0b498): reshaped interpolation
  inside the wrong representation — washboard ripples, worse than C0.
- Fibonacci-lattice sphere field (reverted bbec6969): abandoned the equirect
  grid entirely; traded pole caps for lattice interpolation/moiré issues.

Both treated the *interpolator*; the defect is the *interpolated quantity*.

## Proposed design

### D1 — store the warp as tangent-space angular offsets (fixes R1 + R2)

In `populate_warp_field`, replace the equirect projection with the log map at
`v`: with `c = dot(v, v′)`, `u = v′ − c·v`, the tangent displacement is
`δ = u · (acos(c)/|u|)` (guard `|u|² ≈ 0` → δ = 0). Decompose against the local
frame, both free from the TrigLUT values already fetched for
`pixel_to_vector`:

- `east(x,y)  = (−sinθ, 0, cosθ)`
- `north(x,y) = (cosφ·cosθ, −sinφ, cosφ·sinθ)`

Store `δe = dot(δ, east)`, `δn = dot(δ, north)` as int16 at an angle scale
(`32767/π` covers ±π; `melt_warp` max ≈ 1.14 rad). Both components are smooth,
bounded, seam-free, and **non-degenerate at the poles**: at row 0 the frame is
per-column well-defined (limit along the column's meridian), so the pole row
stores a real per-column displacement — R1 vanishes at the source, and δn is
signed, so warps can pull content *across* the pole (paired with D2).

In `composite_previous_frame`, the inner loop is unchanged except for the final
scale: `ddy = lerp(δn) · Sy` with constant `Sy = (H_VIRT−1)/π / SCALE` (same
shape as today), and `ddx = lerp(δe) · Sx[y]` where
`Sx[y] = (W/2π)/sin_phi[y] / SCALE` is hoisted **per row** (one divide per row,
or a small per-flush table). The 1/sin(φ) factor — the thing that cannot be
interpolated — is now applied *exactly, per pixel*, while interpolation happens
only on smooth quantities. Cost: **+1 multiply per pixel**; in exchange the
entire per-cell `WRAP_PERIOD` 3-branch unwrap block (filter.h ~2059–2129)
is deleted — tangent components never wrap, including across the x=0 seam.

Pole rows (`sin φ = 0`: row 0 always; row H−1 host-only): clamp `Sx` to the
adjacent row's value. The approximation is confined to a single row that
rasterizes as a full-width bar for a single physical point anyway
(DistortedRing pole-fill precedent).

### D2 — pole-crossing bilinear taps (fixes R3)

In `sample_bilinear_prev`, replace the black fallback for out-of-range rows
with spherical reflection, applied per tap only in the (rare) out-of-range
branch — the in-range hot path is untouched:

- north: `row < 0` → `row = −row`, `col = fast_wrap(col + W/2, W)`
- south: reflect about the true pole row `S = H_VIRT − 1`; `row′ = 2S − row`;
  if `row′ ≥ H` (device gap, `H_OFFSET = 3`) the tap stays black — correct,
  nothing is ever rendered there. On host (`H_OFFSET = 0`, S = H−1) reflection
  is live.

This makes the sampler continuous across the north pole and consistent between
host and device where data exists.

### D3 — the generic abstraction

Small, opt-in helpers in `core/math/geometry.h` (no global behavior change —
existing pole contract tests elsewhere stay untouched):

- `tangent_frame<W,H>(x, y)` → east/north from TrigLUT.
- `sphere_log(v, v′)` → `(δe, δn)` tangent components (the D1 encoder).
- `equirect_x_scale<W,H>(y)` → the per-row azimuth scale (D1 decoder side).
- `pole_wrap<W,H>(col, row)` → the D2 tap-reflection rule.

Only `pole_wrap` is in the tree (`core/math/spherical_field.h` uses it); the
three tangent primitives are part of the rejected D1 variant and would have to
be written alongside it. Any future filter that resamples an equirect history
buffer or interpolates a sphere field over an equirect grid should reach for the
same primitives instead of re-deriving them.

Deliberately **not** doing: per-pixel exact reprojection (hoist every pixel to
cartesian, rotate, `atan2`+`acos` back). It is the exactness gold standard but
costs ~41k inverse projections per frame inside an 8 ms margin; the tangent
scheme is exact at coarse points, has bounded isotropic interpolation error
(O(cell²) in a smooth bounded field), and costs one multiply per pixel. It can
serve as the *test oracle* instead (see below).

## Execution plan

1. **Pin the artifacts (before touching code).** Native framebuffer-dump
   harness: two scenes — low-freq/high-amp noise (lensing cap), high-freq
   (pinch) — dumped over ~100 frames to PNG. These are the before/after
   evidence and the eyeball check for every later step.
2. **D1 in `populate_warp_field` + `composite_previous_frame`.** Same arrays,
   same `WarpKey`, same `STORAGE_BYTES`; only the stored quantity and the
   per-row scales change. Delete the seam-unwrap block.
3. **D2 in `sample_bilinear_prev`.**
4. **D3 extraction** once 2–3 are proven, moving the math into geometry.h.
5. **Tests** (all wired into tests/CMakeLists + CI runner):
   - Encoder/decoder round trip: at mid-latitude coarse points the new path
     reproduces the old exact offsets within quantization tolerance.
   - Pole-tap continuity: sampling at `(x, −ε)` equals `(x + W/2, +ε)`.
   - Seam continuity across x = 0 (no unwrap logic left to get wrong).
   - **Rigid-rotation oracle**: with a test-only `SpaceFn` that applies a fixed
     small rotation, N feedback iterations must transport a blob near the pole
     rigidly (compare against the analytically rotated blob). This fails today
     and is the regression lock for both artifacts.
   - Host H_OFFSET injection variant (device value 3) for the south-edge rule.
   - `feedback_divergence` stays green (all math uses existing `hs::`/fast_*
     primitives — no new libm calls, host/device parity preserved).
6. **Perf gates.** Expect ~neutral: populate swaps `atan2+acos` for
   `acos + rsqrt + 2 dots` at coarse rate (noise evals dominate regardless);
   composite trades 3 branches/cell for 1 mult/pixel; sampler hot path
   unchanged. Verify: host cycle A/B of `mf_feedback_flush`, then
   `just profile MeshFeedback` on device — strict per-phase GREEN (every frame
   of every phase, every preset in `MeshFeedback::PRESETS`). Confirm the HS_O3
   region still covers the edited code and byte-audit the ITCM delta (ledger is
   2,072 B under ceiling).
7. **Visual sweep** of every preset, special attention to `MeltingHi` and
   `MeltingLo` (`melt_warp` is pole-directed by design; D2 newly lets the drip
   cross the north pole — confirm it reads as intended, not as a new
   artifact).

## Watch items

- **int16 angle quantization near the pole**: `32767/π` scale → ~0.2 px of x
  quantization at row 1 (amplified by 1/sin φ). If visible as shimmer under
  animated noise, halve the covered range (`16384/rad`, still > melt max) for
  ~0.1 px, or special-case the top/bottom coarse row to float. Decide from the
  frame dumps, not preemptively.
- **Bottom-band extrapolation**: pixel rows below the last coarse row
  (y ≥ H−ds) still constant-extrapolate `cy1 = cy0`. Harmless over a smooth
  tangent field; extending the grid one row is a follow-up only if dumps show
  a residual bottom-edge artifact.
- **Warp cache**: key/layout unchanged; stale-cache semantics unaffected.

## Measured outcome — D1 rejected (2026-07-23)

D1 was implemented as specified and measured against the current path. It is a
large regression for 5 of the 8 presets it was measured against, and must not
land.

### The flaw

The decode `ddx = δe · (W/2π)/sin φ`, `ddy = δn · (H_VIRT−1)/π` is the
*differential* of the equirect projection at `v` — a first-order linearization,
not an exact reprojection. The plan's claim that the scheme is "exact at coarse
points" is false: it is exact only in the limit of small `|δ|`.

The linearization holds while `|δ| ≪ φ`. At rows 1–3 of a 288×144 frame
φ ≈ 0.02–0.07 rad, while the strong presets slide by `amplitude · 0.05` =
0.25–0.41 rad. There `|δ|` is 5–20× the distance to the pole, and the decoded
sample lands somewhere unrelated to the true source — in exactly the band the
change was meant to fix.

### Method

Encode each pixel's own sphere position into the previous frame (`r,g,b` =
`v.x,v.y,v.z`), run one flush, and decode each output pixel back to the position
it actually sampled. Compare against `space_fn(pixel_to_vector(x,y))`, which is
the exact oracle. Reported as mean/max angular error in radians per latitude
band. (Probe kept in the saved patch.)

### North cap (rows 1–3), mean angular sample error

Measured against the 8-style roster of the day (`SlowTwist, Churn, Smoke,
Frozen, Shatter, Drift, Melting, Swirling`); the shipped roster is the 12 styles
in `core/engine/styles.h`. The **amplitude** column is what the verdict turns
on, not the names: the linearization fails once `|δ| = amplitude · 0.05`
approaches φ at the cap. On the shipped roster only `Smoke` (0.51) sits in the
regime D1 improves — the other eleven run 1.56 to 11.25, i.e. at or past the
amplitudes that break here.

| preset | amplitude | current | D1 tangent | |
|---|---|---|---|---|
| SlowTwist | 4.16 | 0.0019 | 0.4082 | 215× worse |
| Churn     | 1.90 | 0.0194 | 0.0131 | 1.5× better |
| Smoke     | 0.51 | 0.0335 | 0.0160 | 2.1× better |
| Frozen    | 2.73 | 0.0126 | 0.0318 | 2.5× worse |
| Shatter   | 8.21 | 0.0027 | 0.4358 | 161× worse |
| Drift     | 4.98 | 0.0053 | 0.2499 | 47× worse |
| Melting   | 6.36 | 0.0014 | 0.3723 | 266× worse |
| Swirling  | 6.36 | 0.0014 | 0.3723 | 266× worse |

Mid-latitudes (rows 16–47) are also worse for every strong preset (e.g. Melting
0.0023 → 0.0365). The two presets that improve are the two weakest ones, the
only ones where `|δ| ≲ φ` holds at the cap — which is the predicted boundary of
the linearization, and confirms the mechanism.

Saturating the decoded azimuth at ±W/2 (needed to keep the sampler's `bx`
inside its single-step wrap contract, since `Sx` reaches ~2086 px/rad at row 1)
was ruled out as the cause: replacing it with a full modulo wrap moves
SlowTwist's cap error only 0.4082 → 0.3561.

### What this leaves

- **R1 is real but narrow.** The current path's cap error is worst for the
  *weak* presets (Smoke 0.0335, Churn 0.0194) and near-perfect for the strong
  ones. Any fix should target the low-amplitude regime, not the representation
  as a whole — on the shipped roster that is `Smoke` alone.
- **D2 (pole-crossing taps) is untestable on its own.** The equirect encoder
  folds a pole crossing into a longitude jump plus a positive Δφ, so `by < 0`
  only ever arises from interpolation undershoot. D2 becomes reachable only
  under a signed-`δn` encoding.
- **The tangent encode/decode pair is not carried.** `equirect_x_scale` is the
  decoder half this measurement invalidates, and `tangent_frame`/`sphere_log`
  encode into a representation nothing decodes. Reviving any of them belongs to
  whichever variant needs them.
- **An exact decode needs a per-pixel reprojection** (`exp` map then
  `fast_acos` + `fast_atan2`), the option this plan rejected on cost. That is
  the only route that keeps the tangent representation.
