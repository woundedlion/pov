# MeshFeedback "Smoke" — cost breakdown and optimization levers

Target: get the Smoke preset's render under **59 ms** (one 62.5 ms display
window, with margin) on Teensy 4.0 in the shipping selective-O3 image.

Measured 2026-07-17 (300 s, window 16, `HS_PROFILE_EPOCH_REVS=2900`, full
8-style cycle, wrap to `Preset: 0/8` confirmed, root cyc vs wall exact to
1.1 ppm). Reproduce:

```
HS_PROFILE_DEEP=1 bash tools/profile_one.sh MeshFeedback profile 300 16 \
  "-D HS_PROFILE_EPOCH_REVS=2900"
```

## Method

`feedback_composite` and `feedback_populate` were split into `HS_PROFILE_DEEP`
sub-scopes. These live in shared render code, so they are compiled out unless a
capture asks for them (`HS_PROFILE_DEEP=1`) — otherwise instrumenting one effect
would tax the whole roster's numbers:

| scope | granularity | covers |
|---|---|---|
| `fb_pop_warp` | per coarse cell | `space_fn` (OpenSimplex2 noise) + `pixel_to_vector` |
| `fb_pop_project` | per coarse cell | `Spherical` (atan2/acos) + delta quantization |
| `fb_comp_cell` | per coarse cell | warp-field bilerp + seam unify |
| `fb_comp_sample` | per pixel | `sample_bilinear_prev` |
| `fb_comp_color` | per **lit** pixel | the colour transform, counted only where the `NEAR_BLACK` skip misses |
| `gamut_clip` | per out-of-gamut pixel | the 16-iteration `gamut_clip_preserve_chroma` bisection |
| `fb_comp_write` | per pixel | the canvas store |

`fb_comp_color`'s call count is therefore the lit-pixel count and `gamut_clip`'s
is the gamut-clip incidence — the counter the R5 ledger deferred.

**Instrumentation cost: +1.96 ms/frame (+2.5% on the flush)**, from the
like-for-like `parse_profile.py presets` figure (Smoke flush 78.84 ms
uninstrumented → 80.80 ms here) — ~8.2 cycles per scope entry across
~143,600 entries/frame. Every figure below is on the instrumented basis; scale
down ~2.5% for shipping.

## Smoke steady-state emit — per frame

Mean over the 16 clean-hold Smoke windows (composite calls == 16/window):

```
frame                  120.94 ms  72.57 Mcyc  100%
  mf_feedback_flush     76.12 ms  45.67 Mcyc   63%
    feedback_composite  68.20 ms  40.92 Mcyc   56%
      fb_comp_color     43.25 ms  25.95 Mcyc   36%  x38495   674 cyc/px
        gamut_clip      21.57 ms  12.94 Mcyc   18%  x6626   1953 cyc/px
      fb_comp_sample     9.34 ms   5.61 Mcyc    8%  x41472   135 cyc/px
      fb_comp_cell       2.20 ms   1.32 Mcyc    2%  x10368   127 cyc/cell
      fb_comp_write      1.91 ms   1.14 Mcyc    2%  x41472    28 cyc/px
    feedback_populate    7.92 ms   4.75 Mcyc    7%
      fb_pop_warp        6.59 ms   3.95 Mcyc    5%  x2592   1525 cyc/cell
      fb_pop_project     1.01 ms   0.61 Mcyc    1%  x2592    234 cyc/cell
    feedback_litscan      0.00 ms
  mf_mesh_draw            6.09 ms   3.65 Mcyc    5%
  mf_timeline_step        0.03 ms
  mf_buffer_wait         38.70 ms  23.22 Mcyc   32%
```

Render avg **82.25 ms**, peak **89.91 ms**, spilled **243/286 frames (85%)** —
a hard 8 fps hold. `mf_buffer_wait` is the round-up idle to the second display
window, not work.

Composite residual (composite minus its four children) is **11.51 ms**, ~17% of
the composite: loop scaffolding, the `NEAR_BLACK` test on every pixel, the
per-pixel `fx` fma, the `xc.clipped(x)` test, and ~1.8 ms of the scopes' own
overhead.

### Where Smoke's cost actually sits

- **92.8%** of the 41,472 pixels are lit — at fade 0.90 the `NEAR_BLACK` skip
  saves only 7%. It is not a lever at this fade.
- **17.2%** of lit pixels leave gamut and pay the bisection, at **1,953 cyc
  each** — 6,626 px/frame, 21.57 ms.
- The base OKLab hue transform is the other **21.68 ms** (338 cyc/px over
  38,495 px): `linear_rgb_to_lms` → 3× `fast_cbrt` → 3×3 → 3 cubes → 3×3 →
  gamut test → 3× `float_to_pixel16`.
- Sampling, the warp field, the writes and the mesh together are only
  **19.5 ms** — a quarter of the render.

## The controlling comparison: Swirling

Same effect, same canvas, same sampler, same warp machinery — the only
difference is `plain_fade` instead of `hue_fade`:

| style | lit % | out-of-gamut % | colour ms | gamut ms | sample ms | populate ms | flush ms |
|---|--:|--:|--:|--:|--:|--:|--:|
| Smoke | 92.8 | 17.2 | 43.25 | 21.57 | 9.34 | 7.92 | **76.12** |
| Drift | 91.9 | 12.7 | 35.93 | 15.00 | 9.13 | 7.47 | 67.73 |
| Churn | 94.3 | 9.8 | 33.69 | 11.69 | 9.38 | 7.58 | 66.30 |
| Melting | 90.5 | 6.7 | 28.25 | 7.22 | 9.38 | 9.15 | 62.31 |
| SlowTwist | 95.2 | 9.9 | 34.43 | 11.91 | 9.55 | 0.03 | 59.95 |
| Shatter | 84.3 | 11.3 | 32.14 | 12.91 | 9.10 | 0.00 | 56.23 |
| Frozen | 60.8 | 16.9 | 27.51 | 13.74 | 8.99 | 0.00 | 50.69 |
| Swirling | 95.6 | 0.0 | **0.97** | 0.00 | 10.25 | 9.37 | **36.02** |

Swirling samples *more* pixels than Smoke and populates a *more* expensive warp
field, yet its flush is 40.10 ms cheaper — and its colour transform is 42.28 ms
cheaper. **The entire Smoke-vs-Swirling gap is the colour transform.** Give
Smoke a `plain_fade`-class colour cost and its flush lands at ~34 ms, render
~40 ms, comfortably 16 fps.

That sets the strategy: the target is reachable through colour alone, and is
*not* reachable without touching colour.

## Levers

### L1 — Cubic-form gamut bisection (exact, −16.3 ms) ★ recommended first

`gamut_clip_preserve_chroma` bisects a chroma scale `s` and, on **every one of
its 16 iterations**, re-runs the whole `oklab_to_linear_rgb` chain: the inverse
OKLab matrix, three cubes, and the LMS→RGB matrix.

But the cube roots enter linearly in `s`:

```
l_ = L + A·s,  m_ = L + B·s,  s_ = L + C·s     (A,B,C fixed per pixel)
```

so each of r, g, b is a **cubic polynomial in s** whose 12 coefficients can be
built once per pixel (~90 flops). Each iteration then collapses to three Horner
evaluations plus the six bound tests — about **20 cyc instead of ~120**.

Verified algebraically exact: over 200,000 random `(L, a, b, s)` the polynomial
form and the direct chain agree to **5.3e-15**, versus a 16-bit output quantum
of 1.53e-05. Same bisection, same iteration count, same result — only the
arithmetic is refactored.

Estimated 1,945 → ~460 cyc/px ⇒ **−16.3 ms/frame**, no visual change.

Optional follow-on: with the cheap polynomial form, dropping 16 → 8 iterations
(chroma resolution 1/256, still far below a visible step, and the bisection is
conservative — it can only under-saturate, never leave gamut) takes it to
~230 cyc/px for another ~1.5 ms. Low value once L1 has landed.

Further out: the binding channel's cubic can be solved by 3 Newton steps from
`s = 1` instead of bisected at all (~150 cyc/px), but that needs care where the
cubic is non-monotone on [0,1].

### L2 — Hoist the canvas base pointers (exact, −4 to −6 ms)

`Canvas::prev()` and `operator()` each do a **relaxed atomic load** of
`prev_`/`cur_`, an indirection through `bufs_[]`, and a runtime
`y * effect_.width_ + x` multiply — and the composite calls them **five times
per pixel** (four bilinear taps plus the store). Compilers will not hoist atomic
loads out of the loop.

Hoisting the two base pointers and the row stride out of `flush()` (W is a
template constant there, so the multiply folds) removes ~5 atomic loads,
5 pointer chases and 5 runtime multiplies per pixel. It bites into
`fb_comp_sample` (135 cyc/px), `fb_comp_write` (28 cyc/px) and the 11.5 ms
composite residual.

No behaviour change. Cheapest confidence-per-effort item on the list.

### L3 — Cheaper `fast_cbrt` (near-exact, −2 to −3 ms)

`fast_cbrt` ends in a float divide, and Cortex-M7's `VDIV.F32` is ~14 cycles and
non-pipelined. The composite calls it **three times per lit pixel** — ~42 cycles
of pure divide latency per pixel, ~2.7 ms/frame. Replacing the Halley step with
a division-free reciprocal-cube-root iteration keeps accuracy in the same band.
Shared win: `fast_cbrt` is on several effects' hot paths.

### L4 — Temporal warp-field reuse (minor look change, −3.3 ms)

`fb_pop_warp` is 6.59 ms/frame — one OpenSimplex2 evaluation per coarse cell,
2,592 cells, 1,525 cyc each, and at Smoke's speed 0.46 the field is rebuilt
every frame. Repopulating every *other* frame and lerping the two cached int16
fields (2,592 cheap lerps) halves it. The `WarpKey` cache already exists; this
is a coarser key plus a second field buffer (+5 KB persistent).

### L5 — Replace the per-pixel OKLab chain (look change, −18 to −37 ms) ★ the one that reaches target

Two variants, both attacking the 43.25 ms:

**L5a — 3D colour LUT.** The transform is a pure function of the sampled
(r,g,b) with frame-constant `k[9]`, so a 17³ trilinear LUT (4,913 entries × 6 B
= **29.5 KB**, against MeshFeedback's ~10 KB of a 266 KB persistent budget)
replaces colour *and* gamut clip with ~8 fetches and a trilinear blend,
~100 cyc/px ⇒ colour+gamut 43.25 → **~6 ms**. Build cost is ~5,000 full
transforms (~5.8 ms), but `k` only changes when fade or hue_shift change, so it
rebuilds once per preset, not per frame. Risk: interpolation error compounding
through the feedback loop, and banding — needs an eyeball, and the "Fade"/"Hue
Shift" sliders are registered as animated params, so a drifting slider would
force rebuilds.

**L5b — Small-angle linear-RGB rotation.** Smoke's `hue_shift` is 0.01 turn —
**3.6°**. A 3×3 fitted directly in linear RGB (9 mul, 6 add, clamp ≈ 30 cyc/px)
approximates so small an OKLab rotation closely, dropping the cbrt/cube round
trip and the gamut search entirely (a linear-RGB clamp suffices). Colour+gamut
43.25 → **~4 ms**. No extra memory, trivial code. Risk: it is a genuine change
to the perceptual model, and hue error compounds across feedback generations —
the accumulation is what makes this an eyeball decision rather than a
measurement one.

### L6 — Coarser warp downsample for high-scale presets (look change, −3.5 ms)

`downsample` 4 → 6 cuts coarse cells 2,592 → 1,152: populate 6.59 → 2.9 ms and
`fb_comp_cell` 2.20 → 1.0 ms. Smoke's noise scale is 23 (fine-grained), so this
visibly softens its warp. Low payoff for the look cost.

### Measured dead — do not re-attempt

- **Occupancy sample-skip** (dilated coarse lit-mask gating sample+colour):
  measured net-negative, 54% → 11% cadence. Smoke is 92.8% lit; almost nothing
  skips and the per-pixel mask lookup swamps it.
- **More `-Os`-wall chasing.** The O3 twin puts the shipping image at 0.99× the
  global-O3 ceiling. There is no compiler headroom left here.
- **Raising `NEAR_BLACK`.** At fade 0.90, 92.8% of pixels are lit; the skip is
  already near its ceiling for this preset.

## Ledger against the 59 ms target

Render peak is **89.91 ms** instrumented (~87.9 ms shipping); the target needs
about **−29 ms**.

| lever | Δ ms | output change |
|---|--:|---|
| L1 cubic-form gamut bisection | −16.3 | none (exact) |
| L2 canvas pointer hoist | −5.0 | none |
| L3 division-free `fast_cbrt` | −2.5 | sub-quantum |
| **exact subtotal** | **−23.8** | **peak ~66 ms — misses** |
| L4 temporal warp reuse | −3.3 | minor |
| **+ L4** | **−27.1** | **peak ~63 ms — still misses** |
| L5a/L5b colour replacement | −18 to −37 | needs eyeball |

**Conclusion.** The exact levers are worth landing on their own merits — L1 and
L2 in particular are free — but they land Smoke at ~63–66 ms peak, which still
spills. **Reaching 59 ms requires L5**: the per-pixel OKLab chain has to get
cheaper, not just better-compiled. With L5a or L5b alone the peak drops to
~53 ms; with L5 plus the exact levers, ~48 ms, which clears the window with real
margin and would likely pull Drift, Churn, Melting and SlowTwist under it too
(they share the same colour path — see the per-preset table).

Recommended order: **L2 → L1** (free and exact, land and re-capture), then take
L5b to the owner as a look decision, with L5a as the fallback if the small-angle
approximation drifts visibly across feedback generations.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs) — ISR CPU is ~21% in these
  windows, the point of profiling under the live driver.
- Attribution inside a scope is approximate at -O3: `HS_OS_CYCLES()` is a
  volatile read, which orders volatile accesses only, so the compiler may sink
  pure float work across a scope boundary. Leaf *call counts* are exact; leaf
  *cycles* carry a few percent of slop.
- Per-scope overhead ~8.2 cyc is included in each leaf's per-call figure;
  subtract it before comparing a leaf against hand-counted flops.
- Savings estimates for L2–L6 are static cycle counts, not measurements. Only
  L1's *exactness* is verified; its 16.3 ms is an estimate like the rest.
- The shipping roster report (`docs/profiles/shipping/`) carries the
  uninstrumented numbers and remains the figure of record; this capture reads
  +2.5% on the flush.
