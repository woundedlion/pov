# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-08-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MeshFeedback`).
Raw capture: `C:/work/hs-rc-mf/build/prof/meshfeedback_ship.log` (worktree
`C:/work/hs-rc-mf` detached at `e68a4e37`, board COM4, mtime 2026-08-05 13:12
local).

**Supersedes the earlier same-dated capture at `062c23ba`** (the finding-20
edge-clipping A/B, mtime 00:22), which this file replaces in place. Two commits
have landed since that pass and both move MeshFeedback, so none of its numbers
are carried over — every figure below is from the `e68a4e37` capture.

| tip | what it did to MeshFeedback |
|---|---|
| `fbaf81e2` | `spherical_field`: keep the RGB sampler's direct band inlinable — **−1.67 ms** on the Smoke clean-hold flush |
| `e68a4e37` | `plot`: take `RasterOptions` by value — **+0.15 ms**, an ITCM placement effect, not codegen |

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, board COM4 |
| Image | `profile` env: `-Os` + newlib-nano + DMA LEDs + `HS_PROFILE_ENABLE`; the `Feedback::flush` body (`filter.h:1833–2344`, covering `feedback_litscan`/`feedback_populate`/`feedback_composite` and the `HS_O3_FN` `sample_bilinear_prev`) and `Plot::rasterize` (`plot.h:1614`) run inside `HS_O3` regions. `core/math/spherical_field.h` carries **no** `HS_O3` region — see the flush note below |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `e68a4e37` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture, `HS_PROFILE_EPOCH_REVS=3400` |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16 '-D HS_PROFILE_EPOCH_REVS=3400'` |

Image size: `FLASH: code:83304, data:150380, headers:9004` / `RAM1:
variables:315168, code:62856, padding:2680, free:143584` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 6545–6560 root counter cyc ÷ 600 MHz
(599,942,386 ÷ 600 = 999,903.98 µs) matches the measured wall sum
(999,906 µs) within **2.0 ppm** (`tools/parse_profile.py ... validate`).

Capture validity: `validate` passes all eight checks — 418 windows, every one
named `MeshFeedback`, first frame `f 1` after a `por` reset (so the board really
rebooted onto this image), no epoch reset across 6,688 frames, all 12 styles
visited over 27 preset advances, and the cycle wrapped twice.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `mf_feedback_flush`
37.3–48.9 ms/f across the 12 styles, worst hold **Smoke** at 48.86 ms/f, peak
frame render **57.70 ms** (frames 3905–3920, SlowDust), spilled **0/6688 frames
(0%)**.

MeshFeedback sets `full_frame` (its `Feedback` filter crosses segments), so it
renders the **FULL canvas — 288×144 = 41,472 px, 4× the 10,368 px quadrant every
other effect draws.** A display window is 62.5 ms; every style holds 16 fps, and
the peak frame clears the window by 4.80 ms. The `canvas_buffer_wait` scope is
the round-up idle to the next display flip, by design — it ranges 19.9 ms on the
cheap wormholes down to 8.4 ms on the worst Smoke hold.

## Phase-by-phase readout

Phase schedule: one fixed solid (icosahedron, 30 wireframe edges), built once in
`init()` and never recompiled, drawn into a persistent feedback buffer that is
composited and warped every frame. The 12 styles cycle every 241 frames
(`Preset: i/12`), and each style is its own steady regime — cost tracks the
style's warp/fade parameters, not any transition. There is no crossfade: a style
switch is instantaneous, and the one window of ramp visible after each marker
(e.g. 41.6 → 45.8 ms/f across frames 241–272) is the feedback buffer converging
to the new style's warp, not a two-mesh draw.

### Peak frame of the capture (window frames 3905–3920, SlowDust)

```
frame                    62.74 ms  37.64 Mcyc  100%
  mf_mesh_draw            4.79 ms   2.87 Mcyc    8%  x1
    filter_blend          0.83 ms   0.50 Mcyc    1%  x8107  61.6 cyc/blend
  mf_feedback_flush      47.55 ms  28.53 Mcyc   76%  x1
    feedback_composite   40.69 ms  24.42 Mcyc   65%  x1
    feedback_populate     6.82 ms   4.09 Mcyc   11%  x1
    feedback_litscan      0.13 us     101 cyc    0%  x1
  mf_timeline_step       40.56 us   24363 cyc    0%
  mf_apply_params         0.13 us      90 cyc    0%
  canvas_clear          410.50 us  246328 cyc    1%
  canvas_buffer_wait      9.95 ms   5.97 Mcyc   16%
```

Wall min/avg/max = 61.07/62.74/65.58 ms; render avg/max = 52.79/57.70 ms. This
is the widest render spread of any hold in the pass (52.8 avg against a 57.7
peak) — SlowDust's warp field moves enough frame to frame that the flush is not
flat, which is why it owns the peak frame even though Smoke owns the higher mean
flush. `feedback_composite` at 40.69 ms is 65% of the frame and sets the
cadence; `feedback_populate` (6.82 ms) is the style-dependent seed pass and is
literally zero on the wormholes. `mf_mesh_draw` is a flat 4.79 ms on every
style, since the solid never changes.

### Worst hold (window frames 6545–6560, Smoke)

```
frame                    62.49 ms  37.50 Mcyc  100%
  mf_mesh_draw            4.79 ms   2.88 Mcyc    8%  x1
    filter_blend          0.83 ms   0.50 Mcyc    1%  x8127  61.2 cyc/blend
  mf_feedback_flush      48.86 ms  29.31 Mcyc   78%  x1
    feedback_composite   42.07 ms  25.24 Mcyc   67%  x1
    feedback_populate     6.75 ms   4.05 Mcyc   11%  x1
    feedback_litscan      0.13 us     100 cyc    0%  x1
  mf_timeline_step       41.13 us   24694 cyc    0%
  mf_apply_params         0.13 us      90 cyc    0%
  canvas_clear          411.06 us  246669 cyc    1%
  canvas_buffer_wait      8.39 ms   5.03 Mcyc   13%
```

Wall min/avg/max = 61.93/62.49/63.06 ms; render avg/max = 54.11/54.69 ms. Smoke
is the most expensive *hold* — 48.86 ms of flush against SlowDust's 47.55 — and
the tightest, with only 0.58 ms between its mean and peak render. It is
therefore the style this effect's headroom should be tracked on: it leaves
8.39 ms of idle, the least of any style, while still peaking 6.80 ms under the
62.5 ms window.

### Light style (window frames 2577–2592, TightWormhole)

```
frame                    62.42 ms  37.45 Mcyc  100%
  mf_mesh_draw            4.79 ms   2.87 Mcyc    8%  x1
    filter_blend          0.83 ms   0.50 Mcyc    1%  x8123  61.0 cyc/blend
  mf_feedback_flush      37.32 ms  22.39 Mcyc   60%  x1
    feedback_composite   37.29 ms  22.38 Mcyc   60%  x1
    feedback_populate     0.00 us      32 cyc    0%  x1
    feedback_litscan      0.13 us     103 cyc    0%  x1
  mf_timeline_step       39.50 us   23725 cyc    0%
  mf_apply_params         0.13 us      90 cyc    0%
  canvas_clear          408.75 us  245255 cyc    1%
  canvas_buffer_wait     19.86 ms  11.92 Mcyc   32%
```

Wall min/avg/max = 61.69/62.42/63.07 ms; render avg/max = 42.56/42.98 ms. The
whole 11.5 ms saving over Smoke is `feedback_populate` collapsing 6.75 → 0.00 ms
plus a 4.8 ms lighter composite; the mesh draw and its ~8,120 blends are
unchanged to three digits. Idle more than doubles to 19.9 ms, which is the
headroom the heavy styles spend.

### Per-preset table

Read from each style's clean-hold windows (`parse_profile.py ... presets`); with
one draw call per frame every hold window qualifies, so all 418 windows are
attributed. All 12 styles visited over 27 advances, and the cycle wrapped back
to style 0 **twice** — the marker sequence runs `2/12 … 12/12, 1/12, 2/12 …
12/12, 1/12, 2/12 … 4/12`. Windows per style: 30–46 (ArcingLightning gets 45
across two segments: 15 in the pre-first-marker opening hold and 30 after each
`1/12` wrap). Ranked by dominant-scope cost.

The marker in this capture is **1-based** (`Preset: 1/12 … 12/12`), so marker
value *v* is array index *v−1*; the `#` column below is the array index, as in
the previous report. Smoke is the 4th hold segment of the cycle throughout.

The `peak render` column is that style's worst *frame* including the transition
that follows it (`parse_profile.py ... buckets`) — the stricter figure the
aggregate README ranks on; `mf_feedback_flush` is the clean hold. The `062c23ba`
columns are the previous report's numbers on the same board configuration.

| # | style | mf_feedback_flush ms | peak render ms | `062c23ba` flush | Δ flush | fps |
|---|---|--:|--:|--:|--:|--:|
| 3 | Smoke | 48.86 | 55.70 | 50.34 | −1.48 | 16 |
| 4 | SlowDust | 47.60 | **57.70** | 49.18 | −1.58 | 16 |
| 1 | SlowFire | 45.78 | 53.95 | 47.42 | −1.64 | 16 |
| 2 | EnergeticFire | 44.46 | 53.39 | 46.16 | −1.70 | 16 |
| 8 | Miasma | 43.93 | 54.06 | 45.62 | −1.69 | 16 |
| 6 | MeltingHi | 43.76 | 52.73 | 45.48 | −1.72 | 16 |
| 5 | WavyTrails | 43.23 | 54.79 | 44.95 | −1.72 | 16 |
| 7 | MeltingLo | 42.40 | 49.56 | 44.13 | −1.73 | 16 |
| 9 | LooseWormhole | 41.15 | 47.61 | 42.71 | −1.56 | 16 |
| 0 | ArcingLightning | 40.76 | 51.51 | 42.58 | −1.82 | 16 |
| 11 | WigglingWormhole | 39.92 | 47.59 | 41.69 | −1.77 | 16 |
| 10 | TightWormhole | 37.32 | 46.17 | 39.02 | −1.70 | 16 |

**Every style's flush fell 1.48–1.82 ms and every style's peak render fell
1.25–1.72 ms.** The span stays at 1.3× (37.3–48.9 ms) — far tighter than the
mesh cyclers, because the full-canvas composite is a fixed cost every style pays
and only `feedback_populate` and the warp parameters move. The three Wormholes
and ArcingLightning are the cheap end; Smoke and SlowDust the expensive one. No
style changed cadence: 16 fps and 0 spilled on all 12, before and after.

### Where the 1.5 ms came from, and how stable it is

Pass averages over all 6,688 frames, against the same figures from the
`062c23ba` capture:

| scope | `062c23ba` | `e68a4e37` | delta |
|---|--:|--:|--:|
| `mf_feedback_flush` | 43.344 ms/f | 41.615 ms/f | **−1.729 ms (−4.0%)** |
| `feedback_composite` | 37.800 ms/f | 36.069 ms/f | **−1.731 ms (−4.6%)** |
| `feedback_populate` | 5.512 ms/f | 5.513 ms/f | +0.001 ms |
| `mf_mesh_draw` (total) | 4.699 ms/f | 4.790 ms/f | +0.091 ms (+1.9%) |
| `filter_blend` | 0.820 ms/f | 0.827 ms/f | +0.007 ms |
| blends/frame | 8,116 | 8,116 | 0 |
| cyc/blend | 60.6 | 61.2 | +0.6 |
| `canvas_buffer_wait` | 13.948 ms/f | 15.587 ms/f | +1.639 ms |
| peak frame render | 58.95 ms | 57.70 ms | **−1.25 ms** |
| spilled | 0/6688 (0%) | 0/6688 (0%) | — |

The saving is entirely inside `feedback_composite` (−1.731 ms, which is the
whole of the flush's −1.729 ms), `feedback_populate` is bit-identical, blend
geometry is unchanged, and the freed time falls straight into
`canvas_buffer_wait` (+1.639 ms) — the signature of a real render saving on a
cadence-locked effect.

**Mechanism (`fbaf81e2`).** `Feedback::flush` samples the previous frame through
the `HS_O3_FN` wrapper `sample_bilinear_prev`, which forwards to
`SphereField::sample_bilinear_rgb` in `core/math/spherical_field.h` — a file
with **no** `HS_O3` region, so the callee itself is `-Os`. `7846bc81` had routed
`sample_bilinear_rgb` through the generic `sample_bilinear`, growing it past
GCC's inline budget; the compositor then made a real call, with a full
register-save prologue, once per canvas pixel — ~41,472 calls per frame.
`fbaf81e2` gives the two samplers a shared pole policy (`in_direct_band()` /
`pole_tap()`) instead of one calling the other and moves the pole-crossing
footprint to a `noinline` helper, so the direct band — every row but the pole
rows — is straight-line again and inlines into the caller's pixel loop.
**This recovered a regression `7846bc81` had introduced; it is not new headroom.**

**This metric sits on a 32-byte ITCM alignment cliff — read future numbers with
that in mind.** `e68a4e37` (`plot`: take `RasterOptions` by value) touches no
code on the flush path at all, yet moved the Smoke clean-hold flush **+0.15 ms**
(48.31 → 48.46) purely by shifting where the flush body lands in ITCM. The
isolating A/B runs for the two commits read Smoke at 49.97 → 48.30 → 48.46,
while this pass reads 48.86 — a further ~0.4 ms of the same placement noise from
an independently linked image. Treat sub-0.3 ms movements on
`mf_feedback_flush` as placement, not codegen, unless a disassembly or an
A/B on one image says otherwise. The `mf_mesh_draw` +0.091 ms above is the same
story from the other side: identical blend geometry, 60.6 → 61.2 cyc/blend, and
`e68a4e37`'s own size A/B *freed* 1,376 B of ITCM.

### Per-pixel figures

The canvas is 41,472 px (full frame, 4× a quadrant). `feedback_composite` writes
all of it every frame: 25.24 Mcyc ÷ 41,472 = **609 cyc per canvas pixel** on the
worst Smoke hold, 589 cyc/px on the SlowDust peak window, 539 cyc/px on
TightWormhole. The mesh draw is a separate, style-invariant 8,116 blended
px/frame = **0.196× canvas coverage** at **61.2 cyc/blend** — 1% of the frame.
Its own scan work (mesh_draw minus `filter_blend`) is 2.38 Mcyc ÷ 8,116 =
**293 cyc per blended pixel** (287 at `062c23ba`). This effect is
composite-bound, not draw-bound: the geometry it feeds back costs a tenth of the
feedback itself.

## Column-ISR / DMA marshaling cost

```
isr_wake        x1157/frame  min/avg/max 0.68/1.99/11.95 us  cpu 3.67%
isr_pack        x145/frame   min/avg/max 6.37/7.27/9.75 us   cpu 1.67%
isr_dma_submit  x145/frame   min/avg/max 0.64/0.94/1.06 us   cpu 0.21%
```

(Window frames 3905–3920, the peak. Identical to three digits at the Smoke and
TightWormhole windows — the ISR load does not vary with style.)

- Submit is 13% of pack's CPU cost — the column marshal dominates, the DMA kick
  is nearly free.
- The wire transfer is asynchronous; submit returns in 0.94 µs and the engine
  runs against the render.
- Total ISR share is 5.55%, and it is already *inside* every render figure here
  (CYCCNT free-runs through the handlers). The peak frame renders in 57.70 ms
  against the 62.5 ms window, clearing it by **4.80 ms** (3.55 ms at
  `062c23ba`); the worst hold (Smoke, 48.86 ms flush, 55.70 ms peak render)
  clears by 6.80 ms. No speedup required.

## Summary ranking

1. `feedback_composite` — 65–67% of the frame, 40.7–42.1 ms: the full-canvas
   warp + fade + OKLab-bound blend at ~609 cyc per canvas pixel. This is the
   effect's cadence, it is inside the R5 `HS_O3` region, and it is the scope
   `fbaf81e2` gave back 1.73 ms of.
2. `feedback_populate` — 11%, 6.8 ms: the style-dependent seed pass, the only
   scope that meaningfully separates the 12 styles (0.00 ms on the wormholes).
3. `mf_mesh_draw` — 8%, 4.79 ms: drawing the fixed 30-edge solid into the
   feedback buffer, 8,116 blends at 61.2 cyc/blend. Style-invariant.
4. `canvas_clear` — under 1%, 0.41 ms: 4× the usual because the canvas is
   full-frame.
5. `mf_timeline_step` / `mf_apply_params` / `feedback_litscan` — under 45 µs
   combined; the parameter and timeline plumbing is free.

The remaining lever is `feedback_composite`'s per-pixel cost. At 609 cyc/px on
41,472 px there is no cull to add — every canvas pixel is genuinely written
every frame — so any further win has to come out of the per-pixel sampler and
colour transform, which is exactly where `fbaf81e2` just came from.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs); the 5.55% ISR share is inside
  every render figure above.
- `filter_blend` parents under `mf_mesh_draw`; its subtree is hidden in windows
  where that parent had 0 calls; calls ≈ blended pixels.
- Selective-O3: the flush path crosses two `HS_O3` regions — `Feedback::flush`
  (`filter.h:1833–2344`, so `feedback_litscan`/`feedback_populate`/
  `feedback_composite` and the `HS_O3_FN` `sample_bilinear_prev` are all `-O3`)
  and `Plot::rasterize` (`plot.h:1614`) under `mf_mesh_draw`.
  `core/math/spherical_field.h` is **not** in a region, which is why
  `sample_bilinear_rgb`'s inline budget is load-bearing (see the mechanism note
  above).
- Dwell compression: `HS_PROFILE_EPOCH_REVS=3400` stretches the epoch so one
  effect instance survives the 12-style cycle. It changes how long the pass
  runs, never a style's per-frame cost, so the ms here are the shipping figures.
- Captured from `C:/work/hs-rc-mf` detached at `e68a4e37` with a clean working
  tree (`git status --porcelain` empty) — no uncommitted state.
- Board COM4. A peer session was capturing a different effect on COM3
  concurrently; the per-board device lock kept the runs disjoint, and the
  capture's `effect=` header, all 418 window names, the `f 1` first frame and
  the log mtime were each checked against this run.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=3400`; `just profile
MeshFeedback` = build + flash + capture. This pass ran
`HS_PROFILE_TREE=/c/work/hs-rc-mf HS_TEENSY_PORT=COM4 bash tools/profile_one.sh
MeshFeedback profile 420 16 "-D HS_PROFILE_EPOCH_REVS=3400"`.
