# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_gsreactiondiffusion.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `GSReactionDiffusion<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile GSReactionDiffusion` |

Image size: FLASH code 46,372 B; ITCM (RAM1 code) 26,696 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to ~1 ppm:
window frames 193–256 root `frame` counter = 4,799,737,155 cyc = 7,999,561.9 µs
vs measured `micros()` window sum 7,999,572 µs (Δ ≈ 1.3 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Unlike its sibling
`BZReactionDiffusion`, GSReactionDiffusion carries an explicit `draw_frame()`
override (behaviorally identical to `ReactionDiffusionBase`), so it exposes the
`grd_buffer_wait` / `grd_timeline_step` scopes.

Render is steady at **~110 ms/frame** — well over one 62.5 ms window but under
two, so wall time snaps to a 2-window cadence, **125 ms/frame (8 fps)**:

- **Render** = `grd_simulate` (Gray-Scott float substeps) + `grd_rasterize`
  (4× SSAA shader raster) ≈ 109.8 ms/frame → overflows into a second window.

`grd_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle — ~15 ms
of the 125 ms frame once render fills the first window.

## Phase-by-phase readout (64-frame windows, per-frame averages)

GSReactionDiffusion is **steady** — every window from frames 65 on is within a
percent of the next (render 85–88 % of frame, rasterize 68 % of render, simulate
31 % of render, wall ~125 ms). The first window (frames 1–64) is the epoch-start
ramp — low `min` wall 58.0 ms, `grd_buffer_wait` 30 % while render is still
cheap — and is skipped for the representative below.

### Steady state — 8 fps (representative: frames 193–256)

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame; leaves add
calls/frame and per-call/per-pixel cost (cycles = µs × 600):

```
frame                125.00 ms  75.00 Mcyc  100%
  grd_render         109.84 ms  65.91 Mcyc   88%
    grd_rasterize     75.33 ms  45.20 Mcyc   60%
      grd_shader_draw 71.38 ms  42.83 Mcyc   57%  x1  10368 px  4131 cyc/px
    grd_simulate      34.51 ms  20.71 Mcyc   28%  x1
  grd_buffer_wait     15.13 ms   9.08 Mcyc   12%
  grd_timeline_step     22 us     13 kcyc     0%
```

Wall: min 124.15 / avg 124.99 / max 125.86 ms. `grd_shader_draw` is the
`Scan::Shader::draw<W,H,4>` 4× SSAA per-pixel kernel-walk raster (the single
hottest leaf, 94 % of `grd_rasterize`); `grd_rasterize` self-time (75.33 − 71.38
= 3.95 ms) is the `orient_nodes` + `fill_hot_flags` cull prep. `grd_render`
self-time is ~0 — render is exactly rasterize + simulate.

**Contrast with `BZReactionDiffusion`**: Gray-Scott spends a far larger share on
the physics substep loop. `grd_simulate` is **28 % of the frame / 31 % of render
≈ 34.5 ms/frame** here, versus ~8 % for BZ. The GS float ping-pong (STEPS_PER_
FRAME substeps + the Q16 convert-in / quantize-out around them) is a first-class
cost, not a rounding error — attack it alongside the rasterizer, not after.

### Per-pixel figures

`Scan::Shader::draw` composites in place and does **not** route through the
`filter_blend` per-pixel counter, so per-pixel cost is derived from the fixed
quadrant coverage: **10,368 px/frame**, no cull. At 4× SSAA the kernel walks
~4 samples/pixel, so `grd_shader_draw` = 42.83 Mcyc ⇒ **4,131 cyc (6.89 µs) per
output pixel ≈ 1,033 cyc (1.72 µs) per SSAA sample** — dominated by the per-
sample Gray-Scott concentration lookup and palette shade, not a blend.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is main-loop-
only), representative steady window frames 193–256 (8.0 s). The pack and submit
run inside the flywheel wake, on the 1-in-8 wakes that render a column. Columns:
rate, per-call min/avg/max, CPU share:

```
isr_wake        18300/s  0.63 /  5.21 / 163.6 us  cpu 9.52%
  isr_pack       2304/s  4.37 / 33.99 / 161.2 us  cpu 7.83%
  isr_dma_submit 2304/s  0.62 /  0.97 /   2.0 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,300/s ≈ the COLUMN_US/8 wake grid, pack +
submit nested inside); `isr_pack` = the 72× `packPixel` marshal, once per
column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush + eDMA
kick).

- **The DMA kick itself is ~1 µs of CPU** (submit min 0.62 µs); the pixel-pack
  marshal around it is the real per-column cost (min 4.37 µs). The wire transfer
  runs asynchronously on eDMA and costs no CPU.
- The `avg`/`max` figures are heavily inflated here — pack avg 34 µs, wake max
  164 µs — because the once-per-window serial dump preempts whatever ISRs fire
  inside it, and the 125 ms frames give preemption a long runway. The stable
  per-call `min` floor (pack 4.37 µs, submit 0.62 µs, wake 0.63 µs) is the true
  cost; the reported `cpu` share rides the same inflation.
- Net: taking the min floor, per-column CPU ≈ pack 4.37 + submit 0.62 + the
  seven non-pack wakes ≈ **~9 µs of the 434 µs column period (~2 %)**, i.e. the
  ISR machinery steals a few percent of the chip. Against the 62.5 ms window
  budget that is ~1.5 ms; render (110 ms) needs ~1.8 windows regardless.

## Summary ranking (steady state, share of the 125 ms frame)

1. `Scan::Shader::draw` 4× SSAA per-pixel raster (`grd_shader_draw`) — **57%**
   (71.4 ms) — the single hottest leaf.
2. Gray-Scott physics substeps + Q16 convert (`grd_simulate`) — **28%**
   (34.5 ms) — a much bigger slice than BZ's ~8 %.
3. display-window sync (`grd_buffer_wait`) — 12% (15.1 ms, idle by design).
4. rasterize cull prep (`orient_nodes` + `fill_hot_flags`, rasterize self) —
   ~3% (4.0 ms).
5. everything else (`grd_timeline_step`) — <1%.

Render is split between two real levers here — the SSAA rasterizer (57 %) and
the simulation (28 %) — where BZ is almost pure rasterizer. The RD perf ledger's
host-side GS work (−37 % via world-space nodes + fused kernel + float substeps)
already reshaped `grd_simulate`; the on-device split shows it is still a third
of render at `-Os` with ISRs live.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Scan::Shader::draw` composites in place, so the
  per-pixel blended-pixel counter never fires; the per-pixel figures are derived
  from the fixed 10,368-px quadrant coverage at 4× SSAA.
- **ISR avg/max inflated**: the once-per-window serial dump preempts in-flight
  ISRs and inflates their measured duration (and the `cpu` share); trust the
  `min` floor.
- Epoch-straddle: the first window (frames 1–64) mixes the ramp / low-`min`
  epoch start into its averages and is skipped when picking a representative.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/GSReactionDiffusion.h` (scopes: `grd_buffer_wait`,
  `grd_timeline_step`, `grd_render`, `grd_simulate`, `grd_rasterize`,
  `grd_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=GSReactionDiffusion`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile GSReactionDiffusion [seconds]`.
