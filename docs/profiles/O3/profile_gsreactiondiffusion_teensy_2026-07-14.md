# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (-O3) (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_gsreactiondiffusion_o3.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm profile flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`, `-D HS_PROFILE_ENABLE`) but built at **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` — the -O3 twin of the shipping `-Os` profile |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `GSReactionDiffusion<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 69,212 B; ITCM (RAM1 code) 47,496 B; RAM2 free 4,736 B.
(At `-Os`: 46,372 / 26,696 / 4,736 — see the -Os counterpart in `docs/profiles/Os/`.)

**Exactness cross-check** — the cycle counter and wall clock agree to ~0.1 ppm:
window frames 193–256 root `frame` counter = 4,799,907,560 cyc = 7,999,845.9 µs
vs measured `micros()` window sum 7,999,847 µs (Δ ≈ 0.13 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Like the -Os run,
GSReactionDiffusion carries an explicit `draw_frame()` override (behaviorally
identical to `ReactionDiffusionBase`), so it exposes the `grd_buffer_wait` /
`grd_timeline_step` scopes.

Even at -O3, render is steady at **~84.9 ms/frame** — faster than -Os (109.8 ms)
but still well over one 62.5 ms window, so wall time still snaps to a 2-window
cadence, **125 ms/frame (8 fps)**. The frame does **not** drop below 125 ms:

- **Render** = `grd_simulate` (Gray-Scott float substeps) + `grd_rasterize`
  (4× SSAA shader raster) ≈ 84.9 ms/frame → still overflows into a second window.

`grd_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle — the
~25 ms -O3 saves in render just widens this idle from ~15 ms (-Os) to **~40 ms**;
it does not lift the fps.

## Phase-by-phase readout (64-frame windows, per-frame averages)

GSReactionDiffusion is **steady** — all eight windows from frames 65 on are
within a percent of each other (render 68 % of frame, rasterize 67 % of render,
simulate 33 % of render, wall ~125 ms). The first window (frames 1–64) is the
epoch-start ramp — low `min` wall 47.4 ms, `grd_buffer_wait` 38 % while render is
still cheap — and is skipped for the representative below.

### Steady state — 8 fps (representative: frames 193–256)

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame; leaves add
calls/frame and per-call/per-pixel cost (cycles = µs × 600):

```
frame                125.00 ms  75.00 Mcyc  100%
  grd_render          84.87 ms  50.92 Mcyc   68%
    grd_rasterize     57.01 ms  34.20 Mcyc   46%
      grd_shader_draw 53.56 ms  32.14 Mcyc   43%  x1  10368 px  3100 cyc/px
    grd_simulate      27.86 ms  16.72 Mcyc   22%  x1
  grd_buffer_wait     40.11 ms  24.07 Mcyc   32%
  grd_timeline_step     20 us     12 kcyc     0%
```

Wall: min 124.19 / avg 125.00 / max 125.72 ms. `grd_shader_draw` is the
`Scan::Shader::draw<W,H,4>` 4× SSAA per-pixel kernel-walk raster (still the
single hottest leaf, 94 % of `grd_rasterize`); `grd_rasterize` self-time (57.01 −
53.56 = 3.45 ms) is the `orient_nodes` + `fill_hot_flags` cull prep. `grd_render`
self-time is ~0 — render is exactly rasterize + simulate (57.01 + 27.86 = 84.87).

**Contrast with `BZReactionDiffusion`**: Gray-Scott still spends a far larger
share on the physics substep loop. `grd_simulate` is **22 % of the frame / 33 %
of render ≈ 27.9 ms/frame** here — a much bigger slice than BZ's single digits.
The GS float ping-pong (STEPS_PER_FRAME substeps + the Q16 convert-in / quantize-
out around them) remains a first-class cost at -O3, not a rounding error.

### Per-pixel figures

`Scan::Shader::draw` composites in place and does **not** route through the
`filter_blend` per-pixel counter, so per-pixel cost is derived from the fixed
quadrant coverage: **10,368 px/frame**, no cull. At 4× SSAA the kernel walks
~4 samples/pixel, so `grd_shader_draw` = 32.14 Mcyc ⇒ **3,100 cyc (5.17 µs) per
output pixel ≈ 775 cyc (1.29 µs) per SSAA sample** — dominated by the per-sample
Gray-Scott concentration lookup and palette shade, not a blend. (At -Os this was
4,131 cyc/px; -O3 shaves ~1,030 cyc/px off the SSAA kernel.)

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is main-loop-
only), representative steady window frames 193–256 (8.0 s). The pack and submit
run inside the flywheel wake, on the 1-in-8 wakes that render a column. Columns:
rate, per-call min/avg/max (ns triple → µs), CPU share:

```
isr_wake        18316/s  0.51 /  4.77 / 167.6 us  cpu 8.74%
  isr_pack       2304/s  4.65 / 31.55 / 166.0 us  cpu 7.26%
  isr_dma_submit 2304/s  0.62 /  0.93 /   1.3 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (18,316/s ≈ the COLUMN_US/8 wake grid, pack +
submit nested inside); `isr_pack` = the 72× `packPixel` marshal, once per column;
`isr_dma_submit` = `submitFrame` (overrun check + dcache flush + eDMA kick). These
figures are driver-side and essentially unchanged from the -Os run — the ISR path
does not depend on the effect's optimization level.

- **The DMA kick itself is ~1 µs of CPU** (submit min 0.62 µs); the pixel-pack
  marshal around it is the real per-column cost (min 4.65 µs). The wire transfer
  runs asynchronously on eDMA and costs no CPU.
- The `avg`/`max` figures are heavily inflated here — pack avg 32 µs, wake max
  168 µs — because the once-per-window serial dump preempts whatever ISRs fire
  inside it, and the 125 ms frames give preemption a long runway. The stable
  per-call `min` floor (pack 4.65 µs, submit 0.62 µs, wake 0.51 µs) is the true
  cost; the reported `cpu` share rides the same inflation.
- Net: taking the min floor, per-column CPU ≈ pack 4.65 + submit 0.62 + the seven
  non-pack wakes ≈ **~9 µs of the 434 µs column period (~2 %)**, i.e. the ISR
  machinery steals a few percent of the chip. Against the 62.5 ms window budget
  that is ~1.5 ms; render (84.9 ms) needs ~1.4 windows regardless.

## Summary ranking (steady state, share of the 125 ms frame)

1. `Scan::Shader::draw` 4× SSAA per-pixel raster (`grd_shader_draw`) — **43%**
   (53.6 ms) — still the single hottest leaf.
2. Gray-Scott physics substeps + Q16 convert (`grd_simulate`) — **22%**
   (27.9 ms) — a much bigger slice than BZ.
3. display-window sync (`grd_buffer_wait`) — 32% (40.1 ms, idle by design; grew
   because -O3 shortened render).
4. rasterize cull prep (`orient_nodes` + `fill_hot_flags`, rasterize self) —
   ~3% (3.4 ms).
5. everything else (`grd_timeline_step`) — <1%.

Render is still split between two real levers — the SSAA rasterizer (43 %) and
the simulation (22 %). -O3 speeds both (see below) but the rasterizer more, so
the SSAA kernel stays the top target; the frame stays pinned at 8 fps because
84.9 ms still spills into a second display window.

## -O3 vs -Os

Same board, same driver/ISRs, same effect — the *only* variable is the base
optimization level (`profile_o3` keeps `-O3 -ffast-math -fno-finite-math-only`;
`profile` forces `-Os`). Representative window frames 193–256 in both runs:

| metric | -Os | -O3 | speedup |
|---|---|---|---|
| `grd_shader_draw` (SSAA raster) | 71.38 ms | 53.56 ms | **1.33×** (−17.8 ms, −25%) |
| `grd_simulate` (GS physics) | 34.51 ms | 27.86 ms | **1.24×** (−6.7 ms, −19%) |
| `grd_render` (rasterize + simulate) | 109.84 ms | 84.87 ms | 1.29× (−25.0 ms) |
| `grd_shader_draw` cyc/output px | 4,131 | 3,100 | 1.33× |
| `grd_buffer_wait` (round-up idle) | 15.13 ms | 40.11 ms | +25.0 ms (absorbs the render win) |
| frame wall / cadence | 125 ms (8 fps) | 125 ms (8 fps) | **unchanged** |
| FLASH code | 46,372 B | 69,212 B | +22,840 B (+49%) |
| ITCM (RAM1 code) | 26,696 B | 47,496 B | +20,800 B (+78%) |

**Headline:** -O3 speeds the SSAA rasterizer (`grd_shader_draw` **1.33×**) *more*
than the physics substep loop (`grd_simulate` **1.24×**) — `-ffast-math` plus
tighter FP scheduling helps the per-sample float shade/lookup more than the
Gray-Scott ping-pong, which is already gated by memory traffic (Q16 convert-in /
quantize-out around the float buffers). Render as a whole drops 25.0 ms (1.29×),
but the **frame stays at 125 ms / 8 fps**: render (84.9 ms) still exceeds one
62.5 ms window, so the cadence still rounds up to two windows and the entire
25 ms saving just widens `grd_buffer_wait` idle from 15 to 40 ms. To actually
lift fps, render must fall below **62.5 ms** — another ~22 ms off the top two
leaves combined.

Note -O3 is **not** the shipping config: it costs +49% FLASH and +78% ITCM here,
and the full 26-effect Phantasm image overflows FlexRAM/ITCM at -O3 (that is
exactly why the ship build is -Os). -O3 is only viable as this single-effect
profile image.

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
- **`-O3` build (does NOT ship)**: this is the `profile_o3` env; the shipping
  Phantasm image is `-Os` because the full 26-effect roster overflows ITCM at
  -O3. This report exists only to quantify the optimization-level headroom.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/GSReactionDiffusion.h` (scopes: `grd_buffer_wait`,
  `grd_timeline_step`, `grd_render`, `grd_simulate`, `grd_rasterize`,
  `grd_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=GSReactionDiffusion`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py`
  (the -O3 twin; the -Os report uses `just profile GSReactionDiffusion`).
