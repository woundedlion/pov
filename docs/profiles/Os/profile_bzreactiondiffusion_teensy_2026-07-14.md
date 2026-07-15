# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_bzreactiondiffusion.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `BZReactionDiffusion<288, 144>` only (single-entry playlist) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile BZReactionDiffusion` |

Image size: FLASH code 46,284 B; ITCM (RAM1 code) 26,856 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to ~1 ppm:
window frames 385–448 root `frame` counter = 4,800,451,770 cyc =
8,000,752.95 µs vs measured `micros()` window sum 8,000,763 µs
(Δ ≈ 1.3 ppm).

## Frame cadence (context for every number below)

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window. This effect renders **one quadrant**
(this segment's 72-row band × the window's 144-column half ≈ **10,368 px** —
it is *not* full-frame) with a 4×-SSAA shader kernel. A frame that misses one
window boundary snaps to the next, so wall time quantizes to multiples of
62.5 ms. The render never fits in one window, so it holds a steady
**2-window cadence, ~125 ms/frame (8 fps)** for the entire pass
(wall min 124.4 / avg 125.0 / max 125.8 ms).

Unlike `DisplacementField`, `BZReactionDiffusion`'s `draw_frame()` is inherited
from the un-editable `ReactionDiffusionBase`, so there is **no**
`buffer_wait`/`timeline_step` scope to instrument. Only `render()` is
instrumented; the gap between the root `frame` (100%) and `bz_render` (83%) —
**~21.6 ms / 17%** — is the base class's `Canvas`-ctor buffer-wait (display
window sync idle) plus its timeline step, both uninstrumented.

## Phase-by-phase readout (64-frame windows, per-frame averages)

The reaction-diffusion pattern reaches a stable working set within the first
window and stays there: windows 2–9 all sit at avg 124.98–125.01 ms with
`bz_raster` at 74–75% of frame in every window. The effect is **steady**, so
one representative window (frames 385–448) suffices; window 1 is an
epoch-straddle (first frame partial, avg 124.1 ms) and is excluded.

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame, and for leaves
calls/frame and per-call cost (cycles = µs × 600):

```
frame                125.01 ms  75.01 Mcyc  100%
  bz_render          103.44 ms  62.07 Mcyc   83%
    bz_raster         93.46 ms  56.08 Mcyc   75%  x1  93.46 ms
    bz_physics         8.86 ms   5.31 Mcyc    7%  x1   8.86 ms
    bz_orient          1.12 ms   0.67 Mcyc    1%  x1   1.12 ms
```

`bz_raster` is the `Scan::Shader::draw` 4×-SSAA rasterize (the hot leaf);
`bz_physics` is the reaction-diffusion `advance_substeps` ping-pong (the
simulation); `bz_orient` is the lattice→world node rotation. `bz_render` self
time (physics/orient/raster dispatch, orient setup) is ≈ 0.01 Mcyc/frame —
essentially all of `render()` lives in the three children. The **17% not
covered by `bz_render`** (~21.6 ms/frame) is the base-class display-sync idle
+ timeline step described above.

### Per-pixel figures

`bz_raster` ≈ 56.08 Mcyc/frame over the ~10,368-px quadrant ⇒ **~5,409 cyc
(9.01 µs) per rendered pixel**. The shader kernel-walks ~4 samples/pixel under
4× SSAA, so that is **~1,352 cyc (2.25 µs) per SSAA sample** — the cost is the
per-sample shade + blend, not the reaction-diffusion field (which is the
separate `bz_physics` phase, ~8.9 ms/frame over the full simulation lattice,
just 7% of frame). `bz_orient` rotates the lattice nodes once per frame for
1.1 ms.

## Column-ISR / DMA marshaling cost (per-window ISR counters)

Measured with `HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree
is main-loop-only). Nested as executed — the pack and submit run inside the
flywheel wake, on the 1-in-8 wakes that render a column. Columns: rate,
per-call min / avg / max (µs), CPU share. The **min is the stable per-call
floor**; the avg/max are inflated by serial-dump preemption during the long
125 ms frames:

```
isr_wake         18432/s  0.57 / 4.2  / ~210-230 us  cpu 7.7-7.9%
  isr_pack        2304/s  4.37 / 26   / ~205-228 us  cpu 6.0-6.2%
  isr_dma_submit  2304/s  0.62 / 0.97 / 1.2-4.0 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid; it
**includes the nested `isr_pack` + `isr_dma_submit`**); `isr_pack` = the 72×
`packPixel` marshal, once per column; `isr_dma_submit` = `submitFrame` (overrun
check + dcache flush + eDMA kick).

- **The DMA submit itself is ~1 µs of CPU** (min 0.62, avg 0.97). The
  marshaling around it (pixel pack) has a stable floor of **~4.4 µs/column**;
  the ~26 µs pack *average* and the ~180–230 µs maxima are preemption
  artifacts (higher-priority USB/DMA-completion ISRs and the serial dumps
  themselves landing inside the scope), not real pack work.
- The wake grid runs at 18,432/s regardless of render load; the per-column
  CPU floor (pack + submit ≈ 5.4 µs of the 434 µs column period) is ~1.2%.
  The window-averaged CPU shares read high (wake 7.7–7.9%) because every
  125 ms frame spans many columns while a serial dump preempts a handful of
  them — the min floors, not the averages, are the shipping steady-state cost.

## Summary ranking (steady frame, share of the 125 ms frame)

1. `Scan::Shader::draw` 4×-SSAA rasterize (`bz_raster`) — **75%** (93.5 ms)
2. base-class display-sync idle + timeline (`frame`−`bz_render` gap) — 17% (21.6 ms, idle by design)
3. reaction-diffusion simulation (`bz_physics`) — 7% (8.9 ms)
4. lattice→world orient (`bz_orient`) — 1% (1.1 ms)

The rasterizer dominates on device: the SSAA shade+blend over the quadrant is
~3× the reaction-diffusion simulation it visualizes. The RD perf ledger's host
wins (BZ −14% via SSAA trig tables, world-space nodes, and the fused kernel)
targeted exactly this leaf; on device at `-Os` with ISRs live it is still 75%
of the frame, so the rasterizer remains the only lever worth pulling.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs
  the flywheel/DMA/USB ISRs that fire inside it. That overhead is the real
  shipping condition (the point of segmented-mode profiling), not a pure-CPU
  algorithm cost.
- **No buffer_wait / timeline_step split**: `BZReactionDiffusion` inherits
  `draw_frame()` from the un-editable `ReactionDiffusionBase`, so only
  `render()` is instrumented. The 17% `frame`−`bz_render` gap is that base
  class's display-sync idle + timeline step lumped together, not measured
  separately.
- Per-frame value = window total ÷ 64; per-call = cyc ÷ calls; µs = cyc ÷ 600.
  Each `bz_*` scope runs once per frame (x1), so per-call equals per-frame.
- **`-Os` build**: shipping Phantasm config; an `-O3` profile would differ.
- Run with **uncommitted working-tree `HS_PROFILE` instrumentation** in
  `effects/BZReactionDiffusion.h` (the `bz_render`/`bz_physics`/`bz_orient`/
  `bz_raster` scopes).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=<EffectClass>`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake,
  pixel pack, DMA submit) each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile BZReactionDiffusion [seconds]`.
