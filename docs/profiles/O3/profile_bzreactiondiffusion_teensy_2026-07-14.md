# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode, -O3 (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_bzreactiondiffusion_o3.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: same Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`, `-D HS_PROFILE_ENABLE`) but keeps the base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forcing `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `BZReactionDiffusion<288, 144>` only (single-entry playlist) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 66,532 B; ITCM (RAM1 code) 45,384 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to <1 ppm:
window frames 385–448 root `frame` counter = 4,800,683,994 cyc =
8,001,139.99 µs vs measured `micros()` window sum 8,001,147 µs
(Δ ≈ 0.88 ppm).

## Frame cadence (context for every number below)

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window. This effect renders **one quadrant**
(this segment's 72-row band × the window's 144-column half ≈ **10,368 px** —
it is *not* full-frame) with a 4×-SSAA shader kernel. A frame that misses one
window boundary snaps to the next, so wall time quantizes to multiples of
62.5 ms. At -O3 the render (`bz_render`) drops to **71.4 ms** — faster than at
-Os but **still over one 62.5 ms window** — so the effect holds the same steady
**2-window cadence, ~125 ms/frame (8 fps)** for the entire pass
(window 385–448 wall min 124.5 / avg 125.0 / max 125.6 ms). The -O3 speedup
does **not** change the frame rate: the saved render time becomes display-sync
idle, not extra frames (see `## -O3 vs -Os`).

Unlike `DisplacementField`, `BZReactionDiffusion`'s `draw_frame()` is inherited
from the un-editable `ReactionDiffusionBase`, so there is **no**
`buffer_wait`/`timeline_step` scope to instrument. Only `render()` is
instrumented; the gap between the root `frame` (100%) and `bz_render` (57%) —
**~53.6 ms / 43%** — is the base class's `Canvas`-ctor buffer-wait (display
window sync idle) plus its timeline step, both uninstrumented. That gap is much
larger than the -Os run's 17% precisely because render shrank while the wall
stayed quantized at 125 ms.

## Phase-by-phase readout (64-frame windows, per-frame averages)

The reaction-diffusion pattern reaches a stable working set within the first
window and stays there: windows 2–9 (frames 65–576) all sit at avg
124.98–125.02 ms with `bz_raster` at ~90% of `bz_render` in every window. The
effect is **steady**, so one representative window (frames 385–448) suffices;
window 1 (frames 1–64) is an epoch-straddle (first frame partial, avg 109.5 ms,
min 58.0 ms) and is excluded.

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame, and for leaves
calls/frame and per-call cost (cycles = µs × 600):

```
frame                125.02 ms  75.01 Mcyc  100%
  bz_render           71.38 ms  42.83 Mcyc   57%
    bz_raster         64.33 ms  38.60 Mcyc   51%  x1  64.33 ms
    bz_physics         6.49 ms   3.89 Mcyc    5%  x1   6.49 ms
    bz_orient          0.57 ms   0.34 Mcyc  0.5%  x1   0.57 ms
```

`bz_raster` is the `Scan::Shader::draw` 4×-SSAA rasterize (the hot leaf);
`bz_physics` is the reaction-diffusion `advance_substeps` ping-pong (the
simulation); `bz_orient` is the lattice→world node rotation. `bz_render` self
time (physics/orient/raster dispatch, orient setup) is ≈ 0.00 Mcyc/frame
(75 cyc/frame) — essentially all of `render()` lives in the three children. The
**43% not covered by `bz_render`** (~53.6 ms/frame) is the base-class
display-sync idle + timeline step described above; it grew from the -Os run's
17% because -O3 cut render time while the 125 ms cadence stayed fixed.

### Per-pixel figures

`bz_raster` ≈ 38.60 Mcyc/frame over the ~10,368-px quadrant ⇒ **~3,723 cyc
(6.20 µs) per rendered pixel**. The shader kernel-walks ~4 samples/pixel under
4× SSAA, so that is **~931 cyc (1.55 µs) per SSAA sample** — the cost is the
per-sample shade + blend, not the reaction-diffusion field (which is the
separate `bz_physics` phase, ~6.5 ms/frame over the full simulation lattice,
just 5% of frame). `bz_orient` rotates the lattice nodes once per frame for
0.57 ms. (At -Os these were 5,409 cyc/px and 1,352 cyc/sample — see the
comparison section.)

## Column-ISR / DMA marshaling cost (per-window ISR counters)

Measured with `HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree
is main-loop-only). Nested as executed — the pack and submit run inside the
flywheel wake, on the 1-in-8 wakes that render a column. Columns: rate,
per-call min / avg / max (µs), CPU share. The **min is the stable per-call
floor**; the avg/max are inflated by serial-dump preemption during the long
125 ms frames:

```
isr_wake         18432/s  0.51 / 4.3  / ~197-210 us  cpu 7.9-8.1%
  isr_pack        2304/s  4.65 / 28   / ~195-208 us  cpu 6.4-6.6%
  isr_dma_submit  2304/s  0.59 / 0.93 / 1.2-3.4 us   cpu 0.21%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid; it
**includes the nested `isr_pack` + `isr_dma_submit`**); `isr_pack` = the 72×
`packPixel` marshal, once per column; `isr_dma_submit` = `submitFrame` (overrun
check + dcache flush + eDMA kick).

- **The DMA submit itself is ~1 µs of CPU** (min 0.59, avg 0.93). The
  marshaling around it (pixel pack) has a stable floor of **~4.65 µs/column**;
  the ~28 µs pack *average* and the ~180–210 µs maxima are preemption
  artifacts (higher-priority USB/DMA-completion ISRs and the serial dumps
  themselves landing inside the scope), not real pack work. The ISR layer is
  ordinary marshaling code, essentially unaffected by -O3 (pack floor 4.65 µs
  vs 4.37 µs at -Os — within measurement noise).
- The wake grid runs at 18,432/s regardless of render load; the per-column
  CPU floor (pack + submit ≈ 5.2 µs of the 434 µs column period) is ~1.2%.
  The window-averaged CPU shares read high (wake 7.9–8.1%) because every
  125 ms frame spans many columns while a serial dump preempts a handful of
  them — the min floors, not the averages, are the shipping steady-state cost.

## Summary ranking (steady frame, share of the 125 ms frame)

1. `Scan::Shader::draw` 4×-SSAA rasterize (`bz_raster`) — **51%** (64.3 ms)
2. base-class display-sync idle + timeline (`frame`−`bz_render` gap) — 43% (53.6 ms, idle by design)
3. reaction-diffusion simulation (`bz_physics`) — 5% (6.5 ms)
4. lattice→world orient (`bz_orient`) — 0.5% (0.6 ms)

The rasterizer still dominates the *render* on device: the SSAA shade+blend
over the quadrant is ~10× the reaction-diffusion simulation it visualizes and
90% of `bz_render`. At -O3 it is 1.45× faster than at -Os, but because render
still overruns one 62.5 ms window the extra time is absorbed as display-sync
idle (rank #2 grows from 17% to 43%) rather than becoming frames — the effect
stays at 8 fps. The rasterizer remains the only lever that could move cadence:
`bz_render` is 71.4 ms, so shaving ~9 ms off `bz_raster` would drop it under
62.5 ms and jump the whole effect to a 1-window / 16 fps cadence.

## -O3 vs -Os

Same board, driver, effect, and instrumentation — the **only** variable is the
optimization level (`profile_o3` keeps `-O3`; `profile` forces `-Os`). Compared
to `docs/profiles/Os/profile_bzreactiondiffusion_teensy_2026-07-14.md`:

| Metric | -Os | -O3 | -O3 speedup |
|---|---|---|---|
| `bz_raster` (hot leaf) | 93.46 ms / 56.08 Mcyc | 64.33 ms / 38.60 Mcyc | **1.45×** (−29.1 ms) |
| `bz_render` (all render) | 103.44 ms / 62.07 Mcyc | 71.38 ms / 42.83 Mcyc | 1.45× (−32.1 ms) |
| `bz_physics` (RD sim) | 8.86 ms | 6.49 ms | 1.37× (−2.4 ms) |
| `bz_orient` | 1.12 ms | 0.57 ms | 1.98× (−0.55 ms) |
| cyc/rendered px | 5,409 (9.01 µs) | 3,723 (6.20 µs) | 1.45× |
| Frame wall / cadence | 125.0 ms / 2-window / 8 fps | 125.0 ms / 2-window / 8 fps | **unchanged** |
| FLASH code | 46,284 B | 66,532 B | +20,248 B (+43.7%) |
| ITCM (RAM1 code) | 26,856 B | 45,384 B | +18,528 B (+69.0%) |
| RAM2 free | 4,736 B | 4,736 B | unchanged |

Headline: **-O3 makes the render ~1.45× faster (rasterizer 93.5→64.3 ms), but
the frame rate is identical** (both 125 ms / 8 fps). Render still exceeds one
62.5 ms window at -O3, so the ~32 ms saved converts to display-sync idle, not
frames — the `frame`−`bz_render` gap grows from 17% to 43%. The cost of -O3 is
+20 KB FLASH and +18.5 KB ITCM (+69%). **-O3 is not the shipping config**: the
full 26-effect Phantasm image overflows FlexRAM/ITCM at -O3 (that ITCM blow-up
is exactly why the ship build is `-Os`); -O3 is only viable here as a
single-effect profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs
  the flywheel/DMA/USB ISRs that fire inside it. That overhead is the real
  shipping condition (the point of segmented-mode profiling), not a pure-CPU
  algorithm cost.
- **No buffer_wait / timeline_step split**: `BZReactionDiffusion` inherits
  `draw_frame()` from the un-editable `ReactionDiffusionBase`, so only
  `render()` is instrumented. The 43% `frame`−`bz_render` gap is that base
  class's display-sync idle + timeline step lumped together, not measured
  separately.
- Per-frame value = window total ÷ 64; per-call = cyc ÷ calls; µs = cyc ÷ 600.
  Each `bz_*` scope runs once per frame (x1), so per-call equals per-frame.
- **-O3 build**: this is the `profile_o3` env, which does **not** ship. The
  shipping Phantasm image is `-Os` because the full effect roster overflows
  ITCM at -O3; -O3 is only usable as this single-effect profile image.
- Run with **uncommitted working-tree `HS_PROFILE` instrumentation** in
  `effects/BZReactionDiffusion.h` (the `bz_render`/`bz_physics`/`bz_orient`/
  `bz_raster` scopes).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=<EffectClass>`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake,
  pixel pack, DMA submit) each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py`
  (the -O3 env has no `just` recipe; the -Os twin is `just profile
  BZReactionDiffusion [seconds]`).
