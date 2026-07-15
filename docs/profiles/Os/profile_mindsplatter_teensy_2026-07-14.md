# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_mindsplatter.log` (64-frame windows, ~80 s single pass, captured
from boot). **Re-captured after the 2026-07-14 plot column-cull batch**
(column arc cull 9ac8cebd, planar-edge cull 62450701, deferred trail shading
708d4b9b) — the saturated scan dropped 107.75 → 96.48 ms (−10.5 %) vs the
previous capture of this report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — the shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `MindSplatter<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile MindSplatter` |

Image size: FLASH code 53,092 B; ITCM (RAM1 code) 36,520 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
saturated window frames 129–192 root counter = 4,804,926,606 cyc =
8,008,211.0 µs vs measured `micros()` window sum 8,008,221 µs (Δ ≈ 1.2 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Wall time snaps
up to a whole number of 62.5 ms windows, so render < 62.5 ms → 62.5 ms/frame
(16 fps); 62.5–125 ms → 125 ms/frame (8 fps).

MindSplatter is **particle-rasterizer-bound and ramps** as its particle pool
fills after the cold start: `Plot::ParticleSystem::draw` (`msp_particle_scan`)
grows from ~27 % to ~77 % of the frame, and wall climbs from ~61 ms to a steady
125 ms over the first ~130 frames. The cadence therefore flips once, from 16 fps
to 8 fps, as more particles land on the quadrant:

- **Cold start / low pool** (frames 1–64): render ~18 ms → fits one window →
  **62.5 ms/frame (16 fps)** (avg 61.2 ms wall; first frame `min` 501 µs is the
  empty pool, so this window is not used as an exactness/representative sample).
- **Growing pool** (frames 65–128): render crosses the 62.5 ms boundary
  mid-window as the pool fills; cadence flips 16→8 fps (wall min 57.4 / max
  128.0 ms — the sub-window minima are boundary catch-up frames right at the
  flip, where the prior frame overran a boundary and the next wait targets the
  nearer one).
- **Saturated pool** (frames 129 onward): render ~100 ms → 2 windows →
  **125 ms/frame (8 fps)**, steady for the rest of the pass.

`msp_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle — ~44 ms
when render fits one window, shrinking to ~25 ms once render fills most of two.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Two regimes bound the ramp. Each block is nested as in the counter tree (each
parent includes its children). Columns: time/frame, cycles/frame, % of frame;
leaves add calls/frame and per-call cost.

### Early / growing pool — cadence climbing 16→8 fps (window: frames 65–128)

```
frame                 104.26 ms  62.56 Mcyc  100%
  msp_draw_particles   66.42 ms  39.85 Mcyc   64%
    msp_particle_scan  66.42 ms  39.85 Mcyc   64%  x1  66.42 ms/call
  msp_particle_step     2.90 ms   1.74 Mcyc    3%
  msp_timeline_step     40 us    23.7 kcyc     0%
  msp_buffer_wait      34.91 ms  20.95 Mcyc   33%
```

Wall: min 57.4 / avg 104.3 / max 128.0 ms. This window straddles the cadence
boundary — early frames still clear one 62.5 ms window (16 fps) while later ones
overflow into a second (125 ms) as particles accumulate, which is why the wall
average sits between the two quantized cadences. `msp_particle_scan` is the sole
render cost (64 %); `msp_draw_particles` self time (particle setup above the
scan) is ~1 µs/frame. Physics (`msp_particle_step`) is already flat at ~3 %.

### Saturated pool — steady 8 fps (window: frames 129–192)

```
frame                 125.13 ms  75.08 Mcyc  100%
  msp_draw_particles   96.49 ms  57.89 Mcyc   77%
    msp_particle_scan  96.48 ms  57.89 Mcyc   77%  x1  96.48 ms/call
  msp_particle_step     3.52 ms   2.11 Mcyc    3%
  msp_timeline_step     44 us    26.2 kcyc     0%
  msp_buffer_wait      25.08 ms  15.05 Mcyc   20%
```

Wall: min 124.5 / avg 125.1 / max 127.2 ms — locked to the 2-window (8 fps)
cadence. The particle pool is full: `msp_particle_scan` alone is **77 % of the
frame** (96.48 ms), the render (`frame` − `buffer_wait` ≈ 100 ms) needs ~1.6
windows, and `msp_buffer_wait` is the 25 ms round-up idle. Physics stays
negligible (`msp_particle_step` ~3 %, `msp_timeline_step` <0.1 %) — the whole
cost is the single `ParticleSystem::draw` rasterize call. This holds for the
remainder of the pass (later windows frames 193–640 all read
`msp_particle_scan` 65–77 % of a steady ~125 ms frame, 81.9–95.3 ms; the swings
are pool-count jitter, not a phase change).

`Plot::ParticleSystem::draw` composites each particle sprite in place and does
not route through `filter_blend`, so there is no per-pixel blended-pixel counter
here and no per-pixel subsection — coverage is sparse and particle-count-driven,
not a fixed full-quadrant shade.

## Column-ISR / DMA marshaling cost (`build/profile_mindsplatter.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative saturated window (frames 129–192, 8.0 s). The
pack and submit run inside the flywheel wake, on the 1-in-8 wakes that render a
column. Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.63 / 2.25 / 174 us  cpu 4.15%
  isr_pack        2304/s  4.25 / 10.53 / 172 us  cpu 2.42%
  isr_dma_submit  2304/s  0.65 / 0.96 / 1.28 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, includes
pack + submit nested); `isr_pack` = the 72× `packPixel` marshal, once per
column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush + eDMA
kick).

- **The DMA kick itself is ~1 µs of CPU** (submit min 0.65 µs); the pixel-pack
  marshal around it is the real per-column cost. Whole per-column CPU ≈ **18 µs
  of the 434 µs column period (~4.2 %)**. The wire transfer runs asynchronously
  on eDMA and costs no CPU.
- The `min` floors (pack 4.25 µs / 2551 cyc, submit 0.65 µs) are the stable
  per-call costs; the avg/max inflate (pack avg 10.5 µs, max 172 µs) because the
  long 125 ms render frames leave more room for serial-dump and USB/DMA-
  completion preemption to land inside the scope — not real pack growth.
- Net: the ISR machinery steals **~4.2 % of the chip**. Against the 62.5 ms
  window that is ~2.6 ms, leaving ~60 ms of render budget per window. The
  saturated render (~100 ms) needs **~1.7×** that budget, hence the 2-window
  (125 ms) cadence.

## Summary ranking (saturated pool, share of the 125 ms frame)

1. `Plot::ParticleSystem::draw` (`msp_particle_scan`, particle-sprite rasterize
   + in-place composite) — **77 %** (96.48 ms) — the entire algorithmic cost.
2. display-window sync (`msp_buffer_wait`) — **20 %** (25.08 ms, idle by design).
3. particle physics advance (`msp_particle_step`) — 3 % (3.52 ms).
4. everything else (`msp_timeline_step`, `msp_draw_particles` self) — <1 %.

MindSplatter is ~100 % particle-rasterizer-bound: the physics `step` is small
and flat (1–3 %) across the whole pass, and the frame cost is set entirely by
how many particles the pool holds and how much of the 10,368-px quadrant they
cover. The only levers are the per-particle scan/composite math and the pool
size, not the simulation.

**vs the previous capture** (same setup, pre column-cull batch): the saturated
scan fell 107.75 → 96.48 ms (**−10.5 %**) and render 111.3 → 100.0 ms from the
column arc cull + planar-edge cull + deferred trail shading. The cadence is
unchanged — render still needs ~1.6 windows, so the saving re-emerges as
`msp_buffer_wait` idle (13.8 → 25.1 ms), consistent with the perf ledger's
call that the scan is at its structural floor for cadence purposes at -Os.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Plot::ParticleSystem::draw` composites in
  place, so the per-pixel blended-pixel counter never fires; coverage is sparse
  and particle-count-driven, so no per-pixel figure is derived.
- The frames 1–64 window is the cold-start empty pool (first frame `min` 501 µs)
  and is skipped when picking exactness/representative windows; per-window sums
  there under-count the eventual steady state.
- 80 s pass — under the 120 s effect epoch, so no epoch-boundary straddle
  windows in this capture.
- `-Os` build: the shipping Phantasm config; the `-O3` twin is
  `docs/profiles/O3/profile_mindsplatter_teensy_2026-07-14.md`.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/MindSplatter.h` (five scopes: `msp_buffer_wait`, `msp_timeline_step`,
  `msp_particle_step`, `msp_draw_particles`, `msp_particle_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=MindSplatter`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile MindSplatter [seconds]`.
