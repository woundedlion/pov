# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_mindsplatter_o3.log` (64-frame windows, single ~80 s pass,
captured from boot). This is the **-O3** twin of the shipping `-Os` report; the only
variable between the two is the optimization level. **Re-captured after the
2026-07-14 plot column-cull batch** (column arc cull 9ac8cebd, planar-edge cull
62450701, deferred trail shading 708d4b9b) — the saturated scan dropped
75.76 → 65.53 ms (−13.5 %) vs the previous capture of this report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: same Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but keeps base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forcing `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — the shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `MindSplatter<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR counters via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 78,340 B; ITCM (RAM1 code) 61,064 B; RAM2 free 4,736 B.
(At -Os: 53,092 / 36,520 / 4,736 — see the `-O3 vs -Os` section for the cost.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
saturated window frames 129–192 root counter = 4,798,033,274 cyc =
7,996,722.1 µs vs measured `micros()` window sum 7,996,724 µs (Δ ≈ 0.24 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Wall time snaps
up to a whole number of 62.5 ms windows, so render < 62.5 ms → 62.5 ms/frame
(16 fps); 62.5–125 ms → 125 ms/frame (8 fps).

MindSplatter is **particle-rasterizer-bound and ramps** as its particle pool
fills after the cold start: `Plot::ParticleSystem::draw` (`msp_particle_scan`)
grows from ~19 % to ~52 % of the frame, and wall climbs from ~61 ms to a
~125 ms steady state over the first ~130 frames. At -O3 the saturated render
(~68 ms) sits **just above** the one-window boundary, so the dominant cadence
is still 8 fps but no longer locked:

- **Cold start / low pool** (frames 1–64): render ~13 ms → fits one window →
  **62.5 ms/frame (16 fps)** (avg 61.0 ms wall; first frame `min` 481 µs is the
  empty pool, so this window is not used as an exactness/representative sample).
- **Growing pool** (frames 65–128): mostly still one-window (avg 68.0 ms wall vs
  104.3 ms at -Os in the same window — the faster scan holds 16 fps deeper into
  the ramp); the flip to 8 fps lands near the end of the window (max 127.2 ms).
- **Saturated pool** (frames 129 onward): render ~68 ms → 2 windows →
  **125 ms/frame (8 fps)** dominant, but when pool coverage momentarily dips the
  render slips under 62.5 ms and single-window frames reappear (window
  577–640: wall min 58.0 / avg 113.5 ms — a mix of 62.5 and 125 ms frames).

`msp_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle. At -O3
it is *larger* than at -Os (~57 ms in the saturated pool, vs ~25 ms at -Os):
the faster scan finishes well inside the 125 ms budget, so the reclaimed time
re-emerges as sync idle — see `-O3 vs -Os`.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Two regimes bound the ramp. Each block is nested as in the counter tree (each
parent includes its children). Columns: time/frame, cycles/frame, % of frame;
leaves add calls/frame and per-call cost.

### Early / growing pool — mostly 16 fps, flipping late (window: frames 65–128)

```
frame                  67.99 ms  40.79 Mcyc  100%
  msp_draw_particles   46.07 ms  27.64 Mcyc   68%
    msp_particle_scan  46.07 ms  27.64 Mcyc   68%  x1  46.07 ms/call
  msp_particle_step     1.82 ms   1.09 Mcyc    3%
  msp_timeline_step     27 us    16.3 kcyc     0%
  msp_buffer_wait      20.08 ms  12.05 Mcyc   30%
```

Wall: min 55.8 / avg 68.0 / max 127.2 ms. Unlike -Os (avg 104.3 ms in this same
window), most of these frames still clear one 62.5 ms window — the 1.4× faster
scan delays the cadence flip until late in the window. `msp_particle_scan` is
the sole render cost (68 %); `msp_draw_particles` self time (particle setup
above the scan) is ~1 µs/frame. Physics (`msp_particle_step`) is flat at ~3 %.

### Saturated pool — 8 fps dominant (window: frames 129–192)

```
frame                 124.95 ms  74.97 Mcyc  100%
  msp_draw_particles   65.54 ms  39.32 Mcyc   52%
    msp_particle_scan  65.53 ms  39.32 Mcyc   52%  x1  65.53 ms/call
  msp_particle_step     2.33 ms   1.40 Mcyc    2%
  msp_timeline_step     56 us    33.6 kcyc     0%
  msp_buffer_wait      57.03 ms  34.22 Mcyc   46%
```

Wall: min 123.7 / avg 124.9 / max 125.6 ms — this window is locked to the
2-window (8 fps) cadence. The particle pool is full: `msp_particle_scan` alone
is **52 % of the frame** (65.53 ms), the render (`frame` − `buffer_wait` ≈
67.9 ms) needs ~1.09 windows, and `msp_buffer_wait` absorbs the remaining 57 ms
as round-up idle. Physics stays negligible (`msp_particle_step` 2 %,
`msp_timeline_step` <0.1 %) — the whole cost is the single
`ParticleSystem::draw` rasterize call. Later windows read `msp_particle_scan`
50–63 % (63.7–78.4 ms) of a ~113–126 ms frame; the swings are pool-count
jitter, and in the lightest windows (513–576, 577–640; wall min 58–60 ms) some
frames dip below one window and run at 16 fps — the render is only ~5 ms above
the one-window budget, so -O3 sits right at the cadence edge.

`Plot::ParticleSystem::draw` composites each particle sprite in place and does
not route through `filter_blend`, so there is no per-pixel blended-pixel counter
here and no per-pixel subsection — coverage is sparse and particle-count-driven,
not a fixed full-quadrant shade.

## Column-ISR / DMA marshaling cost (`build/profile_mindsplatter_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative saturated window (frames 129–192, 8.0 s). The
pack and submit run inside the flywheel wake, on the 1-in-8 wakes that render a
column. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake         18433/s  0.49 / 2.21 / 122 us  cpu 4.06%
  isr_pack        2304/s  4.64 / 10.95 / 120 us  cpu 2.52%
  isr_dma_submit  2304/s  0.66 / 0.91 / 1.21 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, includes
pack + submit nested); `isr_pack` = the 72× `packPixel` marshal, once per
column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush + eDMA
kick).

- **The DMA kick itself is ~1 µs of CPU** (submit min 0.66 µs); the pixel-pack
  marshal around it is the real per-column cost. Whole per-column CPU ≈ **18 µs
  of the 434 µs column period (~4.1 %)**. The wire transfer runs asynchronously
  on eDMA and costs no CPU.
- The `min` floors (pack 4.64 µs / 2786 cyc, submit 0.66 µs) are the stable
  per-call costs; the avg/max inflate (pack avg 10.9 µs, max 120 µs) because the
  long 125 ms render frames leave more room for serial-dump and USB/DMA-
  completion preemption to land inside the scope — not real pack growth.
- Net: the ISR machinery steals **~4.1 % of the chip** (unchanged from -Os — it
  is display-driver work, not effect code). Against the 62.5 ms window that is
  ~2.5 ms, leaving ~60 ms of render budget per window. The saturated render
  (~68 ms) needs **~1.1×** that budget — the closest any MindSplatter config has
  come to one window, but still over it, hence the dominant 2-window cadence.

## Summary ranking (saturated pool, share of the 125 ms frame)

1. `Plot::ParticleSystem::draw` (`msp_particle_scan`, particle-sprite rasterize
   + in-place composite) — **52 %** (65.53 ms) — the entire algorithmic cost.
2. display-window sync (`msp_buffer_wait`) — **46 %** (57.03 ms, idle by design).
3. particle physics advance (`msp_particle_step`) — 2 % (2.33 ms).
4. everything else (`msp_timeline_step`, `msp_draw_particles` self) — <1 %.

MindSplatter is ~100 % particle-rasterizer-bound: the physics `step` is small
and flat (1–3 %) across the whole pass, and the frame cost is set entirely by
how many particles the pool holds and how much of the 10,368-px quadrant they
cover. -O3 shrinks the scan ~1.47× but does not change the shape. What HAS
changed since the previous capture: with the column-cull batch on top of -O3,
the saturated render (~68 ms) is now only **~5 ms above the one-window budget**
— light windows already flicker to 16 fps. Another ~8 % off the scan at -O3
would flip MindSplatter to a steady 16 fps; at -Os (render ~100 ms) that tier
remains far out of reach.

## -O3 vs -Os

Both builds keep the same **dominant 125 ms / 8 fps** saturated cadence — but
at -O3 the render (~68 ms) barely exceeds one 62.5 ms window, so light frames
now dip to 16 fps, where -Os (render ~100 ms) is locked at 8. The rest of the
saving only widens `msp_buffer_wait` idle (25.08 → 57.03 ms). Per-scope, -O3 is
~1.47× faster (representative saturated window 129–192):

| metric | -Os | -O3 | speedup |
|---|---|---|---|
| `msp_particle_scan` (ms/frame) | 96.48 | 65.53 | **1.47×** (−30.9 ms) |
| `msp_draw_particles` (ms/frame) | 96.49 | 65.54 | **1.47×** (−30.9 ms) |
| render (`frame` − `buffer_wait`, ms) | 100.0 | 67.9 | **1.47×** (−32.1 ms) |
| `msp_particle_step` (ms/frame) | 3.52 | 2.33 | 1.51× |
| early-pool `msp_particle_scan` (65–128, ms) | 66.42 | 46.07 | 1.44×* |
| `msp_buffer_wait` (ms/frame, idle) | 25.08 | 57.03 | absorbs the saving |
| frame wall / saturated cadence | 125 ms / 8 fps | 125 ms / 8 fps (occasional 16) | edge, not tier |
| FLASH code | 53,092 B | 78,340 B | **+25,248 B (+47.6 %)** |
| ITCM (RAM1 code) | 36,520 B | 61,064 B | **+24,544 B (+67.2 %)** |

\* the early-pool comparison is approximate: the -O3 ramp holds 16 fps deeper
into the window, so the two builds sample slightly different pool-fill states
across frames 65–128.

The single hot loop (`ParticleSystem::draw`) is the whole story and it gains a
flat ~1.47× at -O3 (particle physics gains ~1.5×). Post column-cull batch, the
win is no longer entirely invisible on device: the saturated render is ~5 ms
over the one-window budget, so the lightest pool states already produce
single-window (16 fps) frames. A steady 16 fps at -O3 needs only ~8 % more off
the scan; at -Os the same tier needs ~38 % and stays out of reach.

**-O3 is not the shipping config.** The full 26-effect Phantasm image overflows
FlexRAM/ITCM at -O3 (a single effect already costs +67 % ITCM here), which is
exactly why the ship build is `-Os`. -O3 is only viable as this single-effect
profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Plot::ParticleSystem::draw` composites in
  place, so the per-pixel blended-pixel counter never fires; coverage is sparse
  and particle-count-driven, so no per-pixel figure is derived.
- The frames 1–64 window is the cold-start empty pool (first frame `min` 481 µs)
  and is skipped when picking exactness/representative windows; per-window sums
  there under-count the eventual steady state.
- 80 s pass — under the 120 s effect epoch, so no epoch-boundary straddle
  windows in this capture.
- **-O3 build**: this is the `profile_o3` env, which does **not** ship. The
  shipping Phantasm image is `-Os` because the full 26-effect roster overflows
  ITCM at -O3.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/MindSplatter.h` (five scopes: `msp_buffer_wait`, `msp_timeline_step`,
  `msp_particle_step`, `msp_draw_particles`, `msp_particle_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=MindSplatter`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py`
  (the `-Os` twin uses env `profile`, reproducible via `just profile MindSplatter`).
