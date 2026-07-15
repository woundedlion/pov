# FlowField on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_flowfield_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `FlowField<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 70,796 B; ITCM (RAM1 code) 54,056 B; RAM2 free 4,736 B.
(At -Os: FLASH 47,764 / ITCM 32,248 — +48% / +67%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 129-192 root counter
4,796,553,233 cyc = 7,994,255.4 µs vs measured `micros()` window sum 7,994,262 µs
(Δ ≈ 0.8 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

FlowField is bound by its particle trail rasterizer — particle bursts that spawn, fan out, and decay to near-empty then re-burst -- there is **no fixed steady state**; the frame cost tracks the live particle population. Render spans
**13.1–111.1 ms** across the pass. `ff_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 129-192)

```
frame                 124.91 ms  74.95 Mcyc  100%
  ff_particle_draw      109.33 ms  65.60 Mcyc   87%
  ff_particle_step        1.70 ms  1.02 Mcyc    1%
  ff_timeline_step          26 us  15.9 kcyc    0%
  ff_buffer_wait         13.85 ms  8.31 Mcyc   11%
```

Wall: min 123.8 / avg 124.9 / max 126.0 ms.
Render (`frame` − `ff_buffer_wait`) = **111.06 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 641-704)

```
frame                  87.72 ms  52.63 Mcyc  100%
  ff_particle_draw       46.31 ms  27.79 Mcyc   52%
  ff_particle_step        1.70 ms  1.02 Mcyc    1%
  ff_timeline_step          22 us  13.0 kcyc    0%
  ff_buffer_wait         39.69 ms  23.81 Mcyc   45%
```

Wall: min 26.3 / avg 87.7 / max 126.8 ms.
Render = **48.04 ms**. The spread between this and the peak is
burst-population-driven, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       18420/s  0.50 / 1.77 / 54.1 us  cpu 3.27%
  isr_pack       2302/s  4.76 / 7.19 / 52.5 us  cpu 1.65%
  isr_dma_submit 2302/s  0.63 / 0.94 / 1.3 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`ff_particle_draw`** — 109.3 ms/frame, 88% of the peak frame; `ff_buffer_wait` is the display-sync idle by design.

Particle trails route through the changed `Plot::ParticleSystem::draw`. Because the population is not pinned across flashes, the -O3 ratio here is soft (the -O3 run caught a heavier burst); the per-fragment -O3 gain is ~1.3-1.45x.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/FlowField.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=FlowField`), + `tools/profile_capture.py`.
