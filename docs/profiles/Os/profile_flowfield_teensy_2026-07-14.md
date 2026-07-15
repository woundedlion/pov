# FlowField on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_flowfield.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `FlowField<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile FlowField` |

Image size: FLASH code 47,764 B; ITCM (RAM1 code) 32,248 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 705-768 root counter
4,802,961,694 cyc = 8,004,936.2 µs vs measured `micros()` window sum 8,004,942 µs
(Δ ≈ 0.7 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

FlowField is bound by its particle trail rasterizer — particle bursts that spawn, fan out, and decay to near-empty then re-burst -- there is **no fixed steady state**; the frame cost tracks the live particle population. Render spans
**19.1–120.7 ms** across the pass. `ff_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 705-768)

```
frame                 125.08 ms  75.05 Mcyc  100%
  ff_particle_draw      118.37 ms  71.02 Mcyc   94%
  ff_particle_step        2.25 ms  1.35 Mcyc    1%
  ff_timeline_step          36 us  21.4 kcyc    0%
  ff_buffer_wait          4.42 ms  2.65 Mcyc    3%
```

Wall: min 123.8 / avg 125.1 / max 126.3 ms.
Render (`frame` − `ff_buffer_wait`) = **120.66 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 65-128)

```
frame                  96.00 ms  57.60 Mcyc  100%
  ff_particle_draw       63.45 ms  38.07 Mcyc   66%
  ff_particle_step        2.19 ms  1.32 Mcyc    2%
  ff_timeline_step          27 us  16.0 kcyc    0%
  ff_buffer_wait         30.33 ms  18.20 Mcyc   31%
```

Wall: min 42.0 / avg 96.0 / max 205.5 ms.
Render = **65.67 ms**. The spread between this and the peak is
burst-population-driven, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       18445/s  0.52 / 1.89 / 87.5 us  cpu 3.48%
  isr_pack       2305/s  4.37 / 7.51 / 85.7 us  cpu 1.73%
  isr_dma_submit 2305/s  0.63 / 0.96 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`ff_particle_draw`** — 118.4 ms/frame, 95% of the peak frame; `ff_buffer_wait` is the display-sync idle by design.

Particle trails route through the changed `Plot::ParticleSystem::draw`. Because the population is not pinned across flashes, the -O3 ratio here is soft (the -O3 run caught a heavier burst); the per-fragment -O3 gain is ~1.3-1.45x.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/FlowField.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=FlowField`), + `tools/profile_capture.py`.
