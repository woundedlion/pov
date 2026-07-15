# ChaoticStrings on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_chaoticstrings.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `ChaoticStrings<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile ChaoticStrings` |

Image size: FLASH code 46,860 B; ITCM (RAM1 code) 32,472 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 833-896 root counter
2,399,619,433 cyc = 3,999,365.7 µs vs measured `micros()` window sum 3,999,372 µs
(Δ ≈ 1.6 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

ChaoticStrings is bound by its multiline rasterizer (idle-bound at 16 fps) — string vertices built then drawn as `Plot::Multiline`; render stays well under one 62.5 ms window, so the effect is display-sync-idle-bound at 16 fps. Render spans
**16.1–20.0 ms** across the pass. `cs_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 833-896)

```
frame                  62.49 ms  37.49 Mcyc  100%
  cs_multiline_draw      14.09 ms  8.46 Mcyc   22%
  cs_build_vertices       5.77 ms  3.46 Mcyc    9%
  cs_noise_prepare           0 us     87 cyc    0%
  cs_timeline_step         122 us  73.3 kcyc    0%
  cs_buffer_wait         42.50 ms  25.50 Mcyc   68%
```

Wall: min 59.4 / avg 62.5 / max 65.7 ms.
Render (`frame` − `cs_buffer_wait`) = **19.99 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 129-192)

```
frame                  62.52 ms  37.51 Mcyc  100%
  cs_multiline_draw      13.66 ms  8.20 Mcyc   21%
  cs_build_vertices       5.76 ms  3.45 Mcyc    9%
  cs_noise_prepare           0 us     87 cyc    0%
  cs_timeline_step         125 us  75.0 kcyc    0%
  cs_buffer_wait         42.98 ms  25.79 Mcyc   68%
```

Wall: min 58.7 / avg 62.5 / max 66.6 ms.
Render = **19.55 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9215/s  0.58 / 1.55 / 77.9 us  cpu 2.86%
  isr_pack       1152/s  4.37 / 5.29 / 76.2 us  cpu 1.21%
  isr_dma_submit 1152/s  0.62 / 0.96 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`cs_multiline_draw`** — 14.1 ms/frame, 23% of the peak frame; `cs_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps throughout. The multiline column cull shaves `cs_multiline_draw` modestly (the render stays far under one window), so the saving is invisible on cadence -- it widens `cs_buffer_wait` idle.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/ChaoticStrings.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=ChaoticStrings`), + `tools/profile_capture.py`.
