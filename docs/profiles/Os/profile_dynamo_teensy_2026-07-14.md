# Dynamo on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_dynamo.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `Dynamo<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile Dynamo` |

Image size: FLASH code 48,372 B; ITCM (RAM1 code) 33,208 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 193-256 root counter
7,280,557,854 cyc = 12,134,263.1 µs vs measured `micros()` window sum 12,134,270 µs
(Δ ≈ 0.6 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

Dynamo is bound by its feedback-buffer (Trails) composite — nodes drawn with `Plot::Line` into a Trails feedback buffer whose replay (`dy_filter_flush`) dominates; the plot draw is a minority of the frame. Render spans
**96.6–165.8 ms** across the pass. `dy_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 193-256)

```
frame                 189.60 ms  113.76 Mcyc  100%
  dy_filter_flush       140.33 ms  84.20 Mcyc   74%
  dy_draw_nodes          25.39 ms  15.23 Mcyc   13%
  dy_timeline_step          91 us  54.7 kcyc    0%
  dy_buffer_wait         23.79 ms  14.27 Mcyc   12%
```

Wall: min 117.9 / avg 189.6 / max 258.3 ms.
Render (`frame` − `dy_buffer_wait`) = **165.81 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 65-128)

```
frame                 175.00 ms  105.00 Mcyc  100%
  dy_filter_flush       109.51 ms  65.71 Mcyc   62%
  dy_draw_nodes          17.71 ms  10.63 Mcyc   10%
  dy_timeline_step         106 us  63.5 kcyc    0%
  dy_buffer_wait         47.66 ms  28.60 Mcyc   27%
```

Wall: min 121.4 / avg 175.0 / max 258.5 ms.
Render = **127.33 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       27958/s  0.62 / 2.67 / 112.7 us  cpu 4.91%
  isr_pack       3494/s  4.25 / 13.76 / 110.8 us  cpu 3.16%
  isr_dma_submit 3494/s  0.65 / 0.97 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`dy_filter_flush`** — 140.3 ms/frame, 74% of the peak frame; `dy_buffer_wait` is the display-sync idle by design.

**Marginal.** `dy_filter_flush` is ~74 % of the frame (untouched by the batch); the plot draw (`dy_draw_nodes`, `Plot::Line`) is ~13 %, so the column cull only shaves that slice.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Dynamo.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=Dynamo`), + `tools/profile_capture.py`.
