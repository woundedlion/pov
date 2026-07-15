# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_meshfeedback.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `MeshFeedback<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile MeshFeedback` |

Image size: FLASH code 64,516 B; ITCM (RAM1 code) 42,744 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 385-448 root counter
5,737,746,906 cyc = 9,562,911.5 µs vs measured `micros()` window sum 9,562,921 µs
(Δ ≈ 1.0 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

MeshFeedback is bound by its feedback-buffer composite (NOT the plot draw) — a mesh morph drawn once per frame into a feedback buffer that is then replayed (`mf_feedback_flush`); the flush, not the plot draw, is the whole cost. Render spans
**92.9–124.4 ms** across the pass. `mf_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 385-448)

```
frame                 149.42 ms  89.65 Mcyc  100%
  mf_timeline_step       11.14 ms  6.68 Mcyc    7%
    mf_morph_draw          11.09 ms  6.65 Mcyc   99%  x2  5545 us/call
  mf_feedback_flush     113.28 ms  67.97 Mcyc   75%
  mf_apply_params            2 us   1.2 kcyc    0%
  mf_buffer_wait         25.00 ms  15.00 Mcyc   16%
```

Wall: min 120.6 / avg 149.4 / max 207.4 ms.
Render (`frame` − `mf_buffer_wait`) = **124.42 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 129-192)

```
frame                 141.72 ms  85.03 Mcyc  100%
  mf_timeline_step        7.72 ms  4.63 Mcyc    5%
    mf_morph_draw           7.67 ms  4.60 Mcyc   99%  x2  3834 us/call
  mf_feedback_flush     111.31 ms  66.78 Mcyc   78%
  mf_apply_params            0 us    116 cyc    0%
  mf_buffer_wait         22.68 ms  13.61 Mcyc   16%
```

Wall: min 121.3 / avg 141.7 / max 189.7 ms.
Render = **119.03 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       22023/s  0.54 / 5.09 / 211.9 us  cpu 9.37%
  isr_pack       2754/s  4.70 / 33.04 / 210.3 us  cpu 7.61%
  isr_dma_submit 2754/s  0.64 / 0.98 / 1.3 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`mf_feedback_flush`** — 113.3 ms/frame, 76% of the peak frame; `mf_buffer_wait` is the display-sync idle by design.

**Not a plot-cull beneficiary.** `Pixel::Feedback::flush` is ~75 % of the frame and the batch does not touch it; the actual plot draw (`mf_morph_draw`) is only 5-7 %. The -O3 speedup is the compiler on the flush loop.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/MeshFeedback.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=MeshFeedback`), + `tools/profile_capture.py`.
