# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_meshfeedback_o3.log` (64-frame windows, ~75 s single pass). This is the **-O3** twin of the shipping `-Os` report (`../Os/profile_meshfeedback_teensy_2026-07-14.md`); the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `MeshFeedback<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 95,948 B; ITCM (RAM1 code) 65,304 B; RAM2 free 4,736 B.
(At -Os: FLASH 64,516 / ITCM 42,744 — +48% / +52%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 385-448 root counter
4,801,318,124 cyc = 8,002,196.9 µs vs measured `micros()` window sum 8,002,203 µs
(Δ ≈ 0.8 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

MeshFeedback is bound by its feedback-buffer composite (NOT the plot draw) — a mesh morph drawn once per frame into a feedback buffer that is then replayed (`mf_feedback_flush`); the flush, not the plot draw, is the whole cost. Render spans
**38.7–69.4 ms** across the pass. `mf_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 385-448)

```
frame                 125.03 ms  75.02 Mcyc  100%
  mf_timeline_step        7.04 ms  4.23 Mcyc    5%
    mf_morph_draw           7.01 ms  4.20 Mcyc   99%  x2  3504 us/call
  mf_feedback_flush      62.31 ms  37.39 Mcyc   49%
  mf_apply_params            0 us     38 cyc    0%
  mf_buffer_wait         55.67 ms  33.40 Mcyc   44%
```

Wall: min 117.9 / avg 125.0 / max 130.2 ms.
Render (`frame` − `mf_buffer_wait`) = **69.36 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 129-192)

```
frame                  80.16 ms  48.10 Mcyc  100%
  mf_timeline_step        4.87 ms  2.92 Mcyc    6%
    mf_morph_draw           4.83 ms  2.90 Mcyc   99%  x2  2415 us/call
  mf_feedback_flush      55.83 ms  33.50 Mcyc   69%
  mf_apply_params            0 us     39 cyc    0%
  mf_buffer_wait         19.47 ms  11.68 Mcyc   24%
```

Wall: min 59.7 / avg 80.2 / max 126.0 ms.
Render = **60.69 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       18435/s  0.55 / 4.93 / 124.6 us  cpu 9.08%
  isr_pack       2304/s  4.85 / 32.52 / 123.0 us  cpu 7.49%
  isr_dma_submit 2304/s  0.75 / 0.96 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`mf_feedback_flush`** — 62.3 ms/frame, 50% of the peak frame; `mf_buffer_wait` is the display-sync idle by design.

**Not a plot-cull beneficiary.** `Pixel::Feedback::flush` is ~75 % of the frame and the batch does not touch it; the actual plot draw (`mf_morph_draw`) is only 5-7 %. The -O3 speedup is the compiler on the flush loop.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/MeshFeedback.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=MeshFeedback`), + `tools/profile_capture.py`.
