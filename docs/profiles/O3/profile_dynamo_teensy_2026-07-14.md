# Dynamo on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_dynamo_o3.log` (64-frame windows, ~75 s single pass). This is the **-O3** twin of the shipping `-Os` report (`../Os/profile_dynamo_teensy_2026-07-14.md`); the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `Dynamo<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 73,820 B; ITCM (RAM1 code) 58,136 B; RAM2 free 4,736 B.
(At -Os: FLASH 48,372 / ITCM 33,208 — +52% / +75%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 449-512 root counter
5,866,121,183 cyc = 9,776,868.6 µs vs measured `micros()` window sum 9,776,873 µs
(Δ ≈ 0.4 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

Dynamo is bound by its feedback-buffer (Trails) composite — nodes drawn with `Plot::Line` into a Trails feedback buffer whose replay (`dy_filter_flush`) dominates; the plot draw is a minority of the frame. Render spans
**65.4–112.9 ms** across the pass. `dy_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 449-512)

```
frame                 152.76 ms  91.66 Mcyc  100%
  dy_filter_flush        92.07 ms  55.24 Mcyc   60%
  dy_draw_nodes          20.79 ms  12.48 Mcyc   13%
  dy_timeline_step          28 us  17.0 kcyc    0%
  dy_buffer_wait         39.87 ms  23.92 Mcyc   26%
```

Wall: min 119.8 / avg 152.8 / max 193.3 ms.
Render (`frame` − `dy_buffer_wait`) = **112.89 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 65-128)

```
frame                 111.78 ms  67.07 Mcyc  100%
  dy_filter_flush        74.88 ms  44.93 Mcyc   66%
  dy_draw_nodes          11.74 ms  7.04 Mcyc   10%
  dy_timeline_step          53 us  31.7 kcyc    0%
  dy_buffer_wait         25.11 ms  15.07 Mcyc   22%
```

Wall: min 59.2 / avg 111.8 / max 193.0 ms.
Render = **86.67 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       22527/s  0.49 / 2.98 / 120.7 us  cpu 5.48%
  isr_pack       2815/s  4.77 / 16.92 / 118.8 us  cpu 3.89%
  isr_dma_submit 2815/s  0.62 / 0.94 / 1.2 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`dy_filter_flush`** — 92.1 ms/frame, 60% of the peak frame; `dy_buffer_wait` is the display-sync idle by design.

**Marginal.** `dy_filter_flush` is ~74 % of the frame (untouched by the batch); the plot draw (`dy_draw_nodes`, `Plot::Line`) is ~13 %, so the column cull only shaves that slice.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Dynamo.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=Dynamo`), + `tools/profile_capture.py`.
