# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_hopffibration.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `HopfFibration<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile HopfFibration` |

Image size: FLASH code 39,188 B; ITCM (RAM1 code) 26,344 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 193-256 root counter
4,283,368,041 cyc = 7,138,946.7 µs vs measured `micros()` window sum 7,138,958 µs
(Δ ≈ 1.6 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

HopfFibration is bound by its trail rasterizer — fibre trails rasterized through `Plot::rasterize` (`hf_trail_raster`, fixed x210 trail segments/frame). Render spans
**30.8–70.2 ms** across the pass. `hf_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 193-256)

```
frame                 111.55 ms  66.93 Mcyc  100%
  hf_render_trails       69.99 ms  41.99 Mcyc   62%
    hf_trail_raster        68.34 ms  41.00 Mcyc   97%  x210  325 us/call
  hf_project_record        235 us  140.7 kcyc    0%
  hf_advance_tumble          1 us    363 cyc    0%
  hf_timeline_step          12 us   7.3 kcyc    0%
  hf_buffer_wait         41.31 ms  24.79 Mcyc   37%
```

Wall: min 42.4 / avg 111.5 / max 128.6 ms.
Render (`frame` − `hf_buffer_wait`) = **70.24 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 641-704)

```
frame                  62.58 ms  37.55 Mcyc  100%
  hf_render_trails       49.78 ms  29.87 Mcyc   79%
    hf_trail_raster        48.21 ms  28.92 Mcyc   96%  x210  230 us/call
  hf_project_record        234 us  140.5 kcyc    0%
  hf_advance_tumble          1 us    373 cyc    0%
  hf_timeline_step          11 us   6.5 kcyc    0%
  hf_buffer_wait         12.55 ms  7.53 Mcyc   20%
```

Wall: min 34.9 / avg 62.6 / max 89.4 ms.
Render = **50.03 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       16448/s  0.52 / 3.35 / 137.6 us  cpu 6.16%
  isr_pack       2056/s  4.37 / 19.41 / 135.7 us  cpu 4.47%
  isr_dma_submit 2056/s  0.63 / 0.97 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`hf_trail_raster`** — 68.3 ms/frame, 61% of the peak frame; `hf_buffer_wait` is the display-sync idle by design.

Trails route through `Plot::rasterize`, but they are short and near-meridian, so they hit `edge_col_span`'s near-meridian no-cull fallback: the **median render is unchanged** from the pre-batch capture. -O3 (the compiler, not the cull) still crosses the peak from 8 to 16 fps.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/HopfFibration.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=HopfFibration`), + `tools/profile_capture.py`.
