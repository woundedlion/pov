# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_hopffibration_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `HopfFibration<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 66,444 B; ITCM (RAM1 code) 52,040 B; RAM2 free 4,736 B.
(At -Os: FLASH 39,188 / ITCM 26,344 — +69% / +97%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 193-256 root counter
2,381,816,133 cyc = 3,969,693.6 µs vs measured `micros()` window sum 3,969,700 µs
(Δ ≈ 1.6 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

HopfFibration is bound by its trail rasterizer — fibre trails rasterized through `Plot::rasterize` (`hf_trail_raster`, fixed x210 trail segments/frame). Render spans
**21.3–44.1 ms** across the pass. `hf_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 193-256)

```
frame                  62.03 ms  37.22 Mcyc  100%
  hf_render_trails       43.94 ms  26.36 Mcyc   70%
    hf_trail_raster        43.03 ms  25.82 Mcyc   97%  x210  205 us/call
  hf_project_record        160 us  96.1 kcyc    0%
  hf_advance_tumble          0 us    260 cyc    0%
  hf_timeline_step          11 us   6.7 kcyc    0%
  hf_buffer_wait         17.92 ms  10.75 Mcyc   28%
```

Wall: min 43.1 / avg 62.0 / max 80.6 ms.
Render (`frame` − `hf_buffer_wait`) = **44.11 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 769-832)

```
frame                  63.12 ms  37.87 Mcyc  100%
  hf_render_trails       33.55 ms  20.13 Mcyc   53%
    hf_trail_raster        32.64 ms  19.58 Mcyc   97%  x210  155 us/call
  hf_project_record        162 us  97.4 kcyc    0%
  hf_advance_tumble          1 us    337 cyc    0%
  hf_timeline_step          10 us   6.0 kcyc    0%
  hf_buffer_wait         29.40 ms  17.64 Mcyc   46%
```

Wall: min 35.1 / avg 63.1 / max 91.7 ms.
Render = **33.72 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9147/s  0.41 / 1.56 / 107.0 us  cpu 2.87%
  isr_pack       1143/s  4.65 / 5.96 / 105.4 us  cpu 1.37%
  isr_dma_submit 1143/s  0.60 / 0.92 / 1.3 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`hf_trail_raster`** — 43.0 ms/frame, 69% of the peak frame; `hf_buffer_wait` is the display-sync idle by design.

Trails route through `Plot::rasterize`, but they are short and near-meridian, so they hit `edge_col_span`'s near-meridian no-cull fallback: the **median render is unchanged** from the pre-batch capture. -O3 (the compiler, not the cull) still crosses the peak from 8 to 16 fps.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/HopfFibration.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=HopfFibration`), + `tools/profile_capture.py`.
