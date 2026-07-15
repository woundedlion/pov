# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_shapeshifter.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `ShapeShifter<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile ShapeShifter` |

Image size: FLASH code 58,268 B; ITCM (RAM1 code) 44,568 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 65-128 root counter
3,544,282,734 cyc = 5,907,137.9 µs vs measured `micros()` window sum 5,907,147 µs
(Δ ≈ 1.5 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

ShapeShifter is bound by its SDF scan rasterizer (mostly `Scan::`, minority `Plot::`) — a morphing SDF solid drawn mostly through `Scan::` (`ss_scan_dispatch`, with a per-pixel `filter_blend` leaf) plus a `Plot::` polygon path (`ss_plot_dispatch`); phase-varying as the shape morphs. Render spans
**18.9–67.8 ms** across the pass. `ss_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 65-128)

```
frame                  92.30 ms  55.38 Mcyc  100%
  ss_timeline_step       67.77 ms  40.66 Mcyc   73%
    ss_draw_all            67.72 ms  40.63 Mcyc   99%  x1  67722 us/call
      ss_scan_dispatch       58.27 ms  34.96 Mcyc   86%  x7  8324 us/call
        filter_blend            8.50 ms  5.10 Mcyc   14%  x20972  0 us/call
      ss_plot_dispatch        9.45 ms  5.67 Mcyc   13%  x7  1350 us/call
  ss_buffer_wait         24.52 ms  14.71 Mcyc   26%
```

Wall: min 24.7 / avg 92.3 / max 125.8 ms.
Render (`frame` − `ss_buffer_wait`) = **67.77 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 385-448)

```
frame                  78.50 ms  47.10 Mcyc  100%
  ss_timeline_step       50.66 ms  30.39 Mcyc   64%
    ss_draw_all            50.61 ms  30.36 Mcyc   99%  x1  50606 us/call
      ss_scan_dispatch       42.90 ms  25.74 Mcyc   84%  x7  6129 us/call
        filter_blend            6.67 ms  4.00 Mcyc   15%  x16735  0 us/call
      ss_plot_dispatch        7.70 ms  4.62 Mcyc   15%  x7  1100 us/call
  ss_buffer_wait         27.84 ms  16.71 Mcyc   35%
```

Wall: min 32.4 / avg 78.5 / max 127.2 ms.
Render = **50.66 ms**. The spread between this and the peak is
reveal/morph-driven, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       13611/s  0.63 / 2.04 / 86.3 us  cpu 3.75%
  isr_pack       1701/s  4.37 / 8.67 / 84.4 us  cpu 1.99%
  isr_dma_submit 1701/s  0.65 / 0.96 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`ss_scan_dispatch`** — 58.3 ms/frame, 63% of the peak frame; `ss_buffer_wait` is the display-sync idle by design.

Mostly `Scan::`-bound: `ss_scan_dispatch` is ~86 % of the draw and only `ss_plot_dispatch` (12-18 %) is on the changed plot path, so the batch effect is small. The -O3 gain (incl. `filter_blend` ~1.7x) is the compiler.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/ShapeShifter.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=ShapeShifter`), + `tools/profile_capture.py`.
