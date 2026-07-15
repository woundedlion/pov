# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_shapeshifter_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `ShapeShifter<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 87,212 B; ITCM (RAM1 code) 73,144 B; RAM2 free 4,736 B.
(At -Os: FLASH 58,268 / ITCM 44,568 — +49% / +64%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 65-128 root counter
3,556,427,776 cyc = 5,927,379.6 µs vs measured `micros()` window sum 5,927,388 µs
(Δ ≈ 1.4 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

ShapeShifter is bound by its SDF scan rasterizer (mostly `Scan::`, minority `Plot::`) — a morphing SDF solid drawn mostly through `Scan::` (`ss_scan_dispatch`, with a per-pixel `filter_blend` leaf) plus a `Plot::` polygon path (`ss_plot_dispatch`); phase-varying as the shape morphs. Render spans
**14.5–50.9 ms** across the pass. `ss_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 65-128)

```
frame                  92.62 ms  55.57 Mcyc  100%
  ss_timeline_step       50.87 ms  30.52 Mcyc   54%
    ss_draw_all            50.82 ms  30.49 Mcyc   99%  x1  50816 us/call
      ss_scan_dispatch       44.31 ms  26.58 Mcyc   87%  x7  6330 us/call
        filter_blend            5.10 ms  3.06 Mcyc   11%  x20972  0 us/call
      ss_plot_dispatch        6.50 ms  3.90 Mcyc   12%  x7  929 us/call
  ss_buffer_wait         41.75 ms  25.05 Mcyc   45%
```

Wall: min 48.4 / avg 92.6 / max 125.8 ms.
Render (`frame` − `ss_buffer_wait`) = **50.87 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 129-192)

```
frame                  62.67 ms  37.60 Mcyc  100%
  ss_timeline_step       32.09 ms  19.26 Mcyc   51%
    ss_draw_all            32.06 ms  19.24 Mcyc   99%  x1  32059 us/call
      ss_scan_dispatch       26.27 ms  15.76 Mcyc   81%  x7  3753 us/call
        filter_blend            1.86 ms  1.12 Mcyc    7%  x9133  0 us/call
      ss_plot_dispatch        5.78 ms  3.47 Mcyc   18%  x7  826 us/call
  ss_buffer_wait         30.57 ms  18.34 Mcyc   48%
```

Wall: min 26.8 / avg 62.7 / max 98.0 ms.
Render = **32.10 ms**. The spread between this and the peak is
reveal/morph-driven, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       13658/s  0.55 / 2.00 / 82.0 us  cpu 3.67%
  isr_pack       1707/s  4.76 / 9.15 / 80.3 us  cpu 2.10%
  isr_dma_submit 1707/s  0.61 / 0.94 / 6.4 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`ss_scan_dispatch`** — 44.3 ms/frame, 48% of the peak frame; `ss_buffer_wait` is the display-sync idle by design.

Mostly `Scan::`-bound: `ss_scan_dispatch` is ~86 % of the draw and only `ss_plot_dispatch` (12-18 %) is on the changed plot path, so the batch effect is small. The -O3 gain (incl. `filter_blend` ~1.7x) is the compiler.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/ShapeShifter.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=ShapeShifter`), + `tools/profile_capture.py`.
