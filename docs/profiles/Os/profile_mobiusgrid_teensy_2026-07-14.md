# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_mobiusgrid.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `MobiusGrid<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile MobiusGrid` |

Image size: FLASH code 50,604 B; ITCM (RAM1 code) 33,752 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 257-320 root counter
2,403,612,729 cyc = 4,006,021.2 µs vs measured `micros()` window sum 4,006,029 µs
(Δ ≈ 1.9 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

MobiusGrid is bound by its curve rasterizer (idle-bound at 16 fps) — a Mobius grid of great-circle lines and rings drawn via `Plot::rasterize`; the grid reveals/wipes so render ramps 2->18 ms, always inside one window (16 fps). Render spans
**2.4–17.7 ms** across the pass. `mg_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 257-320)

```
frame                  62.59 ms  37.56 Mcyc  100%
  mg_draw_grid           17.14 ms  10.29 Mcyc   27%
    mg_lines_draw          10.33 ms  6.20 Mcyc   60%  x1  10330 us/call
    mg_rings_draw           6.81 ms  4.09 Mcyc   39%  x1  6813 us/call
  mg_wipe_rebake           485 us  290.9 kcyc    0%
  mg_timeline_step          78 us  46.6 kcyc    0%
  mg_buffer_wait         44.88 ms  26.93 Mcyc   71%
```

Wall: min 56.0 / avg 62.6 / max 69.6 ms.
Render (`frame` − `mg_buffer_wait`) = **17.71 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 513-576)

```
frame                  62.56 ms  37.54 Mcyc  100%
  mg_draw_grid            8.84 ms  5.30 Mcyc   14%
    mg_lines_draw           5.11 ms  3.07 Mcyc   57%  x1  5109 us/call
    mg_rings_draw           3.73 ms  2.24 Mcyc   42%  x1  3729 us/call
  mg_wipe_rebake           311 us  186.8 kcyc    0%
  mg_timeline_step          66 us  39.3 kcyc    0%
  mg_buffer_wait         53.34 ms  32.01 Mcyc   85%
```

Wall: min 52.9 / avg 62.6 / max 72.5 ms.
Render = **9.22 ms**. The spread between this and the peak is
reveal/morph-driven, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9231/s  0.63 / 1.56 / 85.7 us  cpu 2.87%
  isr_pack       1153/s  4.37 / 5.26 / 83.9 us  cpu 1.21%
  isr_dma_submit 1153/s  0.63 / 0.95 / 1.9 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`mg_draw_grid`** — 17.1 ms/frame, 27% of the peak frame; `mg_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps. Line/ring edges route through `Plot::rasterize` and the cull shaves the peak reveal modestly, but the effect never approaches one window, so there is no cadence change.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/MobiusGrid.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=MobiusGrid`), + `tools/profile_capture.py`.
