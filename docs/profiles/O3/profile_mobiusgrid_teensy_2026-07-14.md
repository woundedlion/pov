# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_mobiusgrid_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `MobiusGrid<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 73,724 B; ITCM (RAM1 code) 55,832 B; RAM2 free 4,736 B.
(At -Os: FLASH 50,604 / ITCM 33,752 — +45% / +65%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 257-320 root counter
2,402,530,080 cyc = 4,004,216.8 µs vs measured `micros()` window sum 4,004,220 µs
(Δ ≈ 0.8 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

MobiusGrid is bound by its curve rasterizer (idle-bound at 16 fps) — a Mobius grid of great-circle lines and rings drawn via `Plot::rasterize`; the grid reveals/wipes so render ramps 2->18 ms, always inside one window (16 fps). Render spans
**1.6–12.4 ms** across the pass. `mg_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 257-320)

```
frame                  62.57 ms  37.54 Mcyc  100%
  mg_draw_grid           12.08 ms  7.25 Mcyc   19%
    mg_lines_draw           7.33 ms  4.40 Mcyc   60%  x1  7329 us/call
    mg_rings_draw           4.75 ms  2.85 Mcyc   39%  x1  4753 us/call
  mg_wipe_rebake           271 us  162.4 kcyc    0%
  mg_timeline_step          72 us  42.9 kcyc    0%
  mg_buffer_wait         50.14 ms  30.08 Mcyc   80%
```

Wall: min 57.9 / avg 62.6 / max 67.5 ms.
Render (`frame` − `mg_buffer_wait`) = **12.43 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 833-896)

```
frame                  62.59 ms  37.55 Mcyc  100%
  mg_draw_grid            6.01 ms  3.61 Mcyc    9%
    mg_lines_draw           3.64 ms  2.18 Mcyc   60%  x1  3637 us/call
    mg_rings_draw           2.38 ms  1.43 Mcyc   39%  x1  2376 us/call
  mg_wipe_rebake           348 us  209.1 kcyc    0%
  mg_timeline_step          86 us  51.4 kcyc    0%
  mg_buffer_wait         56.14 ms  33.68 Mcyc   89%
```

Wall: min 58.5 / avg 62.6 / max 66.6 ms.
Render = **6.45 ms**. The spread between this and the peak is
reveal/morph-driven, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9226/s  0.50 / 1.52 / 79.8 us  cpu 2.79%
  isr_pack       1153/s  4.76 / 5.71 / 78.2 us  cpu 1.31%
  isr_dma_submit 1153/s  0.62 / 0.92 / 5.7 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`mg_draw_grid`** — 12.1 ms/frame, 19% of the peak frame; `mg_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps. Line/ring edges route through `Plot::rasterize` and the cull shaves the peak reveal modestly, but the effect never approaches one window, so there is no cadence change.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/MobiusGrid.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=MobiusGrid`), + `tools/profile_capture.py`.
