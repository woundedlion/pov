# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_dreamballs.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `DreamBalls<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile DreamBalls` |

Image size: FLASH code 64,716 B; ITCM (RAM1 code) 36,328 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 257-320 root counter
5,833,341,373 cyc = 9,722,235.6 µs vs measured `micros()` window sum 9,722,244 µs
(Δ ≈ 0.9 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

DreamBalls is bound by its mesh wireframe rasterizer — a phase-cycling wireframe whose cost tracks the number of orbiting shell copies (`db_mesh_plot` calls/frame), stepping through the preset cycle. Render spans
**20.0–108.8 ms** across the pass. `db_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 257-320)

```
frame                 151.91 ms  91.15 Mcyc  100%
  db_timeline_step      108.76 ms  65.26 Mcyc   71%
    db_draw               108.70 ms  65.22 Mcyc   99%  x2  72465 us/call
      db_draw_scene         108.70 ms  65.22 Mcyc   99%  x2  72464 us/call
        db_mesh_plot          107.99 ms  64.79 Mcyc   99%  x21  5142 us/call
        db_warp_orient           333 us  200.0 kcyc    0%  x21  16 us/call
        db_displace              369 us  221.3 kcyc    0%  x21  18 us/call
      db_mesh_copy               1 us    575 cyc    0%  x2  1 us/call
  db_buffer_wait         43.15 ms  25.89 Mcyc   28%
```

Wall: min 122.9 / avg 151.9 / max 199.9 ms.
Render (`frame` − `db_buffer_wait`) = **108.76 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 577-640)

```
frame                  93.35 ms  56.01 Mcyc  100%
  db_timeline_step       63.24 ms  37.94 Mcyc   67%
    db_draw                63.18 ms  37.91 Mcyc   99%  x2  42120 us/call
      db_draw_scene          63.18 ms  37.91 Mcyc   99%  x2  42119 us/call
        db_mesh_plot           62.65 ms  37.59 Mcyc   99%  x9  6961 us/call
        db_warp_orient           251 us  150.8 kcyc    0%  x9  28 us/call
        db_displace              277 us  165.9 kcyc    0%  x9  31 us/call
      db_mesh_copy               1 us    664 cyc    0%  x2  1 us/call
  db_buffer_wait         30.11 ms  18.07 Mcyc   32%
```

Wall: min 50.9 / avg 93.3 / max 126.9 ms.
Render = **63.24 ms**. The spread between this and the peak is
phase/preset-driven (copy count varies), not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       22401/s  0.63 / 2.42 / 123.2 us  cpu 4.45%
  isr_pack       2800/s  4.37 / 11.88 / 121.1 us  cpu 2.73%
  isr_dma_submit 2800/s  0.63 / 0.97 / 1.3 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`db_mesh_plot`** — 108.0 ms/frame, 71% of the peak frame; `db_buffer_wait` is the display-sync idle by design.

**Real plot-cull beneficiary.** The wireframe's long geodesic/planar edges cross the quadrant, so the new per-edge column cull removes fragments in the wrong arm half: peak render fell ~145 -> ~109 ms vs the pre-batch capture, and -O3 pulls the heavy phase from ~5-8 fps to a steady 8 fps.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/DreamBalls.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=DreamBalls`), + `tools/profile_capture.py`.
