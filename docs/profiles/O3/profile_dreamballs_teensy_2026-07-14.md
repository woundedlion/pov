# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_dreamballs_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `DreamBalls<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 92,548 B; ITCM (RAM1 code) 49,752 B; RAM2 free 4,736 B.
(At -Os: FLASH 64,716 / ITCM 36,328 — +43% / +36%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 257-320 root counter
3,585,651,591 cyc = 5,976,086.0 µs vs measured `micros()` window sum 5,976,092 µs
(Δ ≈ 1.0 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

DreamBalls is bound by its mesh wireframe rasterizer — a phase-cycling wireframe whose cost tracks the number of orbiting shell copies (`db_mesh_plot` calls/frame), stepping through the preset cycle. Render spans
**7.9–72.5 ms** across the pass. `db_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 257-320)

```
frame                  93.38 ms  56.03 Mcyc  100%
  db_timeline_step       72.52 ms  43.51 Mcyc   77%
    db_draw                72.47 ms  43.48 Mcyc   99%  x2  48314 us/call
      db_draw_scene          72.47 ms  43.48 Mcyc   99%  x2  48313 us/call
        db_mesh_plot           72.04 ms  43.22 Mcyc   99%  x21  3430 us/call
        db_warp_orient           232 us  139.0 kcyc    0%  x21  11 us/call
        db_displace              195 us  117.1 kcyc    0%  x21  9 us/call
      db_mesh_copy               1 us    464 cyc    0%  x2  1 us/call
  db_buffer_wait         20.86 ms  12.51 Mcyc   22%
```

Wall: min 61.0 / avg 93.4 / max 126.2 ms.
Render (`frame` − `db_buffer_wait`) = **72.52 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 897-960)

```
frame                  62.09 ms  37.25 Mcyc  100%
  db_timeline_step       36.10 ms  21.66 Mcyc   58%
    db_draw                36.06 ms  21.64 Mcyc   99%  x1  36059 us/call
      db_draw_scene          36.06 ms  21.64 Mcyc   99%  x1  36059 us/call
        db_mesh_plot           35.85 ms  21.51 Mcyc   99%  x10  3585 us/call
        db_warp_orient           111 us  66.8 kcyc    0%  x10  11 us/call
        db_displace               94 us  56.3 kcyc    0%  x10  9 us/call
      db_mesh_copy               0 us    298 cyc    0%  x1  0 us/call
  db_buffer_wait         25.99 ms  15.59 Mcyc   41%
```

Wall: min 44.7 / avg 62.1 / max 80.2 ms.
Render = **36.10 ms**. The spread between this and the peak is
phase/preset-driven (copy count varies), not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       13769/s  0.55 / 2.39 / 118.7 us  cpu 4.40%
  isr_pack       1721/s  4.79 / 12.04 / 117.1 us  cpu 2.77%
  isr_dma_submit 1721/s  0.73 / 0.95 / 1.2 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`db_mesh_plot`** — 72.0 ms/frame, 77% of the peak frame; `db_buffer_wait` is the display-sync idle by design.

**Real plot-cull beneficiary.** The wireframe's long geodesic/planar edges cross the quadrant, so the new per-edge column cull removes fragments in the wrong arm half: peak render fell ~145 -> ~109 ms vs the pre-batch capture, and -O3 pulls the heavy phase from ~5-8 fps to a steady 8 fps.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/DreamBalls.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=DreamBalls`), + `tools/profile_capture.py`.
