# Thrusters on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_thrusters.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `Thrusters<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile Thrusters` |

Image size: FLASH code 45,340 B; ITCM (RAM1 code) 31,384 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 513-576 root counter
2,400,976,840 cyc = 4,001,628.1 µs vs measured `micros()` window sum 4,001,640 µs
(Δ ≈ 3.0 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

Thrusters is bound by its idle-bound at 16 fps (lightest tier) — a thruster ring + distorted ring drawn via `Plot::Ring`/`Plot::DistortedRing`; render ~3 ms, far under one 62.5 ms window (16 fps). Render spans
**2.5–3.4 ms** across the pass. `th_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 513-576)

```
frame                  62.53 ms  37.52 Mcyc  100%
  th_thrusters            1.32 ms  791.3 kcyc    2%
    th_thruster_draw        1.32 ms  790.7 kcyc   99%  x1  981 us/call
  th_warp_step               1 us    360 cyc    0%
  th_timeline_step        2.03 ms  1.22 Mcyc    3%
    th_ring_draw            1.99 ms  1.19 Mcyc   97%  x1  1989 us/call
  th_buffer_wait         59.17 ms  35.50 Mcyc   94%
```

Wall: min 59.0 / avg 62.5 / max 66.2 ms.
Render (`frame` − `th_buffer_wait`) = **3.35 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 769-832)

```
frame                  62.49 ms  37.50 Mcyc  100%
  th_thrusters             997 us  597.9 kcyc    1%
    th_thruster_draw         996 us  597.5 kcyc   99%  x1  996 us/call
  th_warp_step               1 us    573 cyc    0%
  th_timeline_step        2.00 ms  1.20 Mcyc    3%
    th_ring_draw            1.96 ms  1.17 Mcyc   97%  x1  1957 us/call
  th_buffer_wait         59.50 ms  35.70 Mcyc   95%
```

Wall: min 59.5 / avg 62.5 / max 64.4 ms.
Render = **3.00 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9220/s  0.62 / 1.47 / 28.6 us  cpu 2.71%
  isr_pack       1152/s  4.25 / 4.71 / 26.2 us  cpu 1.08%
  isr_dma_submit 1152/s  0.62 / 0.96 / 1.1 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`th_ring_draw`** — 2.0 ms/frame, 3% of the peak frame; `th_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps; render ~3 ms. Ring/DistortedRing edges route through `Plot::rasterize`, but the effect is ~5 % of one window, so any column-cull win is sub-ms and invisible.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Thrusters.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=Thrusters`), + `tools/profile_capture.py`.
