# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_petalflow.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `PetalFlow<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile PetalFlow` |

Image size: FLASH code 38,964 B; ITCM (RAM1 code) 25,752 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 961-1024 root counter
2,399,682,394 cyc = 3,999,470.7 µs vs measured `micros()` window sum 3,999,479 µs
(Δ ≈ 2.1 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

PetalFlow is bound by its ring rasterizer (idle-bound at 16 fps) — petal rings rasterized via `Plot::rasterize` (`pf_ring_scan`, ~x23 rings/frame); steady render ~16 ms, always one window (16 fps). Render spans
**15.3–16.7 ms** across the pass. `pf_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 961-1024)

```
frame                  62.49 ms  37.50 Mcyc  100%
  pf_draw_rings          16.66 ms  10.00 Mcyc   26%
    pf_ring_scan           15.41 ms  9.25 Mcyc   92%  x23  662 us/call
    pf_ring_build           1.14 ms  683.4 kcyc    6%  x23  49 us/call
  pf_timeline_step          22 us  13.0 kcyc    0%
  pf_buffer_wait         45.81 ms  27.49 Mcyc   73%
```

Wall: min 61.5 / avg 62.5 / max 63.5 ms.
Render (`frame` − `pf_buffer_wait`) = **16.68 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 513-576)

```
frame                  62.49 ms  37.49 Mcyc  100%
  pf_draw_rings          16.20 ms  9.72 Mcyc   25%
    pf_ring_scan           14.95 ms  8.97 Mcyc   92%  x23  643 us/call
    pf_ring_build           1.14 ms  685.8 kcyc    7%  x23  49 us/call
  pf_timeline_step          26 us  15.5 kcyc    0%
  pf_buffer_wait         46.26 ms  27.76 Mcyc   74%
```

Wall: min 61.4 / avg 62.5 / max 63.4 ms.
Render = **16.23 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9216/s  0.63 / 1.57 / 65.4 us  cpu 2.90%
  isr_pack       1152/s  4.39 / 5.31 / 63.4 us  cpu 1.22%
  isr_dma_submit 1152/s  0.65 / 0.96 / 3.0 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`pf_ring_scan`** — 15.4 ms/frame, 25% of the peak frame; `pf_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps. Ring edges route through `Plot::rasterize`; the cull shaves `pf_ring_scan` slightly but the render is a quarter of one window, so no cadence change.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/PetalFlow.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=PetalFlow`), + `tools/profile_capture.py`.
