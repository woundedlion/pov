# RingShower on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_ringshower.log` (64-frame windows, ~75 s single pass). Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `RingShower<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `just profile RingShower` |

Image size: FLASH code 41,036 B; ITCM (RAM1 code) 27,000 B; RAM2 free 4,736 B.
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 577-640 root counter
2,398,416,763 cyc = 3,997,361.3 µs vs measured `micros()` window sum 3,997,367 µs
(Δ ≈ 1.4 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

RingShower is bound by its idle-bound at 16 fps (lightest effect) — expanding rings drawn via `Plot::Ring` (`rsh_ring_plot`); render ~3 ms, the lightest effect in the roster (16 fps). Render spans
**0.9–3.1 ms** across the pass. `rsh_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 577-640)

```
frame                  62.46 ms  37.48 Mcyc  100%
  rsh_draw_rings          3.09 ms  1.85 Mcyc    4%
    rsh_ring_plot           3.09 ms  1.85 Mcyc   99%  x3  1184 us/call
  rsh_timeline_step         32 us  19.3 kcyc    0%
  rsh_buffer_wait        59.34 ms  35.60 Mcyc   95%
```

Wall: min 59.2 / avg 62.5 / max 65.6 ms.
Render (`frame` − `rsh_buffer_wait`) = **3.12 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 193-256)

```
frame                  62.50 ms  37.50 Mcyc  100%
  rsh_draw_rings          2.34 ms  1.41 Mcyc    3%
    rsh_ring_plot           2.34 ms  1.41 Mcyc   99%  x2  1350 us/call
  rsh_timeline_step         47 us  28.4 kcyc    0%
  rsh_buffer_wait        60.11 ms  36.07 Mcyc   96%
```

Wall: min 60.5 / avg 62.5 / max 65.5 ms.
Render = **2.39 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9211/s  0.62 / 1.46 / 67.0 us  cpu 2.69%
  isr_pack       1151/s  4.25 / 4.75 / 65.2 us  cpu 1.09%
  isr_dma_submit 1151/s  0.61 / 0.96 / 1.1 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`rsh_ring_plot`** — 3.1 ms/frame, 5% of the peak frame; `rsh_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps (lightest effect); render ~3 ms. Same as Thrusters -- the plot column-cull win is sub-ms and invisible.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- `-Os` build: the shipping Phantasm config; the `-O3` twin is under `../O3/`.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/RingShower.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=RingShower`), + `tools/profile_capture.py`.
