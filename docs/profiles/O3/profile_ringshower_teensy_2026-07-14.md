# RingShower on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_ringshower_o3.log` (64-frame windows, ~75 s single pass). This is the **-O3** twin of the shipping `-Os` report (`../Os/profile_ringshower_teensy_2026-07-14.md`); the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `RingShower<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 72,996 B; ITCM (RAM1 code) 57,688 B; RAM2 free 4,736 B.
(At -Os: FLASH 41,036 / ITCM 27,000 — +77% / +113%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 577-640 root counter
2,398,810,787 cyc = 3,998,018.0 µs vs measured `micros()` window sum 3,998,024 µs
(Δ ≈ 1.5 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

RingShower is bound by its idle-bound at 16 fps (lightest effect) — expanding rings drawn via `Plot::Ring` (`rsh_ring_plot`); render ~3 ms, the lightest effect in the roster (16 fps). Render spans
**0.6–2.2 ms** across the pass. `rsh_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 577-640)

```
frame                  62.47 ms  37.48 Mcyc  100%
  rsh_draw_rings          2.18 ms  1.31 Mcyc    3%
    rsh_ring_plot           2.18 ms  1.31 Mcyc   99%  x3  835 us/call
  rsh_timeline_step         19 us  11.1 kcyc    0%
  rsh_buffer_wait        60.27 ms  36.16 Mcyc   96%
```

Wall: min 60.4 / avg 62.5 / max 64.5 ms.
Render (`frame` − `rsh_buffer_wait`) = **2.20 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 321-384)

```
frame                  62.52 ms  37.51 Mcyc  100%
  rsh_draw_rings          1.66 ms  993.7 kcyc    2%
    rsh_ring_plot           1.66 ms  993.3 kcyc   99%  x2  854 us/call
  rsh_timeline_step         22 us  13.2 kcyc    0%
  rsh_buffer_wait        60.84 ms  36.50 Mcyc   97%
```

Wall: min 60.9 / avg 62.5 / max 64.4 ms.
Render = **1.68 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9212/s  0.43 / 1.40 / 57.1 us  cpu 2.58%
  isr_pack       1151/s  4.65 / 5.28 / 55.7 us  cpu 1.21%
  isr_dma_submit 1151/s  0.61 / 0.92 / 1.3 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`rsh_ring_plot`** — 2.2 ms/frame, 3% of the peak frame; `rsh_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps (lightest effect); render ~3 ms. Same as Thrusters -- the plot column-cull win is sub-ms and invisible.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/RingShower.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=RingShower`), + `tools/profile_capture.py`.
