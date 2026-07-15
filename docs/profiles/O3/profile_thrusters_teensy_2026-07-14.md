# Thrusters on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_thrusters_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `Thrusters<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 69,380 B; ITCM (RAM1 code) 54,328 B; RAM2 free 4,736 B.
(At -Os: FLASH 45,340 / ITCM 31,384 — +53% / +73%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 513-576 root counter
2,400,605,399 cyc = 4,001,009.0 µs vs measured `micros()` window sum 4,001,012 µs
(Δ ≈ 0.8 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

Thrusters is bound by its idle-bound at 16 fps (lightest tier) — a thruster ring + distorted ring drawn via `Plot::Ring`/`Plot::DistortedRing`; render ~3 ms, far under one 62.5 ms window (16 fps). Render spans
**1.7–2.4 ms** across the pass. `th_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 513-576)

```
frame                  62.52 ms  37.51 Mcyc  100%
  th_thrusters             929 us  557.2 kcyc    1%
    th_thruster_draw         928 us  556.6 kcyc   99%  x1  690 us/call
  th_warp_step               1 us    696 cyc    0%
  th_timeline_step        1.42 ms  854.5 kcyc    2%
    th_ring_draw            1.39 ms  834.4 kcyc   97%  x1  1391 us/call
  th_buffer_wait         60.16 ms  36.10 Mcyc   96%
```

Wall: min 60.1 / avg 62.5 / max 64.9 ms.
Render (`frame` − `th_buffer_wait`) = **2.35 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 769-832)

```
frame                  62.49 ms  37.50 Mcyc  100%
  th_thrusters             705 us  423.2 kcyc    1%
    th_thruster_draw         705 us  422.7 kcyc   99%  x1  705 us/call
  th_warp_step               1 us    816 cyc    0%
  th_timeline_step        1.41 ms  843.1 kcyc    2%
    th_ring_draw            1.37 ms  823.2 kcyc   97%  x1  1372 us/call
  th_buffer_wait         60.38 ms  36.23 Mcyc   96%
```

Wall: min 60.5 / avg 62.5 / max 63.8 ms.
Render = **2.11 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9219/s  0.49 / 1.43 / 24.2 us  cpu 2.62%
  isr_pack       1152/s  4.64 / 5.22 / 21.8 us  cpu 1.20%
  isr_dma_submit 1152/s  0.58 / 0.90 / 1.2 us  cpu 0.20%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`th_ring_draw`** — 1.4 ms/frame, 2% of the peak frame; `th_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps; render ~3 ms. Ring/DistortedRing edges route through `Plot::rasterize`, but the effect is ~5 % of one window, so any column-cull win is sub-ms and invisible.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Thrusters.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=Thrusters`), + `tools/profile_capture.py`.
