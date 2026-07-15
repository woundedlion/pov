# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_petalflow_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `PetalFlow<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 65,508 B; ITCM (RAM1 code) 51,224 B; RAM2 free 4,736 B.
(At -Os: FLASH 38,964 / ITCM 25,752 — +68% / +98%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 961-1024 root counter
2,399,781,158 cyc = 3,999,635.3 µs vs measured `micros()` window sum 3,999,641 µs
(Δ ≈ 1.4 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

PetalFlow is bound by its ring rasterizer (idle-bound at 16 fps) — petal rings rasterized via `Plot::rasterize` (`pf_ring_scan`, ~x23 rings/frame); steady render ~16 ms, always one window (16 fps). Render spans
**10.5–11.5 ms** across the pass. `pf_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 961-1024)

```
frame                  62.49 ms  37.50 Mcyc  100%
  pf_draw_rings          11.47 ms  6.88 Mcyc   18%
    pf_ring_scan           10.62 ms  6.37 Mcyc   92%  x23  456 us/call
    pf_ring_build            771 us  462.5 kcyc    6%  x23  33 us/call
  pf_timeline_step          11 us   6.6 kcyc    0%
  pf_buffer_wait         51.01 ms  30.61 Mcyc   81%
```

Wall: min 61.7 / avg 62.5 / max 63.2 ms.
Render (`frame` − `pf_buffer_wait`) = **11.49 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 513-576)

```
frame                  62.49 ms  37.49 Mcyc  100%
  pf_draw_rings          11.15 ms  6.69 Mcyc   17%
    pf_ring_scan           10.29 ms  6.17 Mcyc   92%  x23  442 us/call
    pf_ring_build            773 us  463.6 kcyc    6%  x23  33 us/call
  pf_timeline_step          11 us   6.5 kcyc    0%
  pf_buffer_wait         51.33 ms  30.80 Mcyc   82%
```

Wall: min 61.9 / avg 62.5 / max 63.0 ms.
Render = **11.16 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9216/s  0.43 / 1.45 / 61.8 us  cpu 2.67%
  isr_pack       1152/s  4.65 / 5.59 / 60.3 us  cpu 1.28%
  isr_dma_submit 1152/s  0.61 / 0.91 / 6.3 us  cpu 0.20%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`pf_ring_scan`** — 10.6 ms/frame, 17% of the peak frame; `pf_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps. Ring edges route through `Plot::rasterize`; the cull shaves `pf_ring_scan` slightly but the render is a quarter of one window, so no cadence change.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/PetalFlow.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=PetalFlow`), + `tools/profile_capture.py`.
