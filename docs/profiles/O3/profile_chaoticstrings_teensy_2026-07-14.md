# ChaoticStrings on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_chaoticstrings_o3.log` (64-frame windows, ~75 s single pass). This is the global-`-O3` ceiling twin of the shipping profile; the base is `-Os` and this effect carries no `HS_O3` region, so the only variable is the optimization level. Re-captured 2026-07-14 **after the plot column-cull batch** (9ac8cebd, 62450701, 708d4b9b).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm flags (**`-O3`** (`-ffast-math`), newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `ChaoticStrings<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` |

Image size: FLASH code 77,492 B; ITCM (RAM1 code) 61,368 B; RAM2 free 4,736 B.
(At -Os: FLASH 46,860 / ITCM 32,472 — +65% / +88%.)
(Sizes are the pre-batch build; the column-cull batch adds ~1 KB ITCM shared across `Plot::` effects, absorbed by granule padding.)

**Exactness cross-check** — peak window frames 833-896 root counter
2,399,779,311 cyc = 3,999,632.2 µs vs measured `micros()` window sum 3,999,638 µs
(Δ ≈ 1.5 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver releases
one `draw_frame` per window, each rendering this segment's quadrant (72-row band ×
144-column half ≈ **10,368 px**). Wall time snaps up to a whole number of 62.5 ms
windows (16 fps = 1 window, 8 fps = 2, 5.3 fps = 3).

ChaoticStrings is bound by its multiline rasterizer (idle-bound at 16 fps) — string vertices built then drawn as `Plot::Multiline`; render stays well under one 62.5 ms window, so the effect is display-sync-idle-bound at 16 fps. Render spans
**11.3–14.2 ms** across the pass. `cs_buffer_wait` (the `buffer_free()`
spin in the `Canvas` ctor) is timed separately and is exactly the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Each block is nested as in the counter tree (each parent includes its children).
Columns: time/frame, cycles/frame, % of frame; leaves add calls/frame and per-call
cost.

### Peak render window (frames 833-896)

```
frame                  62.49 ms  37.50 Mcyc  100%
  cs_multiline_draw      10.10 ms  6.06 Mcyc   16%
  cs_build_vertices       3.98 ms  2.39 Mcyc    6%
  cs_noise_prepare           0 us     61 cyc    0%
  cs_timeline_step          86 us  51.5 kcyc    0%
  cs_buffer_wait         48.32 ms  28.99 Mcyc   77%
```

Wall: min 59.8 / avg 62.5 / max 65.2 ms.
Render (`frame` − `cs_buffer_wait`) = **14.17 ms**. This is the heaviest observed
frame — the worst-case cadence.

### Median render window (frames 897-960)

```
frame                  62.48 ms  37.49 Mcyc  100%
  cs_multiline_draw       9.58 ms  5.75 Mcyc   15%
  cs_build_vertices       3.99 ms  2.39 Mcyc    6%
  cs_noise_prepare           0 us     58 cyc    0%
  cs_timeline_step          89 us  53.7 kcyc    0%
  cs_buffer_wait         48.82 ms  29.29 Mcyc   78%
```

Wall: min 61.0 / avg 62.5 / max 63.9 ms.
Render = **13.66 ms**. The spread between this and the peak is
pool/coverage jitter, not a per-frame cost change.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (main-loop counters are separate),
peak window. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake       9216/s  0.54 / 1.55 / 75.2 us  cpu 2.85%
  isr_pack       1151/s  4.76 / 5.73 / 73.6 us  cpu 1.32%
  isr_dma_submit 1151/s  0.61 / 0.92 / 1.3 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (COLUMN_US/8 wake grid, includes pack + submit);
`isr_pack` = the 72× `packPixel` marshal per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The wire transfer runs
asynchronously on eDMA and costs no CPU. The ISR machinery steals ~3–4 % of the
chip (~2–2.5 ms of the 62.5 ms window).

## Summary ranking + plot column-cull verdict

Dominant cost is **`cs_multiline_draw`** — 10.1 ms/frame, 16% of the peak frame; `cs_buffer_wait` is the display-sync idle by design.

Idle-bound at 16 fps throughout. The multiline column cull shaves `cs_multiline_draw` modestly (the render stays far under one window), so the saving is invisible on cadence -- it widens `cs_buffer_wait` idle.

The -O3-vs-Os per-effect speedup is the compiler and is orthogonal to the batch;
see [`../README.md`](../README.md) for the head-to-head render-ms table.

## Caveats

- **ISR time is included** (`CYCCNT` free-runs) — the real shipping condition.
- **Peak vs median**: the two representative windows bound the pass; per-call
  averages stay valid across both.
- Image size is the pre-batch build (+~1 KB ITCM from the batch).
- **-O3 build**: `profile_o3` env, which does not ship (the full roster overflows ITCM at -O3).
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/ChaoticStrings.h`.

## Harness

- `targets/Profile/Profile.ino` (`-D HS_PROFILE_TARGET=ChaoticStrings`), + `tools/profile_capture.py`.
