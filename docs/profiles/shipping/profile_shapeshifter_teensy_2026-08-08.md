# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-08-08, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`). Raw
capture: `build/prof/shapeshifter_ship.log`. This replaces the 2026-08-04
capture after the dense planar-star shaded-fragment shader was deduplicated.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` plus shared selective-`O3` Plot raster/AA/sink regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, nine presets, tip `26a834cb` |
| Method | `HS_PROFILE`, 16-frame windows, 155 s capture, epoch stretched to 1,600 revolutions |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 -D HS_PROFILE_EPOCH_REVS=1600` |

Image size: `FLASH: code:73704, data:148316, headers:8380` / `RAM1:
variables:315040, code:48120, padding:17416, free:143712` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 2273–2288 root counter cycles ÷ 600 MHz
match the measured wall sum within **1.9 ppm**.

## Frame cadence

**Pass aggregate:** peak frame render is **58.22 ms** in frames 2273–2288,
with **0/2,448 spilled frames (0%)**. Clean-hold `ss_draw_all` spans
8.98–52.18 ms/frame across the nine presets.

Every preset holds 16 fps against the 62.5 ms display window. ShapeShifter
renders one quadrant, approximately 10,368 pixels. The slowest frame retains
4.28 ms of render margin. `ss_buffer_wait` is intentional display-sync idle
and is excluded from render time.

## Phase-by-phase readout

Each preset holds for 240 frames and changes by a snap. The capture visits all
nine presets and wraps from preset 9 back to preset 1.

### Adaptive 208-count planar star (frames 2273–2288, worst of capture)

```text
frame                  62.24 ms  37.34 Mcyc  100%
  ss_draw_all          52.18 ms  31.31 Mcyc   83%
    ss_plot_dispatch   52.01 ms  31.21 Mcyc   99%  x147.2  354 us/call
  ss_timeline_step     44.9 us   26.9 kcyc     0%
  ss_buffer_wait       10.00 ms   6.00 Mcyc   16%
    canvas_clear       89.0 us   53.4 kcyc     0%
```

Wall min/avg/max = 52.04/62.24/72.36 ms; render averages 52.33 ms and peaks
at 58.22 ms. Orientation-dependent culling changes the visible contour count.

### Per-preset table

The unlabeled startup bucket and later `Preset: 1/9` bucket are the same first
preset. Rows use the maximum clean-hold draw mean and owned-frame peak.

| # | Shape | `ss_draw_all` ms | render peak ms | fps |
|---:|---|--:|--:|--:|
| 1 | Planar star, adaptive count 208, sides 7.745 | 52.18 | 58.22 | 16 |
| 9 | Flower, count 72, amplitude 1.8721 | 37.47 | 38.15 | 16 |
| 4 | Flower, count 70 | 36.17 | 37.62 | 16 |
| 8 | Spherical polygon, count 144, amplitude 7.0696 | 25.70 | 27.25 | 16 |
| 7 | Spherical polygon, count 144, amplitude 2.377 | 22.20 | 24.61 | 16 |
| 6 | Spherical polygon, count 128, opposite | 18.34 | 19.29 | 16 |
| 2 | Spherical polygon, count 74 | 13.30 | 15.72 | 16 |
| 5 | Planar star, uniform count 72, sides 4.417 | 10.94 | 12.08 | 16 |
| 3 | Planar star, uniform count 43, sides 6.562 | 8.98 | 10.25 | 16 |

### Per-pixel figures

The standard timing image does not enable perturbing Plot count probes, so it
cannot supply blended-pixel or cycles-per-blend figures.

## Column-ISR / DMA marshaling cost

For the slow window:

```text
isr_wake        1148/frame  min/avg/max 0.73/2.00/11.62 us  cpu 3.67%
isr_pack         143/frame  min/avg/max 6.31/7.07/ 9.46 us  cpu 1.62%
isr_dma_submit   143/frame  min/avg/max 0.83/0.94/ 1.09 us  cpu 0.21%
```

DMA transfer continues asynchronously after submit. ISR time is already
absorbed by DWT cycle scopes; combined measured ISR CPU share is 5.50%.

## Summary ranking

1. `ss_plot_dispatch` — 52.01 ms/frame and 99.7% of draw time in the slow window.
2. `ss_buffer_wait` — 10.00 ms/frame of intentional display synchronization.
3. ISR wake/pack/submit — 5.50% CPU, already included in foreground scopes.

The shaded-fragment deduplication leaves one ShapeShifter `Plot::rasterize`
specialization in ITCM. The full Phantasm roster uses 194,888 B of ITCM,
5,136 B less than the preceding image, and passes with 1,720 B of headroom.

## Caveats

- DWT scopes include interrupt time because CYCCNT free-runs.
- The direct-AA path does not enter `filter_blend`.
- Shared selective-`O3` regions are active; ShapeShifter has no local `HS_O3` region.
- Epoch stretching changes capture lifetime, not preset dwell or per-frame work.
- The capture ran from committed tip `26a834cb` with no source changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter` and
`HS_PROFILE_WINDOW=16`; reproduce with the command in Setup.
