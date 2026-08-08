# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-08-04, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`). Raw
capture: `build/prof/shapeshifter_ship.log`. This replaces the earlier
2026-08-04 capture after the adaptive planar-star pole and sampling update.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` plus shared selective-`O3` Plot raster/AA/sink regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, nine presets, tip `0695086b` |
| Method | `HS_PROFILE`, 16-frame windows, 150 s capture, epoch stretched to 1,440 revolutions |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 150 16 -D HS_PROFILE_EPOCH_REVS=1440` |

Image size: `FLASH: code:71180, data:540368, headers:8992` / `RAM1:
variables:315040, code:48440, padding:17096, free:143712` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 2273–2288 root counter cycles ÷ 600 MHz
match the measured wall sum within **3.4 ppm**.

## Frame cadence

**Pass aggregate:** peak frame render is **57.66 ms** in frames 2273–2288,
with **0/2,368 spilled frames (0%)**. Clean-hold `ss_draw_all` means span
8.89–51.69 ms/frame across the nine presets.

Every preset holds 16 fps against the 62.5 ms display window. ShapeShifter
renders one quadrant, approximately 10,368 pixels. The slowest frame retains
4.84 ms of render margin. `ss_buffer_wait` is intentional display-sync idle
and is excluded from render time.

## Phase-by-phase readout

Each preset holds for 240 frames and changes by a snap. The capture visits all
nine presets and wraps from preset 9 back to preset 1.

### Adaptive 208-count planar star (frames 2273–2288, worst of capture)

```text
frame                  62.25 ms  37.35 Mcyc  100%
  ss_draw_all          51.69 ms  31.01 Mcyc   83%
    ss_plot_dispatch   51.49 ms  30.89 Mcyc   99%  x147.2  350 us/call
  ss_timeline_step     46.0 us   27.6 kcyc     0%
  ss_buffer_wait       10.51 ms   6.31 Mcyc   16%
    canvas_clear       88.8 us   53.3 kcyc     0%
```

Wall min/avg/max = 51.98/62.25/72.42 ms; render averages 51.83 ms and peaks
at 57.66 ms. Orientation-dependent culling changes the visible contour count.

### Per-preset table

The unlabeled startup bucket and later `Preset: 1/9` bucket are the same first
preset. Rows use the maximum clean-hold draw mean and owned-frame peak.

| # | Shape | `ss_draw_all` ms | render peak ms | fps |
|---:|---|--:|--:|--:|
| 1 | Planar star, adaptive count 208, sides 7.745 | 51.69 | 57.66 | 16 |
| 9 | Flower, count 72, amplitude 1.8721 | 36.48 | 37.14 | 16 |
| 4 | Flower, count 70 | 35.19 | 36.62 | 16 |
| 8 | Spherical polygon, count 144, amplitude 7.0696 | 25.95 | 27.51 | 16 |
| 7 | Spherical polygon, count 144, amplitude 2.377 | 22.44 | 24.87 | 16 |
| 6 | Spherical polygon, count 128, opposite | 18.54 | 19.50 | 16 |
| 2 | Spherical polygon, count 74 | 13.44 | 15.86 | 16 |
| 5 | Planar star, uniform count 72, sides 4.417 | 10.82 | 11.97 | 16 |
| 3 | Planar star, uniform count 43, sides 6.562 | 8.89 | 10.17 | 16 |

### Per-pixel figures

The standard timing image does not enable perturbing Plot count probes, so it
cannot supply blended-pixel or cycles-per-blend figures.

## Column-ISR / DMA marshaling cost

For the slow window:

```text
isr_wake        1148/frame  min/avg/max 0.76/1.99/11.56 us  cpu 3.67%
isr_pack         144/frame  min/avg/max 6.29/7.07/ 9.36 us  cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.74/0.94/ 1.07 us  cpu 0.21%
```

DMA transfer continues asynchronously after submit. ISR time is already
absorbed by DWT cycle scopes; combined measured ISR CPU share is 5.50%.

## Summary ranking

1. `ss_plot_dispatch` — 51.49 ms/frame and 99.6% of draw time in the slow window.
2. `ss_buffer_wait` — 10.51 ms/frame of intentional display synchronization.
3. ISR wake/pack/submit — 5.50% CPU, already included in foreground scopes.

The adaptive envelope provides approximately 290-count density at the equator
and 145-count density at the poles with 208 actual contours. Small planar-star
edges use fewer azimuthal anchors while large edges retain the six-anchor cap.

## Caveats

- DWT scopes include interrupt time because CYCCNT free-runs.
- The direct-AA path does not enter `filter_blend`.
- Shared selective-`O3` regions are active; ShapeShifter has no local `HS_O3` region.
- Epoch stretching changes capture lifetime, not preset dwell or per-frame work.
- Two unrelated untracked documentation files were present and do not enter the image.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter` and
`HS_PROFILE_WINDOW=16`; reproduce with the command in Setup.
