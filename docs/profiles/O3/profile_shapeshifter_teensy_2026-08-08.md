# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-08-08, **-O3**)

Global-`O3` twin of the [shipping capture](../shipping/profile_shapeshifter_teensy_2026-08-08.md).
Raw capture: `build/prof/shapeshifter_o3.log`. This replaces the 2026-08-04
capture after the dense planar-star shaded-fragment shader was deduplicated.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect reference ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, nine presets, tip `26a834cb` |
| Method | `HS_PROFILE`, 16-frame windows, 155 s capture, epoch stretched to 1,600 revolutions |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile_o3 155 16 -D HS_PROFILE_EPOCH_REVS=1600` |

Image size: `FLASH: code:102320, data:148252, headers:8500` / `RAM1:
variables:315072, code:72136, padding:26168, free:110912` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 2337–2352 root counter cycles ÷ 600 MHz
match the measured wall sum within **2.0 ppm**.

## Frame cadence

**Pass aggregate:** peak frame render is **56.72 ms** in frames 2289–2304,
with **0/2,448 spilled frames (0%)**. Clean-hold `ss_draw_all` spans
8.86–50.51 ms/frame across the nine presets.

Every preset holds 16 fps. The pass peak retains 5.78 ms against the 62.5 ms
display window. `ss_buffer_wait` is intentional display-sync idle and is
excluded from render time.

## Phase-by-phase readout

Each preset holds for 240 frames and changes by a snap. The capture visits all
nine presets and wraps from preset 9 back to preset 1.

### Adaptive 208-count planar star (frames 2289–2304, worst of capture)

```text
frame                  62.83 ms  37.70 Mcyc  100%
  ss_draw_all          50.38 ms  30.23 Mcyc   80%
    ss_plot_dispatch   50.18 ms  30.11 Mcyc   99%  x155.7  322 us/call
  ss_timeline_step     38.1 us   22.9 kcyc     0%
  ss_buffer_wait       12.40 ms   7.44 Mcyc   19%
    canvas_clear       88.4 us   53.0 kcyc     0%
```

Wall min/avg/max = 53.99/62.83/71.72 ms; render averages 50.52 ms and peaks
at 56.72 ms. Orientation-dependent culling changes the visible contour count.

### Per-preset table

The unlabeled startup bucket and later `Preset: 1/9` bucket are the same first
preset. Rows use the maximum clean-hold draw mean and owned-frame peak.

| # | Shape | `ss_draw_all` ms | render peak ms | fps |
|---:|---|--:|--:|--:|
| 1 | Planar star, adaptive count 208, sides 7.745 | 50.51 | 56.72 | 16 |
| 9 | Flower, count 72, amplitude 1.8721 | 36.55 | 37.47 | 16 |
| 4 | Flower, count 70 | 36.27 | 37.42 | 16 |
| 8 | Spherical polygon, count 144, amplitude 7.0696 | 25.23 | 25.55 | 16 |
| 7 | Spherical polygon, count 144, amplitude 2.377 | 21.49 | 23.34 | 16 |
| 6 | Spherical polygon, count 128, opposite | 18.28 | 19.20 | 16 |
| 2 | Spherical polygon, count 74 | 13.21 | 15.58 | 16 |
| 5 | Planar star, uniform count 72, sides 4.417 | 10.49 | 12.74 | 16 |
| 3 | Planar star, uniform count 43, sides 6.562 | 8.86 | 10.08 | 16 |

### Per-pixel figures

The standard timing image does not enable perturbing Plot count probes, so it
cannot supply blended-pixel or cycles-per-blend figures.

## Column-ISR / DMA marshaling cost

For the slow window:

```text
isr_wake        1159/frame  min/avg/max 0.59/1.87/11.25 us  cpu 3.44%
isr_pack         145/frame  min/avg/max 6.32/7.07/ 9.26 us  cpu 1.62%
isr_dma_submit   145/frame  min/avg/max 0.61/0.91/ 1.02 us  cpu 0.21%
```

DMA transfer continues asynchronously after submit. ISR time is already
absorbed by DWT cycle scopes; combined measured ISR CPU share is 5.27%.

## Summary ranking

1. `ss_plot_dispatch` — 50.18 ms/frame and 99.6% of draw time in the slow window.
2. `ss_buffer_wait` — 12.40 ms/frame of intentional display synchronization.
3. ISR wake/pack/submit — 5.27% CPU, already included in foreground scopes.

## Caveats

- DWT scopes include interrupt time because CYCCNT free-runs.
- The direct-AA path does not enter `filter_blend`.
- Global `O3` is a single-effect reference image, not the shipping roster build.
- Epoch stretching changes capture lifetime, not preset dwell or per-frame work.
- The capture ran from committed tip `26a834cb` with no source changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter` and
`HS_PROFILE_WINDOW=16`; reproduce with the command in Setup.

## Global `-O3` vs selective `-O3`

Global `O3` lowers the observed peak from **58.22 to 56.72 ms** (1.50 ms,
2.6%) and lowers the adaptive-star maximum clean-hold draw cost from **52.18
to 50.51 ms** (1.67 ms, 3.2%). It adds **28,616 B FLASH code** and **24,016 B
RAM1 code**. Both images remain green.
