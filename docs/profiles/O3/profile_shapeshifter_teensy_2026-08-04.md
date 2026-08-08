# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-08-04, **-O3**)

Global-`O3` twin of the [shipping capture](../shipping/profile_shapeshifter_teensy_2026-08-04.md).
Raw capture: `build/prof/shapeshifter_o3.log`. This replaces the earlier
2026-08-04 capture after the adaptive planar-star pole and sampling update.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect reference ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, nine presets, tip `0695086b` |
| Method | `HS_PROFILE`, 16-frame windows, 150 s capture, epoch stretched to 1,440 revolutions |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile_o3 150 16 -D HS_PROFILE_EPOCH_REVS=1440` |

Image size: `FLASH: code:101268, data:540376, headers:8592` / `RAM1:
variables:315040, code:73896, padding:24408, free:110944` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 2337–2352 root counter cycles ÷ 600 MHz
match the measured wall sum within **1.6 ppm**.

## Frame cadence

**Pass aggregate:** peak frame render is **59.47 ms** in frames 2145–2160,
with **0/2,368 spilled frames (0%)**. Clean-hold `ss_draw_all` means span
9.18–52.25 ms/frame across the nine presets.

Every preset holds 16 fps. The pass peak retains 3.03 ms against the 62.5 ms
display window. `ss_buffer_wait` is intentional display-sync idle and is
excluded from render time.

## Phase-by-phase readout

Each preset holds for 240 frames and changes by a snap. The capture visits all
nine presets and wraps from preset 9 back to preset 1.

### Adaptive 208-count planar star clean hold (frames 2289–2304)

```text
frame                  62.84 ms  37.70 Mcyc  100%
  ss_draw_all          52.16 ms  31.30 Mcyc   83%
    ss_plot_dispatch   51.96 ms  31.18 Mcyc   99%  x155.7  334 us/call
  ss_timeline_step     31.1 us   18.7 kcyc     0%
  ss_buffer_wait       10.64 ms   6.38 Mcyc   16%
    canvas_clear       89.1 us   53.5 kcyc     0%
```

Wall min/avg/max = 52.45/62.84/73.32 ms; render averages 52.29 ms and peaks
at 59.23 ms. The 59.47 ms pass peak occurs in a preset-snap window.

### Per-preset table

The unlabeled startup bucket and later `Preset: 1/9` bucket are the same first
preset. Rows use the maximum clean-hold draw mean and owned-frame peak.

| # | Shape | `ss_draw_all` ms | render peak ms | fps |
|---:|---|--:|--:|--:|
| 1 | Planar star, adaptive count 208, sides 7.745 | 52.25 | 59.47 | 16 |
| 9 | Flower, count 72, amplitude 1.8721 | 35.69 | 36.55 | 16 |
| 4 | Flower, count 70 | 35.39 | 36.50 | 16 |
| 8 | Spherical polygon, count 144, amplitude 7.0696 | 25.24 | 25.57 | 16 |
| 7 | Spherical polygon, count 144, amplitude 2.377 | 21.56 | 23.41 | 16 |
| 6 | Spherical polygon, count 128, opposite | 18.34 | 19.25 | 16 |
| 2 | Spherical polygon, count 74 | 13.24 | 15.60 | 16 |
| 5 | Planar star, uniform count 72, sides 4.417 | 10.86 | 13.24 | 16 |
| 3 | Planar star, uniform count 43, sides 6.562 | 9.18 | 10.51 | 16 |

### Per-pixel figures

The timing image does not enable perturbing Plot count probes, so blended
pixels and cycles per blend are unavailable.

## Column-ISR / DMA marshaling cost

For the representative adaptive-star window:

```text
isr_wake        1159/frame  min/avg/max 0.56/1.86/11.17 us  cpu 3.42%
isr_pack         145/frame  min/avg/max 6.32/7.07/ 9.25 us  cpu 1.62%
isr_dma_submit   145/frame  min/avg/max 0.63/0.91/ 1.02 us  cpu 0.21%
```

DMA transfer is asynchronous. Combined measured ISR CPU share is 5.25% and
is already included in foreground cycle scopes.

## Summary ranking

1. `ss_plot_dispatch` — 51.96 ms/frame and 99.6% of draw time in the clean hold.
2. `ss_buffer_wait` — 10.64 ms/frame of intentional display synchronization.
3. ISR wake/pack/submit — 5.25% CPU, already included in foreground scopes.

## Caveats

- DWT scopes include interrupt time because CYCCNT free-runs.
- The direct-AA path does not enter `filter_blend`.
- Global `O3` is a single-effect reference image, not the shipping roster build.
- Epoch stretching changes capture lifetime, not preset dwell or per-frame work.
- Two unrelated untracked documentation files were present and do not enter the image.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter` and
`HS_PROFILE_WINDOW=16`; reproduce with the command in Setup.

## Global `-O3` vs selective `-O3`

Global `O3` raises the observed peak from **57.66 to 59.47 ms** (1.81 ms,
3.1%) and raises the adaptive-star maximum clean-hold draw cost from **51.69
to 52.25 ms** (0.56 ms, 1.1%). It adds **30,088 B FLASH code** and
**25,456 B RAM1 code**. Both images remain green, so shipping selective `O3`
is the better size/performance tradeoff.
