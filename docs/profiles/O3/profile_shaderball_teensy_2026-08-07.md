# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-07, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_shaderball_teensy_2026-08-07.md).
Raw capture: `build/prof/shaderball_o3.log`. The first O3 attempt was discarded
because it preceded tip `39d72f18`; this replacement matches shipping exactly.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, tip `39d72f18` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:90872, data:148492, headers:8444` / `RAM1:
variables:314976, code:56824, padding:8712, free:143776` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 4593–4608 root counter cycles divided by
600 MHz match measured wall time within **1.5 ppm**
(`tools/parse_profile.py build/prof/shaderball_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `sb_shader_draw` averages
46.14 ms/frame; total render averages 51.42 ms/frame, peaks at **89.21 ms**
in frames 4593–4608, and spills **1437/4928 frames (29.2%)**.

The 72×144 quadrant contains 10,368 pixels. Three preset buckets spill: presets
5, 6, and the wrap into 1. Every other preset holds 16 fps. The
`canvas_buffer_wait` scope is cadence-quantized display-sync idle.

## Phase-by-phase readout

Presets 1–3 dwell then blend over 60 frames; preset 4 dwells before the first
480-frame parallel bridge; presets 5–12 chain 480-frame blends. Thirteen
markers confirm a complete wrap.

### Dwelled liquid regime (frames 273–288)

```
frame                  62.46 ms  37.47 Mcyc  100%
  sb_shader_draw       28.74 ms  17.24 Mcyc   46%
  sb_timeline_step     55.6 us   33.4 kcyc     0%
  canvas_clear         86.3 us   51.8 kcyc     0%
  canvas_buffer_wait   28.43 ms  17.06 Mcyc   45%
```

The inexpensive liquid anchors keep substantial one-window margin.

### Continuous expensive blend (frames 4593–4608, worst of capture)

```
frame                 124.76 ms  74.85 Mcyc  100%
  sb_shader_draw       81.63 ms  48.98 Mcyc   65%
  sb_timeline_step     58.6 us   35.2 kcyc     0%
  canvas_clear         84.6 us   50.8 kcyc     0%
  canvas_buffer_wait   37.96 ms  22.78 Mcyc   30%
```

Wall min/avg/max = 123.85/124.76/125.37 ms. All 16 frames spill. Global O3
cuts 8.92 ms from the shader relative to the matched shipping worst window,
but the transition still requires two display windows.

### Per-preset table

Presets 5–12 have no clean hold, so rows use the highest window nearest each
anchor; cadence is exact per-frame ownership from `buckets`. Initial unlabeled
frames are combined with preset 1, and the wrap is present.

| # | family | shader ms | render ms | cadence |
|---:|---|--:|--:|---:|
| 1 | liquid mild / wrap blend | 81.63 | 86.80 | 8↔16 |
| 6 | grid orbit A | 72.85 | 78.06 | 8 |
| 5 | liquid-to-grid bridge | 65.14 | 70.51 | 8 |
| 12 | grid/liquid B | 58.69 | 64.04 | 16 |
| 4 | liquid direct | 51.11 | 56.38 | 16 |
| 11 | grid/liquid A | 50.77 | 56.04 | 16 |
| 10 | grid orbit E | 49.36 | 54.66 | 16 |
| 7 | grid orbit B | 48.39 | 53.58 | 16 |
| 8 | grid orbit C | 43.70 | 48.98 | 16 |
| 9 | grid orbit D | 37.19 | 42.53 | 16 |
| 2 | liquid deep | 29.04 | 34.38 | 16 |
| 3 | liquid fine | 28.74 | 34.03 | 16 |

### Per-pixel figures

The worst window shades 10,368 pixels/frame. `sb_shader_draw` costs about
7.87 us/pixel, or 4,724 cycles/pixel. No `filter_blend` count is exposed.

## Column-ISR / DMA marshaling cost

```
isr_wake        2300/frame  min/avg/max 0.66/1.88/11.19 us  cpu 3.46%
isr_pack         288/frame  min/avg/max 6.32/7.03/9.20 us   cpu 1.61%
isr_dma_submit   288/frame  min/avg/max 0.58/0.91/1.04 us   cpu 0.20%
```

- DMA submission averages 13% of packing CPU cost.
- LED transfer continues asynchronously.
- Combined ISR share is 5.27%. The worst preset still needs a **1.43×**
  render speedup to fit one display window.

## Summary ranking

1. `sb_shader_draw` — 81.63 ms/frame and 65% of the worst representative frame.
2. `canvas_buffer_wait` — 37.96 ms/frame of display-sync idle.
3. Timeline and clear together remain below 0.2 ms/frame.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` is not instrumented on this closure-shader path.
- Global O3 is a single-effect ceiling, not a shippable roster image.
- `HS_PROFILE_EPOCH_REVS=3300` only extends effect lifetime.
- The capture ran from clean tip `39d72f18` with the mirrored high-contrast liquid palette.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.

## Global O3 vs selective O3

Global O3 lowers the worst-frame peak from 97.93 to 89.21 ms (**1.10×**) and
the aggregate render average from 58.86 to 51.42 ms. It adds **29,376 B FLASH**
and **26,592 B ITCM** to the single-effect image.
