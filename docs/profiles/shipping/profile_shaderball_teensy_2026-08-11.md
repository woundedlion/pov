# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-12, selective O3)

Final post-optimization capture of the current 23-preset bank. Raw capture:
`build/prof/shaderball_final_green_cycle.log`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, segmented POV driver, COM3 |
| Image | `profile`: `-Os` base plus shipping selective-O3 regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, source `dad52e09` |
| Method | 16-frame windows, 110 s fast cycle, two-frame through-clear transitions, epoch stretched to 1,400 revs |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 110 16 -D HS_PROFILE_EFFECT_HEAP_BYTES=8192 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

The 108-window capture has monotonic frame numbers, visits all 23 presets,
confirms the 22-to-0 wrap, and has no epoch reset. Root cycles divided by 600
MHz agree with wall time within **0.7 ppm**.

## Result

Every preset is green. Peak render is **49.52 ms** with **0/1728 spilled
frames**. The matched pre-optimization capture peaked at 88.07 ms with
849/1680 spills and only 6 of 23 green buckets.

| Preset | Peak render, ms | Spilled/owned frames | Clean shader ms/frame |
|---:|---:|---:|---:|
| 22 | 49.52 | 0/68 | 43.47 |
| 0 | 49.15 | 0/99 | 27.46 |
| 19 | 48.02 | 0/68 | 43.02 |
| 20 | 47.63 | 0/68 | 31.70 |
| 21 | 44.56 | 0/68 | 38.82 |
| 18 | 43.35 | 0/68 | 37.13 |
| 15 | 39.29 | 0/68 | 34.61 |
| 16 | 38.95 | 0/68 | 28.49 |
| 7 | 36.02 | 0/68 | 27.70 |
| 6 | 35.88 | 0/68 | 31.07 |
| 17 | 34.90 | 0/68 | 29.88 |
| 3 | 34.82 | 0/102 | 29.60 |
| 4 | 34.75 | 0/99 | 24.32 |
| 12 | 34.49 | 0/68 | 24.80 |
| 2 | 34.47 | 0/102 | 29.58 |
| 11 | 34.18 | 0/68 | 29.37 |
| 1 | 34.10 | 0/102 | 29.50 |
| 5 | 33.76 | 0/68 | 30.02 |
| 8 | 32.89 | 0/68 | 26.82 |
| 10 | 32.31 | 0/68 | 27.53 |
| 14 | 32.24 | 0/68 | 27.44 |
| 9 | 31.94 | 0/68 | 25.85 |
| 13 | 31.87 | 0/68 | 27.07 |

## What closed the holdouts

The retained implementation combines exact dodecahedral and square-Peirce
folds, selected compile-time inverse-pipeline specializations, prepared
surface-noise frame constants, a single-traversal vector simplex field, and a
persistent 64x16 value/hue field. The hue field removes per-pixel OKLab work
from the animated preset-21 path while bilinear interpolation keeps its added
error within a pinned maximum of 5,278/65,535 and mean of 202 channel units.

## Full-roster size and validation

The attested Phantasm image passes every memory gate:

```text
FLASH code     411,112 B
RAM1 variables 314,816 B
ITCM code      195,960 / 196,608 B (648 B padding)
RAM1 local free 12,864 B
RAM2 allocator free 4,224 B
```

Relative to the matched baseline, the final image saves 384 bytes of ITCM and
adds 9,304 bytes of flash code. ShaderBall's host stack watermark is 7,007
bytes. All 68 native test targets and the release WASM smoke roster pass.

## Caveats

- Counter scopes include live flywheel/DMA ISR time.
- Fast-cycle mode compresses dwell only; it does not change per-frame work.
- Preset buckets own their following transitions, so they are stricter than
  clean-hold shader figures.
- Global O3 is not a valid shipping substitute; its matched twin makes preset
  19 exceed 200 ms.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShaderBall`, and
`HS_PROFILE_WINDOW=16`; `tools/profile_one.sh` performs the locked build,
flash, capture, marker checks, and artifact attestation.
