# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-12, global O3)

Source-matched global-O3 twin of the
[shipping capture](../shipping/profile_shaderball_teensy_2026-08-11.md).
Raw capture: `build/prof/shaderball_final_green_o3_cycle.log`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, segmented POV driver, COM3 |
| Image | `profile_o3`: global `-O3 -ffast-math` compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, source `dad52e09` |
| Method | 16-frame windows, 110 s fast cycle, two-frame through-clear transitions, epoch stretched to 1,400 revs |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 110 16 -D HS_PROFILE_EFFECT_HEAP_BYTES=8192 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

The 97-window capture validates all 23 presets and the 22-to-0 wrap with no
epoch reset. Root cycles agree with wall time within **0.8 ppm**.

## Result

Global O3 is rejected. It produces 21 green and two red preset-owned buckets,
peaks at **207.37 ms**, and spills **66/1552 frames**. The shipping image is
23/23 green, peaks at 49.52 ms, and spills none.

| Preset | O3 peak render, ms | O3 spilled/owned | O3 clean shader, ms/f | Shipping peak, ms |
|---:|---:|---:|---:|---:|
| 19 | 205.81 | 64/68 | 198.33 | 48.02 |
| 20 | 207.37 | 2/68 | 30.64 | 47.63 |
| 22 | 50.40 | 0/59 | 44.95 | 49.52 |
| 0 | 49.87 | 0/65 | 29.17 | 49.15 |
| 21 | 44.96 | 0/68 | 39.23 | 44.56 |

Preset 20's two spills are transition ownership from preset 19. The intrinsic
failure is the Peirce/dodecahedral preset 19: global compiler placement makes
its clean shader 155.31 ms slower than shipping. This is a layout regression,
not an argument for broader O3 promotion.

## Single-effect size delta

| Image | FLASH code | ITCM code |
|---|---:|---:|
| Shipping profile | 108,312 B | 28,536 B |
| Global-O3 profile | 124,456 B | 40,440 B |
| O3 minus shipping | **+16,144 B** | **+11,904 B** |

The source-matched Phantasm gate remains the shipping selective-O3 image at
195,960/196,608 ITCM bytes; the O3 profile is a single-effect diagnostic
ceiling and is not a full-roster shipping configuration.

## Caveats

- Counter scopes include live flywheel/DMA ISR time.
- Fast-cycle mode changes dwell only, not per-frame work.
- Preset-owned transition buckets are stricter than clean-hold shader figures.
- Compiler topology is unusually sensitive in preset 19; compare exact source
  and build flags rather than extrapolating from other presets.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShaderBall`, and
`HS_PROFILE_WINDOW=16`; `tools/profile_one.sh` performs the locked build,
flash, capture, marker checks, and artifact attestation.
