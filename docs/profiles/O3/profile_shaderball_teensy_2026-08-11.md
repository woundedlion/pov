# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-11, global O3)

Global-O3 regression twin of the
[shipping capture](../shipping/profile_shaderball_teensy_2026-08-11.md).
Raw capture: `build/prof/shaderball_final_cycle_o3.log`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, segmented POV driver, COM3 |
| Image | `profile_o3`: global `-O3 -ffast-math` ceiling |
| Method | 16-frame windows, 150 s fast cycle, all 29 presets |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 150 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

The 142-window capture has no epoch reset, visits all 29 presets, confirms the
28-to-0 wrap, and agrees with wall time to 0.5 ppm.

## Result

Peak render is **66.05 ms** with **96/2272 spilled frames (4.23%)**. Strict
buckets are 26 green and three red:

| Preset | Peak render, ms | Spilled/owned frames |
|---:|---:|---:|
| 28 | 66.05 | 63/68 |
| 20 | 64.46 | 31/68 |
| 0 | 65.80 | 2/99 |

Global O3 is no longer a useful shipping direction: it spills 30 more frames
than the selective-O3 image and makes preset 20 red. The shipping topology is
the accepted implementation.

This twin covers the same 29-preset optimization corpus as shipping. Two
presets landed concurrently afterward; post-rebase host tests cover all 31,
but they are not represented in this device cycle.
