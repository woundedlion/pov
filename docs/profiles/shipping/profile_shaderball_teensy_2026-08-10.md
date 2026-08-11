# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-10, **selective -O3**)

Point-in-time snapshot after removing the static cost-point admission gate and
expanding the authored bank to 26 presets. Raw capture:
`build/prof/shaderball_ship.log`. Every preset is measured as authored; the
profile path does not fit, replace, or reject a valid configuration.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; shipping `HS_O3` regions remain active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `ffbd4807` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 140 s capture, two-revolution profile-only preset holds, `HS_PROFILE_EPOCH_REVS=1400` |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 140 16 "-D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400"` |

Image size: `FLASH: code:97744, data:169616, headers:9120` / `RAM1:
variables:323232, code:53816, padding:11720, free:135520` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 865-880 root counter cycles divided by
600 MHz match the measured wall sum within **1.2 ppm**. The capture has 136
windows, monotonically increasing frame numbers, all 26 preset indices, and a
confirmed wrap from preset 25 to preset 0.

## Frame cadence

Peak frame render is **67.71 ms** and **31/2176 frames (1.42%)** spill the
62.5 ms display window. The strict preset buckets are **23 green / 3 red**:
presets 0, 13, and 25 spill; the other 23 do not. ShaderBall renders one
72x144 quadrant, or 10,368 pixels.

The worst clean hold is preset 25 at **57.37 ms/frame** in
`sb_shader_draw`. The larger 67.71 ms pass peak belongs to preset 13's owned
transition, which is why cadence is reported from per-frame ownership rather
than clean-hold averages.

## Per-preset results

`shader ms` is the clean-hold `sb_shader_draw` average. `peak render` owns all
frames displayed for that preset, including its outgoing transition. The
profile-only dwell compression changes duration, not per-frame work.

| # | authored regime | shader ms | peak render ms | cadence | spilled/frames |
|---:|---|--:|--:|---:|---:|
| 0 | liquid stereo 0 | 39.73 | 64.15 | red | 2/99 |
| 1 | liquid stereo 1 | 39.70 | 44.95 | green | 0/102 |
| 2 | liquid stereo 2 | 39.67 | 44.97 | green | 0/102 |
| 3 | liquid stereo 3 | 36.33 | 50.51 | green | 0/102 |
| 4 | liquid stereo 4 | 49.38 | 57.64 | green | 0/102 |
| 5 | liquid stereo 5 | 50.20 | 55.26 | green | 0/102 |
| 6 | liquid stereo 6 | 51.50 | 59.47 | green | 0/102 |
| 7 | liquid stereo 7 | 47.83 | 55.58 | green | 0/102 |
| 8 | liquid stereo 8 | 46.31 | 52.26 | green | 0/102 |
| 9 | liquid stereo 9 | 46.62 | 51.44 | green | 0/102 |
| 10 | liquid stereo 10 | 54.05 | 59.75 | green | 0/102 |
| 11 | liquid stereo 11 | 27.95 | 60.37 | green | 0/102 |
| 12 | liquid stereo 12 | 38.64 | 39.50 | green | 0/71 |
| 13 | liquid stereo 13 | 53.74 | 67.71 | red | 1/68 |
| 14 | liquid stereo 14 | 54.08 | 59.82 | green | 0/68 |
| 15 | diagnostic 0 | 50.28 | 59.21 | green | 0/68 |
| 16 | diagnostic 1 | 41.83 | 55.95 | green | 0/68 |
| 17 | diagnostic 2 | 49.85 | 55.43 | green | 0/68 |
| 18 | diagnostic 3 | 50.77 | 57.73 | green | 0/68 |
| 19 | diagnostic 4 | 38.70 | 57.94 | green | 0/68 |
| 20 | wave shear liquid | 25.99 | 46.78 | green | 0/68 |
| 21 | kaleidoscope mirror | 27.23 | 32.61 | green | 0/68 |
| 22 | gnomonic grid + kaleidoscope | 27.69 | 33.04 | green | 0/68 |
| 23 | gnomonic grid + glitch | 21.07 | 32.95 | green | 0/68 |
| 24 | Bonne primitive lattice | 47.27 | 52.27 | green | 0/68 |
| 25 | Peirce primitive lattice | 57.37 | 63.26 | red | 28/68 |

## Interpretation

`sb_shader_draw` is the dominant scope. Preset 25 is the most expensive held
configuration, while preset 13 owns the worst transition frame. ISR/DMA work
continues during these scopes because CYCCNT free-runs. The normal build is not
fully cadence-safe for the expanded bank; unlike the removed point heuristic,
this conclusion comes from every authored preset running on the device.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- The profile image compresses preset dwell only; it does not change per-frame
  effect work.
- The capture ran from clean source commit `ffbd4807`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.


