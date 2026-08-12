# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-11, selective O3)

Post-optimization shipping capture. Raw cycle capture:
`build/prof/shaderball_final_cycle_ship.log`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, segmented POV driver, COM3 |
| Image | `profile`: shipping `-Os` plus selective O3 |
| Driver | `POVSegmented<288, 4, 480>` |
| Method | 16-frame windows, 150 s fast cycle, 32-frame holds, two-frame through-clear transitions |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 150 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

The capture has 144 windows, monotonically increasing frame numbers, all 29
preset indices, and a confirmed 28-to-0 wrap. Root cycles agree with wall time
to 0.7 ppm.

## Cycle result

Peak render is **68.96 ms** with **66/2304 spilled frames (2.86%)**. Strict
transition-owned buckets are 27 green and two red:

| Preset | Peak render, ms | Spilled/owned frames |
|---:|---:|---:|
| 28 | 68.42 | 64/68 |
| 0 | 68.96 | 2/99 |

All other presets have zero spills. Preset 0 is green in its fixed capture;
its two cycle spills are owned transition frames from preset 28.

## Fixed-preset acceptance runs

Each row is a validated 70 s shipping capture with 32-frame windows.

| Preset | Baseline peak, ms | Final peak, ms | Gain, ms | 59 ms gate |
|---:|---:|---:|---:|---|
| 0 | 65.45 | 53.680 | 11.770 | pass |
| 1 | 45.02 | 44.985 | 0.035 | pass |
| 14 | 59.56 | 51.224 | 8.336 | pass |
| 19 | 166.78 | 49.864 | 116.916 | pass |
| 20 | 47.55 | 47.877 | -0.327 | pass |
| 26 | 391.55 | 54.563 | 336.987 | pass |
| 27 | 76.22 | 59.271 | 16.949 | **miss by 0.271** |
| 28 | 73.92 | 68.284 | 5.636 | **miss by 9.284** |

## Full-roster size

After rebasing onto the two concurrently added presets, `pio run -e phantasm`
passes all region gates: FLASH code 400,488 bytes,
RAM1 variables 314,816 bytes, ITCM code **196,152/196,608 bytes**, and 456
bytes of ITCM headroom. The pre-optimization roster used 196,232 ITCM bytes,
so the optimization set saves 80 bytes while adding the renderer fast paths.

Dense kaleidoscope, Peirce coordinate/fade/seam, gamut-error, legacy-noise
seam, and full native tests pass.

The device cycle predates those two concurrent preset additions and therefore
records the 29-preset optimization corpus. Post-rebase host coverage visits all
31 presets; the added presets were not included in this device timing capture.
