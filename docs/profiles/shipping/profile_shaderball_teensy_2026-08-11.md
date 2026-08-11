# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-11, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw capture:
`build/prof/shaderball_ship.log`. This replaces the 2026-08-10 report and
includes all 29 authored presets, including the dodecahedral grid preset at
index 28.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; `Scan::Shader` and ShaderBall color/sample hot regions use selective O3 |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `0f50fe97` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 150 s capture, 32-frame holds, two-frame blends, `HS_PROFILE_EPOCH_REVS=1400` |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 150 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

Image size: `FLASH: code:97296, data:171612, headers:8596` / `RAM1:
variables:323232, code:29448, padding:3320, free:168288` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 897-912 root counter cycles divided by
600 MHz match the measured wall sum within **0.4 ppm**. The capture has 117
windows, monotonically increasing frame numbers, all 29 preset indices, and a
confirmed wrap from preset 28 to preset 0.

## Frame cadence

Pass aggregate: peak frame render is **370.66 ms** and **233/1872 frames
(12.4%)** spill the 62.5 ms display window. The strict preset buckets are 21
green and 8 red. ShaderBall renders one 72x144 quadrant, or 10,368 pixels.

The new preset 28 holds `sb_shader_draw` at **74.24 ms/frame**. Its owned
frames peak at **79.41 ms**, with **33/34** frames spilling. Even before
display-sync idle, the clean hold exceeds the display window by 11.74 ms.
`canvas_buffer_wait` is the round-up idle to the next display flip.

## Phase-by-phase readout

The compressed cycle alternates 32-frame authored holds with two-frame blends.

### Dodecahedral grid hold (window frames 961-976)

```
frame                  124.88 ms  74.93 Mcyc  100%
  sb_shader_draw        74.24 ms  44.55 Mcyc   59%
  sb_timeline_step      59.3 us   35.6 kcyc     0%
  canvas_clear          84.9 us   50.9 kcyc     0%
  canvas_buffer_wait    46.29 ms  27.77 Mcyc   37%
```

Wall min/avg/max = 124.50/124.88/125.20 ms. Render min/avg/max for the window
is 77.90/78.60/78.91 ms, so every steady frame requires two display windows.

### Transition into dodecahedral grid (window frames 945-960)

```
frame                  121.73 ms  73.04 Mcyc  100%
  sb_shader_draw        65.78 ms  39.47 Mcyc   54%  two endpoints
  sb_timeline_step      64.6 us   38.7 kcyc     0%
  canvas_clear          84.8 us   50.9 kcyc     0%
  canvas_buffer_wait    51.51 ms  30.91 Mcyc   42%
```

Wall min/avg/max = 60.97/121.73/137.45 ms. The two endpoint shader scopes are
combined above; the clear midpoint is the one short frame.

### Per-preset table

`shader ms` is the clean-hold `sb_shader_draw` average. `peak render` and
spilled frames include the transition owned by that preset. Rows are ranked by
clean shader cost.

| # | authored regime | shader ms | peak render ms | cadence | spilled/frames |
|---:|---|--:|--:|---:|---:|
| 26 | Peirce primitive lattice | 364.76 | 370.66 | red | 35/39 |
| 19 | diagnostic Peirce | 142.74 | 148.00 | red | 64/68 |
| 28 | dodecahedral grid + mirror/noise | 74.24 | 79.41 | red | 33/34 |
| 27 | kaleidoscope edge-fade liquid | 62.85 | 370.37 | red | 33/34 |
| 0 | liquid stereo 0 | 60.06 | 78.83 | red | 63/65 |
| 11 | liquid stereo 11 | 54.06 | 58.92 | green | 0/68 |
| 18 | diagnostic curl-flow | 53.57 | 58.33 | green | 0/68 |
| 16 | diagnostic Bonne | 53.27 | 58.37 | green | 0/68 |
| 14 | liquid stereo 14 | 53.03 | 68.61 | red | 1/68 |
| 7 | liquid stereo 7 | 52.17 | 57.84 | green | 0/68 |
| 15 | liquid stereo 15 | 51.57 | 57.59 | green | 0/68 |
| 6 | liquid stereo 6 | 51.02 | 54.89 | green | 0/68 |
| 5 | liquid stereo 5 | 49.68 | 53.98 | green | 0/68 |
| 13 | liquid stereo 13 | 47.87 | 38.27 | green | 0/68 |
| 25 | Bonne primitive lattice | 47.31 | 52.43 | green | 0/68 |
| 8 | liquid stereo 8 | 46.90 | 58.17 | green | 0/68 |
| 20 | diagnostic Airocean | 46.70 | 148.45 | red | 2/68 |
| 9 | liquid stereo 9 | 45.77 | 51.40 | green | 0/68 |
| 17 | diagnostic vector-noise lattice | 44.93 | 59.09 | green | 0/68 |
| 4 | liquid stereo 4 | 44.23 | 50.82 | green | 0/68 |
| 10 | liquid stereo 10 | 43.66 | 50.73 | green | 0/68 |
| 1 | liquid stereo 1 | 40.37 | 64.63 | red | 2/68 |
| 3 | liquid stereo 3 | 40.36 | 45.28 | green | 0/68 |
| 2 | liquid stereo 2 | 40.35 | 45.05 | green | 0/68 |
| 23 | gnomonic grid + kaleidoscope | 32.29 | 36.76 | green | 0/68 |
| 22 | kaleidoscope mirror | 29.05 | 33.55 | green | 0/68 |
| 12 | liquid stereo 12 | 27.51 | 58.82 | green | 0/68 |
| 21 | wave-shear liquid | 27.25 | 54.34 | green | 0/68 |
| 24 | gnomonic grid + glitch | 26.77 | 37.27 | green | 0/68 |

The capture visits every index and wraps 28 to 0. Each preset contributes one
to five clean windows depending on cadence and marker alignment.

### Per-pixel figures

Preset 28 shades 10,368 pixels per frame. Its 74.24 ms shader scope is about
**7.16 us or 4,297 cycles per pixel**. The shallow capture has no separate
`filter_blend` counter, so this is end-to-end shader cost rather than a blend
microbenchmark.

## Column-ISR / DMA marshaling cost

For the preset-28 clean window:

```
isr_wake        2303/frame  min/avg/max 0.73/1.96/18.54 us  cpu 3.61%
isr_pack         288/frame  min/avg/max 6.21/7.13/9.43 us   cpu 1.64%
isr_dma_submit   288/frame  min/avg/max 0.63/0.94/10.85 us  cpu 0.21%
```

Pack dominates submit CPU cost; LED wire transfer proceeds asynchronously
after submission. ISR work consumes 5.46% of the display window, leaving about
59.09 ms of CPU render budget. Preset 28's shader needs about a **1.26x**
speedup to fit that adjusted budget.

## Summary ranking

1. `sb_shader_draw`, Peirce primitive lattice - 364.76 ms held and 370.66 ms peak.
2. `sb_shader_draw`, diagnostic Peirce - 142.74 ms held and 148.00 ms peak.
3. `sb_shader_draw`, dodecahedral grid - 74.24 ms held and 79.41 ms peak.

Preset 28 is the third-costliest shipping hold in this capture and does not
hold 16 fps.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` has no visible subtree in this shallow ShaderBall capture.
- Selective O3 covers `Scan::Shader`, `colorize`, and `sample_pattern`; the dodecahedral reflection fold remains in flash-optimized code.
- Fast-cycle dwell compression changes duration, not per-frame work.
- The capture ran from clean source commit `0f50fe97`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`; `just profile ShaderBall` builds, flashes, and captures.
