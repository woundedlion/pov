# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-11, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_shaderball_teensy_2026-08-11.md).
Raw capture: `build/prof/shaderball_o3.log`. It covers the same 29 presets,
compressed cycle, hardware, and source commit.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `0f50fe97` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 150 s capture, 32-frame holds, two-frame blends, `HS_PROFILE_EPOCH_REVS=1400` |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 150 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

Image size: `FLASH: code:125736, data:171992, headers:8448` / `RAM1:
variables:323232, code:54920, padding:10616, free:135520` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 897-912 root counter cycles divided by
600 MHz match the measured wall sum within **0.5 ppm**. The capture has 138
windows, all 29 indices, no epoch reset, and a confirmed 28-to-0 wrap.

## Frame cadence

Pass aggregate: peak frame render is **97.47 ms** and **168/2208 frames
(7.61%)** spill. The strict buckets are 24 green and 5 red.

Preset 28 holds `sb_shader_draw` at **66.33 ms/frame**, peaks at **72.74 ms**,
and spills **65/68** owned frames. Global O3 improves it, but the hold still
exceeds the 62.5 ms display window by 3.83 ms. `canvas_buffer_wait` is the
round-up idle to the next flip.

## Phase-by-phase readout

### Dodecahedral grid hold (window frames 961-976)

```
frame                  124.75 ms  74.85 Mcyc  100%
  sb_shader_draw        66.33 ms  39.80 Mcyc   53%
  sb_timeline_step      52.4 us   31.4 kcyc     0%
  canvas_clear          84.7 us   50.8 kcyc     0%
  canvas_buffer_wait    53.98 ms  32.39 Mcyc   43%
```

Wall min/avg/max = 124.48/124.75/125.02 ms. Render min/avg/max is
69.52/70.77/72.55 ms; the hold remains on the 8 fps cadence tier.

### Transition into dodecahedral grid (window frames 945-960)

```
frame                  121.49 ms  72.89 Mcyc  100%
  sb_shader_draw        61.10 ms  36.66 Mcyc   50%  two endpoints
  sb_timeline_step      63.0 us   37.8 kcyc     0%
  canvas_clear          84.8 us   50.9 kcyc     0%
  canvas_buffer_wait    55.93 ms  33.56 Mcyc   46%
```

Wall min/avg/max = 63.30/121.49/129.81 ms. The short frame is the clear
midpoint; the other transition frames render both endpoint configurations.

### Per-preset table

| # | authored regime | shader ms | peak render ms | cadence | spilled/frames |
|---:|---|--:|--:|---:|---:|
| 26 | Peirce primitive lattice | 92.09 | 97.47 | red | 64/68 |
| 28 | dodecahedral grid + mirror/noise | 66.33 | 72.74 | red | 65/68 |
| 27 | kaleidoscope edge-fade liquid | 59.20 | 97.14 | red | 33/68 |
| 0 | liquid stereo 0 | 57.05 | 68.85 | red | 5/99 |
| 19 | diagnostic Peirce | 56.02 | 61.17 | green | 0/68 |
| 18 | diagnostic curl-flow | 54.96 | 59.54 | green | 0/68 |
| 16 | diagnostic Bonne | 54.03 | 58.95 | green | 0/68 |
| 14 | liquid stereo 14 | 49.91 | 62.64 | red | 1/68 |
| 11 | liquid stereo 11 | 49.90 | 55.08 | green | 0/68 |
| 25 | Bonne primitive lattice | 49.20 | 53.97 | green | 0/68 |
| 15 | liquid stereo 15 | 48.77 | 54.73 | green | 0/68 |
| 7 | liquid stereo 7 | 48.50 | 58.37 | green | 0/69 |
| 6 | liquid stereo 6 | 48.15 | 51.27 | green | 0/102 |
| 5 | liquid stereo 5 | 45.87 | 51.90 | green | 0/102 |
| 13 | liquid stereo 13 | 44.84 | 36.91 | green | 0/68 |
| 8 | liquid stereo 8 | 42.54 | 52.15 | green | 0/68 |
| 17 | diagnostic vector-noise lattice | 42.50 | 59.30 | green | 0/68 |
| 4 | liquid stereo 4 | 42.12 | 48.56 | green | 0/102 |
| 9 | liquid stereo 9 | 42.02 | 47.09 | green | 0/68 |
| 20 | diagnostic Airocean | 41.37 | 60.31 | green | 0/68 |
| 10 | liquid stereo 10 | 40.51 | 46.66 | green | 0/68 |
| 3 | liquid stereo 3 | 37.98 | 42.60 | green | 0/102 |
| 1 | liquid stereo 1 | 37.97 | 62.17 | green | 0/102 |
| 2 | liquid stereo 2 | 37.96 | 42.66 | green | 0/102 |
| 23 | gnomonic grid + kaleidoscope | 32.81 | 37.34 | green | 0/68 |
| 22 | kaleidoscope mirror | 29.82 | 34.36 | green | 0/68 |
| 21 | wave-shear liquid | 28.41 | 44.58 | green | 0/68 |
| 24 | gnomonic grid + glitch | 27.27 | 37.50 | green | 0/68 |
| 12 | liquid stereo 12 | 26.42 | 54.60 | green | 0/68 |

The capture visits every index and wraps 28 to 0. Each preset contributes two
to seven clean windows depending on cadence and marker alignment.

### Per-pixel figures

Preset 28's 66.33 ms shader scope is about **6.40 us or 3,839 cycles per
pixel** across the 10,368-pixel quadrant.

## Column-ISR / DMA marshaling cost

```
isr_wake        2301/frame  min/avg/max 0.61/1.84/14.00 us  cpu 3.39%
isr_pack         288/frame  min/avg/max 6.31/7.05/9.37 us   cpu 1.62%
isr_dma_submit   288/frame  min/avg/max 0.64/0.91/6.56 us   cpu 0.20%
```

The ISR share is 5.22%, leaving about 59.24 ms of the display window for
rendering. Preset 28 needs about a **1.12x** speedup to fit that adjusted
budget.

## Summary ranking

1. `sb_shader_draw`, Peirce primitive lattice - 92.09 ms held and 97.47 ms peak.
2. `sb_shader_draw`, dodecahedral grid - 66.33 ms held and 72.74 ms peak.
3. `sb_shader_draw`, kaleidoscope edge-fade liquid - 59.20 ms held.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` has no visible subtree in this shallow ShaderBall capture.
- Global O3 is a single-effect ceiling, not a shippable full-roster image.
- Fast-cycle dwell compression changes duration, not per-frame work.
- The capture ran from clean source commit `0f50fe97`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.

## Global O3 vs selective O3

Global O3 lowers preset 28's clean shader cost from 74.24 to 66.33 ms
(**1.12x**) and the pass peak from 370.66 to 97.47 ms. It adds **28,440 B
FLASH** and **25,472 B ITCM** to the single-effect profile image. Neither
configuration keeps preset 28 at 16 fps.
