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
| Effect | ShaderBall 288x144, single-entry playlist, tip `4f00b8d6` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 150 s capture, 32-frame holds, two-frame blends, `HS_PROFILE_EPOCH_REVS=1400` |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 150 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

Image size: `FLASH: code:125736, data:171992, headers:8448` / `RAM1:
variables:323232, code:54920, padding:10616, free:135520` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 897-912 root counter cycles divided by
600 MHz match the measured wall sum within **0.2 ppm**. The capture has 138
windows, all 29 indices, no epoch reset, and a confirmed 28-to-0 wrap.

## Frame cadence

Pass aggregate: peak frame render is **97.47 ms** and **168/2208 frames
(7.61%)** spill. The strict buckets are 24 green and 5 red.

Preset 28 holds `sb_shader_draw` at **66.32 ms/frame**, peaks at **72.74 ms**,
and spills **65/68** owned frames. Global O3 improves it, but the hold still
exceeds the 62.5 ms display window by 3.82 ms. `canvas_buffer_wait` is the
round-up idle to the next flip.

## Phase-by-phase readout

### Dodecahedral grid hold (window frames 961-976)

```
frame                  124.74 ms  74.84 Mcyc  100%
  sb_shader_draw        66.32 ms  39.79 Mcyc   53%
  sb_timeline_step      50.4 us   30.3 kcyc     0%
  canvas_clear          84.9 us   51.0 kcyc     0%
  canvas_buffer_wait    53.95 ms  32.37 Mcyc   43%
```

Wall min/avg/max = 124.51/124.74/125.10 ms. Render min/avg/max is
69.50/70.79/72.50 ms; the hold remains on the 8 fps cadence tier.

### Transition into dodecahedral grid (window frames 945-960)

```
frame                  121.50 ms  72.90 Mcyc  100%
  sb_shader_draw        61.10 ms  36.66 Mcyc   50%  two endpoints
  sb_timeline_step      62.5 us   37.5 kcyc     0%
  canvas_clear          84.7 us   50.8 kcyc     0%
  canvas_buffer_wait    55.94 ms  33.57 Mcyc   46%
```

Wall min/avg/max = 63.56/121.50/129.69 ms. The short frame is the clear
midpoint; the other transition frames render both endpoint configurations.

### Per-preset table

| # | authored regime | shader ms | peak render ms | cadence | spilled/frames |
|---:|---|--:|--:|---:|---:|
| 26 | Peirce primitive lattice | 91.84 | 97.47 | red | 64/68 |
| 28 | dodecahedral grid + mirror/noise | 66.32 | 72.74 | red | 65/68 |
| 27 | kaleidoscope edge-fade liquid | 59.20 | 95.93 | red | 33/68 |
| 0 | liquid stereo 0 | 57.04 | 68.95 | red | 5/99 |
| 19 | diagnostic Peirce | 55.81 | 60.99 | green | 0/68 |
| 18 | diagnostic curl-flow | 54.98 | 59.57 | green | 0/68 |
| 16 | diagnostic Bonne | 54.05 | 58.89 | green | 0/68 |
| 14 | liquid stereo 14 | 49.93 | 62.81 | red | 1/68 |
| 11 | liquid stereo 11 | 49.90 | 55.05 | green | 0/68 |
| 25 | Bonne primitive lattice | 49.20 | 53.98 | green | 0/68 |
| 15 | liquid stereo 15 | 48.77 | 54.76 | green | 0/68 |
| 7 | liquid stereo 7 | 48.51 | 58.33 | green | 0/69 |
| 6 | liquid stereo 6 | 48.15 | 51.28 | green | 0/102 |
| 5 | liquid stereo 5 | 45.87 | 51.91 | green | 0/102 |
| 13 | liquid stereo 13 | 44.84 | 36.94 | green | 0/68 |
| 8 | liquid stereo 8 | 42.53 | 52.15 | green | 0/68 |
| 17 | diagnostic vector-noise lattice | 42.50 | 59.29 | green | 0/68 |
| 4 | liquid stereo 4 | 42.13 | 48.53 | green | 0/102 |
| 9 | liquid stereo 9 | 42.02 | 47.06 | green | 0/68 |
| 20 | diagnostic Airocean | 41.45 | 60.53 | green | 0/68 |
| 10 | liquid stereo 10 | 40.51 | 46.58 | green | 0/68 |
| 3 | liquid stereo 3 | 37.97 | 42.60 | green | 0/102 |
| 1 | liquid stereo 1 | 37.97 | 62.25 | green | 0/102 |
| 2 | liquid stereo 2 | 37.96 | 42.59 | green | 0/102 |
| 23 | gnomonic grid + kaleidoscope | 32.79 | 37.34 | green | 0/68 |
| 22 | kaleidoscope mirror | 29.83 | 34.38 | green | 0/68 |
| 21 | wave-shear liquid | 28.41 | 44.53 | green | 0/68 |
| 24 | gnomonic grid + glitch | 27.25 | 37.55 | green | 0/68 |
| 12 | liquid stereo 12 | 26.42 | 54.65 | green | 0/68 |

The capture visits every index and wraps 28 to 0. Each preset contributes two
to seven clean windows depending on cadence and marker alignment.

### Per-pixel figures

Preset 28's 66.32 ms shader scope is about **6.40 us or 3,839 cycles per
pixel** across the 10,368-pixel quadrant.

## Column-ISR / DMA marshaling cost

```
isr_wake        2300/frame  min/avg/max 0.61/1.84/18.29 us  cpu 3.39%
isr_pack         288/frame  min/avg/max 6.31/7.05/9.36 us   cpu 1.62%
isr_dma_submit   288/frame  min/avg/max 0.66/0.91/8.22 us   cpu 0.21%
```

The ISR share is 5.22%, leaving about 59.24 ms of the display window for
rendering. Preset 28 needs about a **1.12x** speedup to fit that adjusted
budget.

## Summary ranking

1. `sb_shader_draw`, Peirce primitive lattice - 91.84 ms held and 97.47 ms peak.
2. `sb_shader_draw`, dodecahedral grid - 66.32 ms held and 72.74 ms peak.
3. `sb_shader_draw`, kaleidoscope edge-fade liquid - 59.20 ms held.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` has no visible subtree in this shallow ShaderBall capture.
- Global O3 is a single-effect ceiling, not a shippable full-roster image.
- Fast-cycle dwell compression changes duration, not per-frame work.
- The capture ran from clean source commit `4f00b8d6`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.

## Global O3 vs selective O3

Global O3 lowers preset 28's clean shader cost from 74.25 to 66.32 ms
(**1.12x**) and the pass peak from 370.11 to 97.47 ms. It adds **28,440 B
FLASH** and **25,472 B ITCM** to the single-effect profile image. Neither
configuration keeps preset 28 at 16 fps.
