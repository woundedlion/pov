# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw capture:
`build/prof/shaderball_ship.log`. This replaces the 12-preset capture: direct
grid presets now skip the unused coupled-pattern sample, the tiny mixed-preset
lens amount is zero, and the roster includes the new fractional-pattern preset.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; the `Scan::Shader::draw` selective-O3 region is active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, capture tip `99dabc17` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3600`; all 13 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 400 16 "-D HS_PROFILE_EPOCH_REVS=3600"` |

Image size: `FLASH: code:62928, data:147636, headers:8572` / `RAM1:
variables:314976, code:30136, padding:2632, free:176544` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 6321-6336 root counter cycles divided by
600 MHz match measured wall time within **4.8 ppm**.

## Frame cadence

`sb_shader_draw` averages **41.60 ms/frame** and total render averages
**43.68 ms/frame**. Peak frame render is **55.97 ms** and **0/6368 frames**
spill. Every preset holds 16 fps; the worst frame retains 6.53 ms of margin
below the 62.5 ms display window.

ShaderBall shades one 72x144 quadrant, or 10,368 pixels. The
`canvas_buffer_wait` scope is cadence-quantized display-sync idle.

## Phase-by-phase readout

### Held liquid regime (frames 321-336)

```
frame                  62.45 ms  37.47 Mcyc  100%
  sb_shader_draw       30.42 ms  18.25 Mcyc   48%
  sb_timeline_step      65 us    38.9 kcyc     0%
  canvas_clear          86 us    51.7 kcyc     0%
  canvas_buffer_wait   29.97 ms  17.98 Mcyc   47%
```

Render averages 32.48 ms and peaks at 32.87 ms in this window.

### Fractional-pattern preset (frames 6305-6320, worst of capture)

```
frame                  62.47 ms  37.48 Mcyc  100%
  sb_shader_draw       52.73 ms  31.64 Mcyc   84%
  sb_timeline_step      67 us    40.0 kcyc     0%
  canvas_clear          89 us    53.1 kcyc     0%
  canvas_buffer_wait    7.62 ms   4.57 Mcyc   12%
```

Wall min/avg/max is 60.97/62.47/64.83 ms. Render averages 54.84 ms and
peaks at 55.97 ms. Preset 6 has `pattern_mix=0.536`, so it evaluates both
pattern carriers; it is the most expensive preset but does not spill.

### Per-preset table

Initial unlabeled frames are folded into preset 1. Peak render uses exact
per-frame ownership from `parse_profile.py ... buckets`.

| # | peak render ms | cadence | spilled/frames |
|---:|--:|---:|---:|
| 6 | 55.97 | 16 fps | 0/657 |
| 7 | 55.25 | 16 fps | 0/480 |
| 1 | 52.31 | 16 fps | 0/545 |
| 13 | 51.44 | 16 fps | 0/480 |
| 12 | 50.32 | 16 fps | 0/480 |
| 5 | 49.11 | 16 fps | 0/960 |
| 8 | 47.59 | 16 fps | 0/480 |
| 11 | 46.55 | 16 fps | 0/480 |
| 10 | 44.95 | 16 fps | 0/480 |
| 9 | 44.42 | 16 fps | 0/480 |
| 4 | 35.61 | 16 fps | 0/286 |
| 3 | 34.69 | 16 fps | 0/240 |
| 2 | 34.68 | 16 fps | 0/232 |

### Per-pixel figures

The worst shader window costs about **5.09 us/pixel**, or **3,052
cycles/pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.69/1.95/11.33 us  cpu 3.60%
isr_pack         144/frame  min/avg/max 6.24/6.98/9.16 us   cpu 1.60%
isr_dma_submit   144/frame  min/avg/max 0.67/0.94/1.09 us   cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. LED transfer continues
asynchronously, and `filter_blend` is not exposed on this closure-shader path.

## Summary ranking

1. `sb_shader_draw` - 52.73 ms/frame in the worst fractional-pattern window.
2. Palette-cycle and frame work outside the shader - about 2.11 ms/frame.
3. Timeline and clear - below 0.2 ms/frame combined.

The direct-pattern endpoint now returns before building the unused coupled
carrier. Full-cycle peak render fell from 63.05 to 55.97 ms (**1.13x**) even
with the added mixed-pattern preset, and spills fell from 3/6368 to zero.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- The profile image compresses the epoch only; it does not change per-frame
  effect work.
- The capture ran from the clean isolated implementation commit shown above.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
