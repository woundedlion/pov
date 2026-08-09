# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw capture:
`build/prof/shaderball_ship.log`. This replaces the 12-preset capture: direct
grid presets skip the unused coupled-pattern sample, the tiny mixed-preset lens
amount is zero, and the hot closure and pattern evaluator use selective O3.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; `Scan::Shader::draw`, the ShaderBall closure/pattern evaluator, and their selective-O3 leaves are active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, capture tip `17fbf726` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3600`; all 13 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 400 16 "-D HS_PROFILE_EPOCH_REVS=3600"` |

Image size: `FLASH: code:63256, data:147648, headers:8232` / `RAM1:
variables:314976, code:30456, padding:2312, free:176544` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 6321-6336 root counter cycles divided by
600 MHz match measured wall time within **4.0 ppm**.

## Frame cadence

`sb_shader_draw` averages **39.64 ms/frame** and total render averages
**41.72 ms/frame**. Peak frame render is **53.99 ms** and **0/6368 frames**
spill. Every preset holds 16 fps; the worst frame retains 8.51 ms of margin
below the 62.5 ms display window.

ShaderBall shades one 72x144 quadrant, or 10,368 pixels. The
`canvas_buffer_wait` scope is cadence-quantized display-sync idle.

## Phase-by-phase readout

### Held liquid regime (frames 321-336)

```
frame                  62.45 ms  37.47 Mcyc  100%
  sb_shader_draw       28.90 ms  17.34 Mcyc   46%
  sb_timeline_step      64 us    38.4 kcyc     0%
  canvas_clear          86 us    51.5 kcyc     0%
  canvas_buffer_wait   31.50 ms  18.90 Mcyc   50%
```

Render averages 30.96 ms and peaks at 31.29 ms in this window.

### Fractional-pattern preset (frames 6305-6320, worst of capture)

```
frame                  62.46 ms  37.48 Mcyc  100%
  sb_shader_draw       50.81 ms  30.49 Mcyc   81%
  sb_timeline_step      60 us    36.0 kcyc     0%
  canvas_clear          88 us    53.1 kcyc     0%
  canvas_buffer_wait    9.55 ms   5.73 Mcyc   15%
```

Wall min/avg/max is 61.01/62.46/64.78 ms. Render averages 52.92 ms and
peaks at 53.99 ms. Preset 6 has `pattern_mix=0.536`, so it evaluates both
pattern carriers; it is the most expensive preset but does not spill.

### Per-preset table

Initial unlabeled frames are folded into preset 1. Peak render uses exact
per-frame ownership from `parse_profile.py ... buckets`.

| # | peak render ms | cadence | spilled/frames |
|---:|--:|---:|---:|
| 6 | 53.99 | 16 fps | 0/657 |
| 7 | 53.39 | 16 fps | 0/480 |
| 1 | 50.31 | 16 fps | 0/545 |
| 13 | 49.73 | 16 fps | 0/480 |
| 12 | 48.58 | 16 fps | 0/480 |
| 5 | 47.28 | 16 fps | 0/960 |
| 8 | 45.80 | 16 fps | 0/480 |
| 11 | 44.77 | 16 fps | 0/480 |
| 10 | 42.91 | 16 fps | 0/480 |
| 9 | 42.27 | 16 fps | 0/480 |
| 4 | 33.71 | 16 fps | 0/286 |
| 2 | 32.50 | 16 fps | 0/232 |
| 3 | 32.41 | 16 fps | 0/240 |

### Per-pixel figures

The worst shader window costs about **4.90 us/pixel**, or **2,941
cycles/pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.69/1.95/11.39 us  cpu 3.60%
isr_pack         144/frame  min/avg/max 6.24/6.98/9.23 us   cpu 1.60%
isr_dma_submit   144/frame  min/avg/max 0.65/0.94/1.09 us   cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. LED transfer continues
asynchronously, and `filter_blend` is not exposed on this closure-shader path.

## Summary ranking

1. `sb_shader_draw` - 50.81 ms/frame in the worst fractional-pattern window.
2. Palette-cycle and frame work outside the shader - about 2.10 ms/frame.
3. Timeline and clear - below 0.2 ms/frame combined.

Selective O3 on the ShaderBall pixel closure and pattern evaluator lowers the
full-cycle peak from 55.97 to 53.99 ms and keeps all 13 presets spill-free.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- The profile image compresses the epoch only; it does not change per-frame
  effect work.
- The capture ran from the clean isolated implementation commit shown above.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
