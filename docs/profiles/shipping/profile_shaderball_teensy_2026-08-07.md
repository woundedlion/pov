# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-07, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw
capture: `build/prof/shaderball_ship.log`. This is the first ShaderBall report,
replacing the separate retired Liquid2D and Flyby roster rows.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; the closure `Scan::Shader::draw` selective-O3 region is active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, tip `39d72f18` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:61496, data:148448, headers:9192` / `RAM1:
variables:314976, code:30232, padding:2536, free:176544` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 4577–4592 root counter cycles divided by
600 MHz match measured wall time within **1.2 ppm**
(`tools/parse_profile.py build/prof/shaderball_ship.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `sb_shader_draw` averages
53.58 ms/frame; total render averages 58.86 ms/frame, peaks at **97.93 ms**
in frames 4577–4592, and spills **1501/4864 frames (30.9%)**.

One frame shades the 72×144 quadrant (10,368 pixels). The cycle crosses between
16 fps and 8 fps: five preset buckets spill at least one frame, and presets 5,
6, and 1 spend effectively their whole expensive transition at 8 fps. The
`canvas_buffer_wait` scope is the round-up idle to the next display flip.

## Phase-by-phase readout

Presets 1–3 dwell then use 60-frame staggered blends; preset 4 dwells before a
480-frame parallel blend; presets 5–12 chain 480-frame blends. Thirteen markers
cover the full 12-preset wrap and the start of a second cycle.

### Dwelled liquid regime (frames 321–336)

```
frame                  62.45 ms  37.47 Mcyc  100%
  sb_shader_draw       32.70 ms  19.62 Mcyc   52%
  sb_timeline_step      63 us    37.8 kcyc     0%
  canvas_clear          86 us    51.6 kcyc     0%
  canvas_buffer_wait   24.47 ms  14.68 Mcyc   39%
```

This held liquid anchor remains comfortably at 16 fps. The closure shader is
the only material render scope; timeline and clear are negligible.

### Continuous expensive blend (frames 4577–4592, worst of capture)

```
frame                 125.35 ms  75.21 Mcyc  100%
  sb_shader_draw       90.55 ms  54.33 Mcyc   72%
  sb_timeline_step     66.4 us   39.8 kcyc     0%
  canvas_clear         85.3 us   51.2 kcyc     0%
  canvas_buffer_wait   29.63 ms  17.78 Mcyc   24%
```

Wall min/avg/max = 124.21/125.35/126.25 ms. Every frame in this window spills
one display window. The shader is 2.77× costlier than the held liquid anchor;
the rest of the frame changes little.

### Per-preset table

Presets 5–12 have no clean hold: their rows use the highest window nearest each
preset anchor, while cadence comes from exact per-frame ownership in
`parse_profile.py ... buckets`. The initial unlabeled frames are combined with
preset 1, and the wrap to preset 1 is present.

| # | family | shader ms | render ms | cadence |
|---:|---|--:|--:|---:|
| 1 | liquid mild / wrap blend | 90.55 | 95.72 | 8↔16 |
| 6 | grid orbit A | 81.76 | 86.92 | 8 |
| 5 | liquid-to-grid bridge | 74.33 | 79.67 | 8 |
| 12 | grid/liquid B | 66.85 | 72.20 | 16 |
| 4 | liquid direct | 58.70 | 63.98 | 16 |
| 7 | grid orbit B | 58.23 | 63.47 | 8↔16 |
| 11 | grid/liquid A | 56.51 | 61.80 | 8↔16 |
| 10 | grid orbit E | 55.02 | 60.31 | 16 |
| 8 | grid orbit C | 51.47 | 56.74 | 16 |
| 9 | grid orbit D | 44.26 | 49.59 | 16 |
| 3 | liquid fine | 32.70 | 37.98 | 16 |
| 2 | liquid deep | 32.04 | 37.37 | 16 |

### Per-pixel figures

The worst window shades 10,368 pixels/frame. `sb_shader_draw` costs about
8.73 us/pixel, or 5,240 cycles/pixel. The closure path does not expose a
`filter_blend` count.

## Column-ISR / DMA marshaling cost

```
isr_wake        2312/frame  min/avg/max 0.67/1.96/11.62 us  cpu 3.61%
isr_pack         289/frame  min/avg/max 6.37/7.09/9.14 us   cpu 1.63%
isr_dma_submit   289/frame  min/avg/max 0.65/0.94/1.13 us   cpu 0.21%
```

- DMA submission averages 13% of packing CPU cost.
- LED transfer continues asynchronously.
- Combined ISR share is 5.45%. The worst preset needs a **1.57×** render
  speedup to fit one 62.5 ms display window.

## Summary ranking

1. `sb_shader_draw` — 90.55 ms/frame and 72% of the worst representative frame.
2. `canvas_buffer_wait` — 29.63 ms/frame of cadence-quantized sync idle.
3. Timeline and clear together remain below 0.2 ms/frame.

The former standalone Flyby profile peaked at 45.05 ms. ShaderBall's blended
choreography is materially heavier because multiple expensive gates are
simultaneously active during cross-family transitions.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` is not instrumented on this closure-shader path.
- The closure `Scan::Shader::draw` and its hot color/noise callees cross the shipping selective-O3 regions.
- `HS_PROFILE_EPOCH_REVS=3300` only extends effect lifetime; it does not change per-frame cost.
- The capture ran from clean tip `39d72f18` with the mirrored high-contrast liquid palette.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`; `just profile ShaderBall` builds, flashes, and captures.
