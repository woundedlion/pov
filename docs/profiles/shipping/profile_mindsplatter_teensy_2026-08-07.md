# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-08-07, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`). Raw
capture: `build/prof/mindsplatter_ship.log`. This replaces the earlier
2026-08-07 four-preset capture with the current eight-preset full cycle.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; the shared Plot raster selective-O3 region is active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `3ec9284c` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 110 s capture; all eight presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile 110 16` |

Image size: `FLASH: code:61728, data:541016, headers:8584` / `RAM1:
variables:315200, code:40952, padding:24584, free:143552` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 193–208 root counter cycles divided by
600 MHz match measured wall time within **1.6 ppm**
(`tools/parse_profile.py build/prof/mindsplatter_ship.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `msp_draw_particles`
averages 20.86 ms/frame; total render averages 23.06 ms/frame, peaks at
**38.95 ms** in frames 193–208, and spills **0/1728 frames (0%)**.

A display window is 62.5 ms; the effect renders one 72×144 quadrant. All eight
presets hold 16 fps, with 23.55 ms remaining at the worst frame. The
`canvas_buffer_wait` scope is display-sync idle.

## Phase-by-phase readout

The preset timer advances every 160 frames and lerps for the first 48 frames
of each interval. Ten markers cover all eight presets, the wrap, and the start
of a second pass.

### Saturated transition/hold regime (frames 193–208, worst of capture)

```
frame                   62.47 ms  37.48 Mcyc  100%
  msp_draw_particles    31.60 ms  18.96 Mcyc   51%
    msp_particle_scan   31.60 ms  18.96 Mcyc   51%
      plot_ps_raster    22.40 ms  13.44 Mcyc   36%  x348  64.3 us/trail
      plot_ps_gate       4.71 ms   2.82 Mcyc    8%  x1082  4.4 us/call
        plot_ps_cart.. 733.4 us  440.1 kcyc     1%
      plot_ps_tween      3.46 ms   2.08 Mcyc    6%
      plot_ps_deferred 440.8 us  264.5 kcyc     1%
  msp_particle_step      2.28 ms   1.37 Mcyc    4%
  msp_timeline_step     72.6 us   43.6 kcyc     0%
  canvas_clear          86.8 us   52.1 kcyc     0%
  canvas_buffer_wait    28.42 ms  17.05 Mcyc   46%
```

Wall min/avg/max = 53.59/62.47/71.35 ms. Direct AA trail rasterization is the
dominant cost. The varying particle population creates more spread inside a
preset than the fixed timeline and clear costs do.

### Per-preset table

The initial unlabeled windows and later `Preset: 1/8` windows are combined as
preset 1. Rows use the highest clean modal-call-count window for each preset
and are ranked by particle-draw cost. The log wraps through preset 8 to 1.

| # | clean windows | particle draw ms | render ms | fps |
|---:|--:|--:|--:|--:|
| 2 | 20 | 31.60 | 34.04 | 16 |
| 1 | 20 | 28.98 | 31.26 | 16 |
| 6 | 10 | 28.16 | 30.48 | 16 |
| 5 | 10 | 26.95 | 29.25 | 16 |
| 8 | 10 | 26.03 | 28.35 | 16 |
| 4 | 10 | 25.88 | 28.17 | 16 |
| 3 | 18 | 25.74 | 28.05 | 16 |
| 7 | 10 | 24.53 | 26.77 | 16 |

### Per-pixel figures

The direct-AA path does not enter `filter_blend`, so blended-pixel figures do
not apply. The worst window rasterizes 348 accepted trail segments per frame
at 64.3 us/trail and gates 1,082 history samples per frame.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.74/1.94/17.83 us  cpu 3.58%
isr_pack         144/frame  min/avg/max 6.24/6.96/9.34 us   cpu 1.60%
isr_dma_submit   144/frame  min/avg/max 0.78/0.95/10.06 us  cpu 0.21%
```

- DMA submission averages 14% of the CPU cost of packing.
- LED transfer continues asynchronously after submission.
- The 5.39% combined ISR share is already absorbed by the render counters;
  the worst frame remains 1.60× inside the 62.5 ms display budget.

## Summary ranking

1. `plot_ps_raster` — 22.40 ms/frame, 36% of the representative frame.
2. `plot_ps_gate` — 4.71 ms/frame across 1,082 history samples.
3. `plot_ps_tween` — 3.46 ms/frame.
4. `msp_particle_step` — 2.28 ms/frame.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` does not appear on the direct-AA path.
- Plot's shared raster selective-O3 region is active; effect-local wrappers remain `-Os`.
- No dwell-compression or epoch-stretch knobs were used.
- The capture ran from clean tip `3ec9284c`; later tip `39d72f18` changes only ShaderBall palette recipes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MindSplatter`,
`HS_PROFILE_WINDOW=16`; `just profile MindSplatter` builds, flashes, and captures.
