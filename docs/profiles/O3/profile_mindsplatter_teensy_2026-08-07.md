# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-08-07, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_mindsplatter_teensy_2026-08-07.md).
Raw capture: `build/prof/mindsplatter_o3.log`. This replaces the earlier
2026-08-07 four-preset capture with the current eight-preset full cycle.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `3ec9284c` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 110 s capture; all eight presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile_o3 110 16` |

Image size: `FLASH: code:83192, data:541024, headers:8616` / `RAM1:
variables:315232, code:59784, padding:5752, free:143520` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 193–208 root counter cycles divided by
600 MHz match measured wall time within **1.2 ppm**
(`tools/parse_profile.py build/prof/mindsplatter_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `msp_draw_particles`
averages 19.75 ms/frame; total render averages 21.82 ms/frame, peaks at
**38.78 ms** in frames 1713–1728, and spills **0/1728 frames (0%)**.

All eight presets hold 16 fps against the 62.5 ms display window, leaving
23.72 ms at the worst frame. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

The preset timer advances every 160 frames and lerps for 48 frames. Ten
markers cover the full eight-preset wrap and the beginning of a second pass.

### Saturated transition/hold regime (frames 1713–1728, worst of capture)

```
frame                   62.40 ms  37.44 Mcyc  100%
  msp_draw_particles    24.39 ms  14.63 Mcyc   39%
    msp_particle_scan   24.38 ms  14.63 Mcyc   39%
      plot_ps_raster    18.68 ms  11.21 Mcyc   30%  x272  68.7 us/trail
      plot_ps_gate       3.29 ms   1.97 Mcyc    5%  x924  3.6 us/call
        plot_ps_cart.. 594.5 us  356.7 kcyc     1%
      plot_ps_tween      1.86 ms   1.12 Mcyc    3%
      plot_ps_deferred 337.4 us  202.5 kcyc     1%
  msp_particle_step      1.88 ms   1.13 Mcyc    3%
  msp_timeline_step     47.4 us   28.5 kcyc     0%
  canvas_clear          86.3 us   51.8 kcyc     0%
  canvas_buffer_wait    36.00 ms  21.60 Mcyc   58%
```

Wall min/avg/max = 44.39/62.40/76.63 ms. Global O3 reduces the aggregate
particle-draw average by 5.3%, while workload variation leaves the pass peak
nearly unchanged.

### Per-preset table

The initial unlabeled windows and `Preset: 1/8` windows are combined. Rows use
the highest clean modal-call-count window and are ranked by particle-draw cost;
the wrap from preset 8 to 1 is present.

| # | clean windows | particle draw ms | render ms | fps |
|---:|--:|--:|--:|--:|
| 2 | 20 | 29.82 | 32.11 | 16 |
| 6 | 10 | 29.24 | 31.38 | 16 |
| 5 | 10 | 29.06 | 31.21 | 16 |
| 3 | 18 | 26.66 | 28.82 | 16 |
| 1 | 20 | 25.73 | 27.87 | 16 |
| 8 | 10 | 24.91 | 27.10 | 16 |
| 4 | 10 | 24.09 | 26.22 | 16 |
| 7 | 10 | 22.90 | 24.99 | 16 |

### Per-pixel figures

The direct-AA path has no `filter_blend` counter. The representative worst
window rasterizes 272 accepted trail segments/frame at 68.7 us/trail.

## Column-ISR / DMA marshaling cost

```
isr_wake        1151/frame  min/avg/max 0.53/1.82/11.41 us  cpu 3.35%
isr_pack         144/frame  min/avg/max 6.32/6.99/9.38 us   cpu 1.61%
isr_dma_submit   144/frame  min/avg/max 0.65/0.91/1.04 us   cpu 0.20%
```

- DMA submission averages 13% of packing CPU cost.
- LED transfer continues asynchronously.
- Combined ISR share is 5.16%; the worst frame still holds one-window cadence.

## Summary ranking

1. `plot_ps_raster` — 18.68 ms/frame, 30% of the representative frame.
2. `plot_ps_gate` — 3.29 ms/frame.
3. `msp_particle_step` — 1.88 ms/frame.
4. `plot_ps_tween` — 1.86 ms/frame.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` does not appear on the direct-AA path.
- Global O3 is a single-effect ceiling, not a shippable roster image.
- No dwell-compression or epoch-stretch knobs were used.
- The capture ran from clean tip `3ec9284c`; later tip `39d72f18` changes only ShaderWorkbench palette recipes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MindSplatter`,
`HS_PROFILE_WINDOW=16`.

## Global O3 vs selective O3

The worst-frame peak is effectively flat: 38.95 ms shipping versus 38.78 ms
global O3 (1.00×). The aggregate render average falls from 23.06 to 21.82 ms.
Global O3 adds **21,464 B FLASH** and **18,832 B ITCM**.
