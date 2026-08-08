# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-03, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).
Raw capture: `build/prof/gsreactiondiffusion_ship.log`. This replaces the
2026-07-30 report and includes exact integer opaque-SSAA accumulation.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM4, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: GCC 15.2.1, `-Os`, newlib-nano, DMA LEDs, with selective O3 on `draw_grid`, `shade_pixel`, `refine_render_center`, and `step_physics` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `77e86d77` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, epoch stretched to 1200 revs |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:53948, data:630408, headers:8884` / `RAM1:
variables:314496, code:30824, padding:1944, free:177024` / `RAM2:
variables:520064, free:4224`.

Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1
20251203`. Profile ELF SHA-256:
`c93e0bb60c5c74ff822fd627a1367bec953ccff2f073789085703daa77fe71e4`.

Exactness cross-check: window frames 1185–1216 root counter cycles divided by
600 MHz match measured wall time within **0.9 ppm**
(`tools/parse_profile.py build/prof/gsreactiondiffusion_ship.log validate`).

## Frame cadence

**Pass aggregate**: `grd_render` averages 48.480 ms/frame across 64 windows.
The worst window mean is 54.75 ms (frames 1185–1216), peak frame render is
**55.708 ms**, and **0/2048 frames (0%)** spill.

A display window is 62.5 ms. The effect renders one approximately
10,368-pixel quadrant with four sub-samples per pixel. Every captured phase
holds 16 fps; the worst frame has 6.792 ms of display-window margin and is
0.708 ms above the 55 ms stretch target. `canvas_buffer_wait` is display-flip
round-up idle.

## Phase-by-phase readout

GSReactionDiffusion grows from seeded B clusters, reaches a dense plateau,
fades dissolving nodes toward rest, clears, and reseeds. The capture contains
multiple full morphology cycles.

### Dense plateau (window frames 833–864, peak frame in the capture)

```text
frame                  62.40 ms  37.44 Mcyc  100%
  grd_render           54.38 ms  32.63 Mcyc   87%
    grd_rasterize      41.88 ms  25.13 Mcyc   77%
      grd_shader_draw  38.24 ms  22.94 Mcyc   91%
    grd_simulate       11.97 ms   7.18 Mcyc   22%
  rd_timeline_step     37.3 us   22.4 kcyc     0%
  canvas_clear         89.7 us   53.8 kcyc     0%
  canvas_buffer_wait    7.89 ms   4.73 Mcyc   13%
```

Wall min/avg/max = 61.145/62.395/63.676 ms. Raster coverage sets the peak:
the shader consumes 38.24 ms, while six scaled physics integrations remain
nearly constant at 11.97 ms. Raster work outside the shader costs 3.64 ms.

### Dissolve/reseed trough (window frames 1761–1792)

```text
frame                  62.50 ms  37.50 Mcyc  100%
  grd_render           27.87 ms  16.72 Mcyc   45%
    grd_rasterize      15.86 ms   9.51 Mcyc   57%
      grd_shader_draw  11.78 ms   7.07 Mcyc   74%
    grd_simulate       11.96 ms   7.17 Mcyc   43%
  rd_timeline_step     31.8 us   19.1 kcyc     0%
  canvas_clear         86.1 us   51.6 kcyc     0%
  canvas_buffer_wait   34.51 ms  20.71 Mcyc   55%
```

Wall min/avg/max = 57.788/62.504/68.167 ms. Shader cost falls by 26.46 ms as
the fading frontier reduces hot coverage; physics changes by 0.01 ms.

### Per-pixel figures

| Regime | Shader Mcyc/frame | Cycles/sub-sample | Cycles/pixel |
|---|---:|---:|---:|
| Dense plateau | 22.94 | 553.2 | 2,212.9 |
| Dissolve/reseed trough | 7.07 | 170.5 | 681.9 |

The shader overwrites every quadrant pixel, so there is no `filter_blend`
subtree. Integer quarter accumulation replaces nine float channel conversions
and keeps the opaque four-sample result exact.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1151/frame  min/avg/max 0.64/1.95/12.45 us  cpu 3.58%
isr_pack         144/frame  min/avg/max 6.24/7.29/10.38 us  cpu 1.67%
isr_dma_submit   144/frame  min/avg/max 0.69/0.95/4.48 us   cpu 0.21%
```

- Submit remains small next to color/interleave packing; wire transfer is
  asynchronous and absent from these CPU timings.
- The three ISR groups consume 5.46% of wall time in the peak window.
- The captured render peak is 6.792 ms below the display deadline.

## Summary ranking

1. `grd_shader_draw` — 38.24 ms, 70.3% of dense render: four shared-stencil
   interpolation samples over the hot reaction field.
2. `grd_simulate` — 11.97 ms, 22.0%: six scaled Gray-Scott integrations plus
   Q16 conversion and stabilization accumulation.
3. Raster prologue — 3.64 ms: lattice orientation and hot-flag construction.

The preceding shipping capture peaked at 58.606 ms. Exact opaque integer
accumulation reduces the peak by **2.898 ms (4.9%)**. Relative to the original
94.830 ms capture, the current peak is down **39.122 ms (41.3%)**. Exhaustive
testing over all 65,536 channel values proves the four-sample conversion and
rounding remain bit-exact for the fixed opaque palette.

## Caveats

- Cycle scopes include ISR time because CYCCNT free-runs.
- `filter_blend` is not used by this shader.
- Selective O3 reaches the typed grid driver, GS pixel/refinement kernels, and
  physics loop; the surrounding image remains size-oriented.
- `HS_PROFILE_EPOCH_REVS=1200` changes capture lifetime only.
- Six scaled integrations preserve the simulated interval, not bit identity
  with the older ten-stage integration schedule.
- The capture came from clean detached worktree tip `77e86d77`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GSReactionDiffusion`,
`HS_PROFILE_WINDOW=32`; the reproduce command builds, locks, flashes, captures,
validates, and archives the exact artifacts.
