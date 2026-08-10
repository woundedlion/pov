# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-09, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).
Raw capture: `build/prof/gsreactiondiffusion_ship.log`. This replaces the
2026-08-03 report and includes the per-reaction generative-palette rebake.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: GCC 15.2.1, `-Os`, newlib-nano, DMA LEDs, with selective O3 on `draw_grid`, `shade_pixel`, `refine_render_center`, and `step_physics` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `fc1f330a` plus the uncommitted palette-cycle patch later landed as `d3cadd40` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, epoch stretched to 1200 revs |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:57200, data:238728, headers:8200` / `RAM1:
variables:315104, code:28152, padding:4616, free:176416` / `RAM2:
variables:520064, free:4224`.

Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1
20251203`. Profile ELF SHA-256:
`f87cd387987388bff3062a8ee60095e7157cb5b1fd9b1ae06dfe9f1be20bbda6`.

Exactness cross-check: window frames 1601–1632 root counter cycles divided by
600 MHz match measured wall time within **2.3 ppm**
(`tools/parse_profile.py build/prof/gsreactiondiffusion_ship.log validate`).

## Frame cadence

**Pass aggregate**: `grd_render` averages 49.365 ms/frame across 64 windows.
The worst window mean is 55.05 ms (frames 1601–1632), peak frame render is
**56.277 ms**, and **0/2048 frames (0%)** spill.

A display window is 62.5 ms. The effect renders one approximately
10,368-pixel quadrant with four sub-samples per pixel. Every captured phase
holds 16 fps; the worst frame retains 6.223 ms of display-window margin.
`canvas_buffer_wait` is display-flip round-up idle.

## Phase-by-phase readout

GSReactionDiffusion grows from seeded B clusters, reaches a dense plateau,
dissolves to rest, generates and bakes a new palette, then reseeds. The capture
contains four palette rebakes across multiple full morphology cycles.

### Dense plateau (window frames 1601–1632, peak frame in the capture)

```text
frame                  62.49 ms  37.49 Mcyc  100%
  grd_render           55.05 ms  33.03 Mcyc   88%
    grd_rasterize      42.78 ms  25.67 Mcyc   77%
      grd_shader_draw  39.31 ms  23.58 Mcyc   91%
    grd_simulate       11.97 ms   7.18 Mcyc   21%
  rd_timeline_step     34.7 us   20.8 kcyc     0%
  canvas_clear         90.2 us   54.1 kcyc     0%
  canvas_buffer_wait    7.31 ms   4.39 Mcyc   11%
```

Wall min/avg/max = 60.900/62.487/64.035 ms. Raster coverage sets the peak:
the shader consumes 39.31 ms, while six scaled physics integrations remain
nearly constant at 11.97 ms. Raster work outside the shader costs 3.48 ms.

### Dissolve/reseed with palette rebake (window frames 449–480)

```text
frame                  62.13 ms  37.28 Mcyc  100%
  grd_render           29.37 ms  17.62 Mcyc   47%
    grd_palette_rebake 73.4 us   44.0 kcyc     0%  x0.031  2348 us/call
    grd_rasterize      17.24 ms  10.34 Mcyc   58%
      grd_shader_draw  13.15 ms   7.89 Mcyc   76%
    grd_simulate       11.96 ms   7.17 Mcyc   40%
  rd_timeline_step     33.5 us   20.1 kcyc     0%
  canvas_clear         86.2 us   51.7 kcyc     0%
  canvas_buffer_wait   32.64 ms  19.58 Mcyc   52%
```

Wall min/avg/max = 56.584/62.130/69.565 ms. The one palette generation and
LUT refill costs **2.348 ms**; its per-frame average is small because it runs
once in the 32-frame window. Across all four captured transitions the rebake
spans 2.348–2.417 ms and never causes a cadence spill.

### Per-pixel figures

| Regime | Shader Mcyc/frame | Cycles/sub-sample | Cycles/pixel |
|---|---:|---:|---:|
| Dense plateau | 23.58 | 568.7 | 2,274.8 |
| Dissolve/reseed | 7.89 | 190.2 | 761.0 |

The shader overwrites every quadrant pixel, so there is no `filter_blend`
subtree. Coverage culling reduces shader work during the dissolve while the
physics floor and one-time palette refill remain visible separately.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1153/frame  min/avg/max 0.67/2.03/13.72 us  cpu 3.74%
isr_pack         144/frame  min/avg/max 6.51/7.69/10.99 us  cpu 1.77%
isr_dma_submit   144/frame  min/avg/max 0.72/0.95/4.77 us   cpu 0.21%
```

- Submit remains small next to color/interleave packing; wire transfer is
  asynchronous and absent from these CPU timings.
- The three ISR groups consume 5.72% of wall time in the peak window.
- The captured render peak is 6.223 ms below the display deadline.

## Summary ranking

1. `grd_shader_draw` — 39.31 ms at dense coverage: four shared-stencil
   interpolation samples over the hot reaction field.
2. `grd_simulate` — 11.97 ms: six scaled Gray-Scott integrations plus Q16
   conversion and stabilization accumulation.
3. Raster prologue — 3.48 ms: lattice orientation and hot-flag construction.
4. `grd_palette_rebake` — 2.35–2.42 ms once per reaction cycle.

The palette rebake remains well inside the available transition-frame budget;
cycling prebaked flash palettes is not needed.

## Caveats

- Cycle scopes include ISR time because CYCCNT free-runs.
- `filter_blend` is not used by this shader.
- Selective O3 reaches the typed grid driver, GS pixel/refinement kernels, and
  physics loop; the surrounding image remains size-oriented.
- `HS_PROFILE_EPOCH_REVS=1200` changes capture lifetime only.
- The capture used the uncommitted palette-cycle patch atop `fc1f330a`; the
  identical source change landed as `d3cadd40` after rebasing onto live master.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GSReactionDiffusion`,
`HS_PROFILE_WINDOW=32`; `just profile GSReactionDiffusion` builds, locks,
flashes, captures, and validates the shipping configuration.
