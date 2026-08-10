# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-09, **-O3**)

Point-in-time global-O3 reference paired with the
[`shipping report`](../shipping/profile_gsreactiondiffusion_teensy_2026-08-09.md).
Raw capture: `build/prof/gsreactiondiffusion_o3.log`. This replaces the
2026-08-03 reference and includes the per-reaction generative-palette rebake.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile_o3` env: GCC 15.2.1, global `-O3 -ffast-math -fno-finite-math-only`, newlib-nano, DMA LEDs |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `fc1f330a` plus the uncommitted palette-cycle patch later landed as `d3cadd40` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, epoch stretched to 1200 revs |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:68832, data:238740, headers:8844` / `RAM1:
variables:315136, code:38776, padding:26760, free:143616` / `RAM2:
variables:520064, free:4224`.

Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1
20251203`. Profile ELF SHA-256:
`4958aacca14da8a3fb3232c9086333390e608a957adfa91294f35d38863c11b7`.

Exactness cross-check: window frames 1185–1216 root counter cycles divided by
600 MHz match measured wall time within **1.7 ppm**
(`tools/parse_profile.py build/prof/gsreactiondiffusion_o3.log validate`).

## Frame cadence

**Pass aggregate**: `grd_render` averages 50.032 ms/frame across 64 windows.
The worst window mean is 56.14 ms (frames 1185–1216), peak frame render is
**56.972 ms** (frames 353–384), and **0/2048 frames (0%)** spill.

The 62.5 ms display window is held in every morphology phase. The peak retains
5.528 ms of display-window margin. The effect renders one approximately
10,368-pixel quadrant with four samples per pixel; `canvas_buffer_wait` is the
round-up idle to the next display flip.

## Phase-by-phase readout

The capture covers growth, dense reaction, dissolve, palette generation and
rebake, clear, and reseed across four complete transition boundaries.

### Dense plateau (window frames 353–384, peak frame in the capture)

```text
frame                  62.46 ms  37.48 Mcyc  100%
  grd_render           55.98 ms  33.59 Mcyc   89%
    grd_rasterize      43.51 ms  26.11 Mcyc   77%
      grd_shader_draw  40.39 ms  24.24 Mcyc   92%
    grd_simulate       11.95 ms   7.17 Mcyc   21%
  rd_timeline_step     31.4 us   18.9 kcyc     0%
  canvas_clear         89.7 us   53.8 kcyc     0%
  canvas_buffer_wait    6.36 ms   3.81 Mcyc   10%
```

Wall min/avg/max = 61.144/62.459/63.854 ms. Global O3 grows the dense shader
to 40.39 ms while physics remains nearly fixed at 11.95 ms. Raster work
outside the shader costs 3.12 ms.

### Dissolve/reseed with palette rebake (window frames 449–480)

```text
frame                  62.13 ms  37.28 Mcyc  100%
  grd_render           29.67 ms  17.80 Mcyc   47%
    grd_palette_rebake 76.4 us   45.9 kcyc     0%  x0.031  2446 us/call
    grd_rasterize      17.55 ms  10.53 Mcyc   59%
      grd_shader_draw  13.42 ms   8.05 Mcyc   76%
    grd_simulate       11.94 ms   7.17 Mcyc   40%
  rd_timeline_step     30.7 us   18.4 kcyc     0%
  canvas_clear         86.1 us   51.6 kcyc     0%
  canvas_buffer_wait   32.35 ms  19.41 Mcyc   52%
```

Wall min/avg/max = 56.434/62.129/69.756 ms. The one palette generation and
LUT refill costs **2.446 ms**. Across all four captured transitions the rebake
spans 2.396–2.446 ms and never causes a cadence spill.

### Per-pixel figures

| Regime | Shader Mcyc/frame | Cycles/sub-sample | Cycles/pixel |
|---|---:|---:|---:|
| Dense plateau | 24.24 | 584.4 | 2,337.6 |
| Dissolve/reseed | 8.05 | 194.2 | 776.9 |

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.57/1.92/14.89 us  cpu 3.54%
isr_pack         144/frame  min/avg/max 6.52/7.58/10.70 us  cpu 1.74%
isr_dma_submit   144/frame  min/avg/max 0.64/0.92/7.13 us   cpu 0.21%
```

- Packing remains the dominant CPU-side DMA preparation cost.
- The three ISR groups consume 5.49% of wall time in the peak window.
- Async wire-transfer time is excluded from the cycle scopes.

## Summary ranking

1. `grd_shader_draw` — 40.39 ms at dense coverage.
2. `grd_simulate` — 11.95 ms, nearly coverage-independent.
3. Raster prologue — 3.12 ms.
4. `grd_palette_rebake` — 2.40–2.45 ms once per reaction cycle.

## Caveats

- Scope time includes live ISR execution.
- This is a single-effect compiler reference, not a shippable full-roster
  configuration.
- `filter_blend` is not used by this shader.
- Epoch stretch changes capture lifetime only.
- The capture used the uncommitted palette-cycle patch atop `fc1f330a`; the
  identical source change landed as `d3cadd40` after rebasing onto live master.

## Harness

`targets/Profile/Profile.ino`, target `GSReactionDiffusion`, window 32; use the
reproduce command above for the locked build, flash, capture, and validation.

## Global O3 vs selective O3

Global O3 peaks at 56.972 ms versus 56.277 ms for shipping, **0.695 ms
(1.2%) slower**. FLASH code grows by **11,632 B** and ITCM code by
**10,624 B**. The palette rebake cost is effectively unchanged, and the
selective-O3 image remains both smaller and faster.
