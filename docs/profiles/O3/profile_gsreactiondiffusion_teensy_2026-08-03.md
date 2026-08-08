# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-03, **-O3**)

Point-in-time global-O3 reference paired with the
[`shipping report`](../shipping/profile_gsreactiondiffusion_teensy_2026-08-03.md).
Raw capture: `build/prof/gsreactiondiffusion_o3.log`. This replaces the
2026-07-30 reference with the final opaque-accumulation path.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM4, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile_o3` env: GCC 15.2.1, global `-O3 -ffast-math -fno-finite-math-only`, newlib-nano, DMA LEDs |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `77e86d77` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, epoch stretched to 1200 revs |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:71476, data:630376, headers:8796` / `RAM1:
variables:314496, code:47256, padding:18280, free:144256` / `RAM2:
variables:520064, free:4224`.

Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1
20251203`. Profile ELF SHA-256:
`3324c5e751e0e5312f1fc617969973522059881f1bc202dfd84ea9473269dfdb`.

Exactness cross-check: window frames 1601–1632 root counter cycles divided by
600 MHz match measured wall time within **2.1 ppm**.

## Frame cadence

**Pass aggregate**: `grd_render` averages 49.801 ms/frame across 64 windows.
The worst window mean is 55.55 ms (frames 1601–1632), peak frame render is
**56.677 ms**, and **0/2048 frames (0%)** spill.

The 62.5 ms display window is held in every morphology phase. The peak has
5.823 ms of display-window margin and is 1.677 ms above the 55 ms stretch
target. The effect renders one approximately 10,368-pixel quadrant with four
samples per pixel.

## Phase-by-phase readout

The capture covers growth, dense reaction, dissolve, clear, and reseed across
multiple full morphology cycles.

### Dense plateau (window frames 353–384, peak frame in the capture)

```text
frame                  62.42 ms  37.45 Mcyc  100%
  grd_render           55.53 ms  33.32 Mcyc   89%
    grd_rasterize      43.09 ms  25.86 Mcyc   78%
      grd_shader_draw  40.17 ms  24.10 Mcyc   93%
    grd_simulate       11.91 ms   7.14 Mcyc   21%
  rd_timeline_step     26.0 us   15.6 kcyc     0%
  canvas_clear         89.1 us   53.5 kcyc     0%
  canvas_buffer_wait    6.78 ms   4.07 Mcyc   11%
```

Wall min/avg/max = 61.153/62.422/63.838 ms. Global O3 makes physics 0.06 ms
faster but grows the dense shader by 1.93 ms versus shipping.

### Dissolve/reseed trough (window frames 1281–1312)

```text
frame                  61.94 ms  37.17 Mcyc  100%
  grd_render           27.34 ms  16.40 Mcyc   44%
    grd_rasterize      15.32 ms   9.19 Mcyc   56%
      grd_shader_draw  11.23 ms   6.74 Mcyc   73%
    grd_simulate       11.90 ms   7.14 Mcyc   44%
  rd_timeline_step     21.3 us   12.8 kcyc     0%
  canvas_clear         88.7 us   53.2 kcyc     0%
  canvas_buffer_wait   34.50 ms  20.70 Mcyc   56%
```

Wall min/avg/max = 55.115/61.941/70.084 ms. The shader falls to 11.23 ms as
coverage clears, while the fixed simulation floor remains 11.90 ms.

### Per-pixel figures

| Regime | Shader Mcyc/frame | Cycles/sub-sample | Cycles/pixel |
|---|---:|---:|---:|
| Dense plateau | 24.10 | 581.2 | 2,324.8 |
| Dissolve/reseed trough | 6.74 | 162.5 | 650.1 |

## Column-ISR / DMA marshaling cost

```text
isr_wake        1151/frame  min/avg/max 0.55/1.91/12.69 us  cpu 3.52%
isr_pack         144/frame  min/avg/max 6.62/7.61/10.58 us  cpu 1.75%
isr_dma_submit   144/frame  min/avg/max 0.77/0.92/3.91 us   cpu 0.21%
```

- Packing remains the dominant CPU-side DMA preparation cost.
- The three ISR groups consume 5.48% of wall time in the peak window.
- Async wire-transfer time is excluded from the cycle scopes.

## Summary ranking

1. `grd_shader_draw` — 40.17 ms at dense coverage.
2. `grd_simulate` — 11.91 ms, nearly coverage-independent.
3. Raster prologue — 2.92 ms.

## Caveats

- Scope time includes live ISR execution.
- This is a single-effect compiler reference, not a shippable full-roster
  configuration.
- `filter_blend` is not used by this shader.
- Epoch stretch changes capture lifetime only.
- Six scaled integrations preserve the simulated interval, not bit identity
  with the older ten-stage integration schedule.
- The capture used clean detached worktree tip `77e86d77`.

## Harness

`targets/Profile/Profile.ino`, target `GSReactionDiffusion`, window 32; use the
reproduce command above for the locked build, flash, capture, and validation.

## Global O3 vs selective O3

Global O3 peaks at 56.677 ms versus 55.708 ms for shipping, **0.969 ms
(1.7%) slower**. FLASH code grows by **17,528 B** and ITCM code by
**16,432 B**. The selective-O3 image is both smaller and faster for the final
opaque-accumulation path.
