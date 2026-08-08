# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-03, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile BZReactionDiffusion`).
Raw capture: `build/prof/bzreactiondiffusion_ship.log`. This replaces the
2026-08-02 report and includes fused compositing, palette-channel hoisting,
and factorized species-color coefficients.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: GCC 15.2.1, `-Os`, newlib-nano, DMA LEDs, with selective O3 on `draw_grid`, `shade_pixel`, the RD kernel, and physics |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | BZReactionDiffusion 288×144, single-entry playlist, tip `77e86d77` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, epoch stretched to 1200 revs |
| Reproduce | `bash tools/profile_one.sh BZReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:53324, data:630284, headers:8608` / `RAM1:
variables:314496, code:30312, padding:2456, free:177024` / `RAM2:
variables:520064, free:4224`.

Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1
20251203`. Profile ELF SHA-256:
`bca19389cbf59f04724c5f01882080de9d7a3a1db0e771f6faf331e449a871b6`.

Exactness cross-check: window frames 801–832 root counter cycles divided by
600 MHz match the measured wall sum within **2.2 ppm**
(`tools/parse_profile.py build/prof/bzreactiondiffusion_ship.log validate`).

## Frame cadence

**Pass aggregate**: `bz_render` averages 49.178 ms/frame across 64 windows.
The worst window mean is 49.88 ms (frames 801–832), peak frame render is
**50.705 ms**, and **0/2048 frames (0%)** spill.

A display window is 62.5 ms; the effect renders one approximately
10,368-pixel quadrant with four sub-samples per pixel. The peak has 11.795 ms
of display margin and 4.295 ms of margin against the 55 ms target.
`canvas_buffer_wait` is display-flip round-up idle.

## Phase-by-phase readout

Kernel occupancy changes slowly with the simulated concentration field; the
peak window represents its dense raster regime.

### Dense raster (window frames 353–384, peak frame in the capture)

```text
frame                62.43 ms  37.46 Mcyc  100%
  bz_render          49.58 ms  29.75 Mcyc   79%
    bz_raster        44.68 ms  26.81 Mcyc   90%
    bz_physics        4.52 ms   2.71 Mcyc    9%
    bz_orient       375.4 us  225.2 kcyc     1%
  rd_timeline_step   33.1 us   19.9 kcyc     0%
  canvas_clear       88.4 us   53.1 kcyc     0%
  canvas_buffer_wait 12.73 ms   7.64 Mcyc   20%
```

Wall min/avg/max = 61.294/62.435/63.608 ms. Raster remains the dominant
cost; physics and orientation are coverage-independent. The final coefficient
factorization moves palette multiplication out of the four-sample loop.

### Per-pixel figures

The raster overwrites all 10,368 quadrant pixels with four SSAA samples each.
`bz_raster` costs **2,585.7 cycles/pixel**, or **646.4
cycles/sub-sample**. There is no `filter_blend` subtree.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1151/frame  min/avg/max 0.57/1.91/13.33 us  cpu 3.52%
isr_pack         144/frame  min/avg/max 6.20/6.89/10.03 us  cpu 1.58%
isr_dma_submit   144/frame  min/avg/max 0.70/0.94/5.87 us   cpu 0.21%
```

- Submit remains small beside color/interleave packing; wire transfer is
  asynchronous and absent from these CPU timings.
- The three ISR groups consume 5.31% of wall time in the peak window.
- The captured peak remains 11.795 ms below the 62.5 ms display window.

## Summary ranking

1. `bz_raster` — 44.68 ms, 90.1% of render: four shared-stencil Wendland
   samples and concentration-weighted palette compositing per pixel.
2. `bz_physics` — 4.52 ms, 9.1%: two Q16 Lotka-Volterra integrations over
   7,680 nodes.
3. `bz_orient` — 0.38 ms: one lattice rotation per frame.

The matched pre-optimization capture peaked at 59.083 ms. Fused
premultiplication, palette conversion hoisting, and coefficient factorization
reduce the peak by **8.378 ms (14.2%)**. The factorization changes only final
floating-point association: a 256-frame production simulation found no
coverage changes, a one-LSB maximum RGB16 error, and only 64 of 31,850,496
displayed sRGB channels changing by one code (131.46 dB PSNR).

## Caveats

- All scopes include ISR execution because CYCCNT free-runs.
- `filter_blend` is not used by this shader.
- Selective O3 reaches the typed grid driver, BZ pixel/kernel path, and physics.
- `HS_PROFILE_EPOCH_REVS=1200` changes capture lifetime only.
- The capture came from clean detached worktree tip `77e86d77`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=BZReactionDiffusion`,
`HS_PROFILE_WINDOW=32`; the reproduce command builds, locks, flashes, captures,
validates, and archives the exact artifacts.
