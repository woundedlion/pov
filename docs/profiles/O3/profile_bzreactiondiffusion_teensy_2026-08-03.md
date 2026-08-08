# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-03, **-O3**)

Point-in-time global-O3 reference paired with the
[`shipping report`](../shipping/profile_bzreactiondiffusion_teensy_2026-08-03.md).
Raw capture: `build/prof/bzreactiondiffusion_o3.log`. This replaces the
2026-08-02 reference with the final landed raster implementation.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile_o3` env: GCC 15.2.1, global `-O3 -ffast-math -fno-finite-math-only`, newlib-nano, DMA LEDs |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | BZReactionDiffusion 288×144, single-entry playlist, tip `77e86d77` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, epoch stretched to 1200 revs |
| Reproduce | `bash tools/profile_one.sh BZReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:71020, data:630280, headers:8324` / `RAM1:
variables:314496, code:46568, padding:18968, free:144256` / `RAM2:
variables:520064, free:4224`.

Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1
20251203`. Profile ELF SHA-256:
`b550d158981235710c09adc76e002b302d72c7c469b1a2dc244ed58656ebad9b`.

Exactness cross-check: window frames 1601–1632 root counter cycles divided by
600 MHz match the measured wall sum within **3.1 ppm**.

## Frame cadence

**Pass aggregate**: `bz_render` averages 49.165 ms/frame across 64 windows.
The worst window mean is 49.83 ms (frames 1601–1632), peak frame render is
**50.904 ms**, and **0/2048 frames (0%)** spill.

The 62.5 ms display window is held throughout the capture. The peak has
11.596 ms of display margin and 4.096 ms of margin against the 55 ms target.

## Phase-by-phase readout

The capture covers the same slowly varying kernel-occupancy regimes as the
shipping pass.

### Dense raster (window frames 993–1024, peak frame in the capture)

```text
frame                62.45 ms  37.47 Mcyc  100%
  bz_render          49.81 ms  29.89 Mcyc   80%
    bz_raster        44.98 ms  26.99 Mcyc   90%
    bz_physics        4.46 ms   2.67 Mcyc    9%
    bz_orient       374.7 us  224.8 kcyc     1%
  rd_timeline_step   24.5 us   14.7 kcyc     0%
  canvas_clear       89.0 us   53.4 kcyc     0%
  canvas_buffer_wait 12.53 ms   7.52 Mcyc   20%
```

Wall min/avg/max = 60.824/62.450/64.050 ms. Physics is slightly faster than
shipping, while raster is slightly slower; the total difference remains below
one quarter millisecond at the capture peak.

### Per-pixel figures

`bz_raster` costs **2,603.0 cycles/pixel**, or **650.7
cycles/sub-sample**, across the 10,368-pixel quadrant.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.55/1.88/12.36 us  cpu 3.46%
isr_pack         144/frame  min/avg/max 6.52/7.28/10.30 us  cpu 1.67%
isr_dma_submit   144/frame  min/avg/max 0.75/0.91/1.65 us   cpu 0.21%
```

- Packing remains the dominant CPU-side DMA preparation cost.
- Wire transfer is asynchronous and excluded from the cycle scopes.
- The three ISR groups consume 5.34% of wall time in the peak window.

## Summary ranking

1. `bz_raster` — 44.98 ms at dense coverage.
2. `bz_physics` — 4.46 ms, coverage-independent.
3. `bz_orient` — 0.37 ms.

## Caveats

- Scope time includes live ISR execution.
- This is a single-effect compiler reference, not a shippable full-roster image.
- `filter_blend` is not used by this shader.
- Epoch stretch changes capture lifetime only.
- The capture used clean detached worktree tip `77e86d77`.

## Harness

`targets/Profile/Profile.ino`, target `BZReactionDiffusion`, window 32; use the
reproduce command above for the locked build, flash, capture, and validation.

## Global O3 vs selective O3

Global O3 peaks at 50.904 ms versus 50.705 ms for shipping, **0.199 ms
(0.39%) slower**. FLASH code grows by **17,696 B** and ITCM code by
**16,256 B**. Selective O3 is both smaller and marginally faster for this
final raster.
