# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).

Raw capture: `C:\work\temp\gs-reclaim-20261005\combined-com4.log`, captured 2026-10-05 12:20 local time on COM4. Clean source `0f3f33b767f729fd4b379429add9de9f364c13b8`. Paired baseline: `C:\work\temp\gs-reclaim-20261005\baseline-com4.log`, source `e2aede3794e0171fb57ccf528d51875bc28b3e41`. The captures retain compiler, environment, source and ELF hashes in their provenance footers; all artifact hashes were checked.

## Setup

| Setting | Value |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel and DMA ISRs live, COM4 |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GSReactionDiffusion 288×144, shipping stable effect seed, default controls |
| Image | `profile`, shipping selective-O3; per-pixel profiling disabled |
| Capture | 130 s, 32-frame windows, epoch extended to 1200 revolutions (150 s) |
| Runtime frames | 2–2033 (2,032 frames) |
| Compiler | Pinned Teensy GCC 15.2.1 for profile and full-roster shipping builds |

Physics, pigment transport, palette cache construction, modified palette colors, the shared shader color/noise tail and the exact stencil fallback execute from cached flash with `HS_HOT_FLASH_MEMBER`. The fast stencil remains in ITCM. Physics unrolling, six chemistry substeps, pigment cadence, coverage tests and color formulas are unchanged. The shared shader tail replaces two inlined copies.

The shipping Phantasm image was built separately with profiling disabled and passed the size/layout gate:

| Metric | Baseline | Candidate | Delta (bytes) |
|---|--:|--:|--:|
| RAM1 code (ITCM) | 180,424 | 169,368 | -11,056 |
| RAM1 variables | 314,784 | 314,784 | +0 |
| FLASH code | 566,932 | 565,820 | -1,112 |
| FLASH data | 845,660 | 845,660 | +0 |
| Stack/local headroom | 12,896 | 12,896 | +0 |

ITCM falls by **11,056 bytes**, exceeding the earlier GS increase of 8,816 bytes. Six ITCM banks remain allocated; stack/local headroom stays unchanged. The profile image separately uses 25,432 bytes of RAM1 code and 315,328 bytes of RAM1 variables.

## Frame cadence

README cells: peak 🟢 39.371, spilled 🟢 0/2032 (0.00%).

**Peak render: 39.371 ms at frame 395; spilled: 0/2032 (0.00%).** Mean render is 32.854 ms and mean wall is 62.418 ms. The paired baseline peaks at 39.240 ms with a 32.748 ms mean. The candidate retains **13.629 ms** of measured margin below the user-approved 53 ms render ceiling; no captured frame reached that ceiling, including frame 1 (21.814 ms).

A display interval is 62.5 ms at 480 RPM. Each frame visits 10,368 quadrant positions, with four geometric coverage samples. Buffer wait is synchronization idle and is excluded from render telemetry.

`parse_profile.py validate` accepts both captures. The highest-render window's frame cycles / 600 MHz agrees with wall sum within 1.4 ppm. Scope summaries use complete post-startup windows 33–2016 (1,984 frames).

| Scope | Mean ms/frame | Worst window ms/frame |
|---|--:|--:|
| `grd_render` | 32.764 | 38.607 |
| `grd_simulate` | 9.980 | 10.235 |
| `grd_physics` | 5.111 | 5.126 |
| `grd_pigment` | 4.123 | 4.368 |
| `grd_rasterize` | 19.903 | 25.160 |
| `grd_shader_draw` | 17.519 | 22.962 |
| `grd_cull_flags` | 0.524 | 0.939 |
| `grd_orient` | 1.860 | 1.870 |
| `grd_color_noise` | 2.702 | 2.720 |
| `canvas_buffer_wait` | 29.418 | 40.065 |

## Phase-by-phase readout

The capture covers growing and denser reactions plus staged runtime reseeding. The profile has no lifecycle markers, so scope activity identifies reseed windows without assigning exact transition boundaries from timing alone.

### Highest mean render window (frames 385–416)

```text
frame                         62.402 ms 100.0%
  pov_preserve_half            0.138 ms   0.2%
  grd_render                  38.607 ms  61.9%
    grd_rasterize             25.160 ms  40.3%
      grd_shader_draw         22.962 ms  36.8%
      grd_cull_flags           0.335 ms   0.5%
      grd_orient               1.862 ms   3.0%
    grd_simulate              10.194 ms  16.3%
      grd_physics              5.120 ms   8.2%
      grd_pigment              4.322 ms   6.9%
    grd_color_noise            2.715 ms   4.4%
  rd_timeline_step             0.025 ms   0.0%
  canvas_clear                 0.084 ms   0.1%
  canvas_buffer_wait          23.547 ms  37.7%
```

Runtime reseed windows:

| Frames | Mean render ms | Peak frame render ms |
|---|--:|--:|
| 449–480 | 23.476 | 30.056 |
| 865–896 | 21.953 | 27.233 |
| 1281–1312 | 25.564 | 33.834 |
| 1793–1824 | 22.939 | 30.659 |

`grd_seed_reaction` and `grd_color_palette` carry MIXED-PARENT tags; their printed totals are attribution diagnostics and must not be added as exclusive siblings. Geometry and nearest-pigment fidelity remain covered by the native effect oracles.

## Column-ISR / DMA marshaling cost

Complete post-startup windows; packing and DMA submission nest inside wake and are not exclusive siblings.

| Scope | Mean us/call | Max us/call | CPU share |
|---|--:|--:|--:|
| `isr_wake` | 1.678 | 19.175 | 3.09% |
| `isr_pack` | 7.115 | 10.648 | 1.64% |
| `isr_dma_submit` | 0.937 | 11.843 | 0.22% |

Foreground render counters include ISR interruptions. DMA submission times CPU setup; SPI transfer proceeds asynchronously.

## Summary ranking

1. `grd_shader_draw`: 22.962 ms/frame in the highest-render window.
2. `grd_physics`: 5.120 ms/frame in the highest-render window.
3. `grd_pigment`: 4.322 ms/frame in the highest-render window.
4. `grd_color_noise`: 2.715 ms/frame in the highest-render window.
5. `grd_orient`: 1.862 ms/frame in the highest-render window.

The recovered space comes primarily from cached-flash execution of complete O3 kernels and sharing the shader tail. This capture measures that placement under the segmented driver, including foreground ISR interruptions.

## Caveats

This establishes the measured 53 ms ceiling for the default controls, shipping stable seed and captured lifecycle on COM4. It does not bound every control setting or seed; larger hue/shimmer settings use the procedural color fallback. Frame 1 is excluded from runtime means, but its render time is checked separately against 53 ms. Initialization counters preceding the first published frame are outside runtime frame telemetry. The first counter window is excluded from scope means.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1200`.

```bash
HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM4 \
  bash tools/profile_one.sh GSReactionDiffusion profile 130 32 \
  '-D HS_PROFILE_EPOCH_REVS=1200'
```

The supported wrapper holds the shared-device lock throughout build, flash and capture. Raw captures, build logs, environment dumps, provenance and both ELF/map pairs remain in `C:/work/temp/gs-reclaim-20261005`; `validated-comparison.json` records the computed statistics and size deltas.
