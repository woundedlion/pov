# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-07, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`, using the full-cycle command below).
Raw capture: `C:/work/temp/hs-islamic50-20261007/final_full.log`, captured 2026-10-07 15:16 on COM4. Baseline: `C:/work/temp/hs-islamic50-20261007/baseline_full.log`. Replaces the earlier October 7 face-distance snapshot.

The optimization starts from the corrected face-distance implementation at `f90aea8e56577ab7430bfba82af1567bf73c984d`. Final full-cycle source is `0bcce1a109036a35979e13af9c5d3934224a2a24`; the final rebase also incorporates the independent cold parameter-attribute change `be9091aa6`. Source, binaries, flags, and hashes are retained in `C:/work/temp/hs-islamic50-20261007/artifacts/final_full_0bcce1a10903_db1dd75e7bc5`. The compiler is `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, live flywheel + DMA ISRs, COM4 |
| Image | Shipping `profile` environment: `-Os`, newlib-nano, selective `HS_O3` face-distance, mesh-raster, transform, and draw regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist |
| Method | 210 seconds, 16-frame windows, full 23-preset cycle and wrap; transition speed 4; epoch 1920 revolutions |
| Reproduce | `HS_PROFILE_TREE=/c/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4"` |

Profile image: `FLASH: code:128988, data:202060, headers:8916   free for files:1691652 / RAM1: variables:315424, code:40744, padding:24792   free for local variables:143328 / RAM2: variables:520064  free for malloc/new:4224`.

Shipping `phantasm` image; baseline and final builds pass all region-size and layout gates:

| Component | Baseline B | Optimization before peer rebase B | Final rebased B | Final delta B |
|---|--:|--:|--:|--:|
| RAM1 code (ITCM) | 182488 | 181496 | 177864 | -4624 |
| RAM1 variables | 314784 | 314784 | 314784 | 0 |
| FLASH data | 846316 | 846316 | 846316 | 0 |
| RAM2 variables | 520064 | 520064 | 520064 | 0 |
| RAM1 stack headroom | 12896 | 12896 | 12896 | 0 |

The optimization itself reduces ITCM by **992 B**; the additional 3632 B reduction follows the independent peer change during rebase. No ITCM increase is required. Final shipping FLASH code/data/headers are 640060/846316/8660 B.

Exactness cross-check: frames 2577–2592, root cycles / 600 MHz versus measured wall sum differ by **3.2 ppm**. Full-cycle validation confirms one effect/resolution, monotonic frames, complete per-frame telemetry, and all 23 presets followed by wrap.

Profile ELF SHA-256: `db1dd75e7bc5fed4d2ac82b454b68619de9d283781594262c43925895cd73159`. Shipping ELF SHA-256: `4f377365b58e55a1328a5f5d848ca0b53bd48acb8a0e5a4e7ebe79b8eae75883`.

## Frame cadence

README cells: peak 🟢 48.510 (23), spilled 🟢 0/3327 (0.00%).

Peak individual-frame render falls **55.819 to 48.510 ms**, a **7.309 ms (13.09%)** reduction. Peak frame 2584 belongs to `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`. The final pass has 3328 frames in complete windows, including setup frame 1; peak/spill comparisons exclude setup, leaving 3327 measured frames. Baseline likewise has zero spills among 3327 non-setup frames.

`is_timeline_step` averages 21.711 ms/frame across the pass, with a worst window mean of 42.794 ms/frame (frames 2577–2592). These means are separate from the measured individual-frame peak.

A display half-revolution is 62.5 ms at 480 RPM, so all captured phases hold 16 fps. One quadrant contains 10,368 pixels. The final peak leaves 13.990 ms before the display deadline and 1.490 ms below the 50 ms target. `canvas_buffer_wait` is the idle round-up to the next display flip.

An independent fixed-preset-10 stress run at the default transition speed measured **48.575 ms**, frame 585, with 0/1407 non-setup spills (0/1408 including setup). This run used the same algorithm before the peer rebase; its immutable source/build provenance is in `final_worst_default.log`. It preserves the longer default build, hold, and ripple evolution. The largest observed standard-profile peak across the full cycle and this stress run is 48.575 ms, leaving 1.425 ms below the target.

## Phase-by-phase readout

Each preset builds through recipe operations, holds and ripples, then transitions to the next solid. The full cycle visits all 23 presets and returns to its initial shape. Counter scopes include interrupts.

### Costliest clean hold (frames 2577–2592)

```
frame                         62.22 ms  37.33 Mcyc  100%
  pov_preserve_half           143.6 us   86.2 kcyc    0%
  is_timeline_step            42.79 ms  25.68 Mcyc   69%
    is_draw_shape             42.74 ms  25.64 Mcyc   69%
      is_mesh_scan            38.72 ms  23.23 Mcyc   62%
        scan_mesh_raster      32.19 ms  19.31 Mcyc   52%
          filter_blend         1.31 ms  783.2 kcyc    2% x18490 42cyc/bl
        scan_face_setup        6.21 ms   3.72 Mcyc   10% x1082 5.7us/c
      is_face_offsets         484.6 us  290.8 kcyc    1%
      is_mesh_transform        3.53 ms   2.12 Mcyc    6%
  is_ripple_prepare             6.9 us    4.2 kcyc    0%
  canvas_clear                 84.8 us   50.9 kcyc    0%
  canvas_buffer_wait          19.19 ms  11.51 Mcyc   31%
```

Wall min/avg/max: 57.03/62.22/68.95 ms. Timeline work averages 42.79 ms/frame; the window's individual-frame render peak is 48.510 ms. Shared raster scopes can accumulate under their first parent across build and sprite paths. The display wait absorbs the remaining cadence interval.

### Build/advance transition (frames 1025–1040)

```
frame                         65.30 ms  39.18 Mcyc  100%
  pov_preserve_half           144.5 us   86.7 kcyc    0%
  is_timeline_step            25.84 ms  15.51 Mcyc   40%
    is_build_draw             20.91 ms  12.55 Mcyc   32%
      is_build_scan           20.88 ms  12.53 Mcyc   32%
    hk_conway_compile         219.2 us  131.5 kcyc    0%
    hk_conway_sweep           566.3 us  339.8 kcyc    1%
    is_draw_shape              2.63 ms   1.58 Mcyc    4%
      is_mesh_scan             2.62 ms   1.57 Mcyc    4%
        scan_mesh_raster      21.75 ms  13.05 Mcyc   33%
          filter_blend         1.08 ms  648.4 kcyc    2% x16000 41cyc/bl
        scan_face_setup        1.62 ms  972.6 kcyc    2% x287 5.6us/c
      is_face_offsets           6.3 us    3.8 kcyc    0%
      is_mesh_transform         2.2 us    1.3 kcyc    0%
  is_ripple_prepare             0.3 us    0.2 kcyc    0%
  canvas_clear                 85.1 us   51.0 kcyc    0%
  canvas_buffer_wait          39.23 ms  23.54 Mcyc   60%
```

Wall min/avg/max: 61.17/65.30/74.87 ms. Timeline work averages 25.84 ms/frame; the window's individual-frame render peak is 48.295 ms. Shared raster scopes can accumulate under their first parent across build and sprite paths. The display wait absorbs the remaining cadence interval.

### Per-preset table

Ownership buckets include held frames and the following transition. Clean holds exclude advance-straddling windows and select the modal `scan_mesh_raster` call count for each preset. Rows rank by clean timeline cost. V/E/F/I is seed geometry from the spawn marker; the heaviest completed preset, index 10, has V=3240, E=4320, F=1082, I=8640. Counts aggregate repeated visits.

| Shape | Seed geometry | Blends/f | Clean timeline ms/f | Clean render ms/f | Peak ms | Spilled/frames | Clean windows | fps |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | V=20 E=30 F=12 I=60 | 18490 | 42.79 | 43.03 | 48.510 | 0/176 | 6/12 | 16 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | V=60 E=90 F=32 I=180 | 16002 | 38.01 | 38.25 | 41.773 | 0/144 | 4/8 | 16 |
| truncatedOctahedron_gyro_kis_hk17 | V=24 E=36 F=14 I=72 | 18253 | 36.10 | 36.33 | 43.007 | 0/184 | 4/12 | 16 |
| truncatedIcosidodecahedron_bevel5_relax_hk77 | V=120 E=180 F=62 I=360 | 18671 | 35.41 | 35.64 | 39.640 | 0/144 | 4/8 | 16 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin73 | V=60 E=90 F=32 I=180 | 17517 | 35.00 | 35.24 | 38.859 | 0/144 | 4/10 | 16 |
| truncatedIcosahedron_hk54_ambo_hk72 | V=60 E=90 F=32 I=180 | 17772 | 31.63 | 31.87 | 35.846 | 0/140 | 4/8 | 16 |
| truncatedIcosahedron_ambo_relax_truncate33_hk64 | V=60 E=90 F=32 I=180 | 16689 | 29.74 | 29.98 | 33.050 | 0/144 | 4/8 | 16 |
| icosahedron_snub_relax_truncate033_hankin62 | V=12 E=30 F=20 I=60 | 16489 | 28.68 | 28.92 | 30.471 | 0/72 | 2/4 | 16 |
| dodecahedron_ambo_bevel33_relax_hk66 | V=20 E=30 F=12 I=60 | 16027 | 26.80 | 27.04 | 29.749 | 0/156 | 6/10 | 16 |
| truncatedIcosahedron_hk58_chamfer63 | V=60 E=90 F=32 I=180 | 16489 | 26.31 | 26.55 | 28.346 | 0/124 | 4/8 | 16 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | V=120 E=180 F=62 I=360 | 16000 | 25.84 | 26.07 | 48.295 | 0/156 | 2/10 | 16 |
| rhombicuboctahedron_hk63_ambo_hk63 | V=24 E=48 F=26 I=96 | 15133 | 24.97 | 25.21 | 27.093 | 0/140 | 4/10 | 16 |
| dodecahedron_hk72_ambo_dual_hk20 | V=20 E=30 F=12 I=60 | 14711 | 24.10 | 24.34 | 28.892 | 0/103 | 3/5 | 16 |
| dodecahedron_hk54_ambo_hk72 | V=20 E=30 F=12 I=60 | 14646 | 22.53 | 22.77 | 23.414 | 0/140 | 6/10 | 16 |
| dodecahedron_hk62_ambo_hk62 | V=20 E=30 F=12 I=60 | 14126 | 22.24 | 22.48 | 23.924 | 0/138 | 4/10 | 16 |
| icosahedron_ambo_truncate033_hankin59 | V=12 E=30 F=20 I=60 | 13943 | 21.58 | 21.81 | 24.244 | 0/136 | 4/8 | 16 |
| octahedron_hk17_ambo_hk73 | V=6 E=12 F=8 I=24 | 13071 | 19.60 | 19.84 | 21.660 | 0/140 | 4/8 | 16 |
| octahedron_hk34_ambo_hk72 | V=6 E=12 F=8 I=24 | 13040 | 18.83 | 19.07 | 20.316 | 0/140 | 4/8 | 16 |
| snubDodecahedron_truncate5d_ambo_dual | V=60 E=150 F=92 I=300 | 16833 | 18.49 | 18.73 | 38.230 | 0/156 | 4/10 | 16 |
| dodecahedron_bevel2_relax_gyro | V=20 E=30 F=12 I=60 | 15903 | 17.73 | 17.96 | 38.015 | 0/176 | 4/12 | 16 |
| truncatedIcosahedron_truncate50d_ambo_dual | V=60 E=90 F=32 I=180 | 16021 | 17.69 | 17.93 | 33.763 | 0/78 | 2/5 | 16 |
| icosidodecahedron_truncate5d_ambo_dual | V=30 E=60 F=32 I=120 | 14571 | 14.24 | 14.48 | 23.257 | 0/156 | 4/10 | 16 |
| icosahedron_kis_gyro | V=12 E=30 F=20 I=60 | 11530 | 9.37 | 9.61 | 31.114 | 0/240 | 2/14 | 16 |

### Per-pixel figures

Costliest clean-hold window: 18,490 blends/frame (1.78× quadrant coverage), 42.4 cycles/blend; mesh raster averages 1044.4 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```
isr_wake        1148/frame 0.65/1.72/11.73 us cpu 3.17%
isr_pack        144/frame 6.53/7.24/9.83 us cpu 1.66%
isr_dma_submit  144/frame 0.76/0.93/1.11 us cpu 0.21%
```

Pack and submit run inside `isr_wake`; the wake row alone is the inclusive measured flywheel CPU share. DMA completion and other interrupts remain unmeasured. Packing averages 7.8 times the CPU submission cost. ISR CPU shares use the separately captured ISR snapshot interval.

Each segment has 72 pixels. `HD107SFrame<72>` is 300 bytes including its 8-byte end frame, or 600 bytes for image plus black. At 24 MHz SPI and 92 functional clocks/byte at 240 MHz, asynchronous wire time is at most **230 us**, outside CPU submission time. The measured wake share leaves approximately 60.52 ms per display interval after that interrupt alone; this is not an exact foreground budget because other interrupts remain. No measured phase requires a speedup to meet cadence.

## Summary ranking

1. `is_timeline_step`: 68.8% of frame, 42.79 ms/frame.
2. `scan_mesh_raster`: 51.7% of frame, 32.19 ms/frame, included in the timeline scope.
3. `scan_face_setup`: 10.0% of frame, 6.21 ms/frame, included in the timeline scope.

Fixed-preset-10 screens used 30 seconds, 8-frame windows, and transition speed 4. These branch experiments are not independently additive speedups. Their peaks are not interchangeable with the longer full-cycle acceptance capture:

| Experiment | Peak render ms | Worst timeline window mean ms/f |
|---|--:|--:|
| Corrected baseline | 52.700 | 50.427 |
| Outline y-walk | 51.555 | 49.068 |
| Row-edge ownership plus certified outside-distance cap | 51.041 | 48.821 |
| Outline sector search (rejected) | 53.605 | 51.034 |
| Outlined sector plus cheaper arc bounds | 52.181 | 49.635 |
| Outlined sector plus central-edge-first search | 52.245 | 49.957 |
| Bounds plus central-sector search | 50.298 | 48.127 |
| Cache unit rays; one inradius square root | 49.181 | 47.225 |
| Raise sector threshold to 20 (rejected) | 51.136 | 48.514 |
| Early hemisphere rejection | 45.641 | 43.961 |

The intermediate full-cycle peaks were 53.117 ms after bounds/central-sector changes and 52.060 ms after normalized rays. Hemisphere rejection supplied the remaining improvement. A certified outside-distance cap avoids work only beyond the raster rejection threshold; central-edge-first sector search expands until its bound certifies the same distance. Unit boundary rays reduce per-probe projection and distance arithmetic, and inradius calculation takes one square root after minimizing squared distance. Arc bounds reject irrelevant extrema before normalization and defer vertex acos evaluation to extrema. Hemisphere rejection uses the face's conservative phi extent against the render clip including its margin, before projection and per-face setup.

ARM code generation explained a separate register-pressure cost: baseline shipping raster code was 6932 B with 180 B locals and duplicated inline y-walk bodies. Outlining only y-walk produced a 1444 B helper and 4240 B raster body, saving 848 B ITCM. The second candidate had a 4096 B raster body with 164 B locals and a 1128 B helper. Final shipping disassembly has a 3956 B primary raster body, 164 B local allocation, and a 1128 B out-of-line y-walk helper; `compute_full_bounds` is 1322 B. The smaller raster local allocation does not imply a smaller maximum call-chain stack. Edge lambdas were already inlined; forcing them did not change the body, and that trial's loader failure supplied no valid timing. Outlining the sector path and raising its threshold were measured and rejected.

Deep fixed-preset-10 captures compare the same 321-328 animation window. Both diagnostic captures ran on COM3 at 600 MHz. These are instrumented diagnostics, not acceptance timings.

| Deep scope | Baseline ms/f | Final ms/f | Baseline / final share of frame |
|---|--:|--:|--:|
| Timeline | 55.095 | 47.066 | 87.89% / 75.11% |
| Mesh raster | 38.684 | 35.639 | 61.71% / 56.87% |
| Face setup | 11.844 | 6.881 | 18.89% / 10.98% |
| Face bounds | 4.673 | 1.859 | 7.45% / 2.97% |
| Face projection | 2.137 | 1.255 | 3.41% / 2.00% |
| Face phi extent | 0.580 | 0.657 | 0.93% / 1.05% |

Across 57 deep windows, the parser's counter-read-adjusted sector event cost falls 556.8 to 477.7 cycles; sector share of the measured probe stages falls 48.4% to 44.5%. Exact-search event cost falls 307.2 to 302.1 cycles (probe-stage share 21.1% to 22.2%). Instrumented peak render changes 57.472 to 48.879 ms. The higher phi-extent cost is offset by rejecting opposite-hemisphere faces before their projection and setup. Deep logs retain their own source and binary provenance in `baseline_worst_deep.log` and `final_worst_deep.log`.

Image validation renders the actual IslamicStars shader with fixed RNG and time through every recipe, transition and ripple frame. Native Clang 23 Release uses `-O3 -DNDEBUG -ffast-math -fno-finite-math-only`, with no LTO or explicit FP-contraction override. Full-canvas captures at 96×48 and 288×144 cover 1776 frames each, 81,838,080 pixels total: zero changed RGB16 channels for the final distance/bounds algorithm before clip rejection. Four independent fixed-quadrant runs at each resolution exercise clip rejection, preserving the engine's one-pixel render margin: 14,208 frames and 327,352,320 full-frame pixel samples, again zero differences including margins. Baseline repeat determinism and baseline quadrant display parity against full-canvas rendering are exact. The evidence is in `visual/candidate-3/comparison.json`, `visual/candidate-4-quadrants/comparison.json`, and `visual/baseline-quadrants/full_canvas_display_parity.json` under the capture root. This validates native rendered output, not physical spinning-display perception or cross-compiler bit identity.

Validation: final native suite 105 passed, one existing skip, zero failures (88.28 seconds). The additional `HS_EFFECTS_FULL=1 HS_SMOKE_FRAMES=120` effect and rendering sweeps passed 2/2 tests in 345.49 seconds. Debug and Release direct checks each passed 4,234,990 checks. The optimized WASM build and the standard 120-frame WASM smoke suite passed. Clip-oracle mutations that remove the hemisphere sign guards or ignore render margins fail as expected; independent arc-bound and distance oracles cover near-polar arcs, winding, finite rejection caps, and zero-radius certificates.

## Caveats

- CYCCNT includes live ISR time in every scope.
- `filter_blend` and shared raster scopes retain their first-parent attribution; mixed-parent and duplicate-name counters are not exclusive phase costs.
- Deep per-probe profiling adds overhead and informs optimization; acceptance uses the standard non-deep image.
- Selective O3 remains the shipping configuration for geometry, distance, raster, transform, and draw regions.
- Transition speed 4 compresses full-cycle scheduling. The independent default-speed stress run samples the longer choreography; neither finite capture establishes a global maximum over all controls and orientations.
- Intermediate captures used uncommitted experiments with retained source diffs and binary hashes. Final full-cycle provenance names the rebased source above.
- The historical global-O3 reference belongs to its own source pair; its timing and size deltas are not matched comparisons to this optimization.

## Harness

`targets/Profile/Profile.ino` supports target, window, epoch and scheduling controls. `just profile IslamicStars` routes through the supported locked wrapper; use the Reproduce command for a full-cycle capture.
