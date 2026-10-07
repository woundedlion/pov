# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-07, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`, using the full-cycle command below).
Raw capture: `C:/work/temp/hs-finding1-20261007/islamicstars_optimized4_com4.log`, captured 2026-10-07 14:12 on COM4. Baseline: `C:/work/temp/hs-finding1-20261007/islamicstars_opt_baseline_com4.log`. Replaces the October 6 shipping snapshot.

The paired captures isolate the face-distance fix and its exact-search optimizations. Baseline source is `fb0e4975853fda6d1ef8d73445a4f180c482037e`; candidate source is `6e0b0c59f26cc3328658286c38f7c106de739109` plus uncommitted changes to `core/render/sdf/face.h` and `core/render/sdf/face_geometry.h`.
The candidate source diff, ELF files, build logs and ABI attestations are retained under `C:/work/temp/hs-finding1-20261007/artifacts/islamicstars_optimized4_com4_6e0b0c59f26c_b63c5216534c`. These immutable snapshots predate ISR interval-accounting correction `fdc7930c5`; later master changes are outside the comparison.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, live flywheel + DMA ISRs, COM4 |
| Image | Shipping `profile` env: `-Os` plus selective `HS_O3`; face distances, mesh raster and effect draw/transform paths use O3 |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, candidate source snapshot described above |
| Method | 210 seconds, 16-frame windows, matched board and flags, complete preset cycle |
| Reproduce | `bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4"` |

Profile image: `Memory Usage on Teensy 4.0:` / `FLASH: code:129756, data:202212, headers:9020   free for files:1690628` / `RAM1: variables:315424, code:42840, padding:22696   free for local variables:143328` / `RAM2: variables:520064  free for malloc/new:4224`.

Shipping `phantasm` image; both pass the region-size and layout gates:

| Component | Baseline B | Candidate B | Delta B |
|---|--:|--:|--:|
| RAM1 code | 177656 | 181928 | +4272 |
| RAM1 variables | 314784 | 314784 | +0 |
| FLASH data | 846452 | 846452 | +0 |
| RAM2 variables | 520064 | 520064 | +0 |
| RAM1 stack headroom | 12896 | 12896 | +0 |

Exactness cross-check: frames 2577–2592, root cycles / 600 MHz versus measured wall sum differ by **0.3 ppm**. Both full-cycle captures pass `parse_profile.py ... validate`.

## Frame cadence

README cells: peak 🟢 55.842 (23), spilled 🟢 0/3327 (0.00%).

Peak render: **49.843 → 55.842 ms**, delta **+5.999 ms (+12.04%)**. Baseline spilled 0/3327; candidate spilled 0/3327.
Candidate peak is frame 2584, `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`. Every captured non-setup frame fits the 62.5 ms display window at 480 RPM, giving 16 fps. One quadrant contains 10,368 pixels. `canvas_buffer_wait` is the idle round-up to the next display flip; frame 1 is excluded from peak and spill figures. Observed peak headroom is 6.658 ms.

The September 28 global-O3 reference and its image deltas remain historical measurements of their own source pair; they are not a matched optimization baseline for this candidate.

## Phase-by-phase readout

The capture covers 23 ownership buckets and wraps to its initial preset. Held geometry alternates with build/advance transitions; all scopes include ISR time.

### Costliest clean hold (frames 2577–2592)

```
frame                       62.19 ms  37.31 Mcyc 100%
  pov_preserve_half         142.9 us   85.7 kcyc   0% x1.0 142.9 us/call
  is_timeline_step          49.43 ms  29.66 Mcyc  79%
    is_draw_shape           49.37 ms  29.62 Mcyc  79%
      is_mesh_scan          45.36 ms  27.21 Mcyc  73%
        scan_mesh_raster    33.77 ms  20.26 Mcyc  54%
          filter_blend       1.30 ms  778.6 kcyc   2% x18490.1 42 cyc/blend
        scan_face_setup     11.22 ms   6.73 Mcyc  18% x1082.0 10.4 us/call
      is_face_offsets       486.2 us  291.7 kcyc   1% x1.0 486.2 us/call
      is_mesh_transform      3.53 ms   2.12 Mcyc   6% x1.0 3527.5 us/call
  is_ripple_prepare           9.0 us    5.4 kcyc   0% x1.0 9.0 us/call
  canvas_clear               84.9 us   50.9 kcyc   0% x1.0 84.9 us/call
  canvas_buffer_wait        12.52 ms   7.51 Mcyc  20% x1.0 12521.3 us/call
```

Wall min/avg/max: 56.41/62.19/69.80 ms. Timeline work averages 49.43 ms; peak render in this window is 55.842 ms. Scope attribution is inclusive; shared raster counters can first-parent under a different draw path. The wait absorbs the remainder of each display window.

### Build/advance transition (frames 1025–1040)

```
frame                       65.54 ms  39.32 Mcyc 100%
  pov_preserve_half         146.9 us   88.2 kcyc   0% x1.0 146.9 us/call
  is_timeline_step          27.78 ms  16.67 Mcyc  42%
    is_build_draw           22.78 ms  13.67 Mcyc  35%
      is_build_scan         22.74 ms  13.65 Mcyc  35% x0.8 30326.6 us/call
      is_mesh_transform       2.2 us    1.3 kcyc   0% x0.2 8.8 us/call
    hk_conway_compile       223.9 us  134.4 kcyc   0% x0.8 298.5 us/call
    hk_conway_sweep         577.0 us  346.2 kcyc   1% x0.8 769.3 us/call
    is_draw_shape            2.68 ms   1.61 Mcyc   4%
      is_mesh_scan           2.67 ms   1.60 Mcyc   4%
        scan_mesh_raster    22.70 ms  13.62 Mcyc  35%
          filter_blend       1.07 ms  642.5 kcyc   2% x15999.8 40 cyc/blend
        scan_face_setup      2.59 ms   1.56 Mcyc   4% x287.0 9.0 us/call
      is_face_offsets         6.3 us    3.8 kcyc   0% x0.2 25.2 us/call
  is_ripple_prepare           0.1 us    0.1 kcyc   0% x1.0 0.1 us/call
  canvas_clear               87.1 us   52.3 kcyc   0% x1.0 87.1 us/call
  canvas_buffer_wait        37.52 ms  22.51 Mcyc  57% x1.0 37519.7 us/call
```

Wall min/avg/max: 61.17/65.54/74.78 ms. Timeline work averages 27.78 ms; peak render in this window is 52.125 ms. Scope attribution is inclusive; shared raster counters can first-parent under a different draw path. The wait absorbs the remainder of each display window.

### Per-preset table

Each bucket includes its held frames and following transition. Clean holds exclude advance-straddling windows and use the modal `scan_mesh_raster` call count within each preset; rows rank by the clean timeline cost. V/E/F/I describes the IslamicStars seed recorded at spawn, not every intermediate build mesh. Windows aggregate repeated visits to a preset.

| Shape | Geometry | Blends/f | Clean timeline ms/f | Clean render ms/f | Peak ms | Spilled/frames | Clean windows | fps |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | V=20 E=30 F=12 I=60 | 18490 | 49.43 | 49.67 | 55.842 | 0/176 | 6/10 | 16 |
| truncatedOctahedron_gyro_kis_hk17 | V=24 E=36 F=14 I=72 | 18253 | 42.64 | 42.87 | 48.555 | 0/184 | 4/10 | 16 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | V=60 E=90 F=32 I=180 | 16002 | 42.27 | 42.51 | 46.409 | 0/144 | 4/6 | 16 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin73 | V=60 E=90 F=32 I=180 | 17517 | 41.37 | 41.61 | 45.118 | 0/144 | 4/8 | 16 |
| truncatedIcosidodecahedron_bevel5_relax_hk77 | V=120 E=180 F=62 I=360 | 18671 | 39.61 | 39.85 | 43.958 | 0/144 | 4/6 | 16 |
| truncatedIcosahedron_hk54_ambo_hk72 | V=60 E=90 F=32 I=180 | 17772 | 36.33 | 36.57 | 40.743 | 0/140 | 4/6 | 16 |
| truncatedIcosahedron_ambo_relax_truncate33_hk64 | V=60 E=90 F=32 I=180 | 16827 | 33.71 | 33.95 | 36.996 | 0/144 | 4/6 | 16 |
| icosahedron_snub_relax_truncate033_hankin62 | V=12 E=30 F=20 I=60 | 16489 | 33.31 | 33.55 | 35.267 | 0/72 | 2/3 | 16 |
| dodecahedron_ambo_bevel33_relax_hk66 | V=20 E=30 F=12 I=60 | 16027 | 30.37 | 30.62 | 33.026 | 0/156 | 6/8 | 16 |
| truncatedIcosahedron_hk58_chamfer63 | V=60 E=90 F=32 I=180 | 16489 | 28.43 | 28.66 | 30.520 | 0/124 | 4/6 | 16 |
| rhombicuboctahedron_hk63_ambo_hk63 | V=24 E=48 F=26 I=96 | 15133 | 28.33 | 28.56 | 30.487 | 0/140 | 4/8 | 16 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | V=120 E=180 F=62 I=360 | 16035 | 27.82 | 28.05 | 52.125 | 0/156 | 2/8 | 16 |
| dodecahedron_hk72_ambo_dual_hk20 | V=20 E=30 F=12 I=60 | 14711 | 27.12 | 27.35 | 32.395 | 0/103 | 3/4 | 16 |
| dodecahedron_hk54_ambo_hk72 | V=20 E=30 F=12 I=60 | 14646 | 25.34 | 25.58 | 26.281 | 0/140 | 6/8 | 16 |
| dodecahedron_hk62_ambo_hk62 | V=20 E=30 F=12 I=60 | 14126 | 24.73 | 24.97 | 26.423 | 0/138 | 4/7 | 16 |
| icosahedron_ambo_truncate033_hankin59 | V=12 E=30 F=20 I=60 | 13943 | 23.79 | 24.03 | 26.479 | 0/136 | 4/6 | 16 |
| octahedron_hk17_ambo_hk73 | V=6 E=12 F=8 I=24 | 13071 | 21.82 | 22.06 | 23.638 | 0/140 | 4/6 | 16 |
| octahedron_hk34_ambo_hk72 | V=6 E=12 F=8 I=24 | 13040 | 21.08 | 21.32 | 22.797 | 0/140 | 4/6 | 16 |
| snubDodecahedron_truncate5d_ambo_dual | V=60 E=150 F=92 I=300 | 16833 | 19.73 | 19.97 | 40.729 | 0/156 | 4/8 | 16 |
| dodecahedron_bevel2_relax_gyro | V=20 E=30 F=12 I=60 | 15903 | 18.86 | 19.10 | 40.408 | 0/176 | 4/10 | 16 |
| truncatedIcosahedron_truncate50d_ambo_dual | V=60 E=90 F=32 I=180 | 16021 | 18.47 | 18.71 | 35.810 | 0/78 | 2/4 | 16 |
| icosidodecahedron_truncate5d_ambo_dual | V=30 E=60 F=32 I=120 | 14571 | 14.72 | 14.96 | 24.646 | 0/156 | 4/8 | 16 |
| icosahedron_kis_gyro | V=12 E=30 F=20 I=60 | 11530 | 9.49 | 9.73 | 31.988 | 0/240 | 2/12 | 16 |

### Per-pixel figures

Peak-window blends: 18,490/frame (1.78× quadrant coverage), 42.1 cycles/blend. Raster scope: 1095.8 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```
isr_wake         1147/frame  0.59/1.73/11.85 us  cpu 3.18%
isr_pack         143/frame  6.54/7.25/9.98 us  cpu 1.67%
isr_dma_submit   143/frame  0.66/0.94/1.03 us  cpu 0.21%
```

Pack and submit run inside `isr_wake`; the wake row alone is the inclusive measured flywheel CPU share. DMA completion and other interrupts are unmeasured. These captures predate ISR snapshot correction `fdc7930c5`, so the printed CPU-share percentages carry interval-accounting bias and cannot establish an exact foreground render budget.

Average packing costs 7.8 times submission. Each segment contains 72 pixels: `HD107SFrame<72>` uses 300 bytes per image, including its 8-byte end frame. The image-plus-black composite is 600 bytes. At 24 MHz SPI, 92 functional clocks per byte at 240 MHz gives a maximum **230 us** asynchronous wire transfer, excluded from CPU submission time.

The observed 55.842 ms render peak leaves 6.658 ms within the 62.5 ms display window. No captured frame needs a speedup to meet cadence; this is observed margin, not an exact CPU budget.

## Summary ranking

1. `is_timeline_step`: 79% of frame, 49.43 ms/frame.
2. `pov_preserve_half`: 0% of frame, 0.14 ms/frame.
3. `canvas_clear`: 0% of frame, 0.08 ms/frame.
4. `is_ripple_prepare`: 0% of frame, 0.01 ms/frame.

The initial certified-sector candidate peaked at 61.729 ms. Removing angular division/modulo work, strengthening omitted-edge bounds, caching boundary rays, and using an exact row-mask search for angular-backtracking faces reduced the best full-cycle peak to 55.842 ms. The row search is limited to rejected large nonconvex faces; applying it to all stars and a separate bounding-box hierarchy were slower in paired host screens and were not adopted.

The final native direct-raster benchmark covers all 23 completed recipes: 88.235 → 86.610 ms summed at 288×144 (−1.84%) and 24.614 → 24.032 ms at 96×48 (−2.36%). Native tests overlapped these host runs, and this benchmark excludes morph intermediates, so it cannot substitute for the device peak.

Actual IslamicStars shader renders show small black corner notches removed by the corrected sign and subtle palette-band shifts from the corrected distance. The optimized and first corrected versions are bit-identical across 224 sampled recipe-2/7 frames at both resolutions. This is rendered-image evidence, not a visual test on the spinning LED display.


Native timings were used to screen optimizations; their workload and host scheduling differ from a physical quadrant frame. Matched on-device peak render is the acceptance evidence.

## Caveats

- All scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parents under the first entrant; its subtree can be hidden under an inactive parent, while its calls approximate blended pixels.
- Per-pixel scope overhead is avoided; this is standard, non-deep profiling.
- Selective O3 is the shipping configuration; geometry, scan, transform and draw regions remain active.
- Dwell compression and ordered traversal alter scheduling, not per-frame computation; stretched epochs avoid effect teardown.
- Candidate headers were uncommitted at capture; their exact diff, source status and build hashes are retained with the binaries.
- Peaks are observed maxima for these deterministic captures, not global upper bounds.

## Harness

`targets/Profile/Profile.ino` supports target, window, epoch and scheduling knobs. `just profile IslamicStars` routes through the supported locked wrapper; use the Reproduce command for this complete-cycle capture.
