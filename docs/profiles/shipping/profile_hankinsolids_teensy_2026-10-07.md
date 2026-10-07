# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-10-07, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HankinSolids`, using the full-cycle command below).
Raw capture: `C:/work/temp/hs-finding1-20261007/hankinsolids_optimized4_com3.log`, captured 2026-10-07 14:18 on COM3. Baseline: `C:/work/temp/hs-finding1-20261007/hankinsolids_opt_baseline_com3.log`. Replaces the October 6 shipping snapshot.

The paired captures isolate the face-distance fix and its exact-search optimizations. Baseline source is `fb0e4975853fda6d1ef8d73445a4f180c482037e`; candidate source is `6e0b0c59f26cc3328658286c38f7c106de739109` plus uncommitted changes to `core/render/sdf/face.h` and `core/render/sdf/face_geometry.h`.
The candidate source diff, ELF files, build logs and ABI attestations are retained under `C:/work/temp/hs-finding1-20261007/artifacts/hankinsolids_optimized4_com3_6e0b0c59f26c_4cb9f7c05c92`. These immutable snapshots predate ISR interval-accounting correction `fdc7930c5`; later master changes are outside the comparison.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, live flywheel + DMA ISRs, COM3 |
| Image | Shipping `profile` env: `-Os` plus selective `HS_O3`; face distances, mesh raster and effect draw/transform paths use O3 |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, candidate source snapshot described above |
| Method | 280 seconds, 16-frame windows, matched board and flags, complete preset cycle |
| Reproduce | `bash tools/profile_one.sh HankinSolids profile 280 16 "-D HS_PROFILE_EPOCH_REVS=2400 -D HS_PROFILE_ORDERED_CYCLE"` |

Profile image: `Memory Usage on Teensy 4.0:` / `FLASH: code:112172, data:172208, headers:8480   free for files:1738756` / `RAM1: variables:315392, code:40376, padding:25160   free for local variables:143360` / `RAM2: variables:520064  free for malloc/new:4224`.

Shipping `phantasm` image; both pass the region-size and layout gates:

| Component | Baseline B | Candidate B | Delta B |
|---|--:|--:|--:|
| RAM1 code | 177656 | 181928 | +4272 |
| RAM1 variables | 314784 | 314784 | +0 |
| FLASH data | 846452 | 846452 | +0 |
| RAM2 variables | 520064 | 520064 | +0 |
| RAM1 stack headroom | 12896 | 12896 | +0 |

Exactness cross-check: frames 2241–2256, root cycles / 600 MHz versus measured wall sum differ by **1.9 ppm**. Both full-cycle captures pass `parse_profile.py ... validate`.

## Frame cadence

README cells: peak 🟢 34.681 (19), spilled 🟢 0/4447 (0.00%).

Peak render: **34.619 → 34.681 ms**, delta **+0.062 ms (+0.18%)**. Baseline spilled 0/4447; candidate spilled 0/4447.
Candidate peak is frame 2252, `truncatedIcosidodecahedron`. Every captured non-setup frame fits the 62.5 ms display window at 480 RPM, giving 16 fps. One quadrant contains 10,368 pixels. `canvas_buffer_wait` is the idle round-up to the next display flip; frame 1 is excluded from peak and spill figures. Observed peak headroom is 27.819 ms.

The September 28 global-O3 reference and its image deltas remain historical measurements of their own source pair; they are not a matched optimization baseline for this candidate.

## Phase-by-phase readout

The capture covers 19 ownership buckets and wraps to its initial preset. Held geometry alternates with build/advance transitions; all scopes include ISR time.

### Costliest clean hold (frames 2241–2256)

```
frame                       62.80 ms  37.68 Mcyc 100%
  pov_preserve_half         148.9 us   89.4 kcyc   0% x1.0 148.9 us/call
  hk_timeline_step          26.11 ms  15.66 Mcyc  42%
    hk_draw_mesh            25.40 ms  15.24 Mcyc  40%
      hk_mesh_scan          25.37 ms  15.22 Mcyc  40%
        scan_mesh_raster    23.27 ms  13.96 Mcyc  37%
          filter_blend       1.08 ms  648.4 kcyc   2% x15913.6 41 cyc/blend
        scan_face_setup      1.98 ms   1.19 Mcyc   3% x182.0 10.9 us/call
      hk_mesh_transform      31.7 us   19.0 kcyc   0% x1.0 31.7 us/call
    hk_update_hankin        657.8 us  394.7 kcyc   1% x1.0 657.8 us/call
  canvas_clear               84.1 us   50.5 kcyc   0% x1.0 84.1 us/call
  canvas_buffer_wait        36.46 ms  21.88 Mcyc  58% x1.0 36462.6 us/call
```

Wall min/avg/max: 56.06/62.80/68.78 ms. Timeline work averages 26.11 ms; peak render in this window is 34.681 ms. Scope attribution is inclusive; shared raster counters can first-parent under a different draw path. The wait absorbs the remainder of each display window.

### Build/advance transition (frames 1153–1168)

```
frame                       63.32 ms  37.99 Mcyc 100%
  pov_preserve_half         149.9 us   90.0 kcyc   0% x1.0 149.9 us/call
  hk_timeline_step          10.86 ms   6.52 Mcyc  17%
    hk_draw_mesh             2.25 ms   1.35 Mcyc   4%
      hk_mesh_scan           2.25 ms   1.35 Mcyc   4% x0.1 17973.5 us/call
      hk_mesh_transform       1.3 us    0.8 kcyc   0% x0.1 10.5 us/call
    hk_conway_compile        23.9 us   14.3 kcyc   0% x0.9 27.3 us/call
    hk_conway_sweep          60.5 us   36.3 kcyc   0%
        scan_mesh_raster    10.18 ms   6.11 Mcyc  16%
          filter_blend      754.6 us  452.7 kcyc   1% x12169.2 37 cyc/blend
        scan_face_setup     281.6 us  168.9 kcyc   0% x41.0 6.9 us/call
    hk_update_hankin         19.3 us   11.6 kcyc   0% x0.1 154.5 us/call
  canvas_clear               84.2 us   50.5 kcyc   0% x1.0 84.2 us/call
  canvas_buffer_wait        52.22 ms  31.33 Mcyc  82% x1.0 52224.9 us/call
```

Wall min/avg/max: 60.49/63.32/75.55 ms. Timeline work averages 10.86 ms; peak render in this window is 24.981 ms. Scope attribution is inclusive; shared raster counters can first-parent under a different draw path. The wait absorbs the remainder of each display window.

### Per-preset table

Each bucket includes its held frames and following transition. Clean holds exclude advance-straddling windows and use the modal `scan_mesh_raster` call count within each preset; rows rank by the clean timeline cost. V/E/F/I describes the IslamicStars seed recorded at spawn, not every intermediate build mesh. Windows aggregate repeated visits to a preset.

| Shape | Geometry | Blends/f | Clean timeline ms/f | Clean render ms/f | Peak ms | Spilled/frames | Clean windows | fps |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| truncatedIcosidodecahedron | Registry solid | 15914 | 26.11 | 26.34 | 34.681 | 0/125 | 3/7 | 16 |
| snubDodecahedron | Registry solid | 14331 | 23.04 | 23.27 | 28.924 | 0/125 | 4/7 | 16 |
| rhombicosidodecahedron | Registry solid | 13895 | 21.40 | 21.64 | 23.453 | 0/125 | 3/7 | 16 |
| truncatedIcosahedron | Registry solid | 13816 | 21.38 | 21.62 | 24.432 | 0/113 | 3/6 | 16 |
| truncatedDodecahedron | Registry solid | 13812 | 21.07 | 21.30 | 24.865 | 0/113 | 3/6 | 16 |
| truncatedCuboctahedron | Registry solid | 13498 | 20.83 | 21.06 | 24.769 | 0/250 | 6/13 | 16 |
| truncatedCube | Registry solid | 12527 | 19.86 | 20.09 | 24.011 | 0/113 | 3/6 | 16 |
| icosidodecahedron | Registry solid | 12875 | 19.60 | 19.83 | 22.950 | 0/351 | 14/19 | 16 |
| snubCube | Registry solid | 12857 | 19.56 | 19.78 | 24.981 | 0/125 | 3/7 | 16 |
| truncatedOctahedron | Registry solid | 12543 | 18.70 | 18.93 | 21.510 | 0/226 | 6/12 | 16 |
| dodecahedron | Registry solid | 12342 | 17.97 | 18.20 | 20.458 | 0/363 | 12/19 | 16 |
| rhombicuboctahedron | Registry solid | 12393 | 17.36 | 17.59 | 19.852 | 0/113 | 4/7 | 16 |
| truncatedTetrahedron | Registry solid | 11706 | 17.24 | 17.48 | 21.343 | 0/226 | 6/13 | 16 |
| icosahedron | Registry solid | 12156 | 16.81 | 17.05 | 20.207 | 0/226 | 10/13 | 16 |
| cuboctahedron | Registry solid | 11812 | 15.72 | 15.95 | 17.800 | 0/589 | 24/33 | 16 |
| cube | Registry solid | 11493 | 15.44 | 15.67 | 17.432 | 0/463 | 17/26 | 16 |
| octahedron | Registry solid | 11369 | 14.54 | 14.77 | 18.056 | 0/452 | 18/25 | 16 |
| tetrahedron | Registry solid | 11073 | 12.73 | 12.96 | 15.015 | 0/238 | 10/13 | 16 |
| startup | Registry solid | 11071 | 12.23 | 12.46 | 14.569 | 0/111 | 6/6 | 16 |

### Per-pixel figures

Peak-window blends: 15,914/frame (1.53× quadrant coverage), 40.7 cycles/blend. Raster scope: 877.2 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```
isr_wake         1159/frame  0.60/1.70/11.37 us  cpu 3.13%
isr_pack         145/frame  6.51/7.18/9.56 us  cpu 1.65%
isr_dma_submit   145/frame  0.63/0.94/1.02 us  cpu 0.21%
```

Pack and submit run inside `isr_wake`; the wake row alone is the inclusive measured flywheel CPU share. DMA completion and other interrupts are unmeasured. These captures predate ISR snapshot correction `fdc7930c5`, so the printed CPU-share percentages carry interval-accounting bias and cannot establish an exact foreground render budget.

Average packing costs 7.7 times submission. Each segment contains 72 pixels: `HD107SFrame<72>` uses 300 bytes per image, including its 8-byte end frame. The image-plus-black composite is 600 bytes. At 24 MHz SPI, 92 functional clocks per byte at 240 MHz gives a maximum **230 us** asynchronous wire transfer, excluded from CPU submission time.

The observed 34.681 ms render peak leaves 27.819 ms within the 62.5 ms display window. No captured frame needs a speedup to meet cadence; this is observed margin, not an exact CPU budget.

## Summary ranking

1. `hk_timeline_step`: 42% of frame, 26.11 ms/frame.
2. `pov_preserve_half`: 0% of frame, 0.15 ms/frame.
3. `canvas_clear`: 0% of frame, 0.08 ms/frame.

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

`targets/Profile/Profile.ino` supports target, window, epoch and scheduling knobs. `just profile HankinSolids` routes through the supported locked wrapper; use the Reproduce command for this complete-cycle capture.
