# Finding 2 IslamicStars performance comparison — 2026-09-24

Unlanded candidate fa92838f9 vs baseline b4a9848ff, which reverts only finding 2 in face.h and test_sdf.h. Same Teensy 4.0 COM4, 600 MHz, TS=4, 23-shape cycle, 210 seconds per capture. Setup frame 1 excluded.

The comparison measures the widened conservative face bound. TS=4 changes choreography sampling as well as hold duration; these captures are a matched TS=4 workload.

| Config | Matched frames | Before mean ms | After mean ms | Δ mean ms | Before peak ms | After peak ms | Before/after spills |
|---|---:|---:|---:|---:|---:|---:|---|
| ship | 3327 | 23.069846 | 25.262962 | +2.193115 | 50.859 | 60.695 | 0/3327 → 0/3327 |
| o3 | 3327 | 22.799689 | 24.816086 | +2.016397 | 50.353 | 55.674 | 0/3327 → 0/3327 |

Shipping mean render increases 9.51%. Observed peak headroom falls from 11.641 to **1.805 ms**. No spills were observed. The largest shape mean increase is `truncatedIcosidodecahedron_truncate50d_ambo_dual`: 9.432 ms (32.16%).

## ship: scope and workload changes

| Inclusive scope | Before ms/frame | After ms/frame | Δ ms/frame |
|---|---:|---:|---:|
| is_timeline_step | 22.878568 | 25.081302 | +2.202733 |
| is_build_draw | 7.242338 | 7.943362 | +0.701024 |
| is_draw_shape | 14.942907 | 16.441681 | +1.498773 |
| is_mesh_scan | 14.515534 | 16.014202 | +1.498667 |
| scan_face_setup | 3.601600 | 3.883013 | +0.281413 |
| scan_mesh_raster | 17.986360 | 19.900809 | +1.914449 |
| filter_blend | 1.033089 | 1.049728 | +0.016638 |

Blended pixels/frame 14895.963 → 14911.995 (+0.1076%). Inclusive raster cycles/blend 724.479 → 800.730. Face setup and raster are distinct scopes; raster covers build paths outside is_mesh_scan. Nested/mixed-parent totals must not be added as an exclusive cost breakdown.

Candidate minus baseline image code: FLASH +144 B; ITCM +256 B.


## ship: matched shape deltas

Exact paired telemetry: [CSV](matched-ship.csv). Geometry and shape ownership match on all 3327 paired frames.

| Shape | Mean before ms | Mean after ms | Δ mean ms | Peak before ms | Peak after ms |
|---|---:|---:|---:|---:|---:|
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 29.326 | 38.758 | +9.432 | 50.859 | 60.695 |
| snubDodecahedron_truncate5d_ambo_dual | 22.471 | 27.739 | +5.268 | 41.767 | 46.558 |
| truncatedIcosahedron_truncate50d_ambo_dual | 20.682 | 24.981 | +4.299 | 39.858 | 41.540 |
| dodecahedron_bevel2_relax_gyro | 19.733 | 23.460 | +3.727 | 41.431 | 48.242 |
| icosidodecahedron_truncate5d_ambo_dual | 15.360 | 18.495 | +3.134 | 23.484 | 26.548 |
| icosahedron_kis_gyro | 15.448 | 18.307 | +2.860 | 29.161 | 52.233 |
| truncatedOctahedron_gyro_kis_hk17 | 29.053 | 31.756 | +2.704 | 45.656 | 51.070 |
| truncatedIcosidodecahedron_bevel5_relax_hk77 | 32.369 | 34.457 | +2.087 | 44.044 | 46.521 |
| truncatedIcosahedron_hk54_ambo_hk72 | 28.054 | 29.957 | +1.903 | 39.146 | 41.988 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 28.652 | 30.515 | +1.864 | 40.371 | 43.295 |
| icosahedron_snub_relax_truncate033_hankin62 | 23.822 | 25.430 | +1.608 | 31.725 | 33.308 |
| truncatedIcosahedron_hk58_chamfer63 | 24.538 | 26.087 | +1.550 | 29.612 | 31.470 |
| truncatedIcosahedron_ambo_relax_truncate33_hk64 | 25.661 | 27.090 | +1.429 | 33.926 | 35.820 |
| dodecahedron_ambo_bevel33_relax_hk66 | 21.213 | 22.303 | +1.090 | 30.325 | 32.002 |
| dodecahedron_hk72_ambo_dual_hk20 | 19.117 | 20.116 | +0.999 | 29.312 | 32.578 |
| dodecahedron_hk54_ambo_hk72 | 19.664 | 20.515 | +0.850 | 27.638 | 28.816 |
| rhombicuboctahedron_hk63_ambo_hk63 | 22.412 | 23.239 | +0.827 | 29.790 | 31.109 |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 32.715 | 33.536 | +0.821 | 49.581 | 51.735 |
| icosahedron_ambo_truncate033_hankin59 | 18.005 | 18.629 | +0.624 | 25.134 | 26.100 |
| dodecahedron_hk62_ambo_hk62 | 19.151 | 19.756 | +0.605 | 23.378 | 24.200 |
| octahedron_hk17_ambo_hk73 | 17.186 | 17.690 | +0.504 | 21.155 | 21.434 |
| octahedron_hk34_ambo_hk72 | 16.297 | 16.763 | +0.466 | 19.419 | 19.991 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 29.610 | 29.845 | +0.235 | 42.014 | 42.011 |

First complete shape cycle (1774 runtime frames): mean 23.007499 → 25.217662 ms; peak 50.246 → 60.201 ms.

## o3: scope and workload changes

| Inclusive scope | Before ms/frame | After ms/frame | Δ ms/frame |
|---|---:|---:|---:|
| is_timeline_step | 22.607654 | 24.631925 | +2.024271 |
| is_build_draw | 7.206591 | 7.850572 | +0.643982 |
| is_draw_shape | 14.811039 | 16.189293 | +1.378255 |
| is_mesh_scan | 14.384819 | 15.763543 | +1.378724 |
| scan_face_setup | 3.519980 | 3.795960 | +0.275980 |
| scan_mesh_raster | 17.900470 | 19.651979 | +1.751509 |
| filter_blend | 1.141767 | 1.127049 | -0.014718 |

Blended pixels/frame 14878.598 → 14894.526 (+0.1070%). Inclusive raster cycles/blend 721.861 → 791.646. Face setup and raster are distinct scopes; raster covers build paths outside is_mesh_scan. Nested/mixed-parent totals must not be added as an exclusive cost breakdown.

Candidate minus baseline image code: FLASH +32 B; ITCM +32 B.


## o3: matched shape deltas

Exact paired telemetry: [CSV](matched-o3.csv). Geometry and shape ownership match on all 3327 paired frames.

| Shape | Mean before ms | Mean after ms | Δ mean ms | Peak before ms | Peak after ms |
|---|---:|---:|---:|---:|---:|
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 28.587 | 37.341 | +8.754 | 50.353 | 55.674 |
| snubDodecahedron_truncate5d_ambo_dual | 21.996 | 26.881 | +4.885 | 37.570 | 41.971 |
| truncatedIcosahedron_truncate50d_ambo_dual | 20.662 | 24.628 | +3.966 | 38.633 | 40.065 |
| dodecahedron_bevel2_relax_gyro | 19.456 | 22.909 | +3.453 | 37.979 | 44.120 |
| icosidodecahedron_truncate5d_ambo_dual | 15.294 | 18.224 | +2.929 | 23.625 | 26.179 |
| icosahedron_kis_gyro | 15.230 | 17.855 | +2.625 | 28.164 | 51.096 |
| truncatedOctahedron_gyro_kis_hk17 | 28.460 | 30.931 | +2.471 | 44.725 | 49.723 |
| truncatedIcosidodecahedron_bevel5_relax_hk77 | 31.509 | 33.445 | +1.936 | 43.557 | 45.852 |
| truncatedIcosahedron_hk54_ambo_hk72 | 27.676 | 29.431 | +1.755 | 38.196 | 40.870 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 28.234 | 29.967 | +1.733 | 39.992 | 42.736 |
| icosahedron_snub_relax_truncate033_hankin62 | 23.968 | 25.427 | +1.459 | 34.210 | 35.677 |
| truncatedIcosahedron_hk58_chamfer63 | 24.167 | 25.568 | +1.400 | 28.817 | 30.638 |
| truncatedIcosahedron_ambo_relax_truncate33_hk64 | 25.483 | 26.810 | +1.327 | 34.198 | 36.262 |
| dodecahedron_ambo_bevel33_relax_hk66 | 20.852 | 21.830 | +0.978 | 30.158 | 31.670 |
| dodecahedron_hk72_ambo_dual_hk20 | 19.123 | 20.021 | +0.898 | 29.753 | 31.446 |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 32.409 | 33.187 | +0.779 | 48.867 | 50.986 |
| dodecahedron_hk54_ambo_hk72 | 19.404 | 20.166 | +0.762 | 26.719 | 27.667 |
| rhombicuboctahedron_hk63_ambo_hk63 | 22.094 | 22.791 | +0.697 | 26.769 | 27.803 |
| dodecahedron_hk62_ambo_hk62 | 19.208 | 19.748 | +0.541 | 23.974 | 24.760 |
| icosahedron_ambo_truncate033_hankin59 | 17.986 | 18.527 | +0.541 | 23.838 | 24.384 |
| octahedron_hk17_ambo_hk73 | 17.496 | 17.917 | +0.421 | 23.077 | 23.556 |
| octahedron_hk34_ambo_hk72 | 16.215 | 16.594 | +0.379 | 19.419 | 19.806 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 29.358 | 29.579 | +0.221 | 40.670 | 40.691 |

First complete shape cycle (1774 runtime frames): mean 22.774051 → 24.800521 ms; peak 48.867 → 55.674 ms.

## Limits

One sequential run per image/configuration. Interrupt phase and timing noise remain; this is not a repeated-run confidence interval. Scope costs include interrupts and shared parents. Peak values can come from different matched frame numbers. Full per-frame CSVs preserve exact comparisons.

Raw captures and provenance are published in [the capture manifest](capture_manifest.json), with untouched capture text under `data/`. Local ELF binaries and build intermediates are not published. Full baseline and candidate reports carry build sizes, compiler identity, hashes, phase trees, geometry and ISR measurements.

## Full reports and evidence

- [Shipping baseline](baseline_shipping.md) and [candidate](../shipping/profile_islamicstars_teensy_2026-09-24.md).
- [O3 baseline](baseline_O3.md) and [candidate](../O3/profile_islamicstars_teensy_2026-09-24.md).
- [Analysis data](analysis.json), [shipping paired frames](matched-ship.csv), and [O3 paired frames](matched-o3.csv).

Finding 2 remains deferred; the candidate measurements do not describe shipped code.

Recompute timing metrics from the repository root. ELF hash checks will report unavailable local artifacts; the original validation results and hashes are retained as historical evidence:

```powershell
python docs/profiles/review20260924finding2/analyze_profile2.py --manifest docs/profiles/review20260924finding2/capture_manifest.json --out build/prof/review125_analysis
```

## Visual comparison

[Interactive before/after gallery](visual/index.html) | [Contact sheet](visual/comparison.png) | [Capture method and metrics](visual/README.md)

A separate deterministic full-canvas WASM comparison at default transition speed 1 covered 7,104 frames and all 23 shapes. The gallery selects 29 frames, including each shape’s maximum total absolute sRGB difference. Mean changed pixels: 45.59 / 41,472 (0.11%); maximum: 1,414 (3.41%). These are renderer captures, not hardware photographs. Difference views are amplified 16×. Finding 2 was deferred after visual review; see [the disposition](DISPOSITION.md).
