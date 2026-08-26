# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_islamicstars_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"`). Raw capture:
`build/prof/islamicstars_o3.log`, captured 2026-08-26 07:53 local. This refresh replaces the
earlier 02:14 capture with the final clean cadence-reclaim image.
Its `.provenance`, `_envdump.txt`, and `_build.log` sidecars identify the exact
source, compiler flags, ELF hashes, and image sizes used below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, tip `0df961b818ae08c5f58639edb93d812164a9356a` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, flags `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size: `FLASH: code:152,408, data:193,368, headers:8,528` / `RAM1: variables:315,520, code:56,664, padding:8,872, free:143,232` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 2577–2592 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.0 ppm**
(`tools/parse_profile.py build/prof/islamicstars_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `is_timeline_step` avg
22.01 ms/f, worst window 28.12 ms/f
(frames 2801–2816), peak frame render 50.80 ms,
spilled 0/3328 frames (0.0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
The pass is green: every captured frame stays within one display interval, and
the worst frame retains 11.70 ms of render margin. The
`canvas_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: marker-defined presets own their following transitions; the
first section is the worst exact-render window and the second is the
lowest-cost marked regime.

### Worst exact-render regime (window frames 2801–2816)

```text
frame                         64.82ms 38.89Mcyc 100.0%
  is_timeline_step            28.12ms 16.87Mcyc  43.4%
    is_build_draw             23.40ms 14.04Mcyc  36.1%
      is_build_scan           23.37ms 14.02Mcyc  36.0% x0.8 31155us/c
    hk_conway_compile         181.4us 108.8kcyc   0.3% x0.8 242us/c
    hk_conway_sweep           474.0us 284.4kcyc   0.7% x0.8 632us/c
    is_draw_shape              2.69ms  1.61Mcyc   4.2%
      is_mesh_scan             2.68ms  1.61Mcyc   4.1%
          filter_blend         1.22ms 734.9kcyc   1.9% x16324.2 45cyc/blend
      is_face_offsets           7.0us   4.2kcyc   0.0% x0.2 28us/c
  is_ripple_prepare             0.4us    238cyc   0.0% x1.0 238cyc/c
  canvas_clear                 86.7us  52.0kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait          36.61ms 21.97Mcyc  56.5% x1.0 36611us/c
```

Wall min/avg/max = 47.86/64.82/75.90 ms.
`is_timeline_step` averages 28.12 ms/f here, while the exact frame-render
peak is 50.80 ms. The phase retains 11.70 ms against
one display interval; display-sync wait is idle rather than render work.

### Lowest-cost marked regime (window frames 2689–2704)

```text
  filter_blend                827.0us 496.2kcyc   1.3% x11615.6 43cyc/blend
frame                         62.50ms 37.50Mcyc 100.0%
  is_timeline_step             9.38ms  5.63Mcyc  15.0%
    is_build_draw              8.94ms  5.36Mcyc  14.3%
      is_build_scan            8.93ms  5.36Mcyc  14.3% x1.0 8932us/c
      is_mesh_transform         3.4us   2.0kcyc   0.0% x1.0 3us/c
    hk_conway_compile          20.0us  12.0kcyc   0.0% x1.0 20us/c
    hk_conway_sweep            54.4us  32.7kcyc   0.1% x1.0 54us/c
  is_ripple_prepare             0.4us    240cyc   0.0% x1.0 240cyc/c
  canvas_clear                 84.7us  50.8kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait          53.03ms 31.82Mcyc  84.9% x1.0 53035us/c
```

Wall min/avg/max = 60.48/62.50/64.27 ms.
`is_timeline_step` falls to 9.38 ms/f and the exact render peak is
11.63 ms, leaving 50.87 ms of margin. The spread isolates
the preset-dependent mesh build, Conway compilation, and held-shape raster work from the fixed display cadence.

### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `is_timeline_step` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` | V=20 E=30 F=12 I=60 | 12/12 | 18,621.5 | 44.84 | 48.98 | 16.0 |
| `truncatedIcosidodecahedron_bevel5_relax_hk77` | V=120 E=180 F=62 I=360 | 8/8 | 18,755.4 | 39.34 | 43.15 | 16.0 |
| `truncatedIcosidodecahedron_truncate50d_ambo_dual` | V=120 E=180 F=62 I=360 | 10/10 | 19,515.8 | 36.26 | 50.80 | 16.0 |
| `truncatedOctahedron_gyro_kis_hk17` | V=24 E=36 F=14 I=72 | 12/12 | 18,229.4 | 36.09 | 43.80 | 16.0 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin59` | V=60 E=90 F=32 I=180 | 8/8 | 16,123.9 | 35.94 | 38.10 | 16.0 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin73` | V=60 E=90 F=32 I=180 | 10/10 | 17,489.6 | 34.89 | 37.93 | 16.0 |
| `truncatedIcosahedron_hk54_ambo_hk72` | V=60 E=90 F=32 I=180 | 8/8 | 17,847.2 | 31.92 | 36.74 | 16.0 |
| `truncatedIcosahedron_ambo_relax_truncate33_hk64` | V=60 E=90 F=32 I=180 | 8/8 | 16,976.1 | 30.28 | 33.28 | 16.0 |
| `snubDodecahedron_truncate5d_ambo_dual` | V=60 E=150 F=92 I=300 | 10/10 | 18,414.6 | 29.78 | 37.58 | 16.0 |
| `dodecahedron_bevel2_relax_gyro` | V=20 E=30 F=12 I=60 | 12/12 | 18,367.9 | 29.05 | 37.94 | 16.0 |
| `icosahedron_snub_relax_truncate033_hankin62` | V=12 E=30 F=20 I=60 | 4/4 | 16,389.8 | 28.88 | 32.62 | 16.0 |
| `truncatedIcosahedron_hk58_chamfer63` | V=60 E=90 F=32 I=180 | 8/8 | 16,524.4 | 27.65 | 28.56 | 16.0 |
| `truncatedIcosahedron_truncate50d_ambo_dual` | V=60 E=90 F=32 I=180 | 5/5 | 16,503.3 | 27.48 | 38.83 | 16.0 |
| `dodecahedron_ambo_bevel33_relax_hk66` | V=20 E=30 F=12 I=60 | 10/10 | 15,972.4 | 26.45 | 28.67 | 16.0 |
| `rhombicuboctahedron_hk63_ambo_hk63` | V=24 E=48 F=26 I=96 | 10/10 | 15,156.8 | 23.82 | 26.50 | 16.0 |
| `dodecahedron_hk72_ambo_dual_hk20` | V=20 E=30 F=12 I=60 | 5/5 | 14,702.4 | 23.62 | 28.48 | 16.0 |
| `icosahedron_kis_gyro` | V=12 E=30 F=20 I=60 | 14/14 | 16,359.7 | 23.33 | 28.21 | 16.0 |
| `dodecahedron_hk54_ambo_hk72` | V=20 E=30 F=12 I=60 | 10/10 | 14,587.6 | 22.52 | 25.78 | 16.0 |
| `dodecahedron_hk62_ambo_hk62` | V=20 E=30 F=12 I=60 | 10/10 | 14,171.7 | 20.71 | 23.60 | 16.0 |
| `icosahedron_ambo_truncate033_hankin59` | V=12 E=30 F=20 I=60 | 8/8 | 13,933.4 | 20.10 | 22.62 | 16.0 |
| `icosidodecahedron_truncate5d_ambo_dual` | V=30 E=60 F=32 I=120 | 10/10 | 15,025.4 | 19.81 | 23.68 | 16.0 |
| `octahedron_hk17_ambo_hk73` | V=6 E=12 F=8 I=24 | 8/8 | 13,333.4 | 18.99 | 22.32 | 16.0 |
| `octahedron_hk34_ambo_hk72` | V=6 E=12 F=8 I=24 | 8/8 | 13,040.4 | 17.04 | 18.14 | 16.0 |

All 23 marker-defined shapes are visited and a repeated shape confirms cycle closure.

### Per-pixel figures

`filter_blend` averages 16,324.2 blended px/frame versus the 10,368-px quadrant (1.57× coverage), at 45.0 cyc/blend. `is_timeline_step` contributes 1033.7 cycles per blended pixel in the selected window; this is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1196.1/frame  min/avg/max 0.48/1.61/14.77 us  cpu 2.97%
isr_pack         149.5/frame  min/avg/max 5.99/6.96/9.33 us  cpu 1.60%
isr_dma_submit   149.5/frame  min/avg/max 0.68/0.93/7.38 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are
  reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains
  isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.78% CPU. The worst render
  already fits one display interval with 11.70 ms of margin.

## Summary ranking

1. `is_timeline_step` — 35.3% of measured root time, 22.01 ms/f average.
2. `is_draw_shape` — 22.9% of measured root time, 14.29 ms/f average.
3. `is_mesh_scan` — 22.2% of measured root time, 13.86 ms/f average.

No matched WASM/native datum is asserted here; the final paired on-device
shipping capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `is_mesh_scan`; its subtree is hidden
  in windows where that parent has zero calls, and calls approximate blended
  pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths
  identified by the counter tree; global O3 compiles the entire single-effect
  image.
- The capture uses `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` to keep the full cycle inside one effect epoch; these knobs change dwell/epoch length, not per-frame render cost.
- Provenance records clean O3 source `0df961b818ae08c5f58639edb93d812164a9356a` and clean paired
  shipping source `0df961b818ae08c5f58639edb93d812164a9356a`; the capture artifacts record no
  uncommitted source state.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=IslamicStars`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` builds, flashes, and captures the
global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The final paired shipping capture peaks at 50.65 ms versus
50.80 ms here: global O3 raises the peak by 0.15 ms
(0.30%). O3-image minus shipping-image size deltas are **FLASH code
+24,216 B** and **ITCM code +13,168 B**. Both captures are green:
shipping 0/3328 spills versus O3
0/3328.
