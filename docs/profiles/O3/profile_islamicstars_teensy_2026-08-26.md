# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_islamicstars_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"`). Raw capture:
`build/prof/islamicstars_o3.log`, captured 2026-08-26 02:14 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, flags `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size: `FLASH: code:151,432, data:193,360, headers:8,488` / `RAM1: variables:315,520, code:55,608, padding:9,928, free:143,232` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 2577–2592
root counter cycles ÷ 600 MHz match the measured wall sum within **2.3 ppm**
(`tools/parse_profile.py build/prof/islamicstars_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `is_timeline_step` avg
26.15 ms/f, worst window 57.39 ms/f
(frames 2577–2592),
peak frame render 63.74 ms, spilled 1/3328
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
The pass is yellow: 0% of frames exceed one display interval. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 2801–2816)

```text
frame                          69.17ms  41.50Mcyc 100.0%
  is_timeline_step             34.22ms  20.53Mcyc  49.5%
    is_build_draw              29.12ms  17.47Mcyc  42.1%
      is_build_scan            29.09ms  17.45Mcyc  42.1% x0.8 38785us/c
    hk_conway_compile          182.6us  109.6kcyc   0.3% x0.8 244us/c
    hk_conway_sweep            479.9us  288.0kcyc   0.7% x0.8 640us/c
    is_draw_shape               3.05ms   1.83Mcyc   4.4%
      is_mesh_scan              3.04ms   1.83Mcyc   4.4%
          filter_blend          1.27ms  761.5kcyc   1.8% x16364.6 0us/c
      is_face_offsets            6.9us    4.2kcyc   0.0% x0.2 28us/c
  is_ripple_prepare              0.3us     200cyc   0.0% x1.0 0us/c
  canvas_clear                  87.4us   52.5kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           34.86ms  20.92Mcyc  50.4% x1.0 34862us/c
```

Wall min/avg/max = 43.33/69.17/123.53 ms. The `is_timeline_step` scope averages 34.22 ms/f in this window, while the exact frame-render peak is 63.74 ms. That is 1.24 ms of overrun against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 2689–2704)

```text
  filter_blend                 849.1us  509.5kcyc   1.4% x11615.6 0us/c
frame                          62.52ms  37.51Mcyc 100.0%
  is_timeline_step             10.54ms   6.33Mcyc  16.9%
    is_build_draw              10.08ms   6.05Mcyc  16.1%
      is_build_scan            10.08ms   6.05Mcyc  16.1% x1.0 10077us/c
      is_mesh_transform          3.1us    1.8kcyc   0.0% x1.0 3us/c
    hk_conway_compile           29.0us   17.4kcyc   0.0% x1.0 29us/c
    hk_conway_sweep             59.4us   35.7kcyc   0.1% x1.0 59us/c
  is_ripple_prepare              0.2us     139cyc   0.0% x1.0 0us/c
  canvas_clear                  84.6us   50.8kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           51.89ms  31.14Mcyc  83.0% x1.0 51893us/c
```

Wall min/avg/max = 60.30/62.52/64.45 ms. The `is_timeline_step` scope averages 10.54 ms/f in this window, while the exact frame-render peak is 12.96 ms. That is 49.54 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `is_timeline_step` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` | V=20 E=30 F=12 I=60 | 12/12 | 18,621.5 | 57.39 | 61.36 | 16.0 |
| `truncatedIcosidodecahedron_bevel5_relax_hk77` | V=120 E=180 F=62 I=360 | 8/8 | 18,755.4 | 47.54 | 52.52 | 16.0 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin59` | V=60 E=90 F=32 I=180 | 8/8 | 16,123.9 | 44.21 | 46.55 | 16.0 |
| `truncatedOctahedron_gyro_kis_hk17` | V=24 E=36 F=14 I=72 | 12/12 | 18,229.4 | 41.57 | 51.06 | 16.0 |
| `truncatedIcosidodecahedron_truncate50d_ambo_dual` | V=120 E=180 F=62 I=360 | 10/10 | 19,495.6 | 41.51 | 63.74 | 8.0 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin73` | V=60 E=90 F=32 I=180 | 10/10 | 17,489.6 | 41.17 | 44.30 | 16.0 |
| `truncatedIcosahedron_hk54_ambo_hk72` | V=60 E=90 F=32 I=180 | 8/8 | 17,876.2 | 38.78 | 44.52 | 16.0 |
| `truncatedIcosahedron_ambo_relax_truncate33_hk64` | V=60 E=90 F=32 I=180 | 8/8 | 17,061.5 | 36.87 | 39.91 | 16.0 |
| `icosahedron_snub_relax_truncate033_hankin62` | V=12 E=30 F=20 I=60 | 4/4 | 16,389.8 | 34.39 | 38.47 | 16.0 |
| `snubDodecahedron_truncate5d_ambo_dual` | V=60 E=150 F=92 I=300 | 10/10 | 18,384.8 | 33.51 | 41.96 | 16.0 |
| `truncatedIcosahedron_truncate50d_ambo_dual` | V=60 E=90 F=32 I=180 | 5/5 | 16,503.3 | 33.09 | 47.50 | 16.0 |
| `dodecahedron_bevel2_relax_gyro` | V=20 E=30 F=12 I=60 | 12/12 | 18,367.9 | 32.81 | 42.09 | 16.0 |
| `truncatedIcosahedron_hk58_chamfer63` | V=60 E=90 F=32 I=180 | 8/8 | 16,524.4 | 32.64 | 33.70 | 16.0 |
| `dodecahedron_ambo_bevel33_relax_hk66` | V=20 E=30 F=12 I=60 | 10/10 | 15,972.4 | 31.79 | 33.90 | 16.0 |
| `rhombicuboctahedron_hk63_ambo_hk63` | V=24 E=48 F=26 I=96 | 10/10 | 15,246.4 | 28.50 | 32.90 | 16.0 |
| `dodecahedron_hk72_ambo_dual_hk20` | V=20 E=30 F=12 I=60 | 5/5 | 14,702.4 | 27.40 | 32.69 | 16.0 |
| `icosahedron_kis_gyro` | V=12 E=30 F=20 I=60 | 14/14 | 16,359.7 | 26.68 | 31.88 | 16.0 |
| `dodecahedron_hk54_ambo_hk72` | V=20 E=30 F=12 I=60 | 10/10 | 14,540.1 | 26.29 | 31.13 | 16.0 |
| `dodecahedron_hk62_ambo_hk62` | V=20 E=30 F=12 I=60 | 10/10 | 14,171.7 | 24.34 | 29.70 | 16.0 |
| `icosahedron_ambo_truncate033_hankin59` | V=12 E=30 F=20 I=60 | 8/8 | 14,222.1 | 23.77 | 27.32 | 16.0 |
| `icosidodecahedron_truncate5d_ambo_dual` | V=30 E=60 F=32 I=120 | 10/10 | 15,019.8 | 22.29 | 25.45 | 16.0 |
| `octahedron_hk17_ambo_hk73` | V=6 E=12 F=8 I=24 | 8/8 | 13,333.4 | 21.94 | 26.00 | 16.0 |
| `octahedron_hk34_ambo_hk72` | V=6 E=12 F=8 I=24 | 8/8 | 13,039.3 | 19.73 | 21.29 | 16.0 |

208 shape markers cover 23 distinct entries; cycle closure is confirmed by a repeated shape.

### Per-pixel figures

`filter_blend` averages 16,364.6 blended px/frame versus the 10,368-px quadrant (1.58× coverage), at 46.5 cyc/blend. `is_timeline_step` contributes 1254.8 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1276.2/frame  min/avg/max 0.41/1.57/11.02 us  cpu 2.89%
isr_pack          159.5/frame  min/avg/max 5.99/6.91/9.36 us  cpu 1.59%
isr_dma_submit    159.5/frame  min/avg/max 0.61/0.93/1.12 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.69% CPU. The worst render requires 2.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `is_timeline_step` — 41.9% of measured root time, 26.15 ms/f average.
2. `is_draw_shape` — 27.5% of measured root time, 17.17 ms/f average.
3. `is_mesh_scan` — 26.8% of measured root time, 16.74 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `is_mesh_scan`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=IslamicStars`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 64.18 ms versus
63.74 ms here: global O3 lowers the peak by 0.43 ms (0.7%). O3-image minus shipping-image
size deltas are **FLASH code +23,240 B** and **ITCM code
+12,368 B**. Spill fractions are compared rather than raw counts:
shipping 2/3328 (0.1%)
versus O3 1/3328 (0%).
