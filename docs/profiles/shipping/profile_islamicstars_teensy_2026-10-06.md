# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/islamicstars_ship.log`, captured 2026-10-06 18:23 on COM3.
Replaces `profile_islamicstars_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size (`profile` env, this effect only): `FLASH: code:127124, data:202212, headers:8580` / `RAM1: variables:315424, code:37816, padding:27720, free:143328` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2577–2592 root counter cyc ÷ 600 MHz matches the measured wall sum within **4.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `is_timeline_step` averages 22.45 ms/f; its worst window is 44.34 ms/f (frames 2577–2592). Peak frame render is **49.86 ms** (frame 2811), and **0/3327** frames spilled. Setup frame 1 is excluded from both; it rendered 15.60 ms.

The previous shipping report (2026-09-28 18:40) recorded peak 🟢 50.51 (23) and spilled 🟢 0/3327 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 23 shape entries, each announced by a `Spawning Shape` line; each owns its hold and the transition that follows it. The capture runs all 23 from `dodecahedron_hk62_ambo_hk62`, wraps back to it, and repeats 21 of them before the capture ends. The block below is the window holding the pass's peak frame.

### Peak window (frames 2801–2816)

```
frame                    64.94 ms  38.96 Mcyc   100%
  pov_preserve_half      146.1 us   87.7 kcyc     0%
  is_timeline_step       28.14 ms  16.88 Mcyc    43%
    is_build_draw        23.22 ms  13.93 Mcyc    36%  x0.8  30963 us/c
      is_build_scan      23.19 ms  13.91 Mcyc    36%  x0.8  30919 us/c
      is_mesh_transform    2.1 us    1.3 kcyc     0%  x0.2  5174 cyc/c
    hk_conway_compile    192.1 us  115.3 kcyc     0%  x0.8  256 us/c
    hk_conway_sweep      559.3 us  335.6 kcyc     1%  x0.8  746 us/c
    is_draw_shape         2.65 ms   1.59 Mcyc     4%  x0.2  10614 us/c
      is_mesh_scan        2.64 ms   1.59 Mcyc     4%  x0.2  10579 us/c
        scan_mesh_raster 23.19 ms  13.91 Mcyc    36%  x287  81 us/c
          filter_blend    1.09 ms  656.0 kcyc     2%  x16060  41 cyc/c
        scan_face_setup   2.52 ms   1.51 Mcyc     4%  x287  5263 cyc/c
      is_face_offsets      6.4 us    3.8 kcyc     0%  x0.2  26 us/c
  is_ripple_prepare        0.4 us     242 cyc     0%
  canvas_clear            86.3 us   51.8 kcyc     0%
  canvas_buffer_wait     36.56 ms  21.94 Mcyc    56%
```

Wall min/avg/max = 50.89/64.94/73.51 ms. Per-frame values are window averages; `xN` is calls per frame. `scan_mesh_raster` and `scan_face_setup` are entered from both the build scan and the shape scan, so they print under the first parent but carry both callers' time (validate flags them MIXED-PARENT); the build scan holds most of it.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `is_timeline_step` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each shape, windows straddling an advance excluded). V/E/F/I are from each shape's `Spawning Shape` line. The cycle wraps back to the first shape within the capture (validate: shape markers advance and return to the first).

| Shape | V/E/F/I | Peak render ms | Spilled/frames | Clean `is_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| truncatedIcosidodecahedron_truncate50d_ambo_dual | V=120 E=180 F=62 I=360 | 🟢 49.86 | 0/156 | 36.22 | 8/10 |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | V=20 E=30 F=12 I=60 | 🟢 49.09 | 0/176 | 44.34 | 10/12 |
| truncatedOctahedron_gyro_kis_hk17 | V=24 E=36 F=14 I=72 | 🟢 44.27 | 0/184 | 37.61 | 10/12 |
| truncatedIcosidodecahedron_bevel5_relax_hk77 | V=120 E=180 F=62 I=360 | 🟢 42.76 | 0/144 | 38.60 | 6/8 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | V=60 E=90 F=32 I=180 | 🟢 41.75 | 0/144 | 37.60 | 6/8 |
| snubDodecahedron_truncate5d_ambo_dual | V=60 E=150 F=92 I=300 | 🟢 40.40 | 0/156 | 29.33 | 8/10 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin73 | V=60 E=90 F=32 I=180 | 🟢 40.07 | 0/144 | 36.38 | 8/10 |
| dodecahedron_bevel2_relax_gyro | V=20 E=30 F=12 I=60 | 🟢 40.04 | 0/176 | 28.97 | 10/12 |
| truncatedIcosahedron_truncate50d_ambo_dual | V=60 E=90 F=32 I=180 | 🟢 38.23 | 0/78 | 27.80 | 4/5 |
| truncatedIcosahedron_hk54_ambo_hk72 | V=60 E=90 F=32 I=180 | 🟢 37.81 | 0/140 | 33.00 | 6/8 |
| truncatedIcosahedron_ambo_relax_truncate33_hk64 | V=60 E=90 F=32 I=180 | 🟢 32.92 | 0/144 | 30.67 | 6/8 |
| icosahedron_kis_gyro | V=12 E=30 F=20 I=60 | 🟢 31.54 | 0/240 | 23.04 | 12/14 |
| icosahedron_snub_relax_truncate033_hankin62 | V=12 E=30 F=20 I=60 | 🟢 31.07 | 0/72 | 29.45 | 3/4 |
| dodecahedron_ambo_bevel33_relax_hk66 | V=20 E=30 F=12 I=60 | 🟢 30.00 | 0/156 | 27.39 | 8/10 |
| truncatedIcosahedron_hk58_chamfer63 | V=60 E=90 F=32 I=180 | 🟢 29.68 | 0/124 | 27.70 | 6/8 |
| dodecahedron_hk72_ambo_dual_hk20 | V=20 E=30 F=12 I=60 | 🟢 29.60 | 0/103 | 24.42 | 4/5 |
| rhombicuboctahedron_hk63_ambo_hk63 | V=24 E=48 F=26 I=96 | 🟢 27.50 | 0/140 | 25.18 | 8/10 |
| icosidodecahedron_truncate5d_ambo_dual | V=30 E=60 F=32 I=120 | 🟢 24.41 | 0/156 | 19.86 | 8/10 |
| dodecahedron_hk54_ambo_hk72 | V=20 E=30 F=12 I=60 | 🟢 23.35 | 0/140 | 22.55 | 8/10 |
| icosahedron_ambo_truncate033_hankin59 | V=12 E=30 F=20 I=60 | 🟢 23.24 | 0/136 | 21.20 | 6/8 |
| dodecahedron_hk62_ambo_hk62 | V=20 E=30 F=12 I=60 | 🟢 23.14 | 0/138 | 21.78 | 8/10 |
| octahedron_hk17_ambo_hk73 | V=6 E=12 F=8 I=24 | 🟢 21.61 | 0/140 | 19.31 | 6/8 |
| octahedron_hk34_ambo_hk72 | V=6 E=12 F=8 I=24 | 🟢 19.94 | 0/140 | 18.46 | 6/8 |

### Per-pixel figures

`filter_blend` ran 16,060 times per frame in the peak window at 41 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1198/frame  min/avg/max 0.6/1.7/19.1 us  cpu 3.09%
isr_pack         150/frame  min/avg/max 6.2/7.0/9.8 us  cpu 1.62%
isr_dma_submit   150/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `is_timeline_step` — 43% of the peak window, 28.14 ms/f (`is_build_draw` 23.22).
2. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.09 ms/f.
4. `is_ripple_prepare` — 0% of the peak window, 0.00 ms/f.

README cells: peak 🟢 49.86 (23), spilled 🟢 0/3327 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- `is_mesh_transform`, `scan_face_setup` and `scan_mesh_raster` are not exclusive costs (validate INFO: MIXED-PARENT/DUPLICATE-NAME).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: `transform_shape` and `draw_shape` are `HS_O3_FN`. Shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.
- Dwell-compression knobs change how long a shape holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` builds, flashes and captures under the device lock.
