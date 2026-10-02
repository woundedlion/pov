# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/islamicstars_ship.log`, captured 2026-09-28 18:40 on COM3.
Replaces the earlier 2026-09-28 report of the same name, which predated `97eb0bf78`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size (`profile` env, this effect only): `FLASH: code:125392, data:201084, headers:8372` / `RAM1: variables:315424, code:36632, padding:28904, free:143328` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2577–2592 root counter cyc ÷ 600 MHz matches the measured wall sum within **4.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `is_timeline_step` averages 22.61 ms/f; its worst window is 44.69 ms/f (frames 2577–2592). Peak frame render is **50.51 ms** (frame 2811), and **0/3327** frames spilled. Setup frame 1 is excluded from both; it rendered 15.63 ms.

The previous shipping report (2026-09-28 16:50) recorded peak 🟢 50.500 (23) and spilled 🟢 0/3336 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 23 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2801–2816)

```
frame                     65.06 ms   39.04 Mcyc   100%
  pov_preserve_half       143.7 us    86.2 kcyc     0%
  is_timeline_step        28.54 ms   17.12 Mcyc    44%
    is_build_draw         23.48 ms   14.09 Mcyc    36%
      is_build_scan       23.45 ms   14.07 Mcyc    36%
      is_mesh_transform     2.1 us     1.3 kcyc     0%
    hk_conway_compile     187.0 us   112.2 kcyc     0%
    hk_conway_sweep       553.1 us   331.9 kcyc     1%
    is_draw_shape          2.66 ms    1.60 Mcyc     4%
      is_mesh_scan         2.65 ms    1.59 Mcyc     4%
        scan_mesh_raster  23.37 ms   14.02 Mcyc    36%  x287  48857 cyc/c
          filter_blend     1.08 ms   650.7 kcyc     2%  x16060  41 cyc/c
        scan_face_setup    2.59 ms    1.56 Mcyc     4%  x287  5418 cyc/c
      is_face_offsets       6.8 us     4.1 kcyc     0%
  is_ripple_prepare         0.9 us      583 cyc     0%
  canvas_clear             84.7 us    50.8 kcyc     0%
  canvas_buffer_wait      36.29 ms   21.77 Mcyc    56%
```

Wall min/avg/max = 50.48/65.06/73.90 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `is_timeline_step` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `is_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| truncatedIcosidodecahedron_truncate50d_ambo_dual | V=120 E=180 F=62 I=360 | 🟢 50.51 | 0/156 | 36.58 | 10/10 |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | V=20 E=30 F=12 I=60 | 🟢 49.44 | 0/176 | 44.69 | 12/12 |
| truncatedOctahedron_gyro_kis_hk17 | V=24 E=36 F=14 I=72 | 🟢 44.58 | 0/184 | 37.73 | 12/12 |
| truncatedIcosidodecahedron_bevel5_relax_hk77 | V=120 E=180 F=62 I=360 | 🟢 43.09 | 0/144 | 38.91 | 8/8 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | V=60 E=90 F=32 I=180 | 🟢 41.87 | 0/144 | 37.69 | 8/8 |
| snubDodecahedron_truncate5d_ambo_dual | V=60 E=150 F=92 I=300 | 🟢 40.62 | 0/156 | 29.62 | 10/10 |
| dodecahedron_bevel2_relax_gyro | V=20 E=30 F=12 I=60 | 🟢 40.39 | 0/176 | 29.37 | 12/12 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin73 | V=60 E=90 F=32 I=180 | 🟢 40.24 | 0/144 | 36.51 | 10/10 |
| truncatedIcosahedron_truncate50d_ambo_dual | V=60 E=90 F=32 I=180 | 🟢 38.61 | 0/78 | 28.07 | 5/5 |
| truncatedIcosahedron_hk54_ambo_hk72 | V=60 E=90 F=32 I=180 | 🟢 38.06 | 0/140 | 33.16 | 8/8 |
| truncatedIcosahedron_ambo_relax_truncate33_hk64 | V=60 E=90 F=32 I=180 | 🟢 33.08 | 0/144 | 30.84 | 8/8 |
| icosahedron_kis_gyro | V=12 E=30 F=20 I=60 | 🟢 31.82 | 0/240 | 23.26 | 14/14 |
| icosahedron_snub_relax_truncate033_hankin62 | V=12 E=30 F=20 I=60 | 🟢 31.15 | 0/72 | 29.54 | 4/4 |
| dodecahedron_ambo_bevel33_relax_hk66 | V=20 E=30 F=12 I=60 | 🟢 30.08 | 0/156 | 27.47 | 10/10 |
| truncatedIcosahedron_hk58_chamfer63 | V=60 E=90 F=32 I=180 | 🟢 29.84 | 0/124 | 27.85 | 8/8 |
| dodecahedron_hk72_ambo_dual_hk20 | V=20 E=30 F=12 I=60 | 🟢 29.68 | 0/103 | 24.44 | 5/5 |
| rhombicuboctahedron_hk63_ambo_hk63 | V=24 E=48 F=26 I=96 | 🟢 27.56 | 0/140 | 25.23 | 10/10 |
| icosidodecahedron_truncate5d_ambo_dual | V=30 E=60 F=32 I=120 | 🟢 24.55 | 0/156 | 20.02 | 10/10 |
| dodecahedron_hk54_ambo_hk72 | V=20 E=30 F=12 I=60 | 🟢 23.38 | 0/140 | 22.56 | 10/10 |
| icosahedron_ambo_truncate033_hankin59 | V=12 E=30 F=20 I=60 | 🟢 23.26 | 0/136 | 21.22 | 8/8 |
| dodecahedron_hk62_ambo_hk62 | V=20 E=30 F=12 I=60 | 🟢 23.20 | 0/138 | 21.80 | 10/10 |
| octahedron_hk17_ambo_hk73 | V=6 E=12 F=8 I=24 | 🟢 21.64 | 0/140 | 19.30 | 8/8 |
| octahedron_hk34_ambo_hk72 | V=6 E=12 F=8 I=24 | 🟢 19.92 | 0/140 | 18.44 | 8/8 |

### Per-pixel figures

`filter_blend` ran 16,060 times per frame in the peak window at 41 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1201/frame  min/avg/max 0.5/1.7/11.9 us  cpu 3.12%
isr_pack          150/frame  min/avg/max 6.2/7.0/9.8 us  cpu 1.62%
isr_dma_submit    150/frame  min/avg/max 0.7/0.9/3.7 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `is_timeline_step` — 44% of the peak window, 28.54 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `is_ripple_prepare` — 0% of the peak window, 0.00 ms/f.

README cells: peak 🟢 50.51 (23), spilled 🟢 0/3327 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` builds, flashes and captures under the device lock.
