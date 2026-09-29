# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HankinSolids`).
Raw capture: `build/prof/hankinsolids_ship.log`, captured 2026-09-28 19:27 on COM3.
Replaces the earlier 2026-09-28 report of the same name, which predated `97eb0bf78`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 280 s capture, `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh HankinSolids profile 280 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:105880, data:170820, headers:8996` / `RAM1: variables:315392, code:34952, padding:30584, free:143360` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2241–2256 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `hk_timeline_step` averages 12.45 ms/f; its worst window is 24.73 ms/f (frames 2241–2256). Peak frame render is **34.81 ms** (frame 2252), and **0/4447** frames spilled. Setup frame 1 is excluded from both; it rendered 14.52 ms.

The previous shipping report (2026-09-28 16:45) recorded peak 🟢 34.837 (18) and spilled 🟢 0/4137 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 19 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2241–2256)

```
frame                     62.77 ms   37.66 Mcyc   100%
  pov_preserve_half       148.0 us    88.8 kcyc     0%
  hk_timeline_step        24.73 ms   14.84 Mcyc    39%
    hk_draw_mesh          23.97 ms   14.38 Mcyc    38%
      hk_mesh_scan        23.93 ms   14.36 Mcyc    38%
        scan_mesh_raster  21.84 ms   13.10 Mcyc    35%  x182  71995 cyc/c
          filter_blend     1.08 ms   646.5 kcyc     2%  x15914  41 cyc/c
        scan_face_setup    1.96 ms    1.18 Mcyc     3%  x182  6466 cyc/c
      hk_mesh_transform    39.1 us    23.5 kcyc     0%
    hk_update_hankin      710.4 us   426.2 kcyc     1%
  canvas_clear             84.2 us    50.5 kcyc     0%
  canvas_buffer_wait      37.80 ms   22.68 Mcyc    60%
```

Wall min/avg/max = 54.37/62.77/70.49 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `hk_timeline_step` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `hk_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| truncatedIcosidodecahedron | truncatedIcosidodecahedron | 🟢 34.81 | 0/125 | 24.73 | 8/8 |
| snubDodecahedron | snubDodecahedron | 🟢 26.59 | 0/125 | 21.35 | 8/8 |
| truncatedDodecahedron | truncatedDodecahedron | 🟢 23.22 | 0/113 | 19.69 | 7/7 |
| truncatedCuboctahedron | truncatedCuboctahedron | 🟢 23.21 | 0/250 | 19.54 | 15/15 |
| rhombicosidodecahedron | rhombicosidodecahedron | 🟢 23.19 | 0/125 | 20.90 | 8/8 |
| snubCube | snubCube | 🟢 23.18 | 0/125 | 17.98 | 8/8 |
| truncatedIcosahedron | truncatedIcosahedron | 🟢 22.57 | 0/113 | 19.74 | 7/7 |
| truncatedCube | truncatedCube | 🟢 22.32 | 0/113 | 18.43 | 7/7 |
| icosidodecahedron | icosidodecahedron | 🟢 21.39 | 0/351 | 18.40 | 22/22 |
| truncatedOctahedron | truncatedOctahedron | 🟢 19.99 | 0/226 | 17.42 | 14/14 |
| rhombicuboctahedron | rhombicuboctahedron | 🟢 19.86 | 0/113 | 17.37 | 7/7 |
| truncatedTetrahedron | truncatedTetrahedron | 🟢 19.66 | 0/226 | 15.96 | 15/15 |
| dodecahedron | dodecahedron | 🟢 18.73 | 0/363 | 16.73 | 22/22 |
| icosahedron | icosahedron | 🟢 18.45 | 0/226 | 15.49 | 14/14 |
| octahedron | octahedron | 🟢 18.03 | 0/452 | 14.53 | 29/29 |
| cuboctahedron | cuboctahedron | 🟢 17.76 | 0/589 | 15.71 | 36/36 |
| cube | cube | 🟢 17.42 | 0/463 | 15.42 | 29/29 |
| tetrahedron | tetrahedron | 🟢 15.01 | 0/238 | 12.71 | 15/15 |
| 0 | — | 🟢 14.52 | 0/111 | 12.90 | 7/7 |

### Per-pixel figures

`filter_blend` ran 15,914 times per frame in the peak window at 41 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1158/frame  min/avg/max 0.6/1.7/11.4 us  cpu 3.08%
isr_pack          145/frame  min/avg/max 6.5/7.2/9.6 us  cpu 1.65%
isr_dma_submit    145/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `hk_timeline_step` — 39% of the peak window, 24.73 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 34.81 (19), spilled 🟢 0/4447 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HankinSolids`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh HankinSolids profile 280 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
