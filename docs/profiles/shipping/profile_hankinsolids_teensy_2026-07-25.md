# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HankinSolids`).
Raw capture: `build/prof/hankinsolids_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The per-face mesh scan (`scan_mesh_raster`/`raster_scan`) and its `classify_faces_impl` are the hot path; `ConwayMorph` drives the transitions. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, epoch stretched to 1920 revs, ordered cycle (`HS_PROFILE_ORDERED_CYCLE`); this pass walked 19 solids. |
| Reproduce | `bash tools/profile_one.sh HankinSolids profile 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size: `FLASH: code:97804, data:548976, headers:8572` / `RAM1:
variables:315104, code:49160, padding:16376, free:143648` / `RAM2:
variables:519552, free:4736`. (Re-measured 2026-07-25: the hankin parallel-gate
landing grew this from the 2026-07-20 figure.)

Exactness cross-check: window frames 1153–1168 root counter cyc ÷ 600 MHz matches
the measured wall sum within **1.6 ppm**.

## Frame cadence

**Pass aggregate** (`buckets`): all 19 walked solids green — worst peak frame
render **45.22 ms** (truncatedIcosidodecahedron, on a Conway-morph frame),
spilled **0/3328 frames (0%)** — 🟢. The held solids are cheaper (~28 ms for the
62-face truncatedIcosidodecahedron); the peak lives on the Conway sweeps that
morph one solid into the next, which draw two meshes and can carry more faces
than either endpoint.

## Phase-by-phase readout

Phase schedule: HankinSolids loads a solid, holds it under Hankin choreography,
then Conway-morphs to the next. The costliest regime is the morph transition; the
window below straddles one (note the two `hk_draw_mesh` calls at 0.5/frame).

### Conway-morph window (frames 2337–2352, truncatedIcosidodecahedron leg)

```
frame                 63.03 ms  37.82 Mcyc  100%
  hk_timeline_step    21.38 ms  12.83 Mcyc   33%
    hk_draw_mesh (×2)  20.95 ms  12.57 Mcyc   —    two meshes (morph straddle)
      scan_mesh_raster 20.37 ms  12.22 Mcyc   —    x62.0 faces
        raster_scan    20.28 ms  12.17 Mcyc   99%
          raster_shade  7.90 ms   4.74 Mcyc   38%  x12496
            filter_blend 2.36 ms  1.41 Mcyc   29%  x12496  113 cyc/blend
        scan_face_setup 0.57 ms   0.34 Mcyc    4%  x62.0
    hk_conway_op        0.06 ms   0.03 Mcyc    0%
  hk_buffer_wait       41.65 ms  24.99 Mcyc   66%
```

Wall min/avg/max = 60.69/63.03/73.69 ms. `scan_mesh_raster` runs once per face
(x62 = the 62 faces of truncatedIcosidodecahedron); `raster_shade`/`filter_blend`
cover 12,496 blended px = 1.20× the quadrant. `scan_face_setup` (0.57 ms) splits
across `face_bounds` (39%), `face_project`, `face_azimuth`/`face_thetas` — all
sub-ms. `hk_conway_op`/`hk_conway_compile` (the morph itself) are ~0.08 ms.

### Per-preset table

Cycle wrap confirmed (a shape repeats; 19 distinct solids this pass). Worst three
by bucket peak (all green):

| Solid | F | peak ms | spilled | fps |
|---|--:|--:|--:|--:|
| truncatedIcosidodecahedron | 62 | 🟢 45.22 | 0/125 (0%) | 16 |
| truncatedOctahedron | 14 | 🟢 28.85 | 0/113 (0%) | 16 |
| truncatedTetrahedron | 8 | 🟢 25.70 | 0/113 (0%) | 16 |

### Per-pixel figures

`filter_blend` = 12,496 blended px/frame ⇒ **1.20× coverage**; 113 cyc/blend.
`raster_scan` 12.17 Mcyc ÷ 62 faces = 196 kcyc per face on this leg.

## Column-ISR / DMA marshaling cost

```
isr_wake        1163/frame  min/avg/max 0.67/1.77/11.46 us  cpu 3.26%
isr_pack         145/frame  min/avg/max 6.11/6.59/9.38 us  cpu 1.51%
isr_dma_submit   145/frame  min/avg/max 0.61/0.93/1.09 us  cpu 0.21%
```

(window frames 2337–2352.) Total ISR CPU share **4.98%**.

## Summary ranking

1. `raster_scan` — the per-face mesh scan, essentially the whole render.
2. `raster_shade` / `filter_blend` — 7.9 / 2.4 ms, the blended shade.
3. `scan_face_setup` — 0.57 ms, per-face classification (`classify_faces_impl`).

HankinSolids holds a hard 16 fps across all 19 walked solids; the Conway morphs
(the real peaks) stay under 62.5 ms. The global-`-O3` ceiling buys the worst leg
45.22 → 45.37 ms — essentially a tie, because `raster_scan` and
`classify_faces_impl` are already in HS_O3 regions.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `raster_shade`; its calls ≈ blended px.
- Selective -O3: the `-O3` mesh-scan driver + `classify_faces_impl always_inline`
  cover the hot path (dropping the latter crosses phantasm's ITCM granule cliff).
- The peak window straddles a Conway morph (two meshes drawn); held-solid frames
  are cheaper. `scan_mesh_raster` call count = F (62 faces here).
- Sizes re-measured 2026-07-25 (hankin parallel-gate changed the image since
  07-20). Dwell knobs stretch dwell, not per-frame cost. Working tree tip
  `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HankinSolids`,
`HS_PROFILE_WINDOW=16`; `just profile HankinSolids` = build + flash + capture.
