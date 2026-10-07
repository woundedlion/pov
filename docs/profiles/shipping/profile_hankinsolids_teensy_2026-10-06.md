# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/hankinsolids_ship.log`, captured 2026-10-06 19:12 on COM4.
Replaces `profile_hankinsolids_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 280 s capture, `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh HankinSolids profile 280 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:107100, data:172208, headers:8432` / `RAM1: variables:315392, code:35416, padding:30120, free:143360` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2241–2256 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

Method note: the sweep script's budget (210 s, `HS_PROFILE_EPOCH_REVS=1920`) did not close the ordered cycle (28 shape loads, no return to the first), so this re-run uses the previous report's 280 s and 2400 revs.

## Frame cadence

**Pass aggregate**: `hk_timeline_step` averages 12.38 ms/f; its worst window is 24.67 ms/f (frames 2241–2256). Peak frame render is **34.63 ms** (frame 2252), and **0/4447** frames spilled. Setup frame 1 is excluded from both; it rendered 14.49 ms.

The previous shipping report (2026-09-28 19:27) recorded peak 🟢 34.81 (19) and spilled 🟢 0/4447 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 19 shape entries in ordered-cycle order, each announced by a `Loading shape` line after the unannounced initial shape; each owns its hold and the transition that follows it. The capture logs 38 shape loads: the full 30-load cycle, then a return to `truncatedTetrahedron` (the first announced shape) and 7 more. The block below is the window holding the pass's peak frame.

### Peak window (frames 2241–2256)

```
frame                     62.76 ms   37.66 Mcyc   100%
  pov_preserve_half       145.0 us    87.0 kcyc     0%
  hk_timeline_step        24.67 ms   14.80 Mcyc    39%
    hk_draw_mesh          23.95 ms   14.37 Mcyc    38%
      hk_mesh_scan        23.92 ms   14.35 Mcyc    38%
        scan_mesh_raster  21.86 ms   13.11 Mcyc    35%  x182  72052 cyc/c
          filter_blend     1.08 ms   649.9 kcyc     2%  x15914  41 cyc/c
        scan_face_setup    1.94 ms    1.16 Mcyc     3%  x182  6398 cyc/c
      hk_mesh_transform    31.7 us    19.0 kcyc     0%
    hk_update_hankin      666.5 us   399.9 kcyc     1%
  canvas_clear             84.0 us    50.4 kcyc     0%
  canvas_buffer_wait      37.86 ms   22.72 Mcyc    60%
```

Wall min/avg/max = 54.49/62.76/70.30 ms. Per-frame values are window averages; `xN` is calls per frame. The window is `truncatedIcosidodecahedron`'s hold.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `hk_timeline_step` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). Entry `0` is the initial shape, which no `Loading shape` line names. The cycle wraps back to the first announced shape within the capture (validate: shape markers advance and return to the first).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `hk_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| truncatedIcosidodecahedron | truncatedIcosidodecahedron | 🟢 34.63 | 0/125 | 24.67 | 7/8 |
| snubDodecahedron | snubDodecahedron | 🟢 26.46 | 0/125 | 21.30 | 7/8 |
| truncatedDodecahedron | truncatedDodecahedron | 🟢 23.18 | 0/113 | 19.67 | 6/7 |
| truncatedCuboctahedron | truncatedCuboctahedron | 🟢 23.14 | 0/250 | 19.48 | 13/15 |
| rhombicosidodecahedron | rhombicosidodecahedron | 🟢 23.09 | 0/125 | 20.81 | 7/8 |
| snubCube | snubCube | 🟢 23.05 | 0/125 | 17.95 | 7/8 |
| truncatedIcosahedron | truncatedIcosahedron | 🟢 22.52 | 0/113 | 19.70 | 6/7 |
| truncatedCube | truncatedCube | 🟢 22.25 | 0/113 | 18.38 | 6/7 |
| icosidodecahedron | icosidodecahedron | 🟢 21.30 | 0/351 | 18.35 | 19/22 |
| truncatedOctahedron | truncatedOctahedron | 🟢 19.94 | 0/226 | 17.38 | 12/14 |
| rhombicuboctahedron | rhombicuboctahedron | 🟢 19.80 | 0/113 | 17.29 | 7/7 |
| truncatedTetrahedron | truncatedTetrahedron | 🟢 19.61 | 0/226 | 15.93 | 13/15 |
| dodecahedron | dodecahedron | 🟢 18.71 | 0/363 | 16.41 | 19/22 |
| icosahedron | icosahedron | 🟢 18.45 | 0/226 | 15.47 | 13/14 |
| octahedron | octahedron | 🟢 17.96 | 0/452 | 14.45 | 25/29 |
| cuboctahedron | cuboctahedron | 🟢 17.71 | 0/589 | 15.64 | 33/36 |
| cube | cube | 🟢 17.33 | 0/463 | 15.36 | 26/29 |
| tetrahedron | tetrahedron | 🟢 14.94 | 0/238 | 12.66 | 13/15 |
| 0 | — | 🟢 14.49 | 0/111 | 12.85 | 7/7 |

Entry `0`'s 14.49 ms is frame 64 (14.488 ms); setup frame 1 (14.487 ms) is out of its bucket. Every entry is within 0.18 ms of the previous report's peak, and the frame counts match it entry for entry.

### Per-pixel figures

`filter_blend` ran 15,914 times per frame in the peak window at 41 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1158/frame  min/avg/max 0.6/1.7/11.5 us  cpu 3.13%
isr_pack          145/frame  min/avg/max 6.5/7.2/9.8 us  cpu 1.65%
isr_dma_submit    145/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `hk_timeline_step` — 39% of the peak window, 24.67 ms/f (`hk_draw_mesh` 23.95).
2. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 34.63 (19), spilled 🟢 0/4447 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- `hk_draw_mesh`, `hk_mesh_scan`, `hk_mesh_transform`, `scan_face_setup` and `scan_mesh_raster` are not exclusive costs (validate INFO: MIXED-PARENT/DUPLICATE-NAME).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78` on COM3, and the delta between them (peak −0.18 ms) is not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HankinSolids`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh HankinSolids profile 280 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
