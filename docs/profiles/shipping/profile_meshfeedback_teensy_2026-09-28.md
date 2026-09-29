# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MeshFeedback`).
Raw capture: `build/prof/meshfeedback_ship.log`, captured 2026-09-28 23:43 on COM3.
Replaces the 2026-09-28 21:42 report, which predated the polar cap-plane warp reconstruction `a6eb6e3fa`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `a6eb6e3fa` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16` |

Image size (`profile` env, this effect only): `FLASH: code:124328, data:191036, headers:8220   free for files:1708032 / RAM1: variables:315168, code:49528, padding:16008   free for local variables:143584 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: window frames 513–528 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `mf_feedback_flush` averages 33.21 ms/f; its worst window is 37.01 ms/f (frames 513–528). Peak frame render is **48.09 ms** (frames 5537–5552 window), and **0/6688** frames spilled. Setup frame 1 is excluded from both; it rendered 6.70 ms.

The previous shipping report (2026-09-28 21:42) recorded peak 🟢 45.11 (12) and spilled 🟢 0/6688 (0.00%); before the flush work the effect peaked at 61.87.

A display window is 62.5 ms, so render at or under it holds 16 fps. The feedback pipeline renders the full 288×144 canvas (`crosses_segments`). The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 5537–5552)

```
frame                      62.62 ms   37.57 Mcyc  100%
  mf_mesh_draw              9.52 ms    5.71 Mcyc   15%
    filter_blend            1.90 ms    1.14 Mcyc    3%  x14199  80 cyc/c
  mf_feedback_flush        29.61 ms   17.77 Mcyc   47%
    feedback_composite     29.06 ms   17.43 Mcyc   46%
    feedback_populate      481.4 us   288.8 kcyc    1%
    feedback_litscan         1.2 us      747 cyc    0%
  mf_timeline_step          44.8 us    26.9 kcyc    0%
  mf_apply_params           35.1 us    21.1 kcyc    0%
  canvas_clear             406.9 us   244.1 kcyc    1%
  canvas_buffer_wait       23.00 ms   13.80 Mcyc   37%
```

Wall min/avg/max = 54.20/62.62/72.77 ms. Per-frame values are window averages; `xN` is calls per frame. The composite now chroma-scales out-of-gamut trail pixels onto the gamut grid instead of bisecting, blends its taps with Q15 integer weights, skips the colour transform for all-black pairs, and composites the rows under latitude sine ≈ 0.5 at half column resolution; populate reads its lattice origins from the cache built in `init_storage`. Rows within about 30 degrees of a pole reconstruct the warp from cap-plane offsets in `composite_polar_row`, which accounts for the rise from 45.11 ms: the polar rows now land on their true warp targets.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `mf_feedback_flush` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `mf_feedback_flush` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 12 | — | 🟢 48.09 | 0/482 | 32.72 | 30/30 |
| 11 | — | 🟢 45.58 | 0/482 | 29.61 | 30/30 |
| 6 | — | 🟢 45.24 | 0/482 | 36.15 | 31/31 |
| 9 | — | 🟢 44.91 | 0/482 | 36.63 | 30/30 |
| 3 | — | 🟢 44.77 | 0/723 | 37.01 | 45/45 |
| 5 | — | 🟢 44.32 | 0/482 | 36.92 | 30/30 |
| 8 | — | 🟢 43.91 | 0/482 | 35.85 | 30/30 |
| 1 | — | 🟢 43.66 | 0/722 | 33.11 | 45/45 |
| 10 | — | 🟢 43.31 | 0/482 | 29.93 | 30/30 |
| 2 | — | 🟢 40.96 | 0/723 | 36.26 | 46/46 |
| 7 | — | 🟢 40.94 | 0/482 | 34.08 | 30/30 |
| 4 | — | 🟢 40.29 | 0/664 | 36.56 | 41/41 |

### Per-pixel figures

`filter_blend` ran 14199 times per frame in the peak window at 80 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake          1155/frame  min/avg/max 0.6/1.7/18.4 us  cpu 3.21%
isr_pack           144/frame  min/avg/max 6.4/7.3/9.6 us  cpu 1.68%
isr_dma_submit     144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `mf_feedback_flush` — 47% of the peak window, 29.61 ms/f (composite 29.06, populate 0.48).
2. `mf_mesh_draw` — 15% of the peak window, 9.52 ms/f.
3. `canvas_clear` — 1% of the peak window, 0.41 ms/f.
4. `mf_timeline_step` — 0% of the peak window, 0.04 ms/f.

README cells: peak 🟢 48.09 (12), spilled 🟢 0/6688 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the flush (`HS_O3_BEGIN` region in pixel_feedback.h) and its colour helpers are `HS_O3_FN`; `lms_cbrt_scale_to_gamut_lut` joins them.
- The gamut clip lands at most one grid cell's chroma inside the exact boundary, and the rows under latitude sine ≈ 0.5 carry a one-column blur that stays under one row pitch in angle (Pole Half-Res slider, default 1).
- The image was built from source `a6eb6e3fa` (provenance file); master later gained `f14c71a22`, a HyperLattice-only change.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MeshFeedback profile 420 16` builds, flashes and captures under the device lock.
