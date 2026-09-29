# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MeshFeedback`).
Raw capture: `build/prof/meshfeedback_ship.log`, captured 2026-09-28 21:42 on COM3.
Replaces the 2026-09-28 18:36 report, which predated the flush changes landed in `b2dc8ecac..5bcd697fd`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `5bcd697fd` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16` |

Image size (`profile` env, this effect only): `FLASH: code:115968, data:190964, headers:8460   free for files:1716224 / RAM1: variables:315168, code:45272, padding:20264   free for local variables:143584 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: window frames 513–528 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `mf_feedback_flush` averages 30.67 ms/f; its worst window is 34.31 ms/f (frames 513–528). Peak frame render is **45.11 ms** (frames 5537–5552 window), and **0/6688** frames spilled. Setup frame 1 is excluded from both; it rendered 6.48 ms.

The previous shipping report (2026-09-28 18:36) recorded peak 🟢 61.87 (12) and spilled 🟢 0/6687 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The feedback pipeline renders the full 288×144 canvas (`crosses_segments`). The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 5537–5552)

```
frame                      62.62 ms   37.57 Mcyc  100%
  mf_mesh_draw              9.25 ms    5.55 Mcyc   15%
    filter_blend            1.55 ms   927.2 kcyc    2%  x14199  65 cyc/c
  mf_feedback_flush        27.62 ms   16.57 Mcyc   44%
    feedback_composite     27.15 ms   16.29 Mcyc   43%
    feedback_populate      436.4 us   261.9 kcyc    1%
    feedback_litscan         3.6 us     2.2 kcyc    0%
  mf_timeline_step          45.6 us    27.4 kcyc    0%
  mf_apply_params           34.1 us    20.5 kcyc    0%
  canvas_clear             406.9 us   244.2 kcyc    1%
  canvas_buffer_wait       25.26 ms   15.15 Mcyc   40%
```

Wall min/avg/max = 54.98/62.62/71.96 ms. Per-frame values are window averages; `xN` is calls per frame. The composite now chroma-scales out-of-gamut trail pixels onto the gamut grid instead of bisecting, blends its taps with Q15 integer weights, skips the colour transform for all-black pairs, and composites the rows under latitude sine ≈ 0.5 at half column resolution; populate reads its lattice origins from the cache built in `init_storage`.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `mf_feedback_flush` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `mf_feedback_flush` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 12 | — | 🟢 45.11 | 0/482 | 30.40 | 30/30 |
| 11 | — | 🟢 42.68 | 0/482 | 27.62 | 30/30 |
| 6 | — | 🟢 42.26 | 0/482 | 33.46 | 31/31 |
| 9 | — | 🟢 41.99 | 0/482 | 33.97 | 30/30 |
| 3 | — | 🟢 41.85 | 0/723 | 34.31 | 45/45 |
| 5 | — | 🟢 41.40 | 0/482 | 34.15 | 30/30 |
| 8 | — | 🟢 40.99 | 0/482 | 33.11 | 30/30 |
| 1 | — | 🟢 40.79 | 0/722 | 30.38 | 45/45 |
| 10 | — | 🟢 40.51 | 0/482 | 27.90 | 30/30 |
| 7 | — | 🟢 38.14 | 0/482 | 31.37 | 30/30 |
| 2 | — | 🟢 37.98 | 0/723 | 33.58 | 46/46 |
| 4 | — | 🟢 37.40 | 0/664 | 33.69 | 41/41 |

### Per-pixel figures

`filter_blend` ran 14199 times per frame in the peak window at 65 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake          1155/frame  min/avg/max 0.6/1.7/14.3 us  cpu 3.20%
isr_pack           144/frame  min/avg/max 6.4/7.4/9.9 us  cpu 1.69%
isr_dma_submit     144/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `mf_feedback_flush` — 44% of the peak window, 27.62 ms/f (composite 27.15, populate 0.44).
2. `mf_mesh_draw` — 15% of the peak window, 9.25 ms/f.
3. `canvas_clear` — 1% of the peak window, 0.41 ms/f.
4. `mf_timeline_step` — 0% of the peak window, 0.05 ms/f.

README cells: peak 🟢 45.11 (12), spilled 🟢 0/6688 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the flush (`HS_O3_BEGIN` region in pixel_feedback.h) and its colour helpers are `HS_O3_FN`; `lms_cbrt_scale_to_gamut_lut` joins them.
- The gamut clip lands at most one grid cell's chroma inside the exact boundary, and the rows under latitude sine ≈ 0.5 carry a one-column blur that stays under one row pitch in angle.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MeshFeedback profile 420 16` builds, flashes and captures under the device lock.
