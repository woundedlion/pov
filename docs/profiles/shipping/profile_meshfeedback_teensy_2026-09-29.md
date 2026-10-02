# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/meshfeedback_ship.log`, captured 2026-09-29 10:03 on COM4.
Replaces the 2026-09-28 23:43 report, which predated the polar-row ITCM trims `05d6cfeb6` and `58f966be9`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `58f966be9` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16` |

Image size (shipping `phantasm` image): `FLASH: code:512264, data:735628, headers:8556   free for files:775168 / RAM1: variables:314784, code:168856, padding:27752   free for local variables:12896 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: window frames 513–528 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.6 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `mf_feedback_flush` averages 33.77 ms/f; its worst window is 37.65 ms/f (frames 513–528). Peak frame render is **48.61 ms** (frames 5537–5552 window), and **0/6688** frames spilled. Setup frame 1 is excluded from both; it rendered 6.60 ms.

The previous shipping report (2026-09-28 23:43, COM3, before the polar-row ITCM trims) recorded peak 🟢 48.09 (12) and spilled 🟢 0/6688 (0.00%); before the flush work the effect peaked at 61.87.

A display window is 62.5 ms, so render at or under it holds 16 fps. The feedback pipeline renders the full 288×144 canvas (`crosses_segments`). The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 5537–5552)

```
frame                      62.62 ms   37.57 Mcyc  100%
  mf_mesh_draw              9.40 ms    5.64 Mcyc   15%
    filter_blend            1.85 ms    1.11 Mcyc    3%  x14199  78 cyc/c
  mf_feedback_flush        30.26 ms   18.16 Mcyc   48%
    feedback_composite     29.71 ms   17.82 Mcyc   47%
    feedback_populate      480.4 us   288.3 kcyc    1%
    feedback_litscan         2.0 us     1.2 kcyc    0%
  mf_timeline_step          53.8 us    32.3 kcyc    0%
  mf_apply_params           34.3 us    20.6 kcyc    0%
  canvas_clear             406.9 us   244.2 kcyc    1%
  canvas_buffer_wait       22.46 ms   13.47 Mcyc   36%
```

Wall min/avg/max = 54.21/62.62/72.75 ms. Per-frame values are window averages; `xN` is calls per frame. The composite now chroma-scales out-of-gamut trail pixels onto the gamut grid instead of bisecting, blends its taps with Q15 integer weights, skips the colour transform for all-black pairs, and composites the rows under latitude sine ≈ 0.5 at half column resolution; populate reads its lattice origins from the cache built in `init_storage`. Rows within about 30 degrees of a pole reconstruct the warp from cap-plane offsets in `composite_polar_row`, which accounts for the rise from 45.11 ms: the polar rows now land on their true warp targets. Since `58f966be9` those rows run one lane per step, trading about 0.6 ms/frame for 1.5 KB of ITCM.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `mf_feedback_flush` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `mf_feedback_flush` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 12 | — | 🟢 48.61 | 0/482 | 33.35 | 30/30 |
| 11 | — | 🟢 46.10 | 0/482 | 30.26 | 30/30 |
| 6 | — | 🟢 45.73 | 0/482 | 36.75 | 31/31 |
| 9 | — | 🟢 45.46 | 0/482 | 37.28 | 30/30 |
| 3 | — | 🟢 45.33 | 0/723 | 37.65 | 45/45 |
| 5 | — | 🟢 44.86 | 0/482 | 37.55 | 30/30 |
| 8 | — | 🟢 44.43 | 0/482 | 36.42 | 30/30 |
| 1 | — | 🟢 44.23 | 0/722 | 33.56 | 45/45 |
| 10 | — | 🟢 43.86 | 0/482 | 30.54 | 30/30 |
| 2 | — | 🟢 41.52 | 0/723 | 36.87 | 46/46 |
| 7 | — | 🟢 41.47 | 0/482 | 34.57 | 30/30 |
| 4 | — | 🟢 40.88 | 0/664 | 37.17 | 41/41 |

### Per-pixel figures

`filter_blend` ran 14199 times per frame in the peak window at 78 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake          1155/frame  min/avg/max 0.5/1.7/11.6 us  cpu 3.16%
isr_pack           144/frame  min/avg/max 6.4/7.4/9.9 us  cpu 1.69%
isr_dma_submit     144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `mf_feedback_flush` — 48% of the peak window, 30.26 ms/f (composite 29.71, populate 0.48).
2. `mf_mesh_draw` — 15% of the peak window, 9.40 ms/f.
3. `canvas_clear` — 1% of the peak window, 0.41 ms/f.
4. `mf_timeline_step` — 0% of the peak window, 0.05 ms/f.

README cells: peak 🟢 48.61 (12), spilled 🟢 0/6688 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the flush (`HS_O3_BEGIN` region in pixel_feedback.h) and its colour helpers are `HS_O3_FN`; `lms_cbrt_scale_to_gamut_lut` joins them.
- The gamut clip scales chroma to a guarded sampled cell minimum, with no certified deficit bound against the exact boundary. The rows under latitude sine ≈ 0.5 carry a one-column blur that stays under one row pitch in angle (Pole Half-Res slider, default 1).
- The image was built from source `58f966be9` (provenance file).

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MeshFeedback profile 420 16` builds, flashes and captures under the device lock.
