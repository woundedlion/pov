# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/meshfeedback_ship.log`, captured 2026-10-06 18:18 on COM3.
Replaces `profile_meshfeedback_teensy_2026-09-29.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16` |

Image size (`profile` env, this effect only): `FLASH: code:128236, data:193948, headers:8564` / `RAM1: variables:315168, code:48696, padding:16840, free:143584` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 513–528 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.6 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `mf_feedback_flush` averages 33.58 ms/f; its worst window is 37.51 ms/f (frames 513–528). Peak frame render is **48.53 ms** (frame 5543), and **0/6687** frames spilled. Setup frame 1 is excluded from both; it rendered 6.50 ms.

The previous shipping report (2026-09-29 10:03, COM4) recorded peak 🟢 48.61 (12) and spilled 🟢 0/6688 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The feedback pipeline renders the full 288×144 canvas (`crosses_segments`). The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 style presets on a single fixed solid, so there is no per-solid table; each preset owns its hold and the transition that follows it. The capture runs entry 1 through 12 twice, wraps back to entry 1 twice and ends on entry 4. The block below is the window holding the pass's peak frame.

### Peak window (frames 5537–5552)

```
frame                    62.61 ms  37.57 Mcyc   100%
  mf_mesh_draw            9.30 ms   5.58 Mcyc    15%
    filter_blend          1.89 ms   1.14 Mcyc     3%  x14198  80 cyc/c
  mf_feedback_flush      30.06 ms  18.04 Mcyc    48%
    feedback_composite   29.49 ms  17.69 Mcyc    47%
    feedback_populate    488.9 us  293.3 kcyc     1%
    feedback_litscan      18.3 us   11.0 kcyc     0%
  mf_timeline_step        54.6 us   32.8 kcyc     0%
  mf_apply_params         41.3 us   24.8 kcyc     0%
  canvas_clear           406.7 us  244.0 kcyc     1%
  canvas_buffer_wait     22.74 ms  13.64 Mcyc    36%
```

Wall min/avg/max = 54.00/62.61/72.97 ms. Per-frame values are window averages; `xN` is calls per frame. The composite is 98% of the flush in this window; populate is under half a millisecond.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `mf_feedback_flush` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `mf_feedback_flush` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 12 | — | 🟢 48.53 | 0/482 | 29.50 | 28/30 |
| 11 | — | 🟢 45.97 | 0/482 | 29.53 | 28/30 |
| 6 | — | 🟢 45.63 | 0/482 | 36.63 | 29/31 |
| 9 | — | 🟢 45.21 | 0/482 | 37.14 | 28/30 |
| 3 | — | 🟢 45.11 | 0/723 | 37.51 | 42/45 |
| 5 | — | 🟢 44.71 | 0/482 | 37.43 | 29/30 |
| 8 | — | 🟢 44.06 | 0/482 | 36.18 | 28/30 |
| 1 | — | 🟢 44.03 | 0/721 | 32.78 | 43/45 |
| 10 | — | 🟢 43.74 | 0/482 | 30.38 | 28/30 |
| 2 | — | 🟢 41.31 | 0/723 | 36.41 | 43/46 |
| 7 | — | 🟢 41.21 | 0/482 | 33.12 | 28/30 |
| 4 | — | 🟢 40.71 | 0/664 | 37.03 | 39/41 |

### Per-pixel figures

`filter_blend` ran 14,198 times per frame in the peak window at 80 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1155/frame  min/avg/max 0.5/1.7/17.9 us  cpu 3.17%
isr_pack         144/frame  min/avg/max 6.4/7.4/10.1 us  cpu 1.69%
isr_dma_submit   144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `mf_feedback_flush` — 48% of the peak window, 30.06 ms/f (composite 29.49, populate 0.49).
2. `mf_mesh_draw` — 15% of the peak window, 9.30 ms/f.
3. `canvas_clear` — 1% of the peak window, 0.41 ms/f.
4. `mf_timeline_step` — 0% of the peak window, 0.05 ms/f.

README cells: peak 🟢 48.53 (12), spilled 🟢 0/6687 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: the flush (`HS_O3_BEGIN` region in pixel_feedback.h) and its colour helpers are `HS_O3_FN`. Shipping `profile` image at `e2f5b0a3d`; the previous report ran at `58f966be9`, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MeshFeedback profile 420 16` builds, flashes and captures under the device lock.
