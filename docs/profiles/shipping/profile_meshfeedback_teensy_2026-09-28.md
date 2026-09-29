# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MeshFeedback`).
Raw capture: `build/prof/meshfeedback_ship.log`, captured 2026-09-28 18:36 on COM3.
Replaces the earlier 2026-09-28 report of the same name, which predated `97eb0bf78`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16` |

Image size (`profile` env, this effect only): `FLASH: code:116688, data:190948, headers:8780` / `RAM1: variables:315168, code:47240, padding:18296, free:143584` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 3889–3904 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `mf_feedback_flush` averages 43.76 ms/f; its worst window is 50.49 ms/f (frames 3889–3904). Peak frame render is **61.87 ms** (frame 4100), and **0/6687** frames spilled. Setup frame 1 is excluded from both; it rendered 6.46 ms.

The previous shipping report (2026-09-28 14:46) recorded peak 🟢 61.85 (12) and spilled 🟢 0/6688 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 4097–4112)

```
frame                     62.28 ms   37.37 Mcyc   100%
  mf_mesh_draw             6.64 ms    3.99 Mcyc    11%
    filter_blend           1.12 ms   670.5 kcyc     2%  x10396  64 cyc/c
  mf_feedback_flush       46.12 ms   27.67 Mcyc    74%
    feedback_composite    38.24 ms   22.94 Mcyc    61%
    feedback_populate      7.85 ms    4.71 Mcyc    13%
    feedback_litscan        2.6 us     1.6 kcyc     0%
  mf_timeline_step         44.3 us    26.6 kcyc     0%
  mf_apply_params          10.1 us     6.1 kcyc     0%
  canvas_clear            409.6 us   245.8 kcyc     1%
  canvas_buffer_wait       9.05 ms    5.43 Mcyc    15%
```

Wall min/avg/max = 54.69/62.28/64.77 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `mf_feedback_flush` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `mf_feedback_flush` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 6 | — | 🟢 61.87 | 0/482 | 46.12 | 31/31 |
| 5 | — | 🟢 60.96 | 0/482 | 50.49 | 30/30 |
| 3 | — | 🟢 58.90 | 0/723 | 48.52 | 45/45 |
| 11 | — | 🟢 57.99 | 0/482 | 40.67 | 30/30 |
| 9 | — | 🟢 57.76 | 0/482 | 46.47 | 30/30 |
| 12 | — | 🟢 57.24 | 0/482 | 40.38 | 30/30 |
| 4 | — | 🟢 55.84 | 0/664 | 48.12 | 41/41 |
| 2 | — | 🟢 55.02 | 0/723 | 47.00 | 46/46 |
| 1 | — | 🟢 54.59 | 0/721 | 43.86 | 45/45 |
| 8 | — | 🟢 53.79 | 0/482 | 46.53 | 30/30 |
| 10 | — | 🟢 53.01 | 0/482 | 45.35 | 30/30 |
| 7 | — | 🟢 51.55 | 0/482 | 44.57 | 30/30 |

### Per-pixel figures

`filter_blend` ran 10,396 times per frame in the peak window at 64 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1149/frame  min/avg/max 0.6/1.7/19.9 us  cpu 3.18%
isr_pack          144/frame  min/avg/max 6.3/7.2/9.7 us  cpu 1.65%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `mf_feedback_flush` — 74% of the peak window, 46.12 ms/f.
2. `mf_mesh_draw` — 11% of the peak window, 6.64 ms/f.
3. `canvas_clear` — 1% of the peak window, 0.41 ms/f.
4. `mf_timeline_step` — 0% of the peak window, 0.04 ms/f.

README cells: peak 🟢 61.87 (12), spilled 🟢 0/6687 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MeshFeedback profile 420 16` builds, flashes and captures under the device lock.
