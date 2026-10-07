# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/bzreactiondiffusion_ship.log`, captured 2026-10-06 19:05 on COM3.
Replaces `profile_bzreactiondiffusion_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean main checkout at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | BZReactionDiffusion 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, no extra flags |
| Reproduce | `bash tools/profile_one.sh BZReactionDiffusion profile 130 32` |

Image size (`profile` env, this effect only): `FLASH: code:62948, data:243288, headers:9152` / `RAM1: variables:315072, code:24200, padding:8568, free:176448` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1761–1792 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.4 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `bz_render` averages 46.29 ms/f; its worst window is 47.30 ms/f (frames 1761–1792). Peak frame render is **48.34 ms** (frame 1476), and **0/2047** frames spilled. Setup frame 1 is excluded from both; it rendered 77.20 ms, over the 62.5 ms window, and is the one spilled frame the parser footer counts (1/2048).

The previous shipping report (2026-09-28 19:12) recorded peak 🟢 48.84 and spilled 🟢 0/2047 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps; every frame after setup holds 16 fps with at least 14.16 ms of margin. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one regime with a slow oscillation. `bz_render` swings between 45.12 and 47.30 ms/f with a period of roughly 450–600 frames; `bz_physics` stays at 4.51 ms/f, and the swing is entirely `bz_raster`. The blocks below are the crest window holding the pass's peak frame and the trough window.

### Crest (frames 1473–1504, holds the peak frame)

```
frame                    62.44 ms    37.46 Mcyc   100%
  pov_preserve_half      135.4 us     81.2 kcyc     0%
  bz_render              47.04 ms    28.22 Mcyc    75%
    bz_raster            42.15 ms    25.29 Mcyc    68%
    bz_orient            374.5 us    224.7 kcyc     1%
    bz_physics            4.51 ms     2.71 Mcyc     7%
  rd_timeline_step        30.4 us     18.2 kcyc     0%
  canvas_clear            84.4 us     50.7 kcyc     0%
  canvas_buffer_wait     15.15 ms     9.09 Mcyc    24%
```

Wall min/avg/max = 60.87/62.44/64.05 ms. Per-frame values are window averages; `xN` is calls per frame.

### Trough (frames 481–512)

```
frame                    62.46 ms    37.48 Mcyc   100%
  pov_preserve_half      133.6 us     80.2 kcyc     0%
  bz_render              45.12 ms    27.07 Mcyc    72%
    bz_raster            40.23 ms    24.14 Mcyc    64%
    bz_orient            374.4 us    224.6 kcyc     1%
    bz_physics            4.51 ms     2.71 Mcyc     7%
  rd_timeline_step        22.0 us     13.2 kcyc     0%
  canvas_clear            84.6 us     50.8 kcyc     0%
  canvas_buffer_wait     17.10 ms    10.26 Mcyc    27%
```

Wall min/avg/max = 61.95/62.46/63.00 ms. The raster is 1.92 ms/f cheaper than at the crest; physics and orientation are identical, so the swing follows the pattern's raster cost, not the simulation.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/18.7 us  cpu 3.12%
isr_pack         144/frame  min/avg/max 6.3/7.0/10.2 us  cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.8/0.9/2.2 us  cpu 0.21%
```

Crest window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `bz_render` — 75% of the crest window, 47.04 ms/f (`bz_raster` 42.15, `bz_physics` 4.51, `bz_orient` 0.37).
2. `pov_preserve_half` — 0% of the crest window, 0.14 ms/f.
3. `canvas_clear` — 0% of the crest window, 0.08 ms/f.
4. `rd_timeline_step` — 0% of the crest window, 0.03 ms/f.

README cells: peak 🟢 48.34, spilled 🟢 0/2047 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it. Frame 1 overran the display window in the previous report too (77.71 ms).
- The 130 s capture ran without an epoch override; `validate` reports 0 epoch resets across all 2048 frames.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here. The previous report ran on COM4.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=BZReactionDiffusion`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh BZReactionDiffusion profile 130 32` builds, flashes and captures under the device lock.
