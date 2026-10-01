# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile BZReactionDiffusion`).
Raw capture: `build/prof/bzreactiondiffusion_ship.log`, captured 2026-09-28 19:12 on COM4.
Replaces `profile_bzreactiondiffusion_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | BZReactionDiffusion 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture |
| Reproduce | `bash tools/profile_one.sh BZReactionDiffusion profile 130 32` |

Image size (`profile` env, this effect only): `FLASH: code:61448, data:241936, headers:8936` / `RAM1: variables:315072, code:23688, padding:9080, free:176448` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1761–1792 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.9 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `bz_render` averages 46.81 ms/f; its worst window is 47.78 ms/f (frames 1761–1792). Peak frame render is **48.84 ms** (frame 1762), and **0/2047** frames spilled. Setup frame 1 is excluded from both; it rendered 77.71 ms.

The previous shipping report (2026-08-26 01:20) recorded peak 🟢 48.91 and spilled 🟢 0/2048 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 1761–1792)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       134.3 us    80.6 kcyc     0%
  bz_render               47.78 ms   28.67 Mcyc    77%
    bz_raster             42.89 ms   25.73 Mcyc    69%
    bz_orient             374.9 us   225.0 kcyc     1%
    bz_physics             4.52 ms    2.71 Mcyc     7%
  rd_timeline_step         29.3 us    17.6 kcyc     0%
  canvas_clear             84.8 us    50.9 kcyc     0%
  canvas_buffer_wait      14.41 ms    8.65 Mcyc    23%
```

Wall min/avg/max = 60.95/62.44/64.01 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.6/1.7/18.2 us  cpu 3.14%
isr_pack          144/frame  min/avg/max 6.3/7.0/10.2 us  cpu 1.61%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/5.2 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `bz_render` — 77% of the peak window, 47.78 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.13 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `rd_timeline_step` — 0% of the peak window, 0.03 ms/f.

README cells: peak 🟢 48.84, spilled 🟢 0/2047 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=BZReactionDiffusion`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh BZReactionDiffusion profile 130 32` builds, flashes and captures under the device lock.
