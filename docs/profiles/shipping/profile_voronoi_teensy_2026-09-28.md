# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Voronoi`).
Raw capture: `build/prof/voronoi_ship.log`, captured 2026-09-28 19:15 on COM3.
Replaces `profile_voronoi_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Voronoi 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Voronoi profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:38072, data:146800, headers:8664` / `RAM1: variables:314976, code:13384, padding:19384, free:176544` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `vo_shade` averages 7.64 ms/f; its worst window is 7.72 ms/f (frames 353–384). Peak frame render is **8.51 ms** (frame 32), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 15.48 ms.

The previous shipping report (2026-08-26 01:43) recorded peak 🟢 8.96 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 1–32)

```
frame                     59.30 ms   35.58 Mcyc   100%
  pov_preserve_half       146.5 us    87.9 kcyc     0%
  vo_shade                 7.88 ms    4.73 Mcyc    13%
  vo_kdtree               337.7 us   202.6 kcyc     1%
  vo_animate               51.9 us    31.1 kcyc     0%
  canvas_clear             88.4 us    53.0 kcyc     0%
  canvas_buffer_wait      50.79 ms   30.48 Mcyc    86%
```

Wall min/avg/max = 8.26/59.30/62.76 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1373/frame  min/avg/max 0.5/1.4/14.1 us  cpu 2.54%
isr_pack          136/frame  min/avg/max 6.5/6.9/9.7 us  cpu 1.26%
isr_dma_submit    136/frame  min/avg/max 0.6/0.9/1.4 us  cpu 0.16%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `vo_shade` — 13% of the peak window, 7.88 ms/f.
2. `vo_kdtree` — 1% of the peak window, 0.34 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 8.51, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Voronoi`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Voronoi profile 70 32` builds, flashes and captures under the device lock.
