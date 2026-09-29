# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Fishbowl`).
Raw capture: `build/prof/fishbowl_ship.log`, captured 2026-09-28 19:14 on COM4.
Replaces `profile_fishbowl_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Fishbowl 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Fishbowl profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:73840, data:152828, headers:8852` / `RAM1: variables:315072, code:42632, padding:22904, free:143680` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 481–512 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.2 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fish_build_vertices` averages 12.78 ms/f; its worst window is 13.17 ms/f (frames 481–512). Peak frame render is **23.98 ms** (frame 201), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 0.71 ms.

The previous shipping report (2026-08-26 01:22) recorded peak 🟢 23.22 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 193–224)

```
frame                     62.47 ms   37.48 Mcyc   100%
  pov_preserve_half       140.6 us    84.3 kcyc     0%
  fish_multiline_draw      9.35 ms    5.61 Mcyc    15%
    filter_blend          879.7 us   527.8 kcyc     1%  x7658  69 cyc/c
  fish_build_vertices     13.13 ms    7.88 Mcyc    21%
  fish_noise_prepare        0.2 us      159 cyc     0%
  fish_timeline_step      132.7 us    79.6 kcyc     0%
  canvas_clear             84.4 us    50.7 kcyc     0%
  canvas_buffer_wait      39.63 ms   23.78 Mcyc    63%
```

Wall min/avg/max = 60.30/62.47/64.61 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 7,658 times per frame in the peak window at 69 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.6/12.2 us  cpu 2.96%
isr_pack          144/frame  min/avg/max 6.2/6.7/9.3 us  cpu 1.53%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/1.7 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fish_build_vertices` — 21% of the peak window, 13.13 ms/f.
2. `fish_multiline_draw` — 15% of the peak window, 9.35 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `fish_timeline_step` — 0% of the peak window, 0.13 ms/f.

README cells: peak 🟢 23.98, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Fishbowl`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Fishbowl profile 70 32` builds, flashes and captures under the device lock.
