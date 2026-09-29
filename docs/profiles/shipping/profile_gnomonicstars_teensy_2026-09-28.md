# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GnomonicStars`).
Raw capture: `build/prof/gnomonicstars_ship.log`, captured 2026-09-28 19:18 on COM4.
Replaces `profile_gnomonicstars_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GnomonicStars 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh GnomonicStars profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:52752, data:151532, headers:8708` / `RAM1: variables:315072, code:25016, padding:7752, free:176448` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 193–224 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `gn_draw_stars` averages 10.00 ms/f; its worst window is 15.54 ms/f (frames 193–224). Peak frame render is **21.86 ms** (frame 210), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 14.01 ms.

The previous shipping report (2026-08-26 01:25) recorded peak 🟢 29.64 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 193–224)

```
frame                     62.55 ms   37.53 Mcyc   100%
  pov_preserve_half       143.8 us    86.3 kcyc     0%
  gn_draw_stars           15.54 ms    9.32 Mcyc    25%
    gn_star_scan          14.82 ms    8.89 Mcyc    24%  x600  14823 cyc/c
      filter_blend        738.2 us   442.9 kcyc     1%  x7261  61 cyc/c
  gn_timeline_step         40.5 us    24.3 kcyc     0%
  canvas_clear             84.5 us    50.7 kcyc     0%
  canvas_buffer_wait      46.75 ms   28.05 Mcyc    75%
```

Wall min/avg/max = 51.56/62.55/73.35 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 7,261 times per frame in the peak window at 61 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1154/frame  min/avg/max 0.5/1.7/12.3 us  cpu 3.07%
isr_pack          144/frame  min/avg/max 6.2/6.7/9.4 us  cpu 1.55%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `gn_draw_stars` — 25% of the peak window, 15.54 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `gn_timeline_step` — 0% of the peak window, 0.04 ms/f.

README cells: peak 🟢 21.86, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh GnomonicStars profile 70 32` builds, flashes and captures under the device lock.
