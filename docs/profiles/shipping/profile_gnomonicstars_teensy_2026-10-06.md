# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/gnomonicstars_ship.log`, captured 2026-10-06 19:11 on COM3.
Replaces `profile_gnomonicstars_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean main checkout at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GnomonicStars 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh GnomonicStars profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:54844, data:153176, headers:9064` / `RAM1: variables:315072, code:25432, padding:7336, free:176448` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 193–224 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `gn_draw_stars` averages 9.97 ms/f; its worst window is 15.52 ms/f (frames 193–224). Peak frame render is **21.83 ms** (frame 210), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 14.03 ms.

The previous shipping report (2026-09-28 19:18) recorded peak 🟢 21.86 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: a quiet regime near 9.0–9.3 ms/f of `gn_draw_stars`, broken by bursts where it rises to 11.5–15.5 ms/f (frames 161–256, 353–384 and 769–896). The star count is 600 calls per frame in both; the bursts cost more per star.

### Burst (frames 193–224, worst of the capture)

```
frame                    62.56 ms    37.53 Mcyc   100%
  pov_preserve_half      144.0 us     86.4 kcyc     0%
  gn_draw_stars          15.52 ms     9.31 Mcyc    25%
    gn_star_scan         14.80 ms     8.88 Mcyc    24%  x600  25 us/c
      filter_blend       742.2 us    445.3 kcyc     1%  x7261  61 cyc/c
  gn_timeline_step        37.3 us     22.4 kcyc     0%
  canvas_clear            84.2 us     50.5 kcyc     0%
  canvas_buffer_wait     46.77 ms    28.06 Mcyc    75%
```

Wall min/avg/max = 51.53/62.56/73.33 ms. Per-frame values are window averages; `xN` is calls per frame.

### Quiet (frames 641–672)

```
frame                    62.54 ms    37.52 Mcyc   100%
  pov_preserve_half      145.7 us     87.4 kcyc     0%
  gn_draw_stars           9.09 ms     5.45 Mcyc    15%
    gn_star_scan          8.38 ms     5.03 Mcyc    13%  x600  14 us/c
      filter_blend       271.9 us    163.1 kcyc     0%  x2494  65 cyc/c
  gn_timeline_step        38.2 us     22.9 kcyc     0%
  canvas_clear            84.0 us     50.4 kcyc     0%
  canvas_buffer_wait     53.18 ms    31.91 Mcyc    85%
```

Wall min/avg/max = 59.83/62.54/65.20 ms. Blended pixels fall from 7,261 to 2,494 per frame and the per-star scan from 25 to 14 us, so the bursts are star footprint growing on screen, not more stars.

### Per-pixel figures

`filter_blend` ran 7,261 times per frame in the burst window at 61 cycles per blend, and 2,494 times at 65 cycles in the quiet window.

## Column-ISR / DMA marshaling cost

```
isr_wake        1154/frame  min/avg/max 0.5/1.6/19.8 us  cpu 3.02%
isr_pack         144/frame  min/avg/max 6.3/6.8/9.5 us  cpu 1.56%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

Burst window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `gn_draw_stars` — 25% of the burst window, 15.52 ms/f.
2. `pov_preserve_half` — 0% of the burst window, 0.14 ms/f.
3. `canvas_clear` — 0% of the burst window, 0.08 ms/f.
4. `gn_timeline_step` — 0% of the burst window, 0.04 ms/f.

README cells: peak 🟢 21.83, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78` on COM4, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh GnomonicStars profile 70 32` builds, flashes and captures under the device lock.
