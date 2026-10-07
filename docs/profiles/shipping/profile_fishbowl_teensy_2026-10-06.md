# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/fishbowl_ship.log`, captured 2026-10-06 19:07 on COM3.
Replaces `profile_fishbowl_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean main checkout at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Fishbowl 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Fishbowl profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:78428, data:155156, headers:9100` / `RAM1: variables:315072, code:43192, padding:22344, free:143680` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 481–512 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fish_build_vertices` averages 12.44 ms/f; its worst window is 13.15 ms/f (frames 481–512). Peak frame render is **23.78 ms** (frame 203), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 0.70 ms.

The previous shipping report (2026-09-28 19:14) recorded peak 🟢 23.98 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: a ramp over frames 1–128 while `fish_build_vertices` climbs from 1.92 to 12.52 ms/f, then one steady regime at 13.10–13.15 ms/f for the rest of the capture.

### Steady (frames 193–224, holds the peak frame)

```
frame                     62.47 ms    37.48 Mcyc   100%
  pov_preserve_half       141.5 us     84.9 kcyc     0%
  fish_multiline_draw      9.20 ms     5.52 Mcyc    15%
    filter_blend          877.5 us    526.5 kcyc     1%  x7658  69 cyc/c
  fish_build_vertices     13.10 ms     7.86 Mcyc    21%
  fish_noise_prepare        0.4 us      235 cyc      0%
  fish_timeline_step      119.2 us     71.5 kcyc     0%
  canvas_clear             84.6 us     50.7 kcyc     0%
  canvas_buffer_wait      39.82 ms    23.89 Mcyc    64%
```

Wall min/avg/max = 60.31/62.47/64.59 ms. Per-frame values are window averages; `xN` is calls per frame.

### Ramp (frames 33–64)

```
frame                     62.64 ms    37.58 Mcyc   100%
  pov_preserve_half       145.7 us     87.4 kcyc     0%
  fish_multiline_draw      3.65 ms     2.19 Mcyc     6%
    filter_blend          331.7 us    199.0 kcyc     1%  x2882  69 cyc/c
  fish_build_vertices      5.56 ms     3.34 Mcyc     9%
  fish_noise_prepare        0.3 us      182 cyc      0%
  fish_timeline_step      118.6 us     71.2 kcyc     0%
  canvas_clear             84.1 us     50.5 kcyc     0%
  canvas_buffer_wait      53.08 ms    31.85 Mcyc    85%
```

Wall min/avg/max = 61.67/62.64/63.73 ms. Vertex build and the line draw both scale with the growing line set: blends rise from 2,882 to 7,658 per frame at the same 69 cycles each.

### Per-pixel figures

`filter_blend` ran 7,658 times per frame in the steady peak window at 69 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.6/11.6 us  cpu 3.00%
isr_pack         144/frame  min/avg/max 6.2/6.7/9.2 us  cpu 1.53%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/1.1 us  cpu 0.21%
```

Steady peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fish_build_vertices` — 21% of the peak window, 13.10 ms/f.
2. `fish_multiline_draw` — 15% of the peak window, 9.20 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `fish_timeline_step` — 0% of the peak window, 0.12 ms/f.

README cells: peak 🟢 23.78, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78` on COM4, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Fishbowl`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Fishbowl profile 70 32` builds, flashes and captures under the device lock.
