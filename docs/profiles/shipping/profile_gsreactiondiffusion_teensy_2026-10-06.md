# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/gsreactiondiffusion_ship.log`, captured 2026-10-06 19:15 on COM3.
Replaces `profile_gsreactiondiffusion_teensy_2026-10-05.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean main checkout at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, no extra flags |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile 130 32` |

Image size (`profile` env, this effect only): `FLASH: code:70580, data:345028, headers:8324` / `RAM1: variables:315328, code:25704, padding:7064, free:176192` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 385–416 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `grd_render` averages 32.64 ms/f; its worst window is 38.64 ms/f (frames 385–416). Peak frame render is **39.43 ms** (frame 395), and **0/2047** frames spilled. Setup frame 1 is excluded from both; it rendered 21.83 ms.

The previous shipping report (2026-10-05 12:20) recorded peak 🟢 39.371 and spilled 🟢 0/2032 (0.00%). Its peak fell on the same frame 395, and its worst `grd_render` window was the same frames 385–416 at 38.607 ms/f; this capture reads 39.43 ms and 38.64 ms/f, within 0.2% on a different board.

A display window is 62.5 ms, so render at or under it holds 16 fps; the peak leaves 23.07 ms of margin, and 13.57 ms below the 53 ms render ceiling the previous report set. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: the reaction grows from the seed, densifies until `grd_render` crests near 38.3–38.6 ms/f, then a runtime reseed drops it back to 22–26 ms/f and the growth repeats. Reseeds land in the windows carrying `grd_seed_reaction` (frames 449–480, 865–896, 1281–1312, 1793–1824), the same windows as the previous capture. Physics (5.12 ms/f) and colour noise (2.70 ms/f) are flat throughout; the regimes differ in `grd_shader_draw` and, less, `grd_pigment`.

| Crest window | `grd_render` ms/f | Peak frame render ms |
|---|--:|--:|
| 385–416 | 38.64 | 39.43 |
| 801–832 | 38.47 | 39.14 |
| 1217–1248 | 38.27 | 38.96 |
| 1729–1760 | 38.39 | 39.02 |

| Reseed window | `grd_render` ms/f | Peak frame render ms |
|---|--:|--:|
| 449–480 | 23.48 | 30.11 |
| 865–896 | 21.97 | 27.27 |
| 1281–1312 | 25.57 | 33.81 |
| 1793–1824 | 22.95 | 30.70 |

### Dense crest (frames 385–416, worst of the capture)

```
frame                     62.42 ms    37.45 Mcyc   100%
  pov_preserve_half       138.2 us     82.9 kcyc     0%
  grd_render              38.64 ms    23.18 Mcyc    62%
    grd_rasterize         25.19 ms    15.11 Mcyc    40%
      grd_shader_draw     22.98 ms    13.79 Mcyc    37%
      grd_cull_flags      335.0 us    201.0 kcyc     1%
      grd_orient           1.87 ms     1.12 Mcyc     3%
    grd_simulate          10.19 ms     6.12 Mcyc    16%
      grd_physics          5.12 ms     3.07 Mcyc     8%  x6.0  853 us/c
      grd_pigment          4.32 ms     2.59 Mcyc     7%
    grd_color_noise        2.72 ms     1.63 Mcyc     4%
  rd_timeline_step         27.2 us     16.3 kcyc     0%
  canvas_clear             84.2 us     50.5 kcyc     0%
  canvas_buffer_wait      23.53 ms    14.12 Mcyc    38%
```

Wall min/avg/max = 61.45/62.41/63.36 ms. Per-frame values are window averages; `xN` is calls per frame. The shader draw is 59% of the render here, with six physics substeps per frame behind it.

### Reseed (frames 449–480)

```
frame                     62.27 ms    37.36 Mcyc   100%
  pov_preserve_half       143.8 us     86.3 kcyc     0%
  grd_render              23.48 ms    14.09 Mcyc    38%
    grd_rasterize         10.74 ms     6.44 Mcyc    17%  x0.97  11087 us/c
      grd_shader_draw      8.01 ms     4.81 Mcyc    13%  x0.97  8268 us/c
      grd_cull_flags      914.3 us    548.6 kcyc     1%  x0.97  944 us/c
      grd_orient           1.82 ms     1.09 Mcyc     3%  x0.97  1875 us/c
    grd_simulate           8.81 ms     5.29 Mcyc    14%  x0.97  9093 us/c
      grd_physics          4.96 ms     2.98 Mcyc     8%  x5.8  853 us/c
      grd_pigment          3.12 ms     1.87 Mcyc     5%  x0.97  3225 us/c
    grd_color_noise        2.70 ms     1.62 Mcyc     4%
  rd_timeline_step         32.4 us     19.4 kcyc     0%
  canvas_clear             84.4 us     50.6 kcyc     0%
  canvas_buffer_wait      38.53 ms    23.12 Mcyc    62%
grd_color_palette         472.0 us    283.2 kcyc     1%  x0.19  2517 us/c
grd_seed_reaction         546.5 us    327.9 kcyc     1%  x0.06  8744 us/c
```

Wall min/avg/max = 59.97/62.27/67.06 ms. One frame of 32 skips simulate and rasterize while `grd_seed_reaction` (2 calls, 8.74 ms each) and `grd_color_palette` (6 calls) run; those two carry the parser's MIXED-PARENT tag and are attribution diagnostics, not exclusive siblings. The fresh pattern is sparse, so the shader draw falls to 8.01 ms/f and cull flags rise to 0.91 ms/f.

### Regrowth (frames 1825–1856)

```
frame                     62.61 ms    37.56 Mcyc   100%
  pov_preserve_half       141.8 us     85.1 kcyc     0%
  grd_render              22.45 ms    13.47 Mcyc    36%
    grd_rasterize         10.29 ms     6.17 Mcyc    16%
      grd_shader_draw      7.56 ms     4.53 Mcyc    12%
      grd_cull_flags      857.8 us    514.7 kcyc     1%
      grd_orient           1.87 ms     1.12 Mcyc     3%
    grd_simulate           9.45 ms     5.67 Mcyc    15%
      grd_physics          5.12 ms     3.07 Mcyc     8%  x6.0  853 us/c
      grd_pigment          3.58 ms     2.15 Mcyc     6%
    grd_color_noise        2.70 ms     1.62 Mcyc     4%
  rd_timeline_step         25.3 us     15.2 kcyc     0%
  canvas_clear             84.0 us     50.4 kcyc     0%
  canvas_buffer_wait      39.91 ms    23.94 Mcyc    64%
```

Wall min/avg/max = 59.92/62.60/65.07 ms. The first full window after the last reseed is the cheapest of the capture: shader draw is a third of its crest cost, while simulation is within 0.74 ms/f of the crest.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.7/13.5 us  cpu 3.16%
isr_pack         144/frame  min/avg/max 6.6/7.3/10.2 us  cpu 1.68%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/1.1 us  cpu 0.21%
```

Dense crest window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `grd_render` — 62% of the crest window, 38.64 ms/f: `grd_shader_draw` 22.98, `grd_physics` 5.12, `grd_pigment` 4.32, `grd_color_noise` 2.72, `grd_orient` 1.87, `grd_cull_flags` 0.34.
2. `pov_preserve_half` — 0% of the crest window, 0.14 ms/f.
3. `canvas_clear` — 0% of the crest window, 0.08 ms/f.
4. `rd_timeline_step` — 0% of the crest window, 0.03 ms/f.

README cells: peak 🟢 39.43, spilled 🟢 0/2047 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `grd_seed_reaction` and `grd_color_palette` are MIXED-PARENT; their totals are not exclusive costs.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- The 130 s capture ran without an epoch override; `validate` reports 0 epoch resets across all 2048 frames. The previous capture used `-D HS_PROFILE_EPOCH_REVS=1200` and covered 2,032 runtime frames.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`, a descendant of the previous report's `0f3f33b76`; the previous report ran on COM4, this one on COM3.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh GSReactionDiffusion profile 130 32` builds, flashes and captures under the device lock.
