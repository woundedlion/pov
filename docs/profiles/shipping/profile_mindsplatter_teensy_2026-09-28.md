# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`).
Raw capture: `build/prof/mindsplatter_ship.log`, captured 2026-09-28 18:51 on COM3.
Replaces `profile_mindsplatter_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile 110 16` |

Image size (`profile` env, this effect only): `FLASH: code:66832, data:548756, headers:9052` / `RAM1: variables:315200, code:37960, padding:27576, free:143552` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 993–1008 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `msp_draw_particles` averages 21.98 ms/f; its worst window is 48.84 ms/f (frames 993–1008). Peak frame render is **56.11 ms** (frame 1000), and **0/1727** frames spilled. Setup frame 1 is excluded from both; it rendered 0.28 ms.

The previous shipping report (2026-08-26 07:40) recorded peak 🟢 52.77 (9) and spilled 🟢 0/1728 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 8 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 993–1008)

```
frame                     62.14 ms   37.29 Mcyc   100%
  pov_preserve_half       135.3 us    81.2 kcyc     0%
  msp_draw_particles      48.84 ms   29.31 Mcyc    79%
    msp_particle_scan     48.84 ms   29.30 Mcyc    79%
      plot_ps_raster      35.12 ms   21.07 Mcyc    57%  x576  36582 cyc/c
      plot_ps_deferred    579.9 us   348.0 kcyc     1%  x576  604 cyc/c
      plot_ps_gate         7.36 ms    4.42 Mcyc    12%  x1537  2873 cyc/c
        plot_ps_cartesian_gate   1.06 ms   635.7 kcyc     2%  x1537  414 cyc
      plot_ps_tween        4.95 ms    2.97 Mcyc     8%  x1537  1932 cyc/c
  msp_particle_step        2.95 ms    1.77 Mcyc     5%
  msp_timeline_step        54.9 us    32.9 kcyc     0%
  canvas_clear             84.7 us    50.8 kcyc     0%
  canvas_buffer_wait      10.07 ms    6.04 Mcyc    16%
```

Wall min/avg/max = 58.98/62.14/64.85 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `msp_draw_particles` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `msp_draw_particles` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 7 | — | 🟢 56.11 | 0/159 | 48.84 | 10/10 |
| 6 | — | 🟢 52.76 | 0/159 | 41.86 | 10/10 |
| 8 | — | 🟢 45.94 | 0/159 | 30.92 | 10/10 |
| 3 | — | 🟢 37.97 | 0/297 | 28.27 | 18/18 |
| 2 | — | 🟢 36.21 | 0/318 | 30.93 | 20/20 |
| 1 | — | 🟢 35.16 | 0/317 | 28.65 | 20/20 |
| 4 | — | 🟢 34.47 | 0/159 | 28.18 | 10/10 |
| 5 | — | 🟢 22.79 | 0/159 | 17.85 | 10/10 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1147/frame  min/avg/max 0.6/1.7/12.0 us  cpu 3.18%
isr_pack          143/frame  min/avg/max 6.3/7.2/9.5 us  cpu 1.66%
isr_dma_submit    143/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `msp_draw_particles` — 79% of the peak window, 48.84 ms/f.
2. `msp_particle_step` — 5% of the peak window, 2.95 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 56.11 (8), spilled 🟢 0/1727 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MindSplatter profile 110 16` builds, flashes and captures under the device lock.
