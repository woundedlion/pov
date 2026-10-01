# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`).
Raw capture: `build/prof/mindsplatter_ship.log`, captured 2026-09-29 14:48 on COM3.
Replaces `profile_mindsplatter_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `d5ca81403` (landed as `c40c1def8`) |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile 110 16` |

Image size (`profile` env, this effect only): `FLASH: code:66000, data:548756, headers:8860` / `RAM1: variables:315200, code:37128, padding:28408, free:143552` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 993–1008 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `msp_draw_particles` averages 21.53 ms/f; its worst window is 47.77 ms/f (frames 993–1008). Peak frame render is **54.94 ms** (frame 1000), and **0/1727** frames spilled. Setup frame 1 is excluded from both; it rendered 0.28 ms.

The previous shipping report (2026-09-28 18:51) recorded peak 🟢 56.11 (8) and spilled 🟢 0/1727 (0.00%). A same-board baseline of the parent commit `a2c6d2be0` (COM3, 2026-09-29 14:45) recorded peak 56.13 ms; the float min/max sweep in `c40c1def8` cut every entry by 1.2–2.1%.

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 8 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 993–1008)

```
frame                        62.14 ms 37.29 Mcyc  100%
  pov_preserve_half          137.6 us  82.6 kcyc    0%
  msp_draw_particles         47.77 ms 28.66 Mcyc   77%
    msp_particle_scan        47.76 ms 28.66 Mcyc   77%
      plot_ps_raster         34.16 ms 20.49 Mcyc   55%  x576  35580 cyc/c
      plot_ps_deferred       580.9 us 348.6 kcyc    1%  x576  605 cyc/c
      plot_ps_gate            7.27 ms  4.36 Mcyc   12%  x1537  2838 cyc/c
        plot_ps_cartesian_gate  1.06 ms 638.4 kcyc  2%  x1537  415 cyc/c
      plot_ps_tween           4.93 ms  2.96 Mcyc    8%  x1537  1925 cyc/c
  msp_particle_step           2.95 ms  1.77 Mcyc    5%
  msp_timeline_step           54.8 us  32.9 kcyc    0%
  canvas_clear                85.1 us  51.1 kcyc    0%
  canvas_buffer_wait         11.15 ms  6.69 Mcyc   18%
```

Wall min/avg/max = 59.07/62.14/64.80 ms. Per-frame values are window averages; `xN` is calls per frame. Against the 2026-09-28 capture of the same window, `plot_ps_raster` fell from 35.12 to 34.16 ms/f and `plot_ps_gate` from 7.36 to 7.27 ms/f; the rest is unchanged.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `msp_draw_particles` is the costliest modal-call-count window of each entry. Baseline is the same-board `a2c6d2be0` capture.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `msp_draw_particles` ms/f | Clean windows | Baseline peak ms |
|---|---|--:|--:|--:|--:|--:|
| 7 | — | 🟢 54.94 | 0/159 | 47.77 | 10/10 | 56.13 |
| 6 | — | 🟢 51.88 | 0/159 | 41.05 | 10/10 | 52.76 |
| 8 | — | 🟢 44.85 | 0/159 | 30.23 | 10/10 | 45.96 |
| 3 | — | 🟢 37.07 | 0/297 | 27.66 | 18/18 | 37.93 |
| 2 | — | 🟢 35.51 | 0/318 | 30.27 | 20/20 | 36.20 |
| 1 | — | 🟢 34.45 | 0/317 | 28.03 | 20/20 | 35.14 |
| 4 | — | 🟢 33.82 | 0/159 | 27.56 | 10/10 | 34.50 |
| 5 | — | 🟢 22.36 | 0/159 | 17.45 | 10/10 | 22.77 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1147/frame  min/avg/max 0.6/1.7/12.2 us  cpu 3.18%
isr_pack          143/frame  min/avg/max 6.3/7.2/9.7 us  cpu 1.66%
isr_dma_submit    143/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `msp_draw_particles` — 77% of the peak window, 47.77 ms/f; `plot_ps_raster` is 55% of the frame on its own.
2. `msp_particle_step` — 5% of the peak window, 2.95 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 54.94 (8), spilled 🟢 0/1727 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: this effect's `HS_O3` regions are unchanged; `c40c1def8` changed only the Plot stroke and trail float extremes.
- The captured tip `d5ca81403` differs from the landed `c40c1def8` only by an earlier HyperLattice commit, which is not in this single-effect image.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MindSplatter profile 110 16` builds, flashes and captures under the device lock.
