# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/mindsplatter_ship.log`, captured 2026-10-06 18:32 on COM4.
Replaces `profile_mindsplatter_teensy_2026-09-29.md`; the architecture snapshot [profile_mindsplatter_architecture_teensy_2026-10-01.md](profile_mindsplatter_architecture_teensy_2026-10-01.md) is kept.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile 110 16` |

Image size (`profile` env, this effect only): `FLASH: code:70508, data:551440, headers:8832` / `RAM1: variables:315200, code:37960, padding:27576, free:143552` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 993–1008 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `msp_draw_particles` averages 21.89 ms/f; its worst window is 49.27 ms/f (frames 993–1008). Peak frame render is **56.60 ms** (frame 1000), and **0/1727** frames spilled. Setup frame 1 is excluded from both; it rendered 0.28 ms.

The most recent prior capture, the architecture snapshot (2026-10-01 09:15), recorded peak 🟢 54.457 (8) and spilled 🟢 0/1736 (0.00%). The previous un-suffixed shipping report (2026-09-29 14:48) recorded peak 🟢 54.94 (8) and spilled 🟢 0/1727 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 8 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 8, wraps back to entry 1 and ends on entry 3. The block below is the window holding the pass's peak frame.

### Peak window (frames 993–1008)

```
frame                          62.13 ms  37.28 Mcyc   100%
  pov_preserve_half            136.7 us   82.0 kcyc     0%
  msp_draw_particles           49.27 ms  29.56 Mcyc    79%
    msp_particle_scan          49.27 ms  29.56 Mcyc    79%
      plot_ps_raster           36.41 ms  21.85 Mcyc    59%  x576  37932 cyc/c
      plot_ps_deferred         559.9 us  336.0 kcyc     1%  x576  583 cyc/c
      plot_ps_gate              6.57 ms   3.94 Mcyc    11%  x1537  2564 cyc/c
        plot_ps_cartesian_gate 720.9 us  432.5 kcyc     1%  x1537  281 cyc/c
      plot_ps_tween             4.93 ms   2.96 Mcyc     8%  x1537  1923 cyc/c
  msp_particle_step             2.97 ms   1.78 Mcyc     5%
  msp_timeline_step             60.3 us   36.2 kcyc     0%
  canvas_clear                  85.2 us   51.2 kcyc     0%
  canvas_buffer_wait            9.60 ms   5.76 Mcyc    15%
```

Wall min/avg/max = 58.94/62.13/64.88 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `msp_draw_particles` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `msp_draw_particles` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 7 | — | 🟢 56.60 | 0/159 | 49.27 | 9/10 |
| 6 | — | 🟢 53.09 | 0/159 | 41.64 | 9/10 |
| 8 | — | 🟢 45.37 | 0/159 | 30.78 | 9/10 |
| 3 | — | 🟢 37.57 | 0/297 | 28.43 | 17/18 |
| 2 | — | 🟢 36.51 | 0/318 | 31.13 | 18/20 |
| 1 | — | 🟢 35.73 | 0/317 | 28.92 | 18/20 |
| 4 | — | 🟢 34.78 | 0/159 | 28.29 | 9/10 |
| 5 | — | 🟢 22.99 | 0/159 | 17.37 | 9/10 |

Entry 7 sets the peak, 5.9 ms inside the 62.5 ms window.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1146/frame  min/avg/max 0.6/1.7/12.1 us  cpu 3.18%
isr_pack         143/frame  min/avg/max 6.3/7.2/9.5 us  cpu 1.66%
isr_dma_submit   143/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `msp_draw_particles` — 79% of the peak window, 49.27 ms/f; `plot_ps_raster` is 59% of the frame on its own.
2. `msp_particle_step` — 5% of the peak window, 2.97 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 56.60 (8), spilled 🟢 0/1727 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the architecture snapshot ran at `ac4870adb` and the previous un-suffixed report at `d5ca81403` (landed as `c40c1def8`), and the deltas (peak +2.14 ms and +1.66 ms) are not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MindSplatter profile 110 16` builds, flashes and captures under the device lock.
