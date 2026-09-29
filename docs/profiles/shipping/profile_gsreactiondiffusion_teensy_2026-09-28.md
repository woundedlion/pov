# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).
Raw capture: `build/prof/gsreactiondiffusion_ship.log`, captured 2026-09-28 19:21 on COM4.
Replaces `profile_gsreactiondiffusion_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile 130 32` |

Image size (`profile` env, this effect only): `FLASH: code:60928, data:242112, headers:8256` / `RAM1: variables:315168, code:23880, padding:8888, free:176352` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 449–480 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `grd_render` averages 48.57 ms/f; its worst window is 53.71 ms/f (frames 449–480). Peak frame render is **55.53 ms** (frame 489), and **0/2047** frames spilled. Setup frame 1 is excluded from both; it rendered 34.63 ms.

The previous shipping report (2026-08-26 01:28) recorded peak 🟢 55.26 and spilled 🟢 0/2048 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 481–512)

```
frame                     62.40 ms   37.44 Mcyc   100%
  pov_preserve_half       135.6 us    81.4 kcyc     0%
  grd_render              53.67 ms   32.20 Mcyc    86%
    grd_rasterize         41.10 ms   24.66 Mcyc    66%
      grd_shader_draw     37.43 ms   22.46 Mcyc    60%
      grd_cull_flags       3.30 ms    1.98 Mcyc     5%
      grd_orient          373.6 us   224.1 kcyc     1%
    grd_simulate          12.05 ms    7.23 Mcyc    19%
  rd_timeline_step         24.8 us    14.9 kcyc     0%
  canvas_clear             85.0 us    51.0 kcyc     0%
  canvas_buffer_wait       8.48 ms    5.09 Mcyc    14%
```

Wall min/avg/max = 60.24/62.40/64.33 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1151/frame  min/avg/max 0.6/1.8/18.6 us  cpu 3.25%
isr_pack          144/frame  min/avg/max 6.4/7.5/10.6 us  cpu 1.72%
isr_dma_submit    144/frame  min/avg/max 0.6/0.9/3.4 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `grd_render` — 86% of the peak window, 53.67 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.09 ms/f.
4. `rd_timeline_step` — 0% of the peak window, 0.02 ms/f.

README cells: peak 🟢 55.53, spilled 🟢 0/2047 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh GSReactionDiffusion profile 130 32` builds, flashes and captures under the device lock.
