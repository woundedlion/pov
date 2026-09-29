# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Raymarch`).
Raw capture: `build/prof/raymarch_ship.log`, captured 2026-09-28 19:09 on COM3.
Replaces `profile_raymarch_teensy_2026-09-20.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Raymarch 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Raymarch profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:104712, data:187604, headers:8740` / `RAM1: variables:315008, code:30440, padding:2328, free:176512` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.6 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `rm_shader_draw` averages 45.21 ms/f; its worst window is 46.76 ms/f (frames 833–864). Peak frame render is **54.87 ms** (frame 8), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 99.83 ms.

The previous shipping report (2026-09-20 22:51) recorded peak 🟢 56.07 and spilled 🟢 0/1736 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 1–32)

```
frame                     63.08 ms   37.85 Mcyc   100%
  pov_preserve_half       128.1 us    76.9 kcyc     0%
  rm_shader_draw          47.31 ms   28.39 Mcyc    75%
    filter_blend          481.6 us   289.0 kcyc     1%  x5692  51 cyc/c
  rm_timeline_step        375.6 us   225.4 kcyc     1%
  canvas_clear             89.1 us    53.5 kcyc     0%
  canvas_buffer_wait      12.06 ms    7.23 Mcyc    19%
```

Wall min/avg/max = 53.37/63.08/99.83 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 5,692 times per frame in the peak window at 51 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1394/frame  min/avg/max 0.5/1.4/13.9 us  cpu 2.65%
isr_pack          138/frame  min/avg/max 6.2/6.8/9.4 us  cpu 1.25%
isr_dma_submit    138/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.17%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `rm_shader_draw` — 75% of the peak window, 47.31 ms/f.
2. `rm_timeline_step` — 1% of the peak window, 0.38 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.13 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 54.87, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Raymarch`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Raymarch profile 70 32` builds, flashes and captures under the device lock.
