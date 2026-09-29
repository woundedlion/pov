# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HopfFibration`).
Raw capture: `build/prof/hopffibration_ship.log`, captured 2026-09-28 19:02 on COM3.
Replaces `profile_hopffibration_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh HopfFibration profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:63800, data:148516, headers:8868` / `RAM1: variables:315136, code:37512, padding:28024, free:143616` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 193–224 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `hf_render_trails` averages 26.46 ms/f; its worst window is 43.29 ms/f (frames 193–224). Peak frame render is **51.52 ms** (frame 216), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 0.42 ms.

The previous shipping report (2026-08-26 01:30) recorded peak 🟢 48.50 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 193–224)

```
frame                     62.85 ms   37.71 Mcyc   100%
  pov_preserve_half       140.6 us    84.4 kcyc     0%
  hf_render_trails        43.29 ms   25.98 Mcyc    69%
    hf_trail_raster       33.87 ms   20.32 Mcyc    54%  x149  136275 cyc/c
      filter_blend         5.56 ms    3.34 Mcyc     9%  x48901  68 cyc/c
    hf_trail_gate          8.60 ms    5.16 Mcyc    14%  x210  24574 cyc/c
  hf_project_record       184.3 us   110.6 kcyc     0%
  hf_advance_tumble         0.5 us      331 cyc     0%
  hf_timeline_step         17.2 us    10.3 kcyc     0%
  canvas_clear             84.9 us    50.9 kcyc     0%
  canvas_buffer_wait      19.13 ms   11.48 Mcyc    30%
```

Wall min/avg/max = 42.27/62.85/81.88 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 48,901 times per frame in the peak window at 68 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1160/frame  min/avg/max 0.6/1.7/17.5 us  cpu 3.10%
isr_pack          145/frame  min/avg/max 6.3/7.1/9.8 us  cpu 1.63%
isr_dma_submit    145/frame  min/avg/max 0.6/0.9/1.2 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `hf_render_trails` — 69% of the peak window, 43.29 ms/f.
2. `hf_project_record` — 0% of the peak window, 0.18 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 51.52, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HopfFibration`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh HopfFibration profile 70 32` builds, flashes and captures under the device lock.
