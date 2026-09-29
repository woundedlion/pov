# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingSpin`).
Raw capture: `build/prof/ringspin_ship.log`, captured 2026-09-28 19:13 on COM3.
Replaces `profile_ringspin_teensy_2026-09-24.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingSpin profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:54288, data:149908, headers:8796` / `RAM1: variables:315040, code:29496, padding:3272, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 545–576 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.0 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `rs_draw_rings` averages 30.97 ms/f; its worst window is 38.06 ms/f (frames 545–576). Peak frame render is **47.04 ms** (frame 981), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 1.82 ms.

The previous shipping report (2026-09-24 20:15) recorded peak 🟢 49.920 and spilled 🟢 0/1087 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 961–992)

```
frame                     62.39 ms   37.44 Mcyc   100%
  pov_preserve_half       142.7 us    85.6 kcyc     0%
  rs_draw_rings           34.67 ms   20.80 Mcyc    56%
    rs_ring_scan          34.07 ms   20.44 Mcyc    55%  x76.0  448 us/c
      filter_blend         4.46 ms    2.68 Mcyc     7%  x45185  59 cyc/c
  rs_timeline_step         74.3 us    44.6 kcyc     0%
  canvas_clear             84.4 us    50.7 kcyc     0%
  canvas_buffer_wait      27.43 ms   16.46 Mcyc    44%
```

Wall min/avg/max = 40.88/62.39/84.14 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 45,185 times per frame in the peak window at 59 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1151/frame  min/avg/max 0.6/1.7/11.9 us  cpu 3.06%
isr_pack          144/frame  min/avg/max 6.2/6.8/9.4 us  cpu 1.57%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/1.1 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `rs_draw_rings` — 56% of the peak window, 34.67 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `rs_timeline_step` — 0% of the peak window, 0.07 ms/f.

README cells: peak 🟢 47.04, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh RingSpin profile 70 32` builds, flashes and captures under the device lock.
