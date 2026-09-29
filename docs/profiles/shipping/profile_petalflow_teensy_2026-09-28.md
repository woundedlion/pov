# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile PetalFlow`).
Raw capture: `build/prof/petalflow_ship.log`, captured 2026-09-28 19:07 on COM3.
Replaces `profile_petalflow_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh PetalFlow profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:60000, data:148988, headers:9124` / `RAM1: variables:315072, code:32792, padding:32744, free:143680` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 929–960 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.6 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `pf_draw_rings` averages 12.05 ms/f; its worst window is 12.57 ms/f (frames 929–960). Peak frame render is **13.61 ms** (frame 946), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 15.97 ms.

The previous shipping report (2026-08-26 01:35) recorded peak 🟢 11.85 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 929–960)

```
frame                     62.48 ms   37.49 Mcyc   100%
  pov_preserve_half       143.2 us    85.9 kcyc     0%
  pf_draw_rings           12.57 ms    7.54 Mcyc    20%
    pf_ring_scan          11.45 ms    6.87 Mcyc    18%  x23.3  492 us/c
      filter_blend        798.3 us   479.0 kcyc     1%  x6956  69 cyc/c
    pf_ring_build          1.00 ms   601.7 kcyc     2%  x23.3  43 us/c
  pf_timeline_step         23.5 us    14.1 kcyc     0%
  canvas_clear             84.1 us    50.4 kcyc     0%
  canvas_buffer_wait      49.65 ms   29.79 Mcyc    79%
```

Wall min/avg/max = 61.26/62.48/63.67 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 6,956 times per frame in the peak window at 69 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.6/1.6/20.8 us  cpu 3.02%
isr_pack          144/frame  min/avg/max 6.2/6.6/9.4 us  cpu 1.52%
isr_dma_submit    144/frame  min/avg/max 0.6/0.9/10.5 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `pf_draw_rings` — 20% of the peak window, 12.57 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `pf_timeline_step` — 0% of the peak window, 0.02 ms/f.

README cells: peak 🟢 13.61, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=PetalFlow`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh PetalFlow profile 70 32` builds, flashes and captures under the device lock.
