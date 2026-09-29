# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `build/prof/displacementfield_ship.log`, captured 2026-09-28 19:16 on COM4.
Replaces `profile_displacementfield_teensy_2026-09-19.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:73984, data:153336, headers:8200` / `RAM1: variables:315200, code:33032, padding:32504, free:143552` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 609–640 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `df_timeline_step` averages 48.58 ms/f; its worst window is 56.66 ms/f (frames 609–640). Peak frame render is **59.42 ms** (frame 605), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 21.27 ms.

The previous shipping report (2026-09-19 22:17) recorded peak 🟢 58.18 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 577–608)

```
frame                     62.54 ms   37.53 Mcyc   100%
  pov_preserve_half       138.0 us    82.8 kcyc     0%
  df_timeline_step        55.91 ms   33.54 Mcyc    89%
    df_draw_rings         55.82 ms   33.49 Mcyc    89%
      df_hue_table_prep    1.70 ms    1.02 Mcyc     3%  x20.8  82 us/c
      df_lut_bake         10.61 ms    6.37 Mcyc    17%  x41.6  255 us/c
      df_chunk_cull        1.06 ms   635.8 kcyc     2%  x43.4  24 us/c
      df_fused_scan       40.96 ms   24.57 Mcyc    65%
        filter_blend       1.06 ms   638.1 kcyc     2%  x10179  63 cyc/c
  canvas_clear             85.2 us    51.2 kcyc     0%
  canvas_buffer_wait       6.41 ms    3.85 Mcyc    10%
  df_prepare_fields         0.3 us      186 cyc     0%
```

Wall min/avg/max = 58.64/62.54/66.84 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 10,179 times per frame in the peak window at 63 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1154/frame  min/avg/max 0.6/1.7/20.2 us  cpu 3.15%
isr_pack          144/frame  min/avg/max 6.4/7.1/10.0 us  cpu 1.63%
isr_dma_submit    144/frame  min/avg/max 0.6/0.9/2.8 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `df_timeline_step` — 89% of the peak window, 55.91 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.09 ms/f.
4. `df_prepare_fields` — 0% of the peak window, 0.00 ms/f.

README cells: peak 🟢 59.42, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh DisplacementField profile 70 32` builds, flashes and captures under the device lock.
