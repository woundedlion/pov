# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile CosmicEyeball`).
Raw capture: `build/prof/cosmiceyeball_ship.log`, captured 2026-09-28 19:09 on COM4.
Replaces `profile_cosmiceyeball_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | CosmicEyeball 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh CosmicEyeball profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:68376, data:154216, headers:8832` / `RAM1: variables:315040, code:19544, padding:13224, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.0 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 20.81 ms/f; its worst window is 20.83 ms/f (frames 545–576). Peak frame render is **24.11 ms** (frame 464), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 42.70 ms.

The previous shipping report (2026-08-26 03:07) recorded peak 🟢 24.01 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 449–480)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       148.3 us    89.0 kcyc     0%
  fx_shader_draw          20.80 ms   12.48 Mcyc    33%
  fx_prepare_frame        589.8 us   353.9 kcyc     1%
  fx_advance               2.28 ms    1.37 Mcyc     4%
  fx_timeline_step         35.2 us    21.1 kcyc     0%
  canvas_clear             88.2 us    52.9 kcyc     0%
  canvas_buffer_wait      38.48 ms   23.09 Mcyc    62%
```

Wall min/avg/max = 62.22/62.44/62.64 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/22.8 us  cpu 3.07%
isr_pack          144/frame  min/avg/max 6.3/7.1/9.6 us  cpu 1.63%
isr_dma_submit    144/frame  min/avg/max 0.7/1.0/9.9 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 33% of the peak window, 20.80 ms/f.
2. `fx_advance` — 4% of the peak window, 2.28 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.59 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.

README cells: peak 🟢 24.11, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=CosmicEyeball`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh CosmicEyeball profile 70 32` builds, flashes and captures under the device lock.
