# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AshCloud`).
Raw capture: `build/prof/ashcloud_ship.log`, captured 2026-09-28 17:33 on COM4.
Replaces `profile_ashcloud_teensy_2026-08-26.md`.

Captured as the V6 arm of the finding-10 placement A/B: tip `0c02f3912` built with `-D HS_PLACEMENT_VARIANT=6`, whose kernel attributes are exactly what landed in `97eb0bf78`. The commits between that tip and `97eb0bf78` do not touch this effect's render path.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; built with `-D HS_PLACEMENT_VARIANT=6` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, tip `0c02f3912` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AshCloud profile 70 32` |

Image size (`profile` env, this effect only): `.text.itcm 19,728 B`, `.text.code 51,416 B`, `.bss 310,240 B` (section sizes from the profile ELF; no `teensy_size` line).

Exactness cross-check: window frames 545–576 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.4 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 35.41 ms/f; its worst window is 37.66 ms/f (frames 545–576). Peak frame render is **44.57 ms** (frame 872), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 70.92 ms.

The previous shipping report (2026-08-26 02:47) recorded peak 🟢 50.09 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 865–896)

```
frame                     62.41 ms   37.45 Mcyc   100%
  pov_preserve_half       141.5 us    84.9 kcyc     0%
  fx_shader_draw          36.98 ms   22.19 Mcyc    59%
  fx_prepare_frame         3.81 ms    2.29 Mcyc     6%
  fx_advance               2.32 ms    1.39 Mcyc     4%
  fx_timeline_step         72.1 us    43.3 kcyc     0%
  canvas_clear             86.8 us    52.1 kcyc     0%
  canvas_buffer_wait      18.96 ms   11.38 Mcyc    30%
```

Wall min/avg/max = 60.34/62.41/64.51 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/18.0 us  cpu 3.14%
isr_pack          144/frame  min/avg/max 6.4/7.1/9.7 us  cpu 1.64%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/9.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 59% of the peak window, 36.98 ms/f.
2. `fx_prepare_frame` — 6% of the peak window, 3.81 ms/f.
3. `fx_advance` — 4% of the peak window, 2.32 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 44.57, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AshCloud profile 70 32` builds, flashes and captures under the device lock.
