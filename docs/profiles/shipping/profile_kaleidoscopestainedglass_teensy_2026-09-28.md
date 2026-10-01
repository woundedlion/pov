# KaleidoscopeStainedGlass on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeStainedGlass`).
Raw capture: `build/prof/kaleidoscopestainedglass_ship.log`, captured 2026-09-28 17:35 on COM4.
Replaces `profile_kaleidoscopestainedglass_teensy_2026-08-26.md`.

Captured as the V6 arm of the kernel-placement A/B: tip `0c02f3912` built with `-D HS_PLACEMENT_VARIANT=6`, whose kernel attributes are exactly what landed in `97eb0bf78`. The commits between that tip and `97eb0bf78` do not touch this effect's render path.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; built with `-D HS_PLACEMENT_VARIANT=6` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeStainedGlass 288×144, single-entry playlist, tip `0c02f3912` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeStainedGlass profile 70 32` |

Image size (`profile` env, this effect only): `.text.itcm 19,600 B`, `.text.code 51,584 B`, `.bss 310,240 B` (section sizes from the profile ELF; no `teensy_size` line).

Exactness cross-check: window frames 545–576 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 36.99 ms/f; its worst window is 39.06 ms/f (frames 545–576). Peak frame render is **43.50 ms** (frame 598), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 73.61 ms.

The previous shipping report (2026-08-26 03:23) recorded peak 🟢 47.91 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 577–608)

```
frame                     62.43 ms   37.46 Mcyc   100%
  pov_preserve_half       137.6 us    82.6 kcyc     0%
  fx_shader_draw          38.92 ms   23.35 Mcyc    62%
  fx_prepare_frame         1.09 ms   656.7 kcyc     2%
  fx_advance               2.18 ms    1.31 Mcyc     3%
  fx_timeline_step         53.2 us    31.9 kcyc     0%
  canvas_clear             86.1 us    51.6 kcyc     0%
  canvas_buffer_wait      19.93 ms   11.96 Mcyc    32%
```

Wall min/avg/max = 60.54/62.43/64.23 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/21.9 us  cpu 3.19%
isr_pack          144/frame  min/avg/max 6.3/7.3/9.9 us  cpu 1.67%
isr_dma_submit    144/frame  min/avg/max 0.6/1.0/6.6 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 62% of the peak window, 38.92 ms/f.
2. `fx_advance` — 3% of the peak window, 2.18 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.09 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 43.50, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeStainedGlass`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeStainedGlass profile 70 32` builds, flashes and captures under the device lock.
