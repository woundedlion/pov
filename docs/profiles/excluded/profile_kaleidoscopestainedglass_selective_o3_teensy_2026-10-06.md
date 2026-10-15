# KaleidoscopeStainedGlass on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopestainedglass_ship.log`, captured 2026-10-06 18:50 on COM4.
Replaces `profile_kaleidoscopestainedglass_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeStainedGlass 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeStainedGlass profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:70796, data:159724, headers:9092` / `RAM1: variables:315040, code:20216, padding:12552, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 545–576 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.6 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 36.51 ms/f; its worst window is 38.58 ms/f (frames 545–576). Peak frame render is **43.00 ms** (frame 604), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 73.08 ms, over one display window.

The previous shipping report (2026-09-28 17:35) recorded peak 🟢 43.50 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 577–608)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       138.0 us    82.8 kcyc     0%
  fx_shader_draw          38.43 ms   23.06 Mcyc    62%
  fx_prepare_frame         1.10 ms   660.0 kcyc     2%
  fx_advance               2.17 ms    1.30 Mcyc     3%
  fx_timeline_step         51.4 us    30.8 kcyc     0%
  canvas_clear             85.7 us    51.4 kcyc     0%
  canvas_buffer_wait      20.42 ms   12.25 Mcyc    33%
```

Wall min/avg/max = 60.57/62.44/64.30 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.7/21.9 us  cpu 3.17%
isr_pack         144/frame  min/avg/max 6.3/7.3/9.9 us  cpu 1.67%
isr_dma_submit   144/frame  min/avg/max 0.7/1.0/5.3 us  cpu 0.22%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 62% of the peak window, 38.43 ms/f.
2. `fx_advance` — 3% of the peak window, 2.17 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.10 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 43.00, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `0c02f3912` with `-D HS_PLACEMENT_VARIANT=6` (the kernel placement that landed in `97eb0bf78`), and the delta between them (peak -0.50 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeStainedGlass`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeStainedGlass profile 70 32` builds, flashes and captures under the device lock.
