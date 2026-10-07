# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/ashcloud_ship.log`, captured 2026-10-06 19:01 on COM3.
Replaces `profile_ashcloud_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AshCloud profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:70236, data:155452, headers:8804` / `RAM1: variables:315040, code:20280, padding:12488, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 545–576 root counter cyc ÷ 600 MHz matches the measured wall sum within **4.0 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 34.14 ms/f; its worst window is 36.06 ms/f (frames 545–576). Peak frame render is **42.67 ms** (frame 540), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 69.73 ms, over one display window.

The previous shipping report (2026-09-28 17:33) recorded peak 🟢 44.57 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 513–544)

```
frame                     62.43 ms   37.46 Mcyc   100%
  pov_preserve_half       141.8 us    85.1 kcyc     0%
  fx_shader_draw          35.77 ms   21.46 Mcyc    57%
  fx_prepare_frame         3.59 ms    2.15 Mcyc     6%
  fx_advance               2.22 ms    1.33 Mcyc     4%
  fx_timeline_step         72.1 us    43.2 kcyc     0%
  canvas_clear             86.0 us    51.6 kcyc     0%
  canvas_buffer_wait      20.52 ms   12.31 Mcyc    33%
```

Wall min/avg/max = 60.60/62.43/64.17 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/22.6 us  cpu 3.14%
isr_pack         144/frame  min/avg/max 6.3/7.1/9.7 us  cpu 1.64%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/9.2 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 57% of the peak window, 35.77 ms/f.
2. `fx_prepare_frame` — 6% of the peak window, 3.59 ms/f.
3. `fx_advance` — 4% of the peak window, 2.22 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 42.67, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `0c02f3912` with `-D HS_PLACEMENT_VARIANT=6` (the kernel placement that landed in `97eb0bf78`), and the delta between them (peak -1.90 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AshCloud profile 70 32` builds, flashes and captures under the device lock.
