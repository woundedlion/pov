# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/gridspace_ship.log`, captured 2026-10-06 18:48 on COM3.
Replaces `profile_gridspace_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GridSpace 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh GridSpace profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:69068, data:156324, headers:9100` / `RAM1: variables:315040, code:20280, padding:12488, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 17.76 ms/f; its worst window is 17.80 ms/f (frames 577–608). Peak frame render is **24.57 ms** (frame 460), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 40.32 ms.

The previous shipping report (2026-09-28 17:38) recorded peak 🟢 25.31 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 449–480)

```
frame                     62.45 ms   37.47 Mcyc   100%
  pov_preserve_half       147.4 us    88.5 kcyc     0%
  fx_shader_draw          17.76 ms   10.66 Mcyc    28%
  fx_prepare_frame         3.97 ms    2.38 Mcyc     6%
  fx_advance               2.29 ms    1.37 Mcyc     4%
  fx_timeline_step         42.2 us    25.3 kcyc     0%
  canvas_clear             85.8 us    51.5 kcyc     0%
  canvas_buffer_wait      38.10 ms   22.86 Mcyc    61%
```

Wall min/avg/max = 62.15/62.45/62.76 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.7/23.0 us  cpu 3.08%
isr_pack         144/frame  min/avg/max 6.3/6.9/9.5 us  cpu 1.59%
isr_dma_submit   144/frame  min/avg/max 0.7/1.0/5.4 us  cpu 0.22%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 28% of the peak window, 17.76 ms/f.
2. `fx_prepare_frame` — 6% of the peak window, 3.97 ms/f.
3. `fx_advance` — 4% of the peak window, 2.29 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.

README cells: peak 🟢 24.57, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `0c02f3912` with `-D HS_PLACEMENT_VARIANT=6` (the kernel placement that landed in `97eb0bf78`), and the delta between them (peak -0.74 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GridSpace`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh GridSpace profile 70 32` builds, flashes and captures under the device lock.
