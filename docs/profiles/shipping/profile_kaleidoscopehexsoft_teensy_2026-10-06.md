# KaleidoscopeHexSoft on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopehexsoft_ship.log`, captured 2026-10-06 18:37 on COM3.
Replaces `profile_kaleidoscopehexsoft_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexSoft 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexSoft profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:68260, data:155436, headers:8748` / `RAM1: variables:315040, code:20280, padding:12488, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 17.57 ms/f; its worst window is 17.65 ms/f (frames 65–96). Peak frame render is **23.88 ms** (frame 33), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 39.84 ms.

The previous shipping report (2026-09-28 18:32) recorded peak 🟢 24.52 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 33–64)

```
frame                     62.45 ms   37.47 Mcyc   100%
  pov_preserve_half       142.7 us    85.6 kcyc     0%
  fx_shader_draw          17.64 ms   10.58 Mcyc    28%
  fx_prepare_frame         3.55 ms    2.13 Mcyc     6%
  fx_advance               2.20 ms    1.32 Mcyc     4%
  fx_timeline_step         58.9 us    35.4 kcyc     0%
  canvas_clear             86.4 us    51.8 kcyc     0%
  canvas_buffer_wait      38.77 ms   23.26 Mcyc    62%
```

Wall min/avg/max = 62.12/62.45/62.67 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.7/27.2 us  cpu 3.07%
isr_pack         144/frame  min/avg/max 6.3/7.0/13.1 us  cpu 1.60%
isr_dma_submit   144/frame  min/avg/max 0.8/0.9/6.6 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 28% of the peak window, 17.64 ms/f.
2. `fx_prepare_frame` — 6% of the peak window, 3.55 ms/f.
3. `fx_advance` — 4% of the peak window, 2.20 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 23.88, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.64 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexSoft`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeHexSoft profile 70 32` builds, flashes and captures under the device lock.
