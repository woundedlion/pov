# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/mermaidskin_ship.log`, captured 2026-10-06 18:59 on COM3.
Replaces `profile_mermaidskin_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MermaidSkin 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh MermaidSkin profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:70636, data:155300, headers:8556` / `RAM1: variables:315040, code:20264, padding:12504, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.9 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 32.53 ms/f; its worst window is 32.73 ms/f (frames 545–576). Peak frame render is **39.18 ms** (frame 226), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 68.38 ms, over one display window.

The previous shipping report (2026-09-28 18:49) recorded peak 🟢 39.57 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 225–256)

```
frame                     62.42 ms   37.45 Mcyc   100%
  pov_preserve_half       139.0 us    83.4 kcyc     0%
  fx_shader_draw          32.48 ms   19.49 Mcyc    52%
  fx_prepare_frame         4.23 ms    2.54 Mcyc     7%
  fx_advance               1.92 ms    1.15 Mcyc     3%
  fx_timeline_step         66.8 us    40.1 kcyc     0%
  canvas_clear             85.4 us    51.2 kcyc     0%
  canvas_buffer_wait      23.48 ms   14.09 Mcyc    38%
```

Wall min/avg/max = 62.08/62.42/62.77 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.7/18.4 us  cpu 3.21%
isr_pack         144/frame  min/avg/max 6.6/7.5/10.4 us  cpu 1.72%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/4.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 52% of the peak window, 32.48 ms/f.
2. `fx_prepare_frame` — 7% of the peak window, 4.23 ms/f.
3. `fx_advance` — 3% of the peak window, 1.92 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 39.18, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.39 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MermaidSkin`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh MermaidSkin profile 70 32` builds, flashes and captures under the device lock.
