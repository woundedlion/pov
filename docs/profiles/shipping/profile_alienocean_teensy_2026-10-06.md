# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/alienocean_ship.log`, captured 2026-10-06 18:39 on COM3.
Replaces `profile_alienocean_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienOcean 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienOcean profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:68396, data:155184, headers:8864` / `RAM1: variables:315040, code:20328, padding:12440, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.2 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 19.03 ms/f; its worst window is 19.09 ms/f (frames 129–160). Peak frame render is **22.37 ms** (frame 822), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 42.55 ms.

The previous shipping report (2026-09-28 18:34) recorded peak 🟢 22.56 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 801–832)

```
frame                     62.44 ms   37.46 Mcyc   100%
  pov_preserve_half       142.3 us    85.4 kcyc     0%
  fx_shader_draw          19.04 ms   11.42 Mcyc    30%
  fx_prepare_frame        670.9 us   402.5 kcyc     1%
  fx_advance               2.21 ms    1.32 Mcyc     4%
  fx_timeline_step         38.0 us    22.8 kcyc     0%
  canvas_clear             85.1 us    51.1 kcyc     0%
  canvas_buffer_wait      40.23 ms   24.14 Mcyc    64%
```

Wall min/avg/max = 61.97/62.44/62.75 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/22.2 us  cpu 3.10%
isr_pack         144/frame  min/avg/max 6.3/6.9/9.3 us  cpu 1.59%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 30% of the peak window, 19.04 ms/f.
2. `fx_advance` — 4% of the peak window, 2.21 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.67 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 22.37, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.19 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienOcean`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AlienOcean profile 70 32` builds, flashes and captures under the device lock.
