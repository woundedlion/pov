# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/cosmiceyeball_ship.log`, captured 2026-10-06 19:06 on COM4.
Replaces `profile_cosmiceyeball_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | CosmicEyeball 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh CosmicEyeball profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:67676, data:155024, headers:8720` / `RAM1: variables:315040, code:19896, padding:12872, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 20.33 ms/f; its worst window is 20.33 ms/f (frames 545–576). Peak frame render is **23.63 ms** (frame 539), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 41.85 ms.

The previous shipping report (2026-09-28 19:09) recorded peak 🟢 24.11 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 513–544)

```
frame                     62.46 ms   37.48 Mcyc   100%
  pov_preserve_half       149.5 us    89.7 kcyc     0%
  fx_shader_draw          20.32 ms   12.19 Mcyc    33%
  fx_prepare_frame        600.1 us   360.1 kcyc     1%
  fx_advance               2.28 ms    1.37 Mcyc     4%
  fx_timeline_step         41.2 us    24.7 kcyc     0%
  canvas_clear             85.5 us    51.3 kcyc     0%
  canvas_buffer_wait      38.96 ms   23.38 Mcyc    62%
```

Wall min/avg/max = 62.35/62.46/62.58 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/22.5 us  cpu 3.17%
isr_pack         144/frame  min/avg/max 6.4/7.2/9.8 us  cpu 1.65%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/8.3 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 33% of the peak window, 20.32 ms/f.
2. `fx_advance` — 4% of the peak window, 2.28 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.60 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.

README cells: peak 🟢 23.63, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.48 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=CosmicEyeball`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh CosmicEyeball profile 70 32` builds, flashes and captures under the device lock.
