# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/aliencore_ship.log`, captured 2026-10-06 18:42 on COM3.
Replaces `profile_aliencore_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienCore 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienCore profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:68412, data:155184, headers:8848` / `RAM1: variables:315040, code:20328, padding:12440, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 15.24 ms/f; its worst window is 15.23 ms/f (frames 769–800). Peak frame render is **17.86 ms** (frame 978), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 31.66 ms.

The previous shipping report (2026-09-28 18:36) recorded peak 🟢 18.03 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 961–992)

```
frame                     62.43 ms   37.46 Mcyc   100%
  pov_preserve_half       150.9 us    90.6 kcyc     0%
  fx_shader_draw          15.22 ms    9.13 Mcyc    24%
  fx_prepare_frame          2.0 us     1.2 kcyc     0%
  fx_advance               2.30 ms    1.38 Mcyc     4%
  fx_timeline_step         32.4 us    19.5 kcyc     0%
  canvas_clear             85.5 us    51.3 kcyc     0%
  canvas_buffer_wait      44.62 ms   26.77 Mcyc    71%
```

Wall min/avg/max = 62.28/62.43/62.51 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/19.0 us  cpu 3.12%
isr_pack         144/frame  min/avg/max 6.3/7.1/9.8 us  cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/5.8 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 24% of the peak window, 15.22 ms/f.
2. `fx_advance` — 4% of the peak window, 2.30 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 17.86, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.17 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienCore`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AlienCore profile 70 32` builds, flashes and captures under the device lock.
