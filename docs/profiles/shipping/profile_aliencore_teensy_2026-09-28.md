# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienCore`).
Raw capture: `build/prof/aliencore_ship.log`, captured 2026-09-28 18:36 on COM4.
Replaces `profile_aliencore_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienCore 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienCore profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:68768, data:154376, headers:8280` / `RAM1: variables:315040, code:19720, padding:13048, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 15.33 ms/f; its worst window is 15.34 ms/f (frames 33–64). Peak frame render is **18.03 ms** (frame 431), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 31.53 ms.

The previous shipping report (2026-08-26 02:32) recorded peak 🟢 21.09 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 417–448)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       148.3 us    89.0 kcyc     0%
  fx_shader_draw          15.33 ms    9.20 Mcyc    25%
  fx_prepare_frame          3.1 us     1.9 kcyc     0%
  fx_advance               2.29 ms    1.37 Mcyc     4%
  fx_timeline_step         42.3 us    25.4 kcyc     0%
  canvas_clear             86.4 us    51.9 kcyc     0%
  canvas_buffer_wait      44.53 ms   26.72 Mcyc    71%
```

Wall min/avg/max = 62.22/62.44/62.54 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.6/22.8 us  cpu 3.00%
isr_pack          144/frame  min/avg/max 6.3/6.9/9.6 us  cpu 1.59%
isr_dma_submit    144/frame  min/avg/max 0.6/0.9/9.8 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 25% of the peak window, 15.33 ms/f.
2. `fx_advance` — 4% of the peak window, 2.29 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 18.03, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienCore`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AlienCore profile 70 32` builds, flashes and captures under the device lock.
