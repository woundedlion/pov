# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienOcean`).
Raw capture: `build/prof/alienocean_ship.log`, captured 2026-09-28 18:34 on COM4.
Replaces `profile_alienocean_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienOcean 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienOcean profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:68784, data:154376, headers:8264` / `RAM1: variables:315040, code:19720, padding:13048, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 19.12 ms/f; its worst window is 19.20 ms/f (frames 97–128). Peak frame render is **22.56 ms** (frame 810), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 42.60 ms.

The previous shipping report (2026-08-26 02:31) recorded peak 🟢 27.79 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 801–832)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       145.7 us    87.4 kcyc     0%
  fx_shader_draw          19.15 ms   11.49 Mcyc    31%
  fx_prepare_frame        682.4 us   409.4 kcyc     1%
  fx_advance               2.21 ms    1.33 Mcyc     4%
  fx_timeline_step         53.9 us    32.3 kcyc     0%
  canvas_clear             85.0 us    51.0 kcyc     0%
  canvas_buffer_wait      40.10 ms   24.06 Mcyc    64%
```

Wall min/avg/max = 62.03/62.44/62.85 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.6/19.1 us  cpu 3.01%
isr_pack          144/frame  min/avg/max 6.2/6.9/11.6 us  cpu 1.59%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/9.3 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 31% of the peak window, 19.15 ms/f.
2. `fx_advance` — 4% of the peak window, 2.21 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.68 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.

README cells: peak 🟢 22.56, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienOcean`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AlienOcean profile 70 32` builds, flashes and captures under the device lock.
