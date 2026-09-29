# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ChromaticLichen`).
Raw capture: `build/prof/chromaticlichen_ship.log`, captured 2026-09-28 18:47 on COM4.
Replaces `profile_chromaticlichen_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ChromaticLichen 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh ChromaticLichen profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:71064, data:154612, headers:8820` / `RAM1: variables:315040, code:19736, padding:13032, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 29.65 ms/f; its worst window is 29.72 ms/f (frames 801–832). Peak frame render is **36.63 ms** (frame 244), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 63.83 ms.

The previous shipping report (2026-08-26 02:44) recorded peak 🟢 42.87 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 225–256)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       143.1 us    85.9 kcyc     0%
  fx_shader_draw          29.67 ms   17.80 Mcyc    48%
  fx_prepare_frame         4.41 ms    2.65 Mcyc     7%
  fx_advance               2.04 ms    1.22 Mcyc     3%
  fx_timeline_step         68.5 us    41.1 kcyc     0%
  canvas_clear             84.9 us    50.9 kcyc     0%
  canvas_buffer_wait      26.00 ms   15.60 Mcyc    42%
```

Wall min/avg/max = 62.17/62.44/62.77 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/17.4 us  cpu 3.05%
isr_pack          144/frame  min/avg/max 6.2/6.8/9.8 us  cpu 1.57%
isr_dma_submit    144/frame  min/avg/max 0.6/0.9/10.2 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 48% of the peak window, 29.67 ms/f.
2. `fx_prepare_frame` — 7% of the peak window, 4.41 ms/f.
3. `fx_advance` — 3% of the peak window, 2.04 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 36.63, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh ChromaticLichen profile 70 32` builds, flashes and captures under the device lock.
