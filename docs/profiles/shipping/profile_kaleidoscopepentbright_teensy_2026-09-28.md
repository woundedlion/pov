# KaleidoscopePentBright on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopePentBright`).
Raw capture: `build/prof/kaleidoscopepentbright_ship.log`, captured 2026-09-28 18:51 on COM4.
Replaces `profile_kaleidoscopepentbright_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopePentBright 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopePentBright profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:69072, data:154976, headers:8400` / `RAM1: variables:315040, code:19736, padding:13032, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 417–448 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.4 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 22.31 ms/f; its worst window is 23.79 ms/f (frames 417–448). Peak frame render is **27.24 ms** (frame 437), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 49.01 ms.

The previous shipping report (2026-08-26 02:49) recorded peak 🟢 33.37 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 417–448)

```
frame                     62.43 ms   37.46 Mcyc   100%
  pov_preserve_half       144.4 us    86.6 kcyc     0%
  fx_shader_draw          23.79 ms   14.28 Mcyc    38%
  fx_prepare_frame        816.6 us   490.0 kcyc     1%
  fx_advance               1.97 ms    1.18 Mcyc     3%
  fx_timeline_step         45.3 us    27.2 kcyc     0%
  canvas_clear             86.3 us    51.8 kcyc     0%
  canvas_buffer_wait      35.56 ms   21.33 Mcyc    57%
```

Wall min/avg/max = 61.80/62.43/63.10 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/23.5 us  cpu 3.11%
isr_pack          144/frame  min/avg/max 6.3/7.2/9.7 us  cpu 1.66%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/4.4 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 38% of the peak window, 23.79 ms/f.
2. `fx_advance` — 3% of the peak window, 1.97 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.82 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 27.24, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopePentBright`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopePentBright profile 70 32` builds, flashes and captures under the device lock.
