# KaleidoscopeHexBright on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexBright`).
Raw capture: `build/prof/kaleidoscopehexbright_ship.log`, captured 2026-09-28 19:02 on COM4.
Replaces `profile_kaleidoscopehexbright_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexBright 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexBright profile 150 32` |

Image size (`profile` env, this effect only): `FLASH: code:69024, data:154756, headers:8668` / `RAM1: variables:315040, code:19736, padding:13032, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 385–416 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.2 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 22.65 ms/f; its worst window is 24.85 ms/f (frames 385–416). Peak frame render is **31.64 ms** (frame 915), and **0/2367** frames spilled. Setup frame 1 is excluded from both; it rendered 50.65 ms.

The previous shipping report (2026-08-26 03:19) recorded peak 🟢 35.55 (3) and spilled 🟢 0/2368 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 897–928)

```
frame                     62.44 ms   37.46 Mcyc   100%
  pov_preserve_half       141.6 us    85.0 kcyc     0%
  fx_shader_draw          23.89 ms   14.34 Mcyc    38%
  fx_prepare_frame         4.17 ms    2.50 Mcyc     7%
  fx_advance               2.00 ms    1.20 Mcyc     3%
  fx_timeline_step        133.8 us    80.3 kcyc     0%
  canvas_clear             86.4 us    51.8 kcyc     0%
  canvas_buffer_wait      31.99 ms   19.20 Mcyc    51%
```

Wall min/avg/max = 60.18/62.44/64.71 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 🟢 31.64 | 0/1079 | 24.05 | 34/34 |
| 1 | — | 🟢 31.44 | 0/1288 | 24.85 | 40/40 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/15.6 us  cpu 3.13%
isr_pack          144/frame  min/avg/max 6.4/7.2/9.7 us  cpu 1.66%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/2.9 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 38% of the peak window, 23.89 ms/f.
2. `fx_prepare_frame` — 7% of the peak window, 4.17 ms/f.
3. `fx_advance` — 3% of the peak window, 2.00 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 31.64 (2), spilled 🟢 0/2367 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexBright`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeHexBright profile 150 32` builds, flashes and captures under the device lock.
