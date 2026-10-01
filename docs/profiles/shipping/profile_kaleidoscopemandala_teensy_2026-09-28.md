# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeMandala`).
Raw capture: `build/prof/kaleidoscopemandala_ship.log`, captured 2026-09-28 18:39 on COM4.
Replaces `profile_kaleidoscopemandala_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeMandala 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeMandala profile 150 32` |

Image size (`profile` env, this effect only): `FLASH: code:69312, data:154720, headers:8416` / `RAM1: variables:315040, code:19720, padding:13048, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1121–1152 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.0 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 26.22 ms/f; its worst window is 28.41 ms/f (frames 1121–1152). Peak frame render is **35.40 ms** (frame 604), and **0/2367** frames spilled. Setup frame 1 is excluded from both; it rendered 52.33 ms.

The previous shipping report (2026-08-26 03:22) recorded peak 🟢 40.29 (3) and spilled 🟢 0/2368 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 577–608)

```
frame                     62.55 ms   37.53 Mcyc   100%
  pov_preserve_half       143.7 us    86.2 kcyc     0%
  fx_shader_draw          27.38 ms   16.43 Mcyc    44%
  fx_prepare_frame         1.88 ms    1.13 Mcyc     3%
  fx_advance               2.17 ms    1.30 Mcyc     3%
  fx_timeline_step         73.6 us    44.1 kcyc     0%
  canvas_clear             87.0 us    52.2 kcyc     0%
  canvas_buffer_wait      30.79 ms   18.48 Mcyc    49%
```

Wall min/avg/max = 60.56/62.55/64.64 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 🟢 35.40 | 0/1079 | 28.41 | 34/34 |
| 1 | — | 🟢 34.88 | 0/1288 | 27.61 | 40/40 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1154/frame  min/avg/max 0.5/1.7/19.6 us  cpu 3.10%
isr_pack          144/frame  min/avg/max 6.4/7.2/9.7 us  cpu 1.65%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/9.6 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 44% of the peak window, 27.38 ms/f.
2. `fx_advance` — 3% of the peak window, 2.17 ms/f.
3. `fx_prepare_frame` — 3% of the peak window, 1.88 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 35.40 (2), spilled 🟢 0/2367 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeMandala`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeMandala profile 150 32` builds, flashes and captures under the device lock.
