# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopemandala_ship.log`, captured 2026-10-06 18:45 on COM3.
Replaces `profile_kaleidoscopemandala_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeMandala 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeMandala profile 150 32` |

Image size (`profile` env, this effect only): `FLASH: code:72012, data:156860, headers:8692` / `RAM1: variables:315040, code:20536, padding:12232, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 545–576 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 25.64 ms/f; its worst window is 27.36 ms/f (frames 545–576). Peak frame render is **34.97 ms** (frame 604), and **0/2367** frames spilled. Setup frame 1 is excluded from both; it rendered 51.63 ms.

The previous shipping report (2026-09-28 18:39) recorded peak 🟢 35.40 (2) and spilled 🟢 0/2367 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 2 once and wraps back to entry 1. The block below is the window holding the pass's peak frame.

### Peak window (frames 577–608)

```
frame                     62.52 ms   37.51 Mcyc   100%
  pov_preserve_half       140.6 us    84.4 kcyc     0%
  fx_shader_draw          27.10 ms   16.26 Mcyc    43%
  fx_prepare_frame         1.84 ms    1.10 Mcyc     3%
  fx_advance               2.24 ms    1.35 Mcyc     4%
  fx_timeline_step         74.2 us    44.5 kcyc     0%
  canvas_clear             86.9 us    52.1 kcyc     0%
  canvas_buffer_wait      31.00 ms   18.60 Mcyc    50%
```

Wall min/avg/max = 60.49/62.52/64.39 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 🟢 34.97 | 0/1079 | 27.19 | 33/34 |
| 1 | — | 🟢 34.40 | 0/1288 | 27.36 | 39/40 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1154/frame  min/avg/max 0.6/1.7/23.4 us  cpu 3.18%
isr_pack         144/frame  min/avg/max 6.5/7.2/9.8 us  cpu 1.66%
isr_dma_submit   144/frame  min/avg/max 0.8/0.9/8.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 43% of the peak window, 27.10 ms/f.
2. `fx_advance` — 4% of the peak window, 2.24 ms/f.
3. `fx_prepare_frame` — 3% of the peak window, 1.84 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 34.97 (2), spilled 🟢 0/2367 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.43 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeMandala`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeMandala profile 150 32` builds, flashes and captures under the device lock.
