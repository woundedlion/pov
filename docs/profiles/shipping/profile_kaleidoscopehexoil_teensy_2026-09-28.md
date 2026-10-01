# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexOil`).
Raw capture: `build/prof/kaleidoscopehexoil_ship.log`, captured 2026-09-28 18:54 on COM4.
Replaces `profile_kaleidoscopehexoil_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 140 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexOil profile 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:69720, data:158544, headers:8280` / `RAM1: variables:315040, code:19560, padding:13208, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 481–496 root counter cyc ÷ 600 MHz matches the measured wall sum within **4.0 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 32.26 ms/f; its worst window is 34.54 ms/f (frames 481–496). Peak frame render is **38.51 ms** (frame 522), and **0/2207** frames spilled. Setup frame 1 is excluded from both; it rendered 65.26 ms.

The previous shipping report (2026-08-26 02:52) recorded peak 🟢 39.04 (3) and spilled 🟢 0/2208 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 513–528)

```
frame                     62.43 ms   37.46 Mcyc   100%
  pov_preserve_half       138.2 us    82.9 kcyc     0%
  fx_shader_draw          34.23 ms   20.54 Mcyc    55%
  fx_prepare_frame         1.04 ms   621.4 kcyc     2%
  fx_advance               2.21 ms    1.33 Mcyc     4%
  fx_timeline_step         80.8 us    48.5 kcyc     0%
  canvas_clear             86.4 us    51.9 kcyc     0%
  canvas_buffer_wait      24.61 ms   14.77 Mcyc    39%
```

Wall min/avg/max = 61.05/62.43/63.86 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 38.51 | 0/1128 | 34.54 | 71/71 |
| 2 | — | 🟢 38.03 | 0/1079 | 33.66 | 67/67 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.6/1.7/16.9 us  cpu 3.14%
isr_pack          144/frame  min/avg/max 6.4/7.2/10.8 us  cpu 1.66%
isr_dma_submit    144/frame  min/avg/max 0.7/1.0/8.9 us  cpu 0.22%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 55% of the peak window, 34.23 ms/f.
2. `fx_advance` — 4% of the peak window, 2.21 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.04 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 38.51 (2), spilled 🟢 0/2207 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexOil`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeHexOil profile 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
