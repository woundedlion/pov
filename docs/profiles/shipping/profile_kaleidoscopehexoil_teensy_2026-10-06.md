# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopehexoil_ship.log`, captured 2026-10-06 18:47 on COM4.
Replaces `profile_kaleidoscopehexoil_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 140 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexOil profile 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:72132, data:160476, headers:9052` / `RAM1: variables:315040, code:20280, padding:12488, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 497–512 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 30.88 ms/f; its worst window is 33.08 ms/f (frames 497–512). Peak frame render is **37.08 ms** (frame 546), and **0/2207** frames spilled. Setup frame 1 is excluded from both; it rendered 62.82 ms, over one display window.

The previous shipping report (2026-09-28 18:54) recorded peak 🟢 38.51 (2) and spilled 🟢 0/2207 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 2 once and wraps back to entry 1. The block below is the window holding the pass's peak frame.

### Peak window (frames 545–560)

```
frame                     62.40 ms   37.44 Mcyc   100%
  pov_preserve_half       139.6 us    83.7 kcyc     0%
  fx_shader_draw          32.06 ms   19.24 Mcyc    51%
  fx_prepare_frame         1.00 ms   601.6 kcyc     2%
  fx_advance               2.23 ms    1.34 Mcyc     4%
  fx_timeline_step         74.6 us    44.8 kcyc     0%
  canvas_clear             87.6 us    52.6 kcyc     0%
  canvas_buffer_wait      26.77 ms   16.06 Mcyc    43%
```

Wall min/avg/max = 59.74/62.40/65.05 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 37.08 | 0/1128 | 33.08 | 70/71 |
| 2 | — | 🟢 36.83 | 0/1079 | 32.30 | 66/67 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/13.2 us  cpu 3.14%
isr_pack         144/frame  min/avg/max 6.5/7.2/10.1 us  cpu 1.65%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/3.4 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 51% of the peak window, 32.06 ms/f.
2. `fx_advance` — 4% of the peak window, 2.23 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.00 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 37.08 (2), spilled 🟢 0/2207 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -1.43 ms) is not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexOil`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeHexOil profile 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
