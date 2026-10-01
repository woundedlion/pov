# LatticeMelt on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile LatticeMelt`).
Raw capture: `build/prof/latticemelt_ship.log`, captured 2026-09-28 18:45 on COM4.
Replaces `profile_latticemelt_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | LatticeMelt 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture, `-D HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh LatticeMelt profile 110 16 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size (`profile` env, this effect only): `FLASH: code:70696, data:154608, headers:9192` / `RAM1: variables:315040, code:19736, padding:13032, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–16 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 33.18 ms/f; its worst window is 33.33 ms/f (frames 385–400). Peak frame render is **37.18 ms** (frame 601), and **0/1727** frames spilled. Setup frame 1 is excluded from both; it rendered 70.09 ms.

The previous shipping report (2026-08-26 02:42) recorded peak 🟢 43.73 (3) and spilled 🟢 0/1728 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 593–608)

```
frame                     62.41 ms   37.45 Mcyc   100%
  pov_preserve_half       139.1 us    83.5 kcyc     0%
  fx_shader_draw          33.19 ms   19.92 Mcyc    53%
  fx_prepare_frame         1.06 ms   633.2 kcyc     2%
  fx_advance               2.19 ms    1.31 Mcyc     4%
  fx_timeline_step        112.9 us    67.8 kcyc     0%
  canvas_clear             86.3 us    51.8 kcyc     0%
  canvas_buffer_wait      25.61 ms   15.37 Mcyc    41%
```

Wall min/avg/max = 58.85/62.41/65.81 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 🟢 37.18 | 0/1079 | 33.26 | 67/67 |
| 1 | — | 🟢 37.18 | 0/648 | 35.05 | 41/41 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/11.6 us  cpu 3.07%
isr_pack          144/frame  min/avg/max 6.3/6.9/9.8 us  cpu 1.60%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/1.2 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 53% of the peak window, 33.19 ms/f.
2. `fx_advance` — 4% of the peak window, 2.19 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.06 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 37.18 (2), spilled 🟢 0/1727 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=LatticeMelt`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh LatticeMelt profile 110 16 "-D HS_PROFILE_EPOCH_REVS=1200"` builds, flashes and captures under the device lock.
