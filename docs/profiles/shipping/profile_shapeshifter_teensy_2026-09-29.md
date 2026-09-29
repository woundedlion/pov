# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-09-29 14:41 on COM3.
Replaces `profile_shapeshifter_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, tip `d5ca81403` (landed as `c40c1def8`) |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:75152, data:154168, headers:8248` / `RAM1: variables:315040, code:42616, padding:22920, free:143712` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2273–2288 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `ss_draw_all` averages 23.69 ms/f; its worst window is 51.78 ms/f (frames 2273–2288). Peak frame render is **58.16 ms** (frame 2301), and **0/2447** frames spilled. Setup frame 1 is excluded from both; it rendered 50.78 ms.

The previous shipping report (2026-09-28 18:48) recorded peak 🟢 58.24 (9) and spilled 🟢 0/2447 (0.00%). A same-board baseline of the parent commit `a2c6d2be0` (COM3, 2026-09-29 14:55) recorded peak 58.24 ms; the float min/max sweep in `c40c1def8` cut entries 2, 4, 6, 7, 8 and 9 by 2.2–2.5% and left entries 1, 3 and 5 within 0.4%.

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 9 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2289–2304)

```
frame                        62.32 ms 37.39 Mcyc  100%
  pov_preserve_half          138.3 us  83.0 kcyc    0%
  ss_draw_all                51.32 ms 30.79 Mcyc   82%
    ss_plot_dispatch         51.12 ms 30.67 Mcyc   82%  x161  190800 cyc/c
  ss_timeline_step            54.0 us  32.4 kcyc    0%
  ss_buffer_wait             10.79 ms  6.48 Mcyc   17%
    canvas_clear              84.7 us  50.8 kcyc    0%
    canvas_buffer_wait       10.71 ms  6.43 Mcyc   17%
```

Wall min/avg/max = 49.14/62.32/75.34 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `ss_draw_all` is the costliest modal-call-count window of each entry. Baseline is the same-board `a2c6d2be0` capture.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `ss_draw_all` ms/f | Clean windows | Baseline peak ms |
|---|---|--:|--:|--:|--:|--:|
| 1 | — | 🟢 58.16 | 0/478 | 51.78 | 30/30 | 58.24 |
| 9 | — | 🟢 37.47 | 0/240 | 36.47 | 15/15 | 38.32 |
| 4 | — | 🟢 36.58 | 0/240 | 35.35 | 15/15 | 37.46 |
| 8 | — | 🟢 25.80 | 0/240 | 24.24 | 15/15 | 26.45 |
| 7 | — | 🟢 21.89 | 0/240 | 20.34 | 15/15 | 22.44 |
| 6 | — | 🟢 18.56 | 0/240 | 17.69 | 15/15 | 18.97 |
| 2 | — | 🟢 15.39 | 0/289 | 13.08 | 18/18 | 15.77 |
| 5 | — | 🟢 11.34 | 0/240 | 9.71 | 15/15 | 11.39 |
| 3 | — | 🟢 9.50 | 0/240 | 8.16 | 15/15 | 9.48 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1150/frame  min/avg/max 0.6/1.7/11.1 us  cpu 3.07%
isr_pack          144/frame  min/avg/max 6.2/6.9/9.4 us  cpu 1.59%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `ss_draw_all` — 82% of the peak window, 51.32 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `ss_timeline_step` — 0% of the peak window, 0.05 ms/f.

README cells: peak 🟢 58.16 (9), spilled 🟢 0/2447 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: this effect's `HS_O3` regions are unchanged; `c40c1def8` changed only the Plot stroke and trail float extremes.
- The captured tip `d5ca81403` differs from the landed `c40c1def8` only by an earlier HyperLattice commit, which is not in this single-effect image.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
