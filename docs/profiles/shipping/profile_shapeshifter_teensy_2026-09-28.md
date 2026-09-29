# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-09-28 18:48 on COM3.
Replaces `profile_shapeshifter_teensy_2026-09-25.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:76736, data:154168, headers:8712` / `RAM1: variables:315040, code:44456, padding:21080, free:143712` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2273–2288 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `ss_draw_all` averages 24.04 ms/f; its worst window is 51.82 ms/f (frames 2273–2288). Peak frame render is **58.24 ms** (frame 2301), and **0/2447** frames spilled. Setup frame 1 is excluded from both; it rendered 50.95 ms.

The previous shipping report (2026-09-25 07:33) recorded peak 🟢 58.40 (9) and spilled 🟢 0/2457 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 9 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2289–2304)

```
frame                     62.32 ms   37.39 Mcyc   100%
  pov_preserve_half       137.9 us    82.8 kcyc     0%
  ss_draw_all             51.41 ms   30.85 Mcyc    82%
    ss_plot_dispatch      51.19 ms   30.71 Mcyc    82%  x161  191071 cyc/c
  ss_timeline_step         53.9 us    32.4 kcyc     0%
  ss_buffer_wait          10.71 ms    6.43 Mcyc    17%
    canvas_clear           85.1 us    51.1 kcyc     0%
    canvas_buffer_wait    10.62 ms    6.37 Mcyc    17%
```

Wall min/avg/max = 49.21/62.32/75.30 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `ss_draw_all` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `ss_draw_all` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 58.24 | 0/478 | 51.82 | 30/30 |
| 9 | — | 🟢 38.26 | 0/240 | 37.27 | 15/15 |
| 4 | — | 🟢 37.39 | 0/240 | 36.13 | 15/15 |
| 8 | — | 🟢 26.45 | 0/240 | 24.87 | 15/15 |
| 7 | — | 🟢 22.43 | 0/240 | 20.83 | 15/15 |
| 6 | — | 🟢 18.97 | 0/240 | 18.12 | 15/15 |
| 2 | — | 🟢 15.78 | 0/289 | 13.42 | 18/18 |
| 5 | — | 🟢 11.35 | 0/240 | 9.74 | 15/15 |
| 3 | — | 🟢 9.47 | 0/240 | 8.17 | 15/15 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1150/frame  min/avg/max 0.6/1.7/11.5 us  cpu 3.07%
isr_pack          144/frame  min/avg/max 6.2/6.9/9.8 us  cpu 1.59%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/1.1 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `ss_draw_all` — 82% of the peak window, 51.41 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `ss_timeline_step` — 0% of the peak window, 0.05 ms/f.

README cells: peak 🟢 58.24 (9), spilled 🟢 0/2447 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
