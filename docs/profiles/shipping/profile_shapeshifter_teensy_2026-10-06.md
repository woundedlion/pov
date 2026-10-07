# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-10-06 18:29 on COM4.
Replaces `profile_shapeshifter_teensy_2026-09-29.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:89524, data:156720, headers:8728` / `RAM1: variables:315040, code:44312, padding:21224, free:143712` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2257–2272 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `ss_draw_all` averages 19.55 ms/f; its worst window is 43.62 ms/f (frames 2257–2272). Peak frame render is **54.27 ms** (frame 2303), and **0/2447** frames spilled. Setup frame 1 is excluded from both; it rendered 49.19 ms.

The previous shipping report (2026-09-29 22:38) recorded peak 🟢 53.22 (9) and spilled 🟢 0/2448 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 9 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 9, wraps back to entry 1 and ends on entry 2. The block below is the window holding the pass's peak frame.

### Peak window (frames 2289–2304)

```
frame                   62.29 ms  37.38 Mcyc   100%
  pov_preserve_half     137.8 us   82.7 kcyc     0%
  ss_draw_all           42.50 ms  25.50 Mcyc    68%
    ss_plot_dispatch    42.14 ms  25.29 Mcyc    68%  x222  113642 cyc/c
  ss_timeline_step       29.8 us   17.9 kcyc     0%
  ss_buffer_wait        19.62 ms  11.77 Mcyc    31%
    canvas_clear         85.0 us   51.0 kcyc     0%
    canvas_buffer_wait  19.53 ms  11.72 Mcyc    31%
```

Wall min/avg/max = 40.83/62.29/82.86 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `ss_draw_all` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Shape (count) | Peak render ms | Spilled/frames | Clean `ss_draw_all` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | Planar Star (288, screen balanced) | 🟢 54.27 | 0/469 | 43.62 | 28/30 |
| 8 | Spherical Polygon (144) | 🟢 25.26 | 0/239 | 23.86 | 14/14 |
| 9 | Flower (72) | 🟢 24.55 | 0/239 | 24.09 | 14/15 |
| 4 | Flower (70) | 🟢 23.79 | 0/239 | 23.10 | 14/15 |
| 7 | Spherical Polygon (144) | 🟢 21.66 | 0/239 | 20.03 | 14/15 |
| 6 | Spherical Polygon (128) | 🟢 18.02 | 0/239 | 17.34 | 14/15 |
| 2 | Spherical Polygon (74) | 🟢 14.99 | 0/305 | 13.12 | 18/19 |
| 5 | Planar Star (72) | 🟢 8.26 | 0/239 | 5.98 | 14/15 |
| 3 | Planar Star (43) | 🟢 6.28 | 0/239 | 5.04 | 14/15 |

Entry 1 sets the peak, 8.2 ms inside the 62.5 ms window; the next entries sit near 25 ms.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1149/frame  min/avg/max 0.5/1.7/11.2 us  cpu 3.08%
isr_pack         144/frame  min/avg/max 6.2/6.9/9.3 us  cpu 1.59%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `ss_draw_all` — 68% of the peak window, 42.50 ms/f (`ss_plot_dispatch` 42.14).
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `ss_timeline_step` — 0% of the peak window, 0.03 ms/f.

README cells: peak 🟢 54.27 (9), spilled 🟢 0/2447 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `bb27f0c55` plus the pole-run split (landed as `2efde6cee`), and the delta between them (peak +1.05 ms) is not attributed here.
- The cycle samples entry 1 only during its two holds; pin it with `-D HS_PROFILE_PRESET=0` for a worst-case figure.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
