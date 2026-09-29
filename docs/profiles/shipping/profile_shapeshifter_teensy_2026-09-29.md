# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-09-29 16:05 on COM3.
Replaces the 14:41 capture of the same date.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; the per-preset departure branch on `da3084407` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist; each preset dwells 224 frames, then departs through a 16-frame `Segue::Preset::Fade` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:74816, data:153536, headers:8192` / `RAM1: variables:315040, code:42376, padding:23160, free:143712` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2273–2288 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `ss_draw_all` averages 23.62 ms/f; its worst window is 52.05 ms/f (frames 2273–2288). Peak frame render is **58.88 ms** (frame 2299), and **0/2447** frames spilled. Setup frame 1 is excluded from both; it rendered 50.62 ms.

The 14:41 capture of `c40c1def8`, before presets carried their own departures, recorded peak 🟢 58.16 (9) and spilled 🟢 0/2447 (0.00%) on the same board. The fade is now a timeline transition rather than an envelope loop, and the cadence is unchanged at 240 frames per entry. Every entry is within 0.7 ms of that capture; entry 1 is the widest at +0.72 ms (+1.2%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 9 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2289–2304)

```
frame                        62.24 ms 37.35 Mcyc  100%
  pov_preserve_half          135.7 us  81.4 kcyc    0%
  ss_draw_all                51.77 ms 31.06 Mcyc   83%
    ss_plot_dispatch         51.56 ms 30.93 Mcyc   83%  x161  192434 cyc/c
  ss_timeline_step            30.4 us  18.3 kcyc    0%
  ss_buffer_wait             10.30 ms  6.18 Mcyc   17%
    canvas_clear              85.1 us  51.1 kcyc    0%
    canvas_buffer_wait       10.22 ms  6.13 Mcyc   16%
```

Wall min/avg/max = 49.47/62.24/74.89 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded. An entry's marker fires when its parameters are adopted at the fade's dark midpoint, so each bucket holds the second half of the fade into it, its dwell, and the first half of the fade out of it. Clean-hold `ss_draw_all` is the costliest window that neither entry boundary splits. Baseline is the 14:41 capture on the same board.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `ss_draw_all` ms/f | Baseline peak ms |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 58.88 | 0/470 | 52.05 | 58.16 |
| 9 | — | 🟢 37.12 | 0/239 | 36.29 | 37.47 |
| 4 | — | 🟢 35.55 | 0/239 | 34.58 | 36.58 |
| 8 | — | 🟢 25.29 | 0/239 | 23.89 | 25.80 |
| 7 | — | 🟢 21.78 | 0/239 | 20.14 | 21.89 |
| 6 | — | 🟢 18.14 | 0/239 | 17.48 | 18.56 |
| 2 | — | 🟢 14.99 | 0/305 | 13.12 | 15.39 |
| 5 | — | 🟢 11.68 | 0/239 | 9.91 | 11.34 |
| 3 | — | 🟢 9.57 | 0/239 | 8.06 | 9.50 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1148/frame  min/avg/max 0.6/1.7/11.5 us  cpu 3.07%
isr_pack          144/frame  min/avg/max 6.2/6.9/9.7 us  cpu 1.59%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `ss_draw_all` — 83% of the peak window, 51.77 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `ss_timeline_step` — 0% of the peak window, 0.03 ms/f.

README cells: peak 🟢 58.88 (9), spilled 🟢 0/2447 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: this effect's `HS_O3` regions are unchanged; the departure branch touched only the preset choreography.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
