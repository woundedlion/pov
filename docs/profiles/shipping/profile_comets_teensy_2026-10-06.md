# Comets on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/comets_ship.log`, captured 2026-10-06 18:07 on COM3.
Replaces `profile_comets_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Comets 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:69316, data:153896, headers:8208` / `RAM1: variables:315040, code:30168, padding:2600, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1537–1552 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `cm_draw_trail` averages 16.71 ms/f; its worst window is 25.36 ms/f (frames 1537–1552). Peak frame render is **31.92 ms** (frame 645), and **0/4127** frames spilled. Setup frame 1 is excluded from both; it rendered 0.69 ms.

The previous shipping report (2026-09-28 18:28) recorded peak 🟢 30.66 (12) and spilled 🟢 0/4127 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 12 twice and wraps back to entry 1 twice. The block below is the window holding the pass's peak frame.

### Peak window (frames 641–656)

```
frame                    62.22 ms  37.33 Mcyc   100%
  pov_preserve_half      136.8 us   82.1 kcyc     0%
  cm_draw_trail          22.69 ms  13.61 Mcyc    36%
    filter_blend          1.02 ms  610.3 kcyc     2%  x9417  65 cyc/c
  cm_wipe_rebake          2.15 ms   1.29 Mcyc     3%
  cm_timeline_step       113.6 us   68.2 kcyc     0%
  canvas_clear            84.5 us   50.7 kcyc     0%
  canvas_buffer_wait     37.04 ms  22.23 Mcyc    60%
```

Wall min/avg/max = 49.98/62.22/73.73 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `cm_draw_trail` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `cm_draw_trail` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 5 | — | 🟢 31.92 | 0/320 | 22.95 | 18/20 |
| 11 | — | 🟢 30.90 | 0/320 | 24.44 | 18/20 |
| 12 | — | 🟢 30.58 | 0/320 | 23.28 | 18/20 |
| 1 | — | 🟢 29.95 | 0/478 | 21.89 | 27/30 |
| 10 | — | 🟢 29.89 | 0/320 | 25.36 | 18/20 |
| 8 | — | 🟢 29.71 | 0/320 | 24.01 | 18/20 |
| 4 | — | 🟢 29.34 | 0/320 | 24.06 | 18/20 |
| 9 | — | 🟢 28.81 | 0/320 | 24.84 | 18/20 |
| 2 | — | 🟢 27.89 | 0/449 | 22.95 | 26/28 |
| 6 | — | 🟢 22.69 | 0/320 | 17.19 | 18/20 |
| 3 | — | 🟢 22.44 | 0/320 | 17.36 | 18/20 |
| 7 | — | 🟢 17.05 | 0/320 | 12.39 | 18/20 |

### Per-pixel figures

`filter_blend` ran 9,417 times per frame in the peak window at 65 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1148/frame  min/avg/max 0.6/1.6/17.3 us  cpu 3.02%
isr_pack         144/frame  min/avg/max 6.2/6.7/9.4 us  cpu 1.55%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/7.9 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `cm_draw_trail` — 36% of the peak window, 22.69 ms/f.
2. `cm_wipe_rebake` — 3% of the peak window, 2.15 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `cm_timeline_step` — 0% of the peak window, 0.11 ms/f.

README cells: peak 🟢 31.92 (12), spilled 🟢 0/4127 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak +1.26 ms) is not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Comets`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
