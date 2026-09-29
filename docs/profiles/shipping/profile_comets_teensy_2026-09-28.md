# Comets on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Comets`).
Raw capture: `build/prof/comets_ship.log`, captured 2026-09-28 18:28 on COM3.
Replaces `profile_comets_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Comets 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:63920, data:151192, headers:9144` / `RAM1: variables:315040, code:29224, padding:3544, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1537–1552 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.9 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `cm_draw_trail` averages 16.20 ms/f; its worst window is 24.42 ms/f (frames 1537–1552). Peak frame render is **30.66 ms** (frame 645), and **0/4127** frames spilled. Setup frame 1 is excluded from both; it rendered 0.67 ms.

The previous shipping report (2026-08-26 02:05) recorded peak 🟢 33.91 (13) and spilled 🟢 0/4128 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 641–656)

```
frame                     62.24 ms   37.35 Mcyc   100%
  pov_preserve_half       134.9 us    81.0 kcyc     0%
  cm_draw_trail           21.85 ms   13.11 Mcyc    35%
    filter_blend          966.3 us   579.8 kcyc     2%  x9417  62 cyc/c
  cm_wipe_rebake           2.18 ms    1.31 Mcyc     4%
  cm_timeline_step        123.8 us    74.3 kcyc     0%
  canvas_clear             84.1 us    50.5 kcyc     0%
  canvas_buffer_wait      37.87 ms   22.72 Mcyc    61%
```

Wall min/avg/max = 50.51/62.24/72.91 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `cm_draw_trail` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `cm_draw_trail` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 5 | — | 🟢 30.66 | 0/320 | 22.06 | 20/20 |
| 11 | — | 🟢 29.73 | 0/320 | 23.53 | 20/20 |
| 12 | — | 🟢 29.34 | 0/320 | 22.42 | 20/20 |
| 10 | — | 🟢 29.02 | 0/320 | 24.42 | 20/20 |
| 1 | — | 🟢 28.85 | 0/478 | 21.05 | 30/30 |
| 8 | — | 🟢 28.54 | 0/320 | 23.08 | 20/20 |
| 4 | — | 🟢 28.17 | 0/320 | 23.25 | 20/20 |
| 9 | — | 🟢 27.66 | 0/320 | 23.87 | 20/20 |
| 2 | — | 🟢 26.73 | 0/449 | 22.05 | 28/28 |
| 6 | — | 🟢 22.04 | 0/320 | 16.64 | 20/20 |
| 3 | — | 🟢 21.78 | 0/320 | 17.80 | 20/20 |
| 7 | — | 🟢 16.48 | 0/320 | 12.21 | 20/20 |

### Per-pixel figures

`filter_blend` ran 9,417 times per frame in the peak window at 62 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1148/frame  min/avg/max 0.6/1.6/15.7 us  cpu 3.03%
isr_pack          144/frame  min/avg/max 6.2/6.7/9.5 us  cpu 1.54%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/8.4 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `cm_draw_trail` — 35% of the peak window, 21.85 ms/f.
2. `cm_wipe_rebake` — 4% of the peak window, 2.18 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.13 ms/f.
4. `cm_timeline_step` — 0% of the peak window, 0.12 ms/f.

README cells: peak 🟢 30.66 (12), spilled 🟢 0/4127 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Comets`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
