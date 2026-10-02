# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/sphericalharmonics_ship.log`, captured 2026-09-28 19:00 on COM3.
Replaces `profile_sphericalharmonics_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 220 s capture, `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048` |
| Reproduce | `bash tools/profile_one.sh SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` |

Image size (`profile` env, this effect only): `FLASH: code:40856, data:148128, headers:8648` / `RAM1: variables:314944, code:16056, padding:16712, free:176576` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2449–2464 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.0 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `sh_rasterize` averages 9.49 ms/f; its worst window is 12.59 ms/f (frames 2449–2464). Peak frame render is **12.89 ms** (frame 914), and **0/3487** frames spilled. Setup frame 1 is excluded from both; it rendered 17.48 ms.

The previous shipping report (2026-08-26 02:00) recorded peak 🟢 12.64 (24) and spilled 🟢 0/3488 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 24 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 913–928)

```
frame                     62.45 ms   37.47 Mcyc   100%
  pov_preserve_half       144.6 us    86.7 kcyc     0%
  sh_rasterize            12.58 ms    7.55 Mcyc    20%
  sh_timeline_step         23.8 us    14.3 kcyc     0%
  canvas_clear             84.1 us    50.5 kcyc     0%
  canvas_buffer_wait      49.61 ms   29.77 Mcyc    79%
```

Wall min/avg/max = 62.15/62.45/62.60 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `sh_rasterize` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `sh_rasterize` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 21 | — | 🟢 12.89 | 0/128 | 12.59 | 8/8 |
| 20 | — | 🟢 12.89 | 0/128 | 12.59 | 8/8 |
| 19 | — | 🟢 11.74 | 0/128 | 11.41 | 8/8 |
| 22 | — | 🟢 11.72 | 0/128 | 11.40 | 8/8 |
| 12 | — | 🟢 11.23 | 0/192 | 10.92 | 12/12 |
| 13 | — | 🟢 11.23 | 0/161 | 10.92 | 10/10 |
| 18 | — | 🟢 10.33 | 0/128 | 10.02 | 8/8 |
| 23 | — | 🟢 10.30 | 0/128 | 10.01 | 8/8 |
| 11 | — | 🟢 9.93 | 0/192 | 9.61 | 12/12 |
| 14 | — | 🟢 9.93 | 0/128 | 9.62 | 8/8 |
| 24 | — | 🟢 9.57 | 0/128 | 9.27 | 8/8 |
| 17 | — | 🟢 9.56 | 0/128 | 9.26 | 8/8 |
| 7 | — | 🟢 9.44 | 0/190 | 9.49 | 12/12 |
| 6 | — | 🟢 9.44 | 0/128 | 9.12 | 8/8 |
| 16 | — | 🟢 9.27 | 0/128 | 8.94 | 8/8 |
| 15 | — | 🟢 9.18 | 0/128 | 8.86 | 8/8 |
| 10 | — | 🟢 9.17 | 0/192 | 8.86 | 12/12 |
| 9 | — | 🟢 8.87 | 0/192 | 8.57 | 12/12 |
| 1 | — | 🟢 8.84 | 0/128 | 8.52 | 8/8 |
| 8 | — | 🟢 8.79 | 0/192 | 8.47 | 12/12 |
| 5 | — | 🟢 8.76 | 0/128 | 8.45 | 8/8 |
| 4 | — | 🟢 8.46 | 0/128 | 8.16 | 8/8 |
| 3 | — | 🟢 8.32 | 0/128 | 7.98 | 8/8 |
| 2 | — | 🟢 8.30 | 0/128 | 7.97 | 8/8 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.6/1.6/10.9 us  cpu 3.00%
isr_pack          144/frame  min/avg/max 6.3/6.8/9.1 us  cpu 1.56%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `sh_rasterize` — 20% of the peak window, 12.58 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `sh_timeline_step` — 0% of the peak window, 0.02 ms/f.

README cells: peak 🟢 12.89 (24), spilled 🟢 0/3487 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=SphericalHarmonics`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` builds, flashes and captures under the device lock.
