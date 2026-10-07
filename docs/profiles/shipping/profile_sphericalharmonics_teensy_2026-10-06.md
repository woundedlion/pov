# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/sphericalharmonics_ship.log`, captured 2026-10-06 18:42 on COM4.
Replaces `profile_sphericalharmonics_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 220 s capture, `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048`. The modes morph back to back with no clean hold, so each per-mode row is the peak window nearest that mode's anchor |
| Reproduce | `bash tools/profile_one.sh SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` |

Image size (`profile` env, this effect only): `FLASH: code:41100, data:149268, headers:8284` / `RAM1: variables:314944, code:16072, padding:16696, free:176576` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2449–2464 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `sh_rasterize` averages 9.54 ms/f; its worst window is 12.62 ms/f (frames 2449–2464). Peak frame render is **12.92 ms** (frame 2446), and **0/3487** frames spilled. Setup frame 1 is excluded from both; it rendered 17.51 ms.

The previous shipping report (2026-09-28 19:00) recorded peak 🟢 12.89 (24) and spilled 🟢 0/3487 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 24 modes in ordered cycle; each owns its morph and the one that follows it. The capture starts on mode 7, runs through 24, wraps to mode 1 and back to mode 7 twice, and ends on mode 13. The block below is the window holding the pass's peak frame (mode 21).

### Peak window (frames 2433–2448)

```
frame                   62.65 ms  37.59 Mcyc   100%
  pov_preserve_half     146.3 us   87.8 kcyc     0%
  sh_rasterize          12.54 ms   7.52 Mcyc    20%
  sh_timeline_step       22.7 us   13.6 kcyc     0%
  canvas_clear           84.4 us   50.7 kcyc     0%
  canvas_buffer_wait    49.85 ms  29.91 Mcyc    80%
```

Wall min/avg/max = 62.35/62.65/65.19 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded (it belongs to mode 7); the `sh_rasterize` column is from `parse_profile.py ... presets`, the costliest modal-call-count window of each mode with windows straddling an advance excluded, which for these back-to-back morphs is the peak window nearest the mode's anchor. The cycle wraps back to its first mode within the capture (validate: cycle returns to its first index).

| Mode | Meta | Peak render ms | Spilled/frames | Anchor-window `sh_rasterize` ms/f | Windows used |
|---|---|--:|--:|--:|--:|
| 21 | — | 🟢 12.92 | 0/128 | 12.62 | 6/8 |
| 20 | — | 🟢 12.92 | 0/128 | 12.61 | 6/8 |
| 19 | — | 🟢 11.79 | 0/128 | 11.44 | 6/8 |
| 22 | — | 🟢 11.77 | 0/128 | 11.43 | 6/8 |
| 12 | — | 🟢 11.26 | 0/192 | 10.94 | 9/12 |
| 13 | — | 🟢 11.26 | 0/161 | 10.95 | 8/10 |
| 18 | — | 🟢 10.39 | 0/128 | 10.07 | 6/8 |
| 23 | — | 🟢 10.38 | 0/128 | 10.09 | 6/8 |
| 14 | — | 🟢 10.03 | 0/128 | 9.70 | 6/8 |
| 11 | — | 🟢 10.01 | 0/192 | 9.66 | 9/12 |
| 24 | — | 🟢 9.64 | 0/128 | 9.34 | 6/8 |
| 17 | — | 🟢 9.63 | 0/128 | 9.33 | 6/8 |
| 7 | — | 🟢 9.50 | 0/190 | 9.52 | 9/12 |
| 6 | — | 🟢 9.49 | 0/128 | 9.15 | 6/8 |
| 16 | — | 🟢 9.35 | 0/128 | 9.01 | 6/8 |
| 15 | — | 🟢 9.25 | 0/128 | 8.92 | 6/8 |
| 10 | — | 🟢 9.24 | 0/192 | 8.93 | 9/12 |
| 9 | — | 🟢 8.95 | 0/192 | 8.64 | 9/12 |
| 1 | — | 🟢 8.92 | 0/128 | 8.59 | 6/8 |
| 8 | — | 🟢 8.87 | 0/192 | 8.54 | 9/12 |
| 5 | — | 🟢 8.83 | 0/128 | 8.51 | 6/8 |
| 4 | — | 🟢 8.54 | 0/128 | 8.23 | 6/8 |
| 3 | — | 🟢 8.38 | 0/128 | 8.03 | 6/8 |
| 2 | — | 🟢 8.36 | 0/128 | 8.02 | 6/8 |

Mode 7's window column (9.52) sits above its frame-1-excluded peak (9.50) because its first window holds the 17.51 ms setup frame; the window average is not comparable with the per-frame peak there.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1156/frame  min/avg/max 0.6/1.6/11.0 us  cpu 2.99%
isr_pack         144/frame  min/avg/max 6.2/6.8/9.2 us  cpu 1.56%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/1.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `sh_rasterize` — 20% of the peak window, 12.54 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `sh_timeline_step` — 0% of the peak window, 0.02 ms/f.

README cells: peak 🟢 12.92 (24), spilled 🟢 0/3487 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak +0.03 ms) is not attributed here.
- Dwell-compression knobs (ordered cycle, epoch stretch) change how long a mode holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=SphericalHarmonics`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` builds, flashes and captures under the device lock.
