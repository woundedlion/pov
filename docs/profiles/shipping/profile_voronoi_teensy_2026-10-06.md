# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/voronoi_ship.log`, captured 2026-10-06 18:25 on COM4.
Replaces `profile_voronoi_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Voronoi 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Voronoi profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:42716, data:148604, headers:8356` / `RAM1: variables:314976, code:13592, padding:19176, free:176544` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `vo_shade` averages 7.10 ms/f; its worst window is 7.17 ms/f (frames 353–384). Peak frame render is **8.00 ms** (frame 32), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 13.70 ms.

The previous shipping report (2026-09-28 19:15) recorded peak 🟢 8.51 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 1–32)

```
frame                     59.24 ms   35.54 Mcyc   100%
  pov_preserve_half       144.6 us    86.8 kcyc     0%
  vo_shade                 7.28 ms    4.37 Mcyc    12%
  vo_kdtree               354.3 us   212.6 kcyc     1%
  vo_animate               50.6 us    30.4 kcyc     0%
  canvas_clear             88.2 us    53.0 kcyc     0%
  canvas_buffer_wait      51.32 ms   30.79 Mcyc    87%
```

Wall min/avg/max = 7.68/59.24/62.80 ms. This window also holds setup frame 1 and frame 2, the first two frames run before display sync starts (wall = render), so its averages and wall minimum include them. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1373/frame  min/avg/max 0.5/1.4/13.1 us  cpu 2.65%
isr_pack         136/frame  min/avg/max 6.5/6.9/9.5 us  cpu 1.26%
isr_dma_submit   136/frame  min/avg/max 0.8/0.9/1.4 us  cpu 0.16%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `vo_shade` — 12% of the peak window, 7.28 ms/f.
2. `vo_kdtree` — 1% of the peak window, 0.35 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 8.00, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Voronoi`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Voronoi profile 70 32` builds, flashes and captures under the device lock.
