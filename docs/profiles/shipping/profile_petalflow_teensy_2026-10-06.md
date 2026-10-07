# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/petalflow_ship.log`, captured 2026-10-06 18:12 on COM4.
Replaces `profile_petalflow_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh PetalFlow profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:60924, data:150204, headers:9028` / `RAM1: variables:315072, code:33224, padding:32312, free:143680` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 929–960 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `pf_draw_rings` averages 11.89 ms/f; its worst window is 12.26 ms/f (frames 929–960). Peak frame render is **13.29 ms** (frame 307), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 15.58 ms.

The previous shipping report (2026-09-28 19:07) recorded peak 🟢 13.61 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 289–320)

```
frame                     62.48 ms   37.49 Mcyc   100%
  pov_preserve_half       143.3 us    86.0 kcyc     0%
  pf_draw_rings           12.26 ms    7.36 Mcyc    20%
    pf_ring_scan          11.14 ms    6.68 Mcyc    18%  x23.3  478 us/c
      filter_blend        794.5 us   476.7 kcyc     1%  x6959  68 cyc/c
    pf_ring_build          1.00 ms   602.3 kcyc     2%  x23.3  43.1 us/c
  pf_timeline_step         17.6 us    10.6 kcyc     0%
  canvas_clear             84.3 us    50.6 kcyc     0%
  canvas_buffer_wait      49.98 ms   29.99 Mcyc    80%
```

Wall min/avg/max = 61.29/62.48/63.62 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 6,959 times per frame in the peak window at 68 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1153/frame  min/avg/max 0.6/1.6/12.6 us  cpu 3.01%
isr_pack         144/frame  min/avg/max 6.2/6.6/9.1 us  cpu 1.52%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/3.9 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `pf_draw_rings` — 20% of the peak window, 12.26 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `pf_timeline_step` — 0% of the peak window, 0.02 ms/f.

README cells: peak 🟢 13.29, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=PetalFlow`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh PetalFlow profile 70 32` builds, flashes and captures under the device lock.
