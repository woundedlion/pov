# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/ringspin_ship.log`, captured 2026-10-06 18:22 on COM4.
Replaces `profile_ringspin_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingSpin profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:56252, data:151352, headers:8456` / `RAM1: variables:315040, code:30152, padding:2616, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 545–576 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `rs_draw_rings` averages 30.64 ms/f; its worst window is 38.16 ms/f (frames 545–576). Peak frame render is **47.16 ms** (frame 981), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 1.82 ms.

The previous shipping report (2026-09-28 19:13) recorded peak 🟢 47.04 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 961–992)

```
frame                     62.39 ms   37.44 Mcyc   100%
  pov_preserve_half       142.9 us    85.8 kcyc     0%
  rs_draw_rings           34.79 ms   20.88 Mcyc    56%
    rs_ring_scan          34.19 ms   20.52 Mcyc    55%  x76.0  450 us/c
      filter_blend         4.46 ms    2.68 Mcyc     7%  x45185  59 cyc/c
  rs_timeline_step         73.2 us    43.9 kcyc     0%
  canvas_clear             84.4 us    50.7 kcyc     0%
  canvas_buffer_wait      27.30 ms   16.38 Mcyc    44%
```

Wall min/avg/max = 40.96/62.39/84.05 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 45,185 times per frame in the peak window at 59 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1151/frame  min/avg/max 0.6/1.7/19.8 us  cpu 3.08%
isr_pack         144/frame  min/avg/max 6.2/6.8/9.3 us  cpu 1.57%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `rs_draw_rings` — 56% of the peak window, 34.79 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.
4. `rs_timeline_step` — 0% of the peak window, 0.07 ms/f.

README cells: peak 🟢 47.16, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh RingSpin profile 70 32` builds, flashes and captures under the device lock.
