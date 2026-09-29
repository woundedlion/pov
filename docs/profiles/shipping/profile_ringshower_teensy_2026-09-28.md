# RingShower on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingShower`).
Raw capture: `build/prof/ringshower_ship.log`, captured 2026-09-28 19:11 on COM3.
Replaces `profile_ringshower_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingShower 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingShower profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:64792, data:147612, headers:8780` / `RAM1: variables:315040, code:32184, padding:584, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 481–512 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `rsh_draw_rings` averages 1.05 ms/f; its worst window is 2.11 ms/f (frames 481–512). Peak frame render is **4.33 ms** (frame 297), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 4.05 ms.

The previous shipping report (2026-08-26 01:39) recorded peak 🟢 3.98 and spilled 🟢 0/1088 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 289–320)

```
frame                     62.42 ms   37.45 Mcyc   100%
  pov_preserve_half       151.1 us    90.7 kcyc     0%
  rsh_draw_rings           1.66 ms   993.5 kcyc     3%
    rsh_ring_plot          1.65 ms   992.9 kcyc     3%  x2.2  746 us/c
      filter_blend        133.1 us    79.8 kcyc     0%  x1163  69 cyc/c
  rsh_timeline_step       105.7 us    63.4 kcyc     0%
  canvas_clear             84.2 us    50.5 kcyc     0%
  canvas_buffer_wait      60.42 ms   36.25 Mcyc    97%
```

Wall min/avg/max = 61.23/62.41/64.45 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 1,163 times per frame in the peak window at 69 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1151/frame  min/avg/max 0.5/1.6/11.8 us  cpu 2.96%
isr_pack          144/frame  min/avg/max 6.2/6.6/9.1 us  cpu 1.51%
isr_dma_submit    144/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `rsh_draw_rings` — 3% of the peak window, 1.66 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.15 ms/f.
3. `rsh_timeline_step` — 0% of the peak window, 0.11 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 4.33, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingShower`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh RingShower profile 70 32` builds, flashes and captures under the device lock.
