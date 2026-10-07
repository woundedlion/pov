# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/hopffibration_ship.log`, captured 2026-10-06 18:03 on COM4.
Replaces `profile_hopffibration_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh HopfFibration profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:63860, data:149652, headers:8692` / `RAM1: variables:315136, code:37048, padding:28488, free:143616` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 193–224 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.9 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `hf_render_trails` averages 24.36 ms/f; its worst window is 40.82 ms/f (frames 193–224). Peak frame render is **48.61 ms** (frame 216), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 0.44 ms.

The previous shipping report (2026-09-28 19:02) recorded peak 🟢 51.52 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 193–224)

```
frame                     62.84 ms   37.70 Mcyc   100%
  pov_preserve_half       140.0 us    84.0 kcyc     0%
  hf_render_trails        40.82 ms   24.49 Mcyc    65%
    hf_trail_raster       32.37 ms   19.42 Mcyc    52%  x149.1  217 us/c
      filter_blend         5.53 ms    3.32 Mcyc     9%  x48903  68 cyc/c
    hf_trail_gate          7.62 ms    4.57 Mcyc    12%  x210.0  36.3 us/c
  hf_project_record       184.5 us   110.7 kcyc     0%
  hf_advance_tumble         0.5 us      331 cyc     0%
  hf_timeline_step         15.7 us     9.4 kcyc     0%
  canvas_clear             84.6 us    50.8 kcyc     0%
  canvas_buffer_wait      21.59 ms   12.96 Mcyc    34%
```

Wall min/avg/max = 43.32/62.84/80.85 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 48,903 times per frame in the peak window at 68 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1159/frame  min/avg/max 0.6/1.7/15.0 us  cpu 3.19%
isr_pack         145/frame  min/avg/max 6.3/7.1/9.8 us  cpu 1.63%
isr_dma_submit   145/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `hf_render_trails` — 65% of the peak window, 40.82 ms/f.
2. `hf_project_record` — 0% of the peak window, 0.18 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 48.61, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HopfFibration`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh HopfFibration profile 70 32` builds, flashes and captures under the device lock.
