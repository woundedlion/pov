# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/raymarch_ship.log`, captured 2026-10-06 18:17 on COM4.
Replaces `profile_raymarch_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Raymarch 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Raymarch profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:108388, data:188996, headers:8788` / `RAM1: variables:315008, code:31352, padding:1416, free:176512` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.4 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `rm_shader_draw` averages 42.76 ms/f; its worst window is 44.17 ms/f (frames 833–864). Peak frame render is **51.95 ms** (frame 8), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 94.54 ms.

The previous shipping report (2026-09-28 19:09) recorded peak 🟢 54.87 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 1–32)

```
frame                     62.85 ms   37.71 Mcyc   100%
  pov_preserve_half       128.8 us    77.3 kcyc     0%
  rm_shader_draw          44.69 ms   26.81 Mcyc    71%
    filter_blend          456.4 us   273.8 kcyc     1%  x5692  48 cyc/c
  rm_timeline_step        383.2 us   230.0 kcyc     1%
  canvas_clear             89.4 us    53.6 kcyc     0%
  canvas_buffer_wait      14.60 ms    8.76 Mcyc    23%
```

Wall min/avg/max = 50.86/62.85/94.54 ms. This window also holds setup frame 1 (the 94.54 ms wall maximum) and frame 2, both run before display sync starts (wall = render), so its averages include them. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

`filter_blend` ran 5,692 times per frame in the peak window at 48 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1393/frame  min/avg/max 0.5/1.4/13.2 us  cpu 2.63%
isr_pack         138/frame  min/avg/max 6.2/6.8/9.6 us  cpu 1.24%
isr_dma_submit   138/frame  min/avg/max 0.7/0.9/1.7 us  cpu 0.17%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `rm_shader_draw` — 71% of the peak window, 44.69 ms/f.
2. `rm_timeline_step` — 1% of the peak window, 0.38 ms/f.
3. `pov_preserve_half` — 0% of the peak window, 0.13 ms/f.
4. `canvas_clear` — 0% of the peak window, 0.09 ms/f.

README cells: peak 🟢 51.95, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Raymarch`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Raymarch profile 70 32` builds, flashes and captures under the device lock.
