# KaleidoscopePentBright on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopepentbright_ship.log`, captured 2026-10-06 18:44 on COM4.
Replaces `profile_kaleidoscopepentbright_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopePentBright 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopePentBright profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:68492, data:155832, headers:9144` / `RAM1: variables:315040, code:20280, padding:12488, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 417–448 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.9 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 21.81 ms/f; its worst window is 23.50 ms/f (frames 417–448). Peak frame render is **27.04 ms** (frame 881), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 48.30 ms.

The previous shipping report (2026-09-28 18:51) recorded peak 🟢 27.24 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 865–896)

```
frame                     62.43 ms   37.46 Mcyc   100%
  pov_preserve_half       144.2 us    86.5 kcyc     0%
  fx_shader_draw          22.74 ms   13.65 Mcyc    36%
  fx_prepare_frame        867.5 us   520.5 kcyc     1%
  fx_advance               2.00 ms    1.20 Mcyc     3%
  fx_timeline_step         50.7 us    30.4 kcyc     0%
  canvas_clear             84.7 us    50.8 kcyc     0%
  canvas_buffer_wait      36.52 ms   21.91 Mcyc    59%
```

Wall min/avg/max = 60.27/62.43/64.63 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.6/24.5 us  cpu 3.02%
isr_pack         144/frame  min/avg/max 6.2/6.7/9.5 us  cpu 1.54%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/4.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 36% of the peak window, 22.74 ms/f.
2. `fx_advance` — 3% of the peak window, 2.00 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.87 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 27.04, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.20 ms) is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopePentBright`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopePentBright profile 70 32` builds, flashes and captures under the device lock.
