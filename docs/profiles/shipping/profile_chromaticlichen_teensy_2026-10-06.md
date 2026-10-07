# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/chromaticlichen_ship.log`, captured 2026-10-06 18:57 on COM3.
Replaces `profile_chromaticlichen_teensy_2026-09-28.md`; the architecture snapshot [profile_chromaticlichen_architecture_teensy_2026-10-01.md](profile_chromaticlichen_architecture_teensy_2026-10-01.md) is kept.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ChromaticLichen 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh ChromaticLichen profile 70 32` |

Image size (`profile` env, this effect only): `FLASH: code:70580, data:155304, headers:8608` / `RAM1: variables:315040, code:20280, padding:12488, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–32 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 28.92 ms/f; its worst window is 29.23 ms/f (frames 801–832). Peak frame render is **35.53 ms** (frame 924), and **0/1087** frames spilled. Setup frame 1 is excluded from both; it rendered 61.87 ms.

The most recent prior capture, the architecture snapshot (2026-10-01 09:40), recorded peak 🟢 35.773 and spilled 🟢 0/1096 (0.00%). The previous un-suffixed shipping report (2026-09-28 18:47) recorded peak 🟢 36.63 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one steady regime. The block below is the window holding the pass's peak frame.

### Peak window (frames 897–928)

```
frame                     62.43 ms   37.46 Mcyc   100%
  pov_preserve_half       140.6 us    84.4 kcyc     0%
  fx_shader_draw          29.06 ms   17.43 Mcyc    47%
  fx_prepare_frame         3.95 ms    2.37 Mcyc     6%
  fx_advance               2.00 ms    1.20 Mcyc     3%
  fx_timeline_step         82.4 us    49.4 kcyc     0%
  canvas_clear             85.0 us    51.0 kcyc     0%
  canvas_buffer_wait      27.08 ms   16.25 Mcyc    43%
```

Wall min/avg/max = 62.16/62.43/62.73 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.7/20.0 us  cpu 3.07%
isr_pack         144/frame  min/avg/max 6.2/6.9/9.5 us  cpu 1.58%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/8.3 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 47% of the peak window, 29.06 ms/f.
2. `fx_prepare_frame` — 6% of the peak window, 3.95 ms/f.
3. `fx_advance` — 3% of the peak window, 2.00 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 35.53, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the architecture snapshot ran at source `d25dd85de` and the previous un-suffixed report at `97eb0bf78`, and the deltas (peak -0.24 ms and -1.10 ms) are not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh ChromaticLichen profile 70 32` builds, flashes and captures under the device lock.
