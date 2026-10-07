# LatticeMelt on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/latticemelt_ship.log`, captured 2026-10-06 18:55 on COM3.
Replaces `profile_latticemelt_teensy_2026-09-28.md`; the architecture snapshot [profile_latticemelt_architecture_teensy_2026-10-01.md](profile_latticemelt_architecture_teensy_2026-10-01.md) is kept.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | LatticeMelt 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture, `-D HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh LatticeMelt profile 110 16 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size (`profile` env, this effect only): `FLASH: code:73148, data:156540, headers:8900` / `RAM1: variables:315040, code:20472, padding:12296, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–16 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 31.94 ms/f; its worst window is 32.08 ms/f (frames 385–400). Peak frame render is **35.83 ms** (frame 385), and **0/1727** frames spilled. Setup frame 1 is excluded from both; it rendered 67.42 ms, over one display window.

The most recent prior capture, the architecture snapshot (2026-10-01 09:38), recorded peak 🟢 36.583 (2) and spilled 🟢 0/1736 (0.00%). The previous un-suffixed shipping report (2026-09-28 18:45) recorded peak 🟢 37.18 (2) and spilled 🟢 0/1727 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 2 once and wraps back to entry 1. The block below is the window holding the pass's peak frame.

### Peak window (frames 385–400)

```
frame                     62.44 ms   37.46 Mcyc   100%
  pov_preserve_half       140.6 us    84.4 kcyc     0%
  fx_shader_draw          32.08 ms   19.25 Mcyc    51%
  fx_prepare_frame         1.19 ms   713.5 kcyc     2%
  fx_advance               2.11 ms    1.27 Mcyc     3%
  fx_timeline_step         67.2 us    40.4 kcyc     0%
  canvas_clear            106.9 us    64.2 kcyc     0%
  canvas_buffer_wait      26.69 ms   16.02 Mcyc    43%
```

Wall min/avg/max = 62.27/62.44/62.65 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 35.83 | 0/648 | 33.67 | 40/41 |
| 2 | — | 🟢 35.77 | 0/1079 | 31.98 | 66/67 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/12.1 us  cpu 3.09%
isr_pack         144/frame  min/avg/max 6.3/7.0/9.7 us  cpu 1.61%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/4.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 51% of the peak window, 32.08 ms/f.
2. `fx_advance` — 3% of the peak window, 2.11 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.19 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 35.83 (2), spilled 🟢 0/1727 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the architecture snapshot ran at source `d25dd85de` and the previous un-suffixed report at `97eb0bf78`, and the deltas (peak -0.75 ms and -1.35 ms) are not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=LatticeMelt`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh LatticeMelt profile 110 16 "-D HS_PROFILE_EPOCH_REVS=1200"` builds, flashes and captures under the device lock.
