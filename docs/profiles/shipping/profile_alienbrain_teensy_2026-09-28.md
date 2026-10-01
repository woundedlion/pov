# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienBrain`).
Raw capture: `build/prof/alienbrain_ship.log`, captured 2026-09-28 18:30 on COM4.
Replaces `profile_alienbrain_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienBrain 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 300 s capture, `-D HS_PROFILE_EPOCH_REVS=2600` |
| Reproduce | `bash tools/profile_one.sh AlienBrain profile 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"` |

Image size (`profile` env, this effect only): `FLASH: code:69384, data:154708, headers:8356` / `RAM1: variables:315040, code:19736, padding:13032, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–16 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.0 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 21.84 ms/f; its worst window is 22.13 ms/f (frames 3553–3568). Peak frame render is **25.89 ms** (frame 2044), and **0/4767** frames spilled. Setup frame 1 is excluded from both; it rendered 48.74 ms.

The previous shipping report (2026-08-26 02:27) recorded peak 🟢 31.61 (5) and spilled 🟢 0/4768 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 4 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2033–2048)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       144.9 us    87.0 kcyc     0%
  fx_shader_draw          22.10 ms   13.26 Mcyc    35%
  fx_prepare_frame         1.06 ms   634.2 kcyc     2%
  fx_advance               2.20 ms    1.32 Mcyc     4%
  fx_timeline_step        133.1 us    79.9 kcyc     0%
  canvas_clear             86.1 us    51.6 kcyc     0%
  canvas_buffer_wait      36.70 ms   22.02 Mcyc    59%
```

Wall min/avg/max = 62.13/62.44/62.59 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 3 | — | 🟢 25.89 | 0/1079 | 22.13 | 68/68 |
| 2 | — | 🟢 25.89 | 0/1079 | 22.09 | 67/67 |
| 4 | — | 🟢 25.86 | 0/1079 | 22.13 | 67/67 |
| 1 | — | 🟢 25.80 | 0/1530 | 23.34 | 96/96 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.6/16.4 us  cpu 3.02%
isr_pack          144/frame  min/avg/max 6.2/6.8/9.3 us  cpu 1.57%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/6.9 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 35% of the peak window, 22.10 ms/f.
2. `fx_advance` — 4% of the peak window, 2.20 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.06 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 25.89 (4), spilled 🟢 0/4767 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienBrain`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh AlienBrain profile 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"` builds, flashes and captures under the device lock.
