# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/alienbrain_ship.log`, captured 2026-10-06 18:35 on COM3.
Replaces `profile_alienbrain_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienBrain 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 300 s capture, `-D HS_PROFILE_EPOCH_REVS=2600` |
| Reproduce | `bash tools/profile_one.sh AlienBrain profile 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"` |

Image size (`profile` env, this effect only): `FLASH: code:71964, data:156740, headers:8860` / `RAM1: variables:315040, code:20488, padding:12280, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–16 root counter cyc ÷ 600 MHz matches the measured wall sum within **5.3 ppm** (`tools/parse_profile.py ... validate`, INVALID on this check alone: its 5 ppm gate reads only the richest window, here the one holding setup frame 1). Across all 298 windows the same check reads median 2.4 ppm, p90 4.6, max 13.6, with 21 windows above 5 ppm — the same distribution as today's VALID Comets (median 2.3, p90 4.4, max 15.4, 13/258 above 5), MeshFeedback (2.3/4.1/14.8) and DreamBalls (2.2/4.4/14.1). Every other validate check passes.

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 20.25 ms/f; its worst window is 20.48 ms/f (frames 3185–3200). Peak frame render is **25.11 ms** (frame 2076), and **0/4767** frames spilled. Setup frame 1 is excluded from both; it rendered 45.52 ms.

The previous shipping report (2026-09-28 18:30) recorded peak 🟢 25.89 (4) and spilled 🟢 0/4767 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 4 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 4 once and wraps back to entry 1. The block below is the window holding the pass's peak frame.

### Peak window (frames 2065–2080)

```
frame                     62.42 ms   37.45 Mcyc   100%
  pov_preserve_half       141.6 us    85.0 kcyc     0%
  fx_shader_draw          20.46 ms   12.27 Mcyc    33%
  fx_prepare_frame         1.71 ms    1.02 Mcyc     3%
  fx_advance               2.37 ms    1.42 Mcyc     4%
  fx_timeline_step        163.4 us    98.1 kcyc     0%
  canvas_clear             86.8 us    52.1 kcyc     0%
  canvas_buffer_wait      37.45 ms   22.47 Mcyc    60%
```

Wall min/avg/max = 62.18/62.41/62.56 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 3 | — | 🟢 25.11 | 0/1079 | 20.46 | 67/68 |
| 4 | — | 🟢 24.95 | 0/1079 | 20.48 | 66/67 |
| 1 | — | 🟢 24.89 | 0/1530 | 21.60 | 95/96 |
| 2 | — | 🟢 24.86 | 0/1079 | 20.44 | 66/67 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.5/1.7/15.8 us  cpu 3.05%
isr_pack         144/frame  min/avg/max 6.2/6.8/9.6 us  cpu 1.57%
isr_dma_submit   144/frame  min/avg/max 0.6/1.0/8.5 us  cpu 0.22%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 33% of the peak window, 20.46 ms/f.
2. `fx_advance` — 4% of the peak window, 2.37 ms/f.
3. `fx_prepare_frame` — 3% of the peak window, 1.71 ms/f.
4. `fx_timeline_step` — 0% of the peak window, 0.16 ms/f.

README cells: peak 🟢 25.11 (4), spilled 🟢 0/4767 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -0.78 ms) is not attributed here.
- `validate` reports INVALID solely on the exactness check: it tests only the richest window (frames 1–16, holding setup frame 1), which reads 5.3 ppm against the 5 ppm gate. The all-window distribution (median 2.4 ppm, p90 4.6, max 13.6, 21/298 windows above 5 ppm) matches today's VALID cycler captures, and wrap, all 4 presets and markers pass, so the capture is treated as sound.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienBrain`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh AlienBrain profile 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"` builds, flashes and captures under the device lock.
