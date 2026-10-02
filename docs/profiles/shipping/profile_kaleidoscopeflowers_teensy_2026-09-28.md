# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopeflowers_ship.log`, captured 2026-09-28 19:07 on COM4.
Replaces `profile_kaleidoscopeflowers_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeFlowers 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeFlowers profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:69752, data:154860, headers:8860` / `RAM1: variables:315040, code:19736, padding:13032, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 3953–3968 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.9 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 26.57 ms/f; its worst window is 28.56 ms/f (frames 3953–3968). Peak frame render is **32.52 ms** (frame 3965), and **0/4127** frames spilled. Setup frame 1 is excluded from both; it rendered 54.33 ms.

The previous shipping report (2026-08-26 03:06) recorded peak 🟢 36.69 (4) and spilled 🟢 0/4128 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 3 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 3953–3968)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       140.7 us    84.4 kcyc     0%
  fx_shader_draw          28.56 ms   17.14 Mcyc    46%
  fx_prepare_frame        788.2 us   473.0 kcyc     1%
  fx_advance               2.04 ms    1.22 Mcyc     3%
  fx_timeline_step        134.1 us    80.5 kcyc     0%
  canvas_clear             88.6 us    53.2 kcyc     0%
  canvas_buffer_wait      30.67 ms   18.40 Mcyc    49%
```

Wall min/avg/max = 61.14/62.44/63.73 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 🟢 32.52 | 0/1371 | 28.56 | 85/85 |
| 3 | — | 🟢 32.38 | 0/1079 | 28.52 | 68/68 |
| 1 | — | 🟢 32.31 | 0/1677 | 28.47 | 105/105 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.5/1.7/15.8 us  cpu 3.10%
isr_pack          144/frame  min/avg/max 6.5/7.2/9.7 us  cpu 1.65%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/5.7 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 46% of the peak window, 28.56 ms/f.
2. `fx_advance` — 3% of the peak window, 2.04 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.79 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 32.52 (3), spilled 🟢 0/4127 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeFlowers`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeFlowers profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
