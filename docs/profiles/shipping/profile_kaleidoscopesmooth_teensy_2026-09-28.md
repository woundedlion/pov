# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeSmooth`).
Raw capture: `build/prof/kaleidoscopesmooth_ship.log`, captured 2026-09-28 18:59 on COM4.
Replaces `profile_kaleidoscopesmooth_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeSmooth profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:69592, data:154860, headers:9020` / `RAM1: variables:315040, code:19736, padding:13032, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 3969–3984 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 24.20 ms/f; its worst window is 28.89 ms/f (frames 3969–3984). Peak frame render is **32.93 ms** (frame 3979), and **0/4127** frames spilled. Setup frame 1 is excluded from both; it rendered 52.89 ms.

The previous shipping report (2026-08-26 02:59) recorded peak 🟢 35.85 (5) and spilled 🟢 0/4128 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 4 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 3969–3984)

```
frame                     62.42 ms   37.45 Mcyc   100%
  pov_preserve_half       139.1 us    83.5 kcyc     0%
  fx_shader_draw          28.89 ms   17.33 Mcyc    46%
  fx_prepare_frame        834.1 us   500.5 kcyc     1%
  fx_advance               2.09 ms    1.25 Mcyc     3%
  fx_timeline_step        131.8 us    79.1 kcyc     0%
  canvas_clear             86.6 us    51.9 kcyc     0%
  canvas_buffer_wait      30.24 ms   18.14 Mcyc    48%
```

Wall min/avg/max = 60.97/62.42/63.88 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 32.93 | 0/890 | 28.89 | 56/56 |
| 2 | — | 🟢 32.81 | 0/1079 | 28.75 | 67/67 |
| 4 | — | 🟢 29.26 | 0/1079 | 25.13 | 67/67 |
| 3 | — | 🟢 28.62 | 0/1079 | 25.19 | 68/68 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.6/1.6/11.4 us  cpu 3.04%
isr_pack          144/frame  min/avg/max 6.2/6.8/9.5 us  cpu 1.57%
isr_dma_submit    144/frame  min/avg/max 0.7/0.9/2.1 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 46% of the peak window, 28.89 ms/f.
2. `fx_advance` — 3% of the peak window, 2.09 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.83 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 32.93 (4), spilled 🟢 0/4127 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeSmooth`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeSmooth profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
