# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopesmooth_ship.log`, captured 2026-10-06 18:55 on COM4.
Replaces `profile_kaleidoscopesmooth_teensy_2026-09-28.md`; the architecture snapshot [profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md](profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md) is kept.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeSmooth profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:72236, data:156900, headers:8428` / `RAM1: variables:315040, code:20488, padding:12280, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 3953–3968 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 22.86 ms/f; its worst window is 27.79 ms/f (frames 3953–3968). Peak frame render is **32.02 ms** (frame 3983), and **0/4127** frames spilled. Setup frame 1 is excluded from both; it rendered 49.23 ms.

The most recent prior capture, the architecture snapshot (2026-10-01 11:05), recorded peak 🟢 32.230 (4) and spilled 🟢 0/4136 (0.00%). The previous un-suffixed shipping report (2026-09-28 18:59) recorded peak 🟢 32.93 (4) and spilled 🟢 0/4127 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 4 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 4 once and wraps back to entry 1. The block below is the window holding the pass's peak frame.

### Peak window (frames 3969–3984)

```
frame                     62.40 ms   37.44 Mcyc   100%
  pov_preserve_half       142.2 us    85.4 kcyc     0%
  fx_shader_draw          27.78 ms   16.67 Mcyc    45%
  fx_prepare_frame        851.5 us   510.9 kcyc     1%
  fx_advance               2.15 ms    1.29 Mcyc     3%
  fx_timeline_step        165.8 us    99.5 kcyc     0%
  canvas_clear             86.4 us    51.9 kcyc     0%
  canvas_buffer_wait      31.19 ms   18.72 Mcyc    50%
```

Wall min/avg/max = 60.79/62.40/64.03 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 32.02 | 0/890 | 27.79 | 55/56 |
| 2 | — | 🟢 31.23 | 0/1079 | 27.14 | 66/67 |
| 4 | — | 🟢 28.02 | 0/1079 | 23.78 | 66/67 |
| 3 | — | 🟢 27.36 | 0/1079 | 23.86 | 67/68 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/11.8 us  cpu 3.06%
isr_pack         144/frame  min/avg/max 6.2/6.8/9.5 us  cpu 1.57%
isr_dma_submit   144/frame  min/avg/max 0.6/1.0/2.2 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 45% of the peak window, 27.78 ms/f.
2. `fx_advance` — 3% of the peak window, 2.15 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.85 ms/f.
4. `fx_timeline_step` — 0% of the peak window, 0.17 ms/f.

README cells: peak 🟢 32.02 (4), spilled 🟢 0/4127 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the architecture snapshot ran at source `0156d0d74` plus an unretained source patch, the previous un-suffixed report at `97eb0bf78`, and the deltas (peak -0.21 ms and -0.91 ms) are not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeSmooth`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeSmooth profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
