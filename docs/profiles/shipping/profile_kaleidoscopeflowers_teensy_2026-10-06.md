# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopeflowers_ship.log`, captured 2026-10-06 19:04 on COM4.
Replaces `profile_kaleidoscopeflowers_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeFlowers 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeFlowers profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size (`profile` env, this effect only): `FLASH: code:72300, data:156900, headers:8364` / `RAM1: variables:315040, code:20472, padding:12296, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2193–2208 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.2 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 25.22 ms/f; its worst window is 27.28 ms/f (frames 2193–2208). Peak frame render is **31.42 ms** (frame 3115), and **0/4127** frames spilled. Setup frame 1 is excluded from both; it rendered 51.85 ms.

The previous shipping report (2026-09-28 19:07) recorded peak 🟢 32.52 (3) and spilled 🟢 0/4127 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 3 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 3 once, wraps back to entry 1, and advances into entry 2 before it ends. The block below is the window holding the pass's peak frame.

### Peak window (frames 3105–3120)

```
frame                     62.42 ms   37.45 Mcyc   100%
  pov_preserve_half       146.4 us    87.8 kcyc     0%
  fx_shader_draw          27.01 ms   16.21 Mcyc    43%
  fx_prepare_frame         1.08 ms   647.9 kcyc     2%
  fx_advance               2.34 ms    1.40 Mcyc     4%
  fx_timeline_step        168.4 us   101.0 kcyc     0%
  canvas_clear             87.1 us    52.3 kcyc     0%
  canvas_buffer_wait      31.55 ms   18.93 Mcyc    51%
```

Wall min/avg/max = 61.26/62.42/63.54 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 31.42 | 0/1677 | 27.26 | 103/105 |
| 2 | — | 🟢 31.07 | 0/1371 | 26.88 | 84/85 |
| 3 | — | 🟢 30.78 | 0/1079 | 27.28 | 67/68 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/19.1 us  cpu 3.18%
isr_pack         144/frame  min/avg/max 6.4/7.4/10.1 us  cpu 1.69%
isr_dma_submit   144/frame  min/avg/max 0.8/1.0/11.0 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 43% of the peak window, 27.01 ms/f.
2. `fx_advance` — 4% of the peak window, 2.34 ms/f.
3. `fx_prepare_frame` — 2% of the peak window, 1.08 ms/f.
4. `fx_timeline_step` — 0% of the peak window, 0.17 ms/f.

README cells: peak 🟢 31.42 (3), spilled 🟢 0/4127 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them (peak -1.10 ms) is not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeFlowers`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeFlowers profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes and captures under the device lock.
