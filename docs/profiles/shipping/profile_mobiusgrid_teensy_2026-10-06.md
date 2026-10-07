# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/mobiusgrid_ship.log`, captured 2026-10-06 18:08 on COM4.
Replaces `profile_mobiusgrid_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean detached worktree of master at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh MobiusGrid profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:72692, data:157296, headers:8600` / `RAM1: variables:315040, code:20872, padding:11896, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–16 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.2 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 18.18 ms/f; its worst window is 18.22 ms/f (frames 273–288). Peak frame render is **21.77 ms** (frame 2049), and **0/2687** frames spilled. Setup frame 1 is excluded from both; it rendered 38.12 ms.

The previous shipping report (2026-09-28 19:05) recorded peak 🟢 22.05 (2) and spilled 🟢 0/2687 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2049–2064)

```
frame                     62.41 ms   37.45 Mcyc   100%
  pov_preserve_half       143.6 us    86.1 kcyc     0%
  fx_shader_draw          18.21 ms   10.93 Mcyc    29%
  fx_prepare_frame        692.0 us   415.2 kcyc     1%
  fx_advance               2.30 ms    1.38 Mcyc     4%
  fx_timeline_step        179.4 us   107.7 kcyc     0%
  canvas_clear             86.4 us    51.9 kcyc     0%
  canvas_buffer_wait      40.75 ms   24.45 Mcyc    65%
```

Wall min/avg/max = 62.26/62.41/62.53 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is from `parse_profile.py ... presets` (modal-call-count windows of each entry). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

Entry 1's clean figure (19.23 ms/f) is the setup window 1–16, which holds setup frame 1; that is why it exceeds the pass's worst window above.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 21.77 | 0/1608 | 19.23 | 100/101 |
| 2 | — | 🟢 21.68 | 0/1079 | 18.21 | 66/67 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.6/1.7/15.9 us  cpu 3.05%
isr_pack         144/frame  min/avg/max 6.3/6.7/9.5 us  cpu 1.54%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/8.5 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 29% of the peak window, 18.21 ms/f.
2. `fx_advance` — 4% of the peak window, 2.30 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.69 ms/f.
4. `fx_timeline_step` — 0% of the peak window, 0.18 ms/f.

README cells: peak 🟢 21.77 (2), spilled 🟢 0/2687 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MobiusGrid`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MobiusGrid profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
