# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/mobiusgrid_ship.log`, captured 2026-09-28 19:05 on COM3.
Replaces `profile_mobiusgrid_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh MobiusGrid profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:69992, data:155160, headers:8320` / `RAM1: variables:315040, code:19720, padding:13048, free:176480` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1–16 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `fx_shader_draw` averages 18.61 ms/f; its worst window is 18.66 ms/f (frames 289–304). Peak frame render is **22.05 ms** (frame 2030), and **0/2687** frames spilled. Setup frame 1 is excluded from both; it rendered 38.85 ms.

The previous shipping report (2026-08-26 01:34) recorded peak 🟢 22.37 (3) and spilled 🟢 0/2688 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 2 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2017–2032)

```
frame                     62.44 ms   37.47 Mcyc   100%
  pov_preserve_half       141.8 us    85.1 kcyc     0%
  fx_shader_draw          18.63 ms   11.18 Mcyc    30%
  fx_prepare_frame        719.4 us   431.7 kcyc     1%
  fx_advance               2.25 ms    1.35 Mcyc     4%
  fx_timeline_step        136.3 us    81.8 kcyc     0%
  canvas_clear             86.3 us    51.8 kcyc     0%
  canvas_buffer_wait      40.46 ms   24.28 Mcyc    65%
```

Wall min/avg/max = 62.21/62.44/62.60 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `fx_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 🟢 22.05 | 0/1608 | 19.69 | 101/101 |
| 2 | — | 🟢 22.05 | 0/1079 | 18.65 | 67/67 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.6/1.6/11.4 us  cpu 2.97%
isr_pack          144/frame  min/avg/max 6.2/6.8/9.5 us  cpu 1.55%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/2.8 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 30% of the peak window, 18.63 ms/f.
2. `fx_advance` — 4% of the peak window, 2.25 ms/f.
3. `fx_prepare_frame` — 1% of the peak window, 0.72 ms/f.
4. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.

README cells: peak 🟢 22.05 (2), spilled 🟢 0/2687 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MobiusGrid`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MobiusGrid profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
