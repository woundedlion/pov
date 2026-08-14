# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-13, **-O3**)

Global-O3 twin of the [shipping selective-O3 capture](../shipping/profile_shaderball_teensy_2026-08-13.md).
Raw capture: `build/prof/shaderball_o3.log`. Both images were built from the
same clean, pinned source at `6f2e7e77`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3`: global `-O3 -ffast-math` compiler reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, tip `6f2e7e77` |
| Method | `HS_PROFILE`, 16-frame windows, 140 s fast cycle, two-frame transitions, epoch stretched to 1,400 revs |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 140 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

Image size: `FLASH: code:128,608, data:171,952, headers:8,688` / `RAM1:
variables:323,200, code:38,616, padding:26,920, free:135,552` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 545–560 root cycles divided by 600 MHz
match measured wall within **0.1 ppm**. The parser validates 90 windows, all 17
presets, the 16→0 wrap, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate:** peak frame render is **244.46 ms**, with **510/1,440
spilled frames (35.4%)**. Five preset buckets are green and twelve are red.
The global-O3 reference is materially slower than shipping for this bank.

## Phase-by-phase readout

Fast-cycle mode gives each preset held windows separated by a two-frame
through-clear transition.

### Worst clean hold: preset 16 (window frames 545–560)

```text
frame                 249.89 ms  149.93 Mcyc  100%
  sb_shader_draw      236.66 ms  142.00 Mcyc   94%  x1.0
  sb_timeline_step     45.75 us   27.47 kcyc    0%
  canvas_clear         86.00 us   51.60 kcyc    0%
  canvas_buffer_wait    6.51 ms    3.91 Mcyc    2%
```

Wall min/avg/max = 248.83/249.89/250.53 ms. Preset 16 falls to roughly four
frames per second; global O3 makes its dominant shader 1.77× slower than the
shipping image.

### Transition into preset 16 (window frames 529–544)

```text
frame                  74.66 ms   44.79 Mcyc  100%
  sb_shader_draw       49.87 ms   29.92 Mcyc   66%  x0.94
  sb_timeline_step     30.75 us   18.47 kcyc    0%
  canvas_clear         87.31 us   52.42 kcyc    0%
  canvas_buffer_wait   18.24 ms   10.95 Mcyc   24%
```

Wall min/avg/max = 21.01/74.66/299.64 ms. The window crosses the clear into
preset 16; its peak render is 243.28 ms even though the window mean is lower.

### Per-preset table

Rows are ranked by maximum clean-hold shader cost. `Peak render` includes the
transition owned by that preset. The capture contains multiple 16→0 wraps and
2–4 clean windows per preset.

| # | Authored topology | shader ms | peak render ms | spilled/frames |
|---:|---|--:|--:|--:|
| 16 | Gnomonic dodecahedral vector-noise + mirror grid | 236.66 | 244.46 | 64/68 |
| 14 | Sinusoidal curl lattice, coarse field | 167.81 | 175.39 | 68/68 |
| 13 | Sinusoidal curl lattice, fine field | 167.74 | 179.54 | 66/68 |
| 1 | Kaleidoscope twin-wave + inner mirror | 104.00 | 121.66 | 96/102 |
| 8 | Peirce dodecahedral grid | 68.83 | 81.25 | 77/83 |
| 11 | Gnomonic dodecahedral wave/mirror grid | 63.41 | 71.42 | 64/68 |
| 12 | Gnomonic affine lattice contour | 57.51 | 72.38 | 64/68 |
| 7 | Dodecahedral grid + mirror | 53.88 | 61.23 | 0/102 |
| 10 | Dodecahedral noise lattice + mirror | 53.52 | 66.30 | 2/68 |
| 6 | Kaleidoscope noise grid + edge fade | 50.32 | 57.44 | 0/102 |
| 2 | Gnomonic kaleidoscope grid + mirror | 48.53 | 119.00 | 3/102 |
| 9 | Dodecahedral noise grid | 48.33 | 75.31 | 2/68 |
| 0 | Glitch wave-shear grid | 46.48 | 243.97 | 2/99 |
| 3 | Gnomonic glitch grid + mirror | 44.42 | 58.10 | 0/102 |
| 15 | Stereographic prism polar-wave lattice | 40.44 | 175.24 | 2/68 |
| 4 | Bonne kaleidoscope lattice + mirror | 38.03 | 51.41 | 0/102 |
| 5 | Peirce kaleidoscope lattice | 27.64 | 44.63 | 0/102 |

### Per-pixel figures

The worst held shader consumes about **13,695 cycles/sample** over the 10,368
quadrant samples. Per-pixel stage counters are disabled in the timing image.

## Column-ISR / DMA marshaling cost

```text
isr_wake       4608/frame  min/avg/max 0.62/1.91/27.09 us  cpu 3.52%
isr_pack        576/frame  min/avg/max 6.32/7.29/9.71 us   cpu 1.68%
isr_dma_submit  576/frame  min/avg/max 0.62/0.91/9.86 us   cpu 0.21%
```

- Packing dominates submit CPU time; LED transfer remains asynchronous.
- ISR CPU share is 5.41% in the worst held window.
- Preset 16 needs roughly a 3.91× peak-render speedup to fit 62.5 ms.

## Summary ranking

1. `sb_shader_draw` — 94% of the worst held frame, 236.66 ms.
2. Column wake/pack/submit ISRs — 5.41% CPU share under the live driver.
3. `canvas_buffer_wait` — 2%, 6.51 ms of cadence round-up.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` remains parented under the first scope that enters it.
- Per-pixel stage scopes are disabled in this timing image.
- Global O3 is a single-effect compiler reference, not a shippable roster image.
- Fast-cycle mode compresses dwell only; it does not change per-frame work.
- The capture used a clean, pinned worktree at `6f2e7e77`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`; `profile_one.sh` performs locked build, flash, capture,
marker validation, and artifact attestation.

## Global -O3 vs selective -O3

Global O3 peaks at **244.46 ms**, **97.11 ms (65.9%) slower** than shipping's
147.35 ms peak. Its single-effect image adds **15,336 B flash code** and
**11,536 B ITCM** (`128,608/38,616 B` versus `113,272/27,080 B`).

