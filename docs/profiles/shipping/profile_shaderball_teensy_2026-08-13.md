# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-13, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw capture:
`build/prof/shaderball_ship.log`. This replaces the 2026-08-11 report and is
the first capture of the curated 17-preset bank.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os` base with the shipping selective-O3 shader regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, tip `6f2e7e77` |
| Method | `HS_PROFILE`, 16-frame windows, 140 s fast cycle, two-frame transitions, epoch stretched to 1,400 revs |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 140 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

Image size: `FLASH: code:113,272, data:171,564, headers:9,052` / `RAM1:
variables:323,200, code:27,080, padding:5,688, free:168,320` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 1137–1152 root cycles divided by 600 MHz
match measured wall within **0.2 ppm**. The parser validates 92 windows, all 17
presets, the 16→0 wrap, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate:** peak frame render is **147.35 ms**, with **664/1,472
spilled frames (45.1%)**. Five preset buckets are green and twelve are red.
ShaderBall renders one 10,368-sample quadrant per frame; `canvas_buffer_wait`
is intentional cadence round-up idle.

## Phase-by-phase readout

Fast-cycle mode gives each preset held windows separated by a two-frame
through-clear transition.

### Worst clean hold: preset 16 (window frames 1137–1152)

```text
frame                 187.39 ms  112.44 Mcyc  100%
  sb_shader_draw      133.79 ms   80.27 Mcyc   71%  x1.0
  sb_timeline_step     64.63 us   38.81 kcyc    0%
  canvas_clear         85.94 us   51.58 kcyc    0%
  canvas_buffer_wait   41.46 ms   24.87 Mcyc   22%
```

Wall min/avg/max = 186.13/187.39/188.65 ms. The new dodecahedral
vector-noise/mirror preset is the worst held topology and cannot sustain the
62.5 ms display window.

### Transition into preset 16 (window frames 1121–1136)

```text
frame                 179.79 ms  107.87 Mcyc  100%
  sb_shader_draw      126.00 ms   75.60 Mcyc   70%  x0.94
  sb_timeline_step     68.25 us   40.98 kcyc    0%
  canvas_clear         87.06 us   52.26 kcyc    0%
  canvas_buffer_wait   41.64 ms   24.98 Mcyc   23%
```

Wall min/avg/max = 54.68/179.79/197.34 ms. One clear frame omits the shader;
the remaining frames already carry preset 16's cost.

### Per-preset table

Rows are ranked by maximum clean-hold shader cost. `Peak render` includes the
transition owned by that preset. The capture contains multiple 16→0 wraps and
2–4 clean windows per preset.

| # | Authored topology | shader ms | peak render ms | spilled/frames |
|---:|---|--:|--:|--:|
| 16 | Gnomonic dodecahedral vector-noise + mirror grid | 133.79 | 147.35 | 66/68 |
| 13 | Sinusoidal curl lattice, fine field | 85.34 | 107.89 | 64/68 |
| 14 | Sinusoidal curl lattice, coarse field | 84.98 | 97.91 | 68/68 |
| 7 | Dodecahedral grid + mirror | 69.73 | 82.27 | 98/102 |
| 15 | Stereographic prism polar-wave lattice | 68.77 | 97.44 | 66/68 |
| 8 | Peirce dodecahedral grid | 67.64 | 81.58 | 98/102 |
| 11 | Gnomonic dodecahedral wave/mirror grid | 66.34 | 79.08 | 66/68 |
| 6 | Kaleidoscope noise grid + edge fade | 62.68 | 75.36 | 64/102 |
| 10 | Dodecahedral noise lattice + mirror | 56.80 | 78.25 | 64/68 |
| 9 | Dodecahedral noise grid | 51.93 | 79.59 | 6/81 |
| 12 | Gnomonic affine lattice contour | 43.91 | 78.86 | 2/68 |
| 4 | Bonne kaleidoscope lattice + mirror | 37.94 | 49.30 | 0/102 |
| 0 | Glitch wave-shear grid | 35.50 | 146.53 | 2/99 |
| 5 | Peirce kaleidoscope lattice | 30.21 | 49.35 | 0/102 |
| 3 | Gnomonic glitch grid + mirror | 21.33 | 33.36 | 0/102 |
| 2 | Gnomonic kaleidoscope grid + mirror | 21.13 | 32.89 | 0/102 |
| 1 | Kaleidoscope twin-wave + inner mirror | 19.63 | 47.81 | 0/102 |

### Per-pixel figures

The worst held shader consumes about **7,742 cycles/sample** over the 10,368
quadrant samples. Per-pixel stage counters are disabled in the timing image.

## Column-ISR / DMA marshaling cost

```text
isr_wake       3455/frame  min/avg/max 0.58/1.96/20.08 us  cpu 3.60%
isr_pack        432/frame  min/avg/max 6.21/7.03/9.47 us   cpu 1.61%
isr_dma_submit  432/frame  min/avg/max 0.61/0.93/12.92 us  cpu 0.21%
```

- Packing dominates submit CPU time; LED transfer remains asynchronous.
- ISR CPU share is 5.42% in the worst held window.
- Preset 16 needs roughly a 2.36× peak-render speedup to fit 62.5 ms.

## Summary ranking

1. `sb_shader_draw` — 71% of the worst held frame, 133.79 ms.
2. `canvas_buffer_wait` — 22%, 41.46 ms of display cadence round-up.
3. Column wake/pack/submit ISRs — 5.42% CPU share under the live driver.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` remains parented under the first scope that enters it.
- Per-pixel stage scopes are disabled in this timing image.
- Fast-cycle mode compresses dwell only; it does not change per-frame work.
- The capture used a clean, pinned worktree at `6f2e7e77`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`; `profile_one.sh` performs locked build, flash, capture,
marker validation, and artifact attestation.

