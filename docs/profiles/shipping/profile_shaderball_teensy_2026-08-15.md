# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-15, **selective -O3**)

Point-in-time snapshot of the extracted core pullback pipeline (regenerate with
`just profile ShaderBall`). Raw capture: `build/prof/shaderball_ship.log`.
This replaces the earlier 2026-08-15 capture at `910a7fe6`; the current bank
has 13 presets and includes the stereographic dodecahedral program.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os` base; selective O3 on the hot flash-resident FastNoise vector-noise leaves, with the pullback helpers forced inline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, tip `e920b4fa81c1` |
| Method | `HS_PROFILE`, 16-frame windows, 140 s fast cycle, two-frame transitions, epoch stretched to 1,400 revs |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 140 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

Image size: `FLASH: code:111,696, data:169,716, headers:8,380` / `RAM1:
variables:323,200, code:28,392, padding:4,376, free:168,320` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 1025–1040 root counter cycles divided by
600 MHz match the measured wall sum within **3.3 ppm**. Validation found 138
complete windows, all 13 presets, repeated 12→0 wraps, monotonic frames, no
epoch reset, and the `LANDED` pullback arm.

## Frame cadence

**Pass aggregate:** render averages **36.98 ms/frame**, peaks at **59.20 ms**
(frames 1041–1056), and spills **0/2,208 frames**. All 13 preset buckets are
green. The worst frame retains **3.30 ms** against the 62.5 ms display window.

ShaderBall renders one 10,368-sample quadrant per frame. The
`canvas_buffer_wait` scope is the intentional round-up idle before the next
display flip.

## Phase-by-phase readout

Fast-cycle mode gives every preset several held windows separated by a clear
and a two-frame pullback transition. Presets 8 and 9 share one compiled
pipeline, so that parameter-only handoff emits no separate endpoint pair.

### Worst held program: preset 4 (window frames 1025–1040)

```text
frame                  62.41 ms  37.45 Mcyc  100%
  sb_shader_draw       54.64 ms  32.78 Mcyc   87%  x1.00
  sb_timeline_step      61.4 us  36.86 kcyc    0%
  canvas_clear          92.6 us  55.54 kcyc    0%
  canvas_buffer_wait    5.08 ms   3.05 Mcyc    8%
```

Wall min/avg/max = 60.35/62.41/64.80 ms. The Bonne kaleidoscope lattice
program is now the shipping ceiling, but its peak render remains below one
display window. The previous preset-11 bottleneck is no longer dominant.

### Worst transition: preset 4 → 5 (window frames 1041–1056)

```text
frame                  62.12 ms  37.27 Mcyc  100%
  sb_shader_draw       48.99 ms  29.39 Mcyc   78%  x0.94
  sb_timeline_step      66.3 us  39.82 kcyc    0%
  canvas_clear          91.1 us  54.70 kcyc    0%
  canvas_buffer_wait    9.69 ms   5.81 Mcyc   15%
```

Wall min/avg/max = 9.94/62.12/110.99 ms. This window includes the clear frame
and the Bonne-to-Peirce endpoint handoff. Its worst render is **59.20 ms**;
the long wall maximum is display synchronization, not render work.

### Per-preset table

Rows are ranked by maximum clean-hold shader cost. Peak render includes the
transition owned by the preset. The capture contains at least six clean
windows per preset and five complete 12→0 wraps.

| # | Pullback program | shader ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 4 | Bonne kaleidoscope lattice mirror | 54.64 | 59.20 | 0/170 | 16 |
| 11 | Gnomonic dodecahedral vector mirror | 44.71 | 54.99 | 0/170 | 16 |
| 10 | Stereographic prism polar-wave lattice | 44.06 | 54.24 | 0/170 | 16 |
| 5 | Peirce dodecahedral grid | 40.74 | 56.95 | 0/170 | 16 |
| 8 | Sinusoidal curl lattice, variant A | 40.07 | 48.69 | 0/170 | 16 |
| 9 | Sinusoidal curl lattice, variant B | 40.03 | 44.74 | 0/170 | 16 |
| 7 | Gnomonic affine lattice contour | 36.95 | 47.98 | 0/170 | 16 |
| 6 | Gnomonic dodecahedral wave mirror | 35.74 | 51.89 | 0/170 | 16 |
| 12 | Stereographic dodecahedral grid inner mirror | 31.52 | 48.81 | 0/170 | 16 |
| 0 | Glitch noise grid wave-shear | 23.42 | 36.47 | 0/168 | 16 |
| 3 | Gnomonic glitch grid mirror | 19.21 | 23.10 | 0/170 | 16 |
| 1 | Kaleidoscope twin-wave inner mirror | 18.51 | 26.85 | 0/170 | 16 |
| 2 | Gnomonic kaleidoscope grid mirror | 18.30 | 22.38 | 0/170 | 16 |

### Per-pixel figures

The worst held shader consumes about **3,162 cycles/sample** over the 10,368
quadrant samples. ShaderBall writes the scan directly; there is no blended-pixel
population to report. Per-stage and per-pixel counters were disabled to avoid
perturbing this acceptance image.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1152/frame  min/avg/max 0.71/1.98/14.29 us  cpu 3.65%
isr_pack        144/frame  min/avg/max 6.22/7.13/9.43 us   cpu 1.64%
isr_dma_submit  144/frame  min/avg/max 0.69/0.94/6.59 us   cpu 0.21%
```

- Pack plus submit costs about 1.16 ms of CPU per frame in the worst held window.
- LED transfer remains asynchronous.
- ISR time is already absorbed by the render counters; no speedup is required
  for cadence, and the worst observed frame has 3.30 ms of margin.

## Summary ranking

1. Preset 4 `sb_shader_draw` — 54.64 ms, 87% of its held frame.
2. Preset 4→5 transition — 59.20 ms peak render, still cadence-safe.
3. Preset 11 `sb_shader_draw` — 44.71 ms, down from 65.60 ms in the replaced
   pre-specialization capture.

The selective vector-noise specialization removes every prior spill and moves
the shipping peak from 74.42 ms to 59.20 ms. Because the program bank also
gained preset 12, only like-for-like preset comparisons are strictly direct.

## Caveats

- All scopes absorb live ISR time because CYCCNT free-runs.
- `filter_blend` does not appear in this direct scan; if enabled elsewhere its
  subtree remains attached to the first entering parent.
- Per-stage and per-pixel profiling scopes were disabled in this timing image.
- Shipping selective O3 touches only the hot flash-resident FastNoise
  vector-noise leaves; the rest of the image retains the `-Os` base policy.
- `HS_PROFILE_SHADERBALL_FAST_CYCLE` compresses dwell and transitions but does
  not change per-frame program cost.
- The run used a clean detached worktree at `e920b4fa81c1`; report/spec edits in
  the main checkout did not participate in the build.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShaderBall`, and
`HS_PROFILE_WINDOW=16`; `profile_one.sh` performs locked build, flash, capture,
marker validation, and artifact attestation.
