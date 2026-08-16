# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-15, **-O3**)

Global-O3 reference twin of the
[shipping capture](../shipping/profile_shaderball_teensy_2026-08-15.md).
Raw capture: `build/prof/shaderball_o3.log`. This replaces the earlier
2026-08-15 O3 capture at `910a7fe6` and covers the current 13-preset bank.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3`: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, tip `e920b4fa81c1` |
| Method | `HS_PROFILE`, 16-frame windows, 140 s fast cycle, two-frame transitions, epoch stretched to 1,400 revs |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 140 16 -D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400` |

Image size: `FLASH: code:128,960, data:170,084, headers:9,180` / `RAM1:
variables:323,200, code:40,872, padding:24,664, free:135,552` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 817–832 root counter cycles divided by
600 MHz match the measured wall sum within **1.9 ppm**. Validation found 138
complete windows, all 13 presets, repeated 12→0 wraps, monotonic frames, no
epoch reset, and the `LANDED` pullback arm.

## Frame cadence

**Pass aggregate:** render averages **32.83 ms/frame**, peaks at **50.72 ms**
(frames 2129–2144), and spills **0/2,208 frames**. All 13 program buckets hold
16 fps, with 11.79 ms of worst-frame margin against one display window.

ShaderBall renders one 10,368-sample quadrant per frame. The
`canvas_buffer_wait` scope is display-sync idle by design.

## Phase-by-phase readout

Fast-cycle mode gives every preset several held windows separated by a clear
and a two-frame pullback transition. Presets 8 and 9 share one compiled
pipeline, so that parameter-only handoff emits no separate endpoint pair.

### Worst held program: preset 11 (window frames 817–832)

```text
frame                  62.25 ms  37.35 Mcyc  100%
  sb_shader_draw       43.63 ms  26.18 Mcyc   70%  x1.00
  sb_timeline_step      77.1 us  46.25 kcyc    0%
  canvas_clear          87.8 us  52.71 kcyc    0%
  canvas_buffer_wait   14.88 ms   8.93 Mcyc   23%
```

Wall min/avg/max = 60.04/62.25/63.15 ms. Preset 11 remains the costliest O3
hold, but it retains 14.72 ms of render margin in this window.

### Worst transition: preset 10 → 11 (window frames 2129–2144)

```text
frame                  63.51 ms  38.10 Mcyc  100%
  sb_shader_draw       29.28 ms  17.57 Mcyc   46%  x0.94
  sb_timeline_step      72.8 us  43.67 kcyc    0%
  canvas_clear          86.3 us  51.82 kcyc    0%
  canvas_buffer_wait   30.87 ms  18.52 Mcyc   48%
```

Wall min/avg/max = 33.40/63.51/110.43 ms. The clear frame lowers the window
mean, while the first preset-11 endpoint frame sets the **50.72 ms** pass peak.
The wall maximum is synchronization, not a render overrun.

### Per-preset table

Rows are ranked by maximum clean-hold shader cost. Peak render includes the
transition owned by the preset. The capture contains at least six clean
windows per preset and five complete 12→0 wraps.

| # | Pullback program | shader ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 11 | Gnomonic dodecahedral vector mirror | 43.63 | 50.72 | 0/170 | 16 |
| 4 | Bonne kaleidoscope lattice mirror | 41.15 | 44.13 | 0/170 | 16 |
| 5 | Peirce dodecahedral grid | 38.84 | 48.80 | 0/170 | 16 |
| 8 | Sinusoidal curl lattice, variant A | 38.73 | 45.57 | 0/170 | 16 |
| 9 | Sinusoidal curl lattice, variant B | 38.64 | 43.32 | 0/170 | 16 |
| 10 | Stereographic prism polar-wave lattice | 35.45 | 43.56 | 0/170 | 16 |
| 6 | Gnomonic dodecahedral wave mirror | 34.53 | 48.83 | 0/170 | 16 |
| 7 | Gnomonic affine lattice contour | 34.19 | 38.91 | 0/170 | 16 |
| 12 | Stereographic dodecahedral grid inner mirror | 31.17 | 48.01 | 0/170 | 16 |
| 0 | Glitch noise grid wave-shear | 22.42 | 35.30 | 0/168 | 16 |
| 3 | Gnomonic glitch grid mirror | 19.32 | 22.52 | 0/170 | 16 |
| 1 | Kaleidoscope twin-wave inner mirror | 18.10 | 25.79 | 0/170 | 16 |
| 2 | Gnomonic kaleidoscope grid mirror | 18.08 | 21.98 | 0/170 | 16 |

### Per-pixel figures

The worst held shader consumes about **2,524 cycles/sample** over the 10,368
quadrant samples. ShaderBall writes the scan directly; there is no blended-pixel
population to report. Per-stage and per-pixel counters were disabled.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1149/frame  min/avg/max 0.56/1.88/19.06 us  cpu 3.46%
isr_pack        144/frame  min/avg/max 6.32/7.31/9.53 us   cpu 1.68%
isr_dma_submit  144/frame  min/avg/max 0.62/0.92/4.92 us   cpu 0.21%
```

- Pack plus submit costs about 1.18 ms of CPU per frame in the worst held window.
- LED transfer remains asynchronous.
- ISR time is already included in the render counters; the worst frame needs
  no cadence speedup and retains 11.79 ms of margin.

## Summary ranking

1. Preset 11 `sb_shader_draw` — 43.63 ms, 70% of its held frame.
2. Preset 10→11 transition — 50.72 ms peak render.
3. Preset 4 `sb_shader_draw` — 41.15 ms, the largest remaining gap versus the
   selective shipping image.

No current O3 preset misses the 62.5 ms display budget.

## Caveats

- All scopes absorb live ISR time because CYCCNT free-runs.
- `filter_blend` does not appear in this direct scan; if enabled elsewhere its
  subtree remains attached to the first entering parent.
- Per-stage and per-pixel profiling scopes were disabled in this timing image.
- Global O3 is a compiler reference, not a shippable roster image.
- `HS_PROFILE_SHADERBALL_FAST_CYCLE` changes dwell duration, not per-frame cost.
- The run used a clean detached worktree at `e920b4fa81c1`; report/spec edits in
  the main checkout did not participate in the build.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShaderBall`, and
`HS_PROFILE_WINDOW=16`; use the reproduce command above for the matched O3
reference.

## Global O3 versus selective O3

Global O3 lowers the matched peak from **59.20 ms to 50.72 ms** (14.3%) and
the pass-average render from **36.98 ms to 32.83 ms** (11.2%). It adds
**17,264 B flash code** and **12,480 B ITCM**. Selective O3 has nearly closed
the preset-11 gap (44.71 ms versus 43.63 ms); preset 4 is now the useful
codegen target (54.64 ms versus 41.15 ms).
