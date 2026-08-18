# CurlLattice on-device profile — Teensy 4.0, segmented mode (2026-08-17, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile CurlLattice`).
Raw capture: `build/prof/curllattice_ship.log`.
Replaces `profile_curllattice_teensy_2026-08-16.md`: re-captured after the
pullback instance-pipeline refactor landed, confirming the stage-dispatch
inline regression it initially introduced (shader 40.8 ms → 56.4 ms) is fixed
by `4d93e449` and the tip now shades ~11% faster than the 08-16 report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | CurlLattice 288×144, single-entry playlist, tip `4d93e449` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 140 s capture, -D HS_PROFILE_EPOCH_REVS=1400 |
| Reproduce | `bash tools/profile_one.sh CurlLattice profile 140 16 -D HS_PROFILE_EPOCH_REVS=1400` |

Image sections (`arm-none-eabi-size -A`):

```text
.text.code     43,160   .text.progmem  145,488
.text.itcm     25,792   .data            4,800
.bss          310,272   .bss.dma       520,064
```

Exactness cross-check: window frames 2177-2192 root counter cycles ÷ 600 MHz
match the measured wall sum within **2.5 ppm**. Validation reports
138 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**45.47 ms** at frames 1-16, spilled **0/2,208 frames** (0%).

A display window is 62.5 ms. CurlLattice is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **17.03 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns 1,080 frames. The regimes below are
the worst held window and a representative transition window; the spread
across presets is in the per-preset table.

### Worst held window (window frames 2177–2192)

```text
frame                  62.44 ms  37.46 Mcyc  100%
  fx_shader_draw       39.86 ms  23.91 Mcyc   63%
  fx_prepare_frame      1.11 ms  667.02 kcyc    1%
  fx_advance            2.15 ms   1.29 Mcyc    3%
  fx_timeline_step      61.8 us  37.10 kcyc    0%
  canvas_clear          87.7 us  52.64 kcyc    0%
  canvas_buffer_wait   19.15 ms  11.49 Mcyc   30%
```

Wall min/avg/max = 61.88/62.44/62.92 ms. `fx_shader_draw` is 63% of the
frame at 39.86 ms; the worst render is 43.55 ms and the remaining
19.15 ms (30%) is `canvas_buffer_wait` idle. `fx_prepare_frame` costs
1.11 ms — the frame-context snapshot plus the pipeline's per-frame
`prepare_stages` pass, which resolves the curl loop point once per frame —
and `fx_advance` 2.15 ms for the runtime clocks, spatial frames and
generated-palette cycler step. None of that is per-sample work. The wall
spread is about a millisecond, which is the flywheel holding cadence rather
than render jitter.

### Transition window (window frames 641–656)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  fx_shader_draw       39.39 ms  23.64 Mcyc   63%
  fx_prepare_frame      1.10 ms  661.15 kcyc    1%
  fx_advance            2.16 ms   1.29 Mcyc    3%
  fx_timeline_step     135.0 us  81.02 kcyc    0%
  canvas_clear          87.6 us  52.57 kcyc    0%
  canvas_buffer_wait   19.56 ms  11.74 Mcyc   31%
```

Wall min/avg/max = 62.03/62.45/62.75 ms. During a morph the timeline lerp
re-interpolates the parameter block each frame; that lands in
`fx_timeline_step`, 135 µs here against 62 µs held. `fx_shader_draw` reads
39.39 ms (63% of the frame), 0.47 ms below the worst hold: the shade cost
tracks whichever pair of presets it is between, so a transition is bounded by
its two endpoints rather than being a cost peak of its own.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the first
`Preset:` marker; the capture wraps back to preset 1, confirming a full cycle.

| bucket | preset | fx_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 1 | open-curl | 39.86 | 43.55 | 0/530 | 16 |
| 0 | open-curl (initial hold) | 39.80 | 45.47 | 0/599 | 16 |
| 2 | dense-curl | 39.79 | 43.24 | 0/1,079 | 16 |

Span across presets is 1.002× of shader cost — every preset is green.

### Per-pixel figures

The dominant scope spends 2,307 cycles per quadrant sample over the
10,368 samples of a frame (2,598 in the 08-16 report — the per-frame prepared
state and the union-free slot types shaved ~11%). The fixed pipeline writes
the scan directly through `Scan::Shader::draw_cached`, so there is no
blended-pixel population and no `filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.81/2.02/20.32 us  cpu 3.72%
isr_pack         144/frame  min/avg/max 6.36/6.97/9.22 us  cpu 1.60%
isr_dma_submit   144/frame  min/avg/max 0.78/0.94/4.10 us  cpu 0.21%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.53%, already absorbed by the render
  counters. With 17.03 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fx_shader_draw` — 63% of the frame, 39.86 ms: the inlined pullback shade
   over 10,368 samples. A `FoldedSinusoidal` projection over a simplex
   `CurlNoise` surface integrated with Euler steps, feeding a
   `PrimitiveLattice` source; both planar warps are identity. The curl
   loop point is resolved once per frame into the pipeline's prepared state.
2. `fx_advance` — 3% of the frame, 2.15 ms: runtime clocks, spatial frames and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 1.11 ms: the frame-context snapshot,
   the hue LUTs, and the pipeline `prepare_stages` pass.

The perf ledger has no WASM or native baseline for this effect, so there is no
cross-target comparison to make yet.

## Caveats

- All scopes absorb live ISR time because CYCCNT free-runs.
- `filter_blend` does not appear in this direct scan; if it were entered it would
  parent under whichever scope reached it first, and its subtree would be hidden
  in windows where that parent had 0 calls.
- Shipping selective O3: `Scan::Shader::draw_cached` is `HS_O3_FN` with cached-
  flash placement, so the whole inlined shade path compiles at `-O3`; the OKLab
  and gamut helpers the generated palette bakes through are `HS_O3_FN` as well.
  The rest of the image keeps the `-Os` base policy.
- The inlined-pipeline property is inline-budget-sensitive: wrapping the stage
  policy dispatch in lambdas spilled `Surface::curl_noise` out of line and cost
  16 ms/frame until `4d93e449`. A `Pullback::Surface::curl_noise` symbol in the
  profile ELF is the regression tell.
- The curl-noise surface integration is the one stage here that is not a
  closed-form pullback, so its cost scales with the integrator step count
  rather than with coverage.
- The capture used `-D HS_PROFILE_EPOCH_REVS=1400` to fit the cycle inside one
  epoch. Epoch length changes how long a preset holds, never its per-frame
  cost.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=CurlLattice` and
`HS_PROFILE_WINDOW=16`; `just profile CurlLattice` performs the locked build,
flash, capture, marker validation, and artifact attestation.
