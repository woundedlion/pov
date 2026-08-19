# CurlLattice on-device profile — Teensy 4.0, segmented mode (2026-08-18, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile CurlLattice`).
Raw capture: `build/prof/curllattice_ship.log`.
Replaces `profile_curllattice_teensy_2026-08-17.md`: re-captured after the
ranked stage-families cut-over (`7da9afbdc`) landed, as the A/B gate for that
landing. A same-parameters baseline capture of the pre-cut-over master
(`4565d605d`, second bench board) read 39.16 ms dense-curl / 39.10 ms
open-curl clean-hold shader against the tip's 39.17 / 39.15 — parity within
0.1%, so the chain-normalization rewrite is perf-neutral on a displacing
effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | CurlLattice 288×144, single-entry playlist, tip `7da9afbdc` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture, -D HS_PROFILE_EPOCH_REVS=1200 |
| Reproduce | `bash tools/profile_one.sh CurlLattice profile 110 16 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:68984, data:150344, headers:9024` / `RAM1:
variables:315072, code:25752, padding:7016, free:176448` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 1089–1104 root counter cycles ÷ 600 MHz
match the measured wall sum within **3.2 ppm**. Validation reports 108
complete windows, monotonic frame numbers, no epoch reset, and the cycle
wrapping back to preset 1.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**45.18 ms** at frames 1-16, spilled **0/1,728 frames** (0%).

A display window is 62.5 ms. CurlLattice is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()` nor
`persists_pixels()`, so there is no full-canvas multiplier. The worst frame of
the pass keeps **17.32 ms** of the window in hand, so the effect holds 16 fps
for the whole capture. `canvas_buffer_wait` is the round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns ~1,080 frames. The regimes below are
the worst held window and a representative transition window; the spread
across presets is in the per-preset table.

### Worst held window (window frames 1089–1104, dense-curl)

```text
frame                  62.44 ms  37.46 Mcyc  100%
  fx_shader_draw       39.17 ms  23.50 Mcyc   62%
  fx_prepare_frame      0.88 ms  526.59 kcyc    1%
  fx_advance            2.11 ms   1.27 Mcyc    3%
  fx_timeline_step      75.2 us  45.13 kcyc    0%
  canvas_clear          88.3 us  52.99 kcyc    0%
  canvas_buffer_wait   20.10 ms  12.06 Mcyc   32%
```

Wall min/avg/max = 61.96/62.44/62.79 ms. `fx_shader_draw` is 62% of the frame
at 39.17 ms; the window's worst render is 42.53 ms and the remaining 20.10 ms
(32%) is `canvas_buffer_wait` idle. `fx_prepare_frame` costs 0.88 ms — the
frame-context snapshot plus the pipeline's per-frame prepare pass, which
resolves the curl loop point once per frame — and `fx_advance` 2.11 ms for the
runtime clocks, spatial frames and generated-palette cycler step. None of that
is per-sample work. The wall spread is under a millisecond, which is the
flywheel holding cadence rather than render jitter.

### Transition window (window frames 641–656)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  fx_shader_draw       38.66 ms  23.20 Mcyc   61%
  fx_prepare_frame      1.12 ms  670.90 kcyc    1%
  fx_advance            2.18 ms   1.31 Mcyc    3%
  fx_timeline_step     147.3 us  88.40 kcyc    0%
  canvas_clear          88.6 us  53.15 kcyc    0%
  canvas_buffer_wait   20.25 ms  12.15 Mcyc   32%
```

Wall min/avg/max = 62.19/62.45/62.70 ms. During a morph the timeline lerp
re-interpolates the parameter block each frame; that lands in
`fx_timeline_step`, 147 µs here against 75 µs held. `fx_shader_draw` reads
38.66 ms (61% of the frame), below both holds: the shade cost tracks whichever
pair of presets it is between, so a transition is bounded by its endpoints
rather than being a cost peak of its own.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the
first `Preset:` marker; the capture wraps back to preset 1, confirming a full
cycle (bucket 1 is the post-wrap sliver).

| bucket | preset | fx_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 2 | dense-curl | 39.17 | 42.66 | 0/1,079 | 16 |
| 0 | open-curl (initial hold) | 39.15 | 45.18 | 0/599 | 16 |
| 1 | open-curl (wrap) | 38.92 | 42.54 | 0/50 | 16 |

Span across presets is 1.006× of shader cost — every preset is green.

### Per-pixel figures

The dominant scope spends 2,267 cycles per quadrant sample over the 10,368
samples of a frame (2,307 in the 08-17 report). The pipeline writes the scan
directly through `Scan::Shader::draw_cached`, so there is no blended-pixel
population and no `filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.80/2.36/25.72 us  cpu 4.34%
isr_pack         144/frame  min/avg/max 6.26/7.01/9.54 us  cpu 1.61%
isr_dma_submit   144/frame  min/avg/max 0.78/0.94/1.03 us  cpu 0.21%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the
  CPU.
- Total ISR CPU share is 6.16%, already absorbed by the render counters. With
  17.32 ms of margin at the worst frame, no speedup is required for cadence.

## Summary ranking

1. `fx_shader_draw` — 62% of the frame, 39.17 ms: the inlined pullback shade
   over 10,368 samples. A `FoldedSinusoidal` projection over a simplex
   `CurlNoise` surface integrated with Euler steps, feeding a
   `PrimitiveLattice` source, now expressed as a ranked-family chain
   (`Displace → Project` inside the single `OUT_OF_LINE_FLASH` sphere run).
   The curl loop point is resolved once per frame into the prepared state.
2. `fx_advance` — 3% of the frame, 2.11 ms: runtime clocks, spatial frames and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 0.88 ms: the frame-context snapshot,
   the hue LUTs, and the pipeline prepare pass.

The perf ledger has no WASM or native baseline for this effect, so there is no
cross-target comparison to make yet.

## Caveats

- All scopes absorb live ISR time because CYCCNT free-runs.
- `filter_blend` does not appear in this direct scan; if it were entered it
  would parent under whichever scope reached it first, and its subtree would
  be hidden in windows where that parent had 0 calls.
- Shipping selective O3: `Scan::Shader::draw_cached` is `HS_O3_FN` with
  cached-flash placement, so the whole inlined shade path compiles at `-O3`;
  the OKLab and gamut helpers the generated palette bakes through are
  `HS_O3_FN` as well. The rest of the image keeps the `-Os` base policy.
- The inlined-pipeline property is inline-budget-sensitive: the map-file gate
  for it is exactly one `OUT_OF_LINE_FLASH` `BoundPlaced::run` sphere-run
  symbol per displacing effect and no standalone hot stage symbols. A
  `Pullback::Surface::curl_noise` symbol in the profile ELF is the regression
  tell.
- The curl-noise surface integration is the one stage here that is not a
  closed-form pullback, so its cost scales with the integrator step count
  rather than with coverage.
- The capture used `-D HS_PROFILE_EPOCH_REVS=1200` (150 s epoch) to fit the
  110 s pass inside one epoch. Epoch length changes how long a preset holds,
  never its per-frame cost.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=CurlLattice` and
`HS_PROFILE_WINDOW=16`; `tools/profile_one.sh` performs the locked build,
flash, capture, marker validation, and artifact attestation.
