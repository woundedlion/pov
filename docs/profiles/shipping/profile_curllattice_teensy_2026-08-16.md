# CurlLattice on-device profile — Teensy 4.0, segmented mode (2026-08-16, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile CurlLattice`).
Raw capture: `build/prof/curllattice_ship.log`.
First profile of this effect: it entered the roster with the fixed-pipeline
workbench migration (`69d4751c`) and had no prior report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | CurlLattice 288×144, single-entry playlist, tip `69d4751cc077` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 140 s capture, -D HS_PROFILE_EPOCH_REVS=1400 |
| Reproduce | `bash tools/profile_one.sh CurlLattice profile 140 16 -D HS_PROFILE_EPOCH_REVS=1400` |

Image size:

```text
FLASH: code:67,024, data:149,112, headers:9,144
       free for files:1,806,336
RAM1:  variables:315,072, code:25,480, padding:7,288
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2177-2192 root counter cycles ÷ 600 MHz
match the measured wall sum within **1.8 ppm**. Validation reports
138 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**48.67 ms** at frames 2177-2192, spilled **0/2,208 frames** (0%).


A display window is 62.5 ms. CurlLattice is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **13.83 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns 1,080 frames. The regimes below are
the worst held window and the worst transition window; the spread across
presets is in the per-preset table.

### Worst held window (window frames 2177–2192)

```text
frame                  62.42 ms  37.45 Mcyc  100%
  cl_shader_draw       44.89 ms  26.94 Mcyc   71%
  cl_prepare_frame      1.11 ms  666.62 kcyc    1%
  cl_advance            2.21 ms   1.33 Mcyc    3%
  cl_timeline_step      49.2 us  29.59 kcyc    0%
  canvas_clear          88.2 us  52.95 kcyc    0%
  canvas_buffer_wait   14.06 ms   8.44 Mcyc   22%
```

Wall min/avg/max = 61.69/62.42/63.02 ms. `cl_shader_draw` is 71% of the
frame at 44.89 ms; the whole render is 48.36 ms and the remaining
14.06 ms (22%) is `canvas_buffer_wait` idle. `cl_prepare_frame` costs
1.11 ms and `cl_advance` 2.21 ms — the runtime advance, the spatial-frame
update and the generated-palette cycler step, plus preset transition
bookkeeping. None of that is per-sample work. The wall spread is under a
millisecond, which is the flywheel holding cadence rather than render jitter.

### Worst transition window (window frames 593–608)

```text
frame                  62.42 ms  37.45 Mcyc  100%
  cl_shader_draw       44.22 ms  26.53 Mcyc   70%
  cl_prepare_frame      1.13 ms  678.45 kcyc    1%
  cl_advance            2.12 ms   1.27 Mcyc    3%
  cl_timeline_step      48.4 us  29.08 kcyc    0%
  canvas_clear          88.0 us  52.81 kcyc    0%
  canvas_buffer_wait   14.80 ms   8.88 Mcyc   23%
```

Wall min/avg/max = 60.07/62.41/64.47 ms. During a
morph the parameter block is re-interpolated every frame and the palette mapping
is lerped, but that work lands in `cl_advance`, which is 2.12 ms here
against 2.21 ms held. `cl_shader_draw` reads 44.22 ms (70% of the frame),
0.67 ms below the worst hold: the shade cost tracks whichever pair of
presets it is between, so a transition is bounded by its two endpoints rather
than being a cost peak of its own.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the first
`Preset:` marker; the capture wraps back to preset 1, confirming a full cycle.

| bucket | preset | cl_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 1 | open-curl | 44.89 | 48.67 | 0/527 | 16 |
| 0 | open-curl (initial hold) | 44.78 | 48.56 | 0/600 | 16 |
| 2 | dense-curl | 44.66 | 48.08 | 0/1,081 | 16 |

Span across presets is 1.01× of shader cost — every preset is green.

### Per-pixel figures

The dominant scope spends 2,598 cycles per quadrant sample over the
10,368 samples of a frame. The fixed pipeline writes the scan directly through
`Scan::Shader::draw_cached`, so there is no blended-pixel population and no
`filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1151/frame  min/avg/max 0.73/1.99/15.72 us  cpu 3.67%
isr_pack         144/frame  min/avg/max 6.36/6.96/9.16 us  cpu 1.60%
isr_dma_submit   144/frame  min/avg/max 0.68/0.94/7.86 us  cpu 0.21%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.48%, already absorbed by the render
  counters. With 13.83 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `cl_shader_draw` — 71% of the frame, 44.89 ms: the inlined pullback shade
   over 10,368 samples. A `FoldedSinusoidal` projection over a simplex
   `CurlNoise` surface integrated with Euler steps, feeding a
   `PrimitiveLattice` source; both planar warps are identity.
2. `cl_advance` — 3% of the frame, 2.21 ms: preset transition bookkeeping and
   the generated-palette cycler step.
3. `cl_prepare_frame` — 1% of the frame, 1.11 ms: per-frame snapshot of the
   parameter block and the hue LUTs.

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
- The curl-noise surface integration is the one stage here that is not a
  closed-form pullback, so its cost scales with the integrator step count
  rather than with coverage.
- The capture used `-D HS_PROFILE_EPOCH_REVS=1400` to fit the cycle inside one
  epoch. Epoch length changes how long a preset holds, never its per-frame
  cost.
- The capture ran with the `HS_PROFILE` scopes this sweep added to the fixed-
  pipeline draw path, since landed as `ac6fe641`. They compile to nothing
  without `HS_PROFILE_ENABLE`, so the shipping image is unchanged.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=CurlLattice` and
`HS_PROFILE_WINDOW=16`; `just profile CurlLattice` performs the locked build,
flash, capture, marker validation, and artifact attestation.
