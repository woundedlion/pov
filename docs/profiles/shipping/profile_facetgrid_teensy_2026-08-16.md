# FacetGrid on-device profile — Teensy 4.0, segmented mode (2026-08-16, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile FacetGrid`).
Raw capture: `build/prof/facetgrid_ship.log`.
First profile of this effect: it entered the roster with the fixed-pipeline
workbench migration (`69d4751c`) and had no prior report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | FacetGrid 288×144, single-entry playlist, tip `69d4751cc077` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, -D HS_PROFILE_EPOCH_REVS=2400 |
| Reproduce | `bash tools/profile_one.sh FacetGrid profile 260 16 -D HS_PROFILE_EPOCH_REVS=2400` |

Image size:

```text
FLASH: code:66,480, data:149,556, headers:8,220
       free for files:1,807,360
RAM1:  variables:315,072, code:25,480, padding:7,288
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 3921-3936 root counter cycles ÷ 600 MHz
match the measured wall sum within **4.0 ppm**. Validation reports
258 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**35.52 ms** at frames 1057-1072, spilled **0/4,128 frames** (0%).


A display window is 62.5 ms. FacetGrid is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **26.98 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns 1,080 frames. The regimes below are
the worst held window and the worst transition window; the spread across
presets is in the per-preset table.

### Worst held window (window frames 1057–1072)

```text
frame                  62.44 ms  37.46 Mcyc  100%
  fg_shader_draw       31.53 ms  18.92 Mcyc   50%
  fg_prepare_frame     887.2 us  532.37 kcyc    1%
  fg_advance            2.01 ms   1.20 Mcyc    3%
  fg_timeline_step      54.5 us  32.71 kcyc    0%
  canvas_clear          86.4 us  51.86 kcyc    0%
  canvas_buffer_wait   27.86 ms  16.72 Mcyc   44%
```

Wall min/avg/max = 60.34/62.44/64.30 ms. `fg_shader_draw` is 50% of the
frame at 31.53 ms; the whole render is 34.57 ms and the remaining
27.86 ms (44%) is `canvas_buffer_wait` idle. `fg_prepare_frame` costs
887 µs and `fg_advance` 2.01 ms — the runtime advance, the spatial-frame update
and the generated-palette cycler step, plus preset transition bookkeeping. None
of that is per-sample work. The wall spread is under a millisecond, which is
the flywheel holding cadence rather than render jitter.

### Worst transition window (window frames 1041–1056)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  fg_shader_draw       31.27 ms  18.76 Mcyc   50%
  fg_prepare_frame     904.9 us  542.95 kcyc    1%
  fg_advance            2.01 ms   1.21 Mcyc    3%
  fg_timeline_step      63.4 us  38.09 kcyc    0%
  canvas_clear          86.3 us  51.82 kcyc    0%
  canvas_buffer_wait   28.11 ms  16.87 Mcyc   45%
```

Wall min/avg/max = 60.17/62.45/64.63 ms. During a
morph the parameter block is re-interpolated every frame and the palette mapping
is lerped, but that work lands in `fg_advance`, which is 2.01 ms here
against 2.01 ms held. `fg_shader_draw` reads 31.27 ms (50% of the frame),
0.26 ms below the worst hold: the shade cost tracks whichever pair of
presets it is between, so a transition is bounded by its two endpoints rather
than being a cost peak of its own.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the first
`Preset:` marker; the capture wraps back to preset 1, confirming a full cycle.

| bucket | preset | fg_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 1 | coupled-grid | 32.23 | 35.34 | 0/285 | 16 |
| 2 | direct-grid | 32.18 | 35.52 | 0/1,081 | 16 |
| 4 | stretched-grid | 29.58 | 31.38 | 0/1,081 | 16 |
| 0 | coupled-grid (initial hold) | 29.28 | 33.28 | 0/600 | 16 |
| 3 | double-map | 27.10 | 31.62 | 0/1,081 | 16 |

Span across presets is 1.19× of shader cost — every preset is green.

### Per-pixel figures

The dominant scope spends 1,825 cycles per quadrant sample over the
10,368 samples of a frame. The fixed pipeline writes the scan directly through
`Scan::Shader::draw_cached`, so there is no blended-pixel population and no
`filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.73/1.98/17.52 us  cpu 3.65%
isr_pack         144/frame  min/avg/max 6.23/6.95/9.54 us  cpu 1.60%
isr_dma_submit   144/frame  min/avg/max 0.68/0.93/7.65 us  cpu 0.21%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.46%, already absorbed by the render
  counters. With 26.98 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fg_shader_draw` — 50% of the frame, 31.53 ms: the inlined pullback shade
   over 10,368 samples. Stereographic projection through a
   `DodecahedralKaleidoscope` lens, an inner `MirrorTile` warp, a `Grid`
   source, and a `ProjectionSquared` coverage term.
2. `fg_advance` — 3% of the frame, 2.01 ms: preset transition bookkeeping and
   the generated-palette cycler step.
3. `fg_prepare_frame` — 1% of the frame, 887 µs: per-frame snapshot of the
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
- The capture used `-D HS_PROFILE_EPOCH_REVS=2400` to fit the cycle inside one
  epoch. Epoch length changes how long a preset holds, never its per-frame
  cost.
- The capture ran with the `HS_PROFILE` scopes this sweep added to the fixed-
  pipeline draw path, since landed as `ac6fe641`. They compile to nothing
  without `HS_PROFILE_ENABLE`, so the shipping image is unchanged.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=FacetGrid` and
`HS_PROFILE_WINDOW=16`; `just profile FacetGrid` performs the locked build,
flash, capture, marker validation, and artifact attestation.
