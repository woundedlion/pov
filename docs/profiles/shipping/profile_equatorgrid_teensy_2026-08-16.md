# EquatorGrid on-device profile — Teensy 4.0, segmented mode (2026-08-16, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile EquatorGrid`).
Raw capture: `build/prof/equatorgrid_ship.log`.
First profile of this effect: it entered the roster with the fixed-pipeline
workbench migration (`69d4751c`) and had no prior report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | EquatorGrid 288×144, single-entry playlist, tip `69d4751cc077` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, -D HS_PROFILE_EPOCH_REVS=2400 |
| Reproduce | `bash tools/profile_one.sh EquatorGrid profile 260 16 -D HS_PROFILE_EPOCH_REVS=2400` |

Image size:

```text
FLASH: code:68,096, data:149,580, headers:8,628
       free for files:1,805,312
RAM1:  variables:315,072, code:26,488, padding:6,280
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 1105-1120 root counter cycles ÷ 600 MHz
match the measured wall sum within **4.0 ppm**. Validation reports
258 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**35.91 ms** at frames 1905-1920, spilled **0/4,128 frames** (0%).


A display window is 62.5 ms. EquatorGrid is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **26.59 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns 1,080 frames. The regimes below are
the worst held window and the worst transition window; the spread across
presets is in the per-preset table.

### Worst held window (window frames 1345–1360)

```text
frame                  62.43 ms  37.46 Mcyc  100%
  fx_shader_draw       31.90 ms  19.14 Mcyc   51%
  fx_prepare_frame     944.2 us  566.58 kcyc    1%
  fx_advance            1.99 ms   1.20 Mcyc    3%
  fx_timeline_step      45.9 us  27.55 kcyc    0%
  canvas_clear          86.4 us  51.89 kcyc    0%
  canvas_buffer_wait   27.46 ms  16.48 Mcyc   43%
```

Wall min/avg/max = 60.74/62.43/64.23 ms. `fx_shader_draw` is 51% of the
frame at 31.90 ms; the whole render is 34.98 ms and the remaining
27.46 ms (43%) is `canvas_buffer_wait` idle. `fx_prepare_frame` costs
944 µs and `fx_advance` 1.99 ms — the runtime advance, the spatial-frame update
and the generated-palette cycler step, plus preset transition bookkeeping. None
of that is per-sample work. The wall spread is under a millisecond, which is
the flywheel holding cadence rather than render jitter.

### Worst transition window (window frames 1905–1920)

```text
frame                  62.44 ms  37.46 Mcyc  100%
  fx_shader_draw       31.39 ms  18.83 Mcyc   50%
  fx_prepare_frame      1.32 ms  794.49 kcyc    2%
  fx_advance            1.98 ms   1.19 Mcyc    3%
  fx_timeline_step      53.2 us  31.95 kcyc    0%
  canvas_clear          86.2 us  51.72 kcyc    0%
  canvas_buffer_wait   27.60 ms  16.56 Mcyc   44%
```

Wall min/avg/max = 60.33/62.44/64.55 ms. During a
morph the parameter block is re-interpolated every frame and the palette mapping
is lerped, but that work lands in `fx_advance`, which is 1.98 ms here
against 1.99 ms held. `fx_shader_draw` reads 31.39 ms (50% of the frame),
0.51 ms below the worst hold: the shade cost tracks whichever pair of
presets it is between, so a transition is bounded by its two endpoints rather
than being a cost peak of its own.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the first
`Preset:` marker; the capture wraps back to preset 1, confirming a full cycle.

| bucket | preset | fx_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 2 | open-grid | 32.25 | 35.90 | 0/1,369 | 16 |
| 1 | double-map | 32.07 | 35.84 | 0/1,080 | 16 |
| 3 | fine-grid | 31.39 | 35.91 | 0/1,080 | 16 |
| 0 | double-map (initial hold) | 31.32 | 35.43 | 0/599 | 16 |

Span across presets is 1.03× of shader cost — every preset is green.

### Per-pixel figures

The dominant scope spends 1,846 cycles per quadrant sample over the
10,368 samples of a frame. The fixed pipeline writes the scan directly through
`Scan::Shader::draw_cached`, so there is no blended-pixel population and no
`filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.73/2.00/13.66 us  cpu 3.68%
isr_pack         144/frame  min/avg/max 6.33/7.04/9.23 us  cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.65/0.93/5.62 us  cpu 0.21%
```

- Pack plus submit costs 1.15 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.51%, already absorbed by the render
  counters. With 26.59 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fx_shader_draw` — 51% of the frame, 31.90 ms: the inlined pullback shade
   over 10,368 samples. Equirectangular projection through a
   `DodecahedralKaleidoscope` lens, an inner `MirrorTile` warp, and a `Grid`
   source under an analogous noise-hue palette.
2. `fx_advance` — 3% of the frame, 1.99 ms: preset transition bookkeeping and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 944 µs: per-frame snapshot of the
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
- `USES_CENTRAL_MERIDIAN` pins the equirectangular chart to a moving meridian;
  the wander is per-frame and is included in `fx_advance`, not in the shader
  draw.
- The capture used `-D HS_PROFILE_EPOCH_REVS=2400` to fit the cycle inside one
  epoch. Epoch length changes how long a preset holds, never its per-frame
  cost.
- The capture ran with the `HS_PROFILE` scopes this sweep added to the fixed-
  pipeline draw path, since landed as `ac6fe641`. They compile to nothing
  without `HS_PROFILE_ENABLE`, so the shipping image is unchanged.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=EquatorGrid` and
`HS_PROFILE_WINDOW=16`; `just profile EquatorGrid` performs the locked build,
flash, capture, marker validation, and artifact attestation.
