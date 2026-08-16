# PrismLattice on-device profile — Teensy 4.0, segmented mode (2026-08-16, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile PrismLattice`).
Raw capture: `build/prof/prismlattice_ship.log`.
First profile of this effect: it entered the roster with the fixed-pipeline
workbench migration (`69d4751c`) and had no prior report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | PrismLattice 288×144, single-entry playlist, tip `69d4751cc077` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh PrismLattice profile 70 32` |

Image size:

```text
FLASH: code:67,008, data:149,640, headers:8,632
       free for files:1,806,336
RAM1:  variables:315,072, code:26,520, padding:6,248
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 769-800 root counter cycles ÷ 600 MHz
match the measured wall sum within **0.7 ppm**. Validation reports
34 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**30.79 ms** at frames 1-32, spilled **0/1,088 frames** (0%).


The peak frame falls in the first window, the post-init transient where the
palette bake and the first parameter snapshot run; the steady regime below is
the costliest window after it.

A display window is 62.5 ms. PrismLattice is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **31.71 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: PrismLattice has a single preset and one steady regime — the
camera and warp animation vary the shade cost slightly frame to frame but
never change the pipeline. The window below is the costliest of the pass.

### Steady state (window frames 481–512)

```text
frame                  62.47 ms  37.48 Mcyc  100%
  fx_shader_draw       25.95 ms  15.57 Mcyc   41%
  fx_prepare_frame     677.6 us  406.56 kcyc    1%
  fx_advance            1.92 ms   1.15 Mcyc    3%
  fx_timeline_step      49.2 us  29.50 kcyc    0%
  canvas_clear          86.0 us  51.61 kcyc    0%
  canvas_buffer_wait   33.78 ms  20.27 Mcyc   54%
```

Wall min/avg/max = 61.16/62.47/63.73 ms. `fx_shader_draw` is 41% of the
frame at 25.95 ms; the whole render is 28.69 ms and the remaining
33.78 ms (54%) is `canvas_buffer_wait` idle. `fx_prepare_frame` costs
678 µs and `fx_advance` 1.92 ms — the runtime advance, the spatial-frame update
and the generated-palette cycler step. None of that is per-sample work. The
wall spread is under a millisecond, which is the flywheel holding cadence
rather than render jitter.

### Per-pixel figures

The dominant scope spends 1,502 cycles per quadrant sample over the
10,368 samples of a frame. The fixed pipeline writes the scan directly through
`Scan::Shader::draw_cached`, so there is no blended-pixel population and no
`filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.68/2.01/18.49 us  cpu 3.71%
isr_pack         144/frame  min/avg/max 6.35/7.21/9.60 us  cpu 1.66%
isr_dma_submit   144/frame  min/avg/max 0.66/0.93/4.48 us  cpu 0.21%
```

- Pack plus submit costs 1.17 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.58%, already absorbed by the render
  counters. With 31.71 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fx_shader_draw` — 41% of the frame, 25.95 ms: the inlined pullback shade
   over 10,368 samples. Stereographic projection through a
   `TriangularPrismKaleidoscope` lens, an outer `PolarChart` (linear-polar)
   plus inner `WaveShear` warp pair, and a `PrimitiveLattice` source under an
   analogous palette.
2. `fx_advance` — 3% of the frame, 1.92 ms: preset transition bookkeeping and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 678 µs: per-frame snapshot of the
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
- The capture ran with the `HS_PROFILE` scopes this sweep added to the fixed-
  pipeline draw path, since landed as `ac6fe641`. They compile to nothing
  without `HS_PROFILE_ENABLE`, so the shipping image is unchanged.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=PrismLattice` and
`HS_PROFILE_WINDOW=32`; `just profile PrismLattice` performs the locked build,
flash, capture, marker validation, and artifact attestation.
