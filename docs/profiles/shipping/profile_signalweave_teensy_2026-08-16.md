# SignalWeave on-device profile — Teensy 4.0, segmented mode (2026-08-16, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile SignalWeave`).
Raw capture: `build/prof/signalweave_ship.log`.
First profile of this effect: it entered the roster with the fixed-pipeline
workbench migration (`69d4751c`) and had no prior report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | SignalWeave 288×144, single-entry playlist, tip `69d4751cc077` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 300 s capture, -D HS_PROFILE_EPOCH_REVS=2600 |
| Reproduce | `bash tools/profile_one.sh SignalWeave profile 300 16 -D HS_PROFILE_EPOCH_REVS=2600` |

Image size:

```text
FLASH: code:67,544, data:149,528, headers:8,208
       free for files:1,806,336
RAM1:  variables:315,072, code:26,472, padding:6,296
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 3249-3264 root counter cycles ÷ 600 MHz
match the measured wall sum within **3.5 ppm**. Validation reports
298 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**31.95 ms** at frames 1-16, spilled **0/4,768 frames** (0%).


The peak frame falls in the first window, the post-init transient where the
palette bake and the first parameter snapshot run; the steady regime below is
the costliest window after it.

A display window is 62.5 ms. SignalWeave is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **30.55 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns 1,080 frames. The regimes below are
the worst held window and the worst transition window; the spread across
presets is in the per-preset table.

### Worst held window (window frames 657–672)

```text
frame                  62.43 ms  37.46 Mcyc  100%
  fx_shader_draw       26.98 ms  16.19 Mcyc   43%
  fx_prepare_frame     988.2 us  592.96 kcyc    1%
  fx_advance            2.27 ms   1.36 Mcyc    3%
  fx_timeline_step      57.5 us  34.51 kcyc    0%
  canvas_clear          86.0 us  51.60 kcyc    0%
  canvas_buffer_wait   32.04 ms  19.22 Mcyc   51%
```

Wall min/avg/max = 62.05/62.43/62.81 ms. `fx_shader_draw` is 43% of the
frame at 26.98 ms; the whole render is 30.39 ms and the remaining
32.04 ms (51%) is `canvas_buffer_wait` idle. `fx_prepare_frame` costs
988 µs and `fx_advance` 2.27 ms — the runtime advance, the spatial-frame update
and the generated-palette cycler step, plus preset transition bookkeeping. None
of that is per-sample work. The wall spread is under a millisecond, which is
the flywheel holding cadence rather than render jitter.

### Worst transition window (window frames 3345–3360)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  fx_shader_draw       26.86 ms  16.12 Mcyc   43%
  fx_prepare_frame     859.6 us  515.80 kcyc    1%
  fx_advance            2.23 ms   1.34 Mcyc    3%
  fx_timeline_step      44.9 us  26.97 kcyc    0%
  canvas_clear          86.1 us  51.69 kcyc    0%
  canvas_buffer_wait   32.37 ms  19.42 Mcyc   51%
```

Wall min/avg/max = 62.14/62.45/62.72 ms. During a
morph the parameter block is re-interpolated every frame and the palette mapping
is lerped, but that work lands in `fx_advance`, which is 2.23 ms here
against 2.27 ms held. `fx_shader_draw` reads 26.86 ms (43% of the frame),
0.12 ms below the worst hold: the shade cost tracks whichever pair of
presets it is between, so a transition is bounded by its two endpoints rather
than being a cost peak of its own.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the first
`Preset:` marker; the capture wraps back to preset 1, confirming a full cycle.

| bucket | preset | fx_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 4 | signal-weave-4 | 27.00 | 30.28 | 0/1,080 | 16 |
| 2 | signal-weave-2 | 26.98 | 30.61 | 0/1,409 | 16 |
| 1 | signal-weave | 26.93 | 30.47 | 0/1,080 | 16 |
| 0 | signal-weave (initial hold) | 26.87 | 31.95 | 0/119 | 16 |
| 3 | signal-weave-3 | 26.18 | 30.08 | 0/1,080 | 16 |

Span across presets is 1.03× of shader cost — every preset is green.

### Per-pixel figures

The dominant scope spends 1,561 cycles per quadrant sample over the
10,368 samples of a frame. The fixed pipeline writes the scan directly through
`Scan::Shader::draw_cached`, so there is no blended-pixel population and no
`filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.68/1.94/13.50 us  cpu 3.58%
isr_pack         144/frame  min/avg/max 6.23/6.76/9.06 us  cpu 1.55%
isr_dma_submit   144/frame  min/avg/max 0.69/0.93/5.82 us  cpu 0.21%
```

- Pack plus submit costs 1.11 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.34%, already absorbed by the render
  counters. With 30.55 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fx_shader_draw` — 43% of the frame, 26.98 ms: the inlined pullback shade
   over 10,368 samples. Stereographic projection through a Glitch lens, an
   outer `WaveShear` planar warp, a `Grid` source, and a triadic generated
   palette driven by noise hue.
2. `fx_advance` — 3% of the frame, 2.27 ms: preset transition bookkeeping and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 988 µs: per-frame snapshot of the
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
- `ANIMATED_PROJECTION` keeps the camera wander live in every preset, so there
  is no static hold: the reported per-preset cost is the worst clean-hold
  window within that preset's dwell.
- The capture used `-D HS_PROFILE_EPOCH_REVS=2600` to fit the cycle inside one
  epoch. Epoch length changes how long a preset holds, never its per-frame
  cost.
- The capture ran with the `HS_PROFILE` scopes this sweep added to the fixed-
  pipeline draw path, since landed as `ac6fe641`. They compile to nothing
  without `HS_PROFILE_ENABLE`, so the shipping image is unchanged.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=SignalWeave` and
`HS_PROFILE_WINDOW=16`; `just profile SignalWeave` performs the locked build,
flash, capture, marker validation, and artifact attestation.
