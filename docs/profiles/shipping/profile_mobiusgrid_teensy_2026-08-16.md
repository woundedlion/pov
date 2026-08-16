# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-08-16, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MobiusGrid`).
Raw capture: `build/prof/mobiusgrid_ship.log`.
This replaces `profile_mobiusgrid_teensy_2026-07-25.md`, which profiled the
pre-migration curve-raster MobiusGrid; the effect of that name is now a
fixed-pipeline look and shares nothing with the old render path.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, tip `69d4751cc077` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, -D HS_PROFILE_EPOCH_REVS=1600 |
| Reproduce | `bash tools/profile_one.sh MobiusGrid profile 170 16 -D HS_PROFILE_EPOCH_REVS=1600` |

Image size:

```text
FLASH: code:67,688, data:149,992, headers:8,624
       free for files:1,805,312
RAM1:  variables:315,072, code:26,584, padding:6,184
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 785-800 root counter cycles ÷ 600 MHz
match the measured wall sum within **2.8 ppm**. Validation reports
168 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**24.20 ms** at frames 801-816, spilled **0/2,688 frames** (0%).


A display window is 62.5 ms. MobiusGrid is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **38.30 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns 1,080 frames. The regimes below are
the worst held window and the worst transition window; the spread across
presets is in the per-preset table.

### Worst held window (window frames 2145–2160)

```text
frame                  62.41 ms  37.45 Mcyc  100%
  fx_shader_draw       20.86 ms  12.51 Mcyc   33%
  fx_prepare_frame     654.7 us  392.83 kcyc    1%
  fx_advance            2.21 ms   1.33 Mcyc    3%
  fx_timeline_step      63.0 us  37.83 kcyc    0%
  canvas_clear          85.8 us  51.51 kcyc    0%
  canvas_buffer_wait   38.53 ms  23.12 Mcyc   61%
```

Wall min/avg/max = 61.94/62.41/62.68 ms. `fx_shader_draw` is 33% of the
frame at 20.86 ms; the whole render is 23.88 ms and the remaining
38.53 ms (61%) is `canvas_buffer_wait` idle. `fx_prepare_frame` costs
655 µs and `fx_advance` 2.21 ms — the runtime advance, the spatial-frame update
and the generated-palette cycler step, plus preset transition bookkeeping. None
of that is per-sample work. The wall spread is under a millisecond, which is
the flywheel holding cadence rather than render jitter.

### Worst transition window (window frames 801–816)

```text
frame                  62.44 ms  37.47 Mcyc  100%
  fx_shader_draw       20.91 ms  12.55 Mcyc   33%
  fx_prepare_frame     626.8 us  376.07 kcyc    1%
  fx_advance            2.30 ms   1.38 Mcyc    3%
  fx_timeline_step      59.4 us  35.64 kcyc    0%
  canvas_clear          85.8 us  51.50 kcyc    0%
  canvas_buffer_wait   38.45 ms  23.07 Mcyc   61%
```

Wall min/avg/max = 62.17/62.44/62.76 ms. During a
morph the parameter block is re-interpolated every frame and the palette mapping
is lerped, but that work lands in `fx_advance`, which is 2.30 ms here
against 2.21 ms held. `fx_shader_draw` reads 20.91 ms (33% of the frame),
0.06 ms above the worst hold: the shade cost tracks whichever pair of
presets it is between, so a transition is bounded by its two endpoints rather
than being a cost peak of its own.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the first
`Preset:` marker; the capture wraps back to preset 1, confirming a full cycle.

| bucket | preset | fx_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 2 | mobius-grid-2 | 20.92 | 24.20 | 0/1,080 | 16 |
| 1 | mobius-grid | 20.90 | 24.15 | 0/1,009 | 16 |
| 0 | mobius-grid (initial hold) | 20.60 | 23.50 | 0/599 | 16 |

Span across presets is 1.02× of shader cost — every preset is green.

### Per-pixel figures

The dominant scope spends 1,207 cycles per quadrant sample over the
10,368 samples of a frame. The fixed pipeline writes the scan directly through
`Scan::Shader::draw_cached`, so there is no blended-pixel population and no
`filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1151/frame  min/avg/max 0.73/1.97/16.46 us  cpu 3.63%
isr_pack         143/frame  min/avg/max 6.23/6.88/9.32 us  cpu 1.58%
isr_dma_submit   143/frame  min/avg/max 0.77/0.93/9.03 us  cpu 0.21%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.42%, already absorbed by the render
  counters. With 38.30 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fx_shader_draw` — 33% of the frame, 20.86 ms: the inlined pullback shade
   over 10,368 samples. Stereographic projection through a parameterized
   `Mobius` lens, an inner `MirrorTile` warp, and a `TwinWave` source under a
   complementary path-length palette.
2. `fx_advance` — 3% of the frame, 2.21 ms: preset transition bookkeeping and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 655 µs: per-frame snapshot of the
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
- `ANIMATED_MOBIUS` keeps the Mobius warp animating through preset morphs, so
  the lens parameters differ frame to frame within a single held preset.
- The capture used `-D HS_PROFILE_EPOCH_REVS=1600` to fit the cycle inside one
  epoch. Epoch length changes how long a preset holds, never its per-frame
  cost.
- The capture ran with the `HS_PROFILE` scopes this sweep added to the fixed-
  pipeline draw path, since landed as `ac6fe641`. They compile to nothing
  without `HS_PROFILE_ENABLE`, so the shipping image is unchanged.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=MobiusGrid` and
`HS_PROFILE_WINDOW=16`; `just profile MobiusGrid` performs the locked build,
flash, capture, marker validation, and artifact attestation.
