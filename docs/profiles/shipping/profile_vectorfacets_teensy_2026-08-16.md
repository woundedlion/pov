# VectorFacets on-device profile — Teensy 4.0, segmented mode (2026-08-16, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile VectorFacets`).
Raw capture: `build/prof/vectorfacets_ship.log`.
First profile of this effect: it entered the roster with the fixed-pipeline
workbench migration (`69d4751c`) and had no prior report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | VectorFacets 288×144, single-entry playlist, tip `69d4751cc077` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh VectorFacets profile 70 32` |

Image size:

```text
FLASH: code:69,720, data:153,688, headers:9,040
       free for files:1,799,168
RAM1:  variables:315,072, code:26,600, padding:6,168
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 769-800 root counter cycles ÷ 600 MHz
match the measured wall sum within **3.5 ppm**. Validation reports
34 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**47.20 ms** at frames 769-800, spilled **0/1,088 frames** (0%).


A display window is 62.5 ms. VectorFacets is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **15.30 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: VectorFacets has a single preset and one steady regime — the
camera and warp animation vary the shade cost slightly frame to frame but
never change the pipeline. The window below is the costliest of the pass.

### Steady state (window frames 769–800)

```text
frame                  62.44 ms  37.46 Mcyc  100%
  fx_shader_draw       43.36 ms  26.02 Mcyc   69%
  fx_prepare_frame      1.06 ms  636.80 kcyc    1%
  fx_advance            2.20 ms   1.32 Mcyc    3%
  fx_timeline_step      65.8 us  39.47 kcyc    0%
  canvas_clear          87.8 us  52.71 kcyc    0%
  canvas_buffer_wait   15.65 ms   9.39 Mcyc   25%
```

Wall min/avg/max = 61.70/62.44/63.15 ms. `fx_shader_draw` is 69% of the
frame at 43.36 ms; the whole render is 46.79 ms and the remaining
15.65 ms (25%) is `canvas_buffer_wait` idle. `fx_prepare_frame` costs
1.06 ms and `fx_advance` 2.20 ms — the runtime advance, the spatial-frame
update and the generated-palette cycler step. None of that is per-sample work.
The wall spread is under a millisecond, which is the flywheel holding cadence
rather than render jitter.

### Per-pixel figures

The dominant scope spends 2,509 cycles per quadrant sample over the
10,368 samples of a frame. The fixed pipeline writes the scan directly through
`Scan::Shader::draw_cached`, so there is no blended-pixel population and no
`filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1151/frame  min/avg/max 0.68/2.02/17.16 us  cpu 3.73%
isr_pack         143/frame  min/avg/max 6.31/7.22/9.62 us  cpu 1.66%
isr_dma_submit   143/frame  min/avg/max 0.60/0.93/9.30 us  cpu 0.21%
```

- Pack plus submit costs 1.17 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.60%, already absorbed by the render
  counters. With 15.30 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fx_shader_draw` — 69% of the frame, 43.36 ms: the inlined pullback shade
   over 10,368 samples. Folded-hemisphere Gnomonic projection through a
   `DodecahedralKaleidoscope` lens, an outer `VectorNoise` warp with a flat
   envelope, an inner `MirrorTile` warp, and a `Grid` source.
2. `fx_advance` — 3% of the frame, 2.20 ms: preset transition bookkeeping and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 1.06 ms: per-frame snapshot of the
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
- `Warp::VectorNoise` calls `FastNoiseLite::GetVectorNoiseSingle`, which is
  `HS_HOT_FLASH_MEMBER` — `HS_O3_FN` plus hot/noinline flash placement — so
  that leaf is already `-O3` in the shipping image. This is the heaviest effect
  in the fixed-pipeline group and the only one whose warp stage carries a per-
  sample noise lookup.
- The capture ran with the `HS_PROFILE` scopes this sweep added to the fixed-
  pipeline draw path, since landed as `ac6fe641`. They compile to nothing
  without `HS_PROFILE_ENABLE`, so the shipping image is unchanged.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=VectorFacets` and
`HS_PROFILE_WINDOW=32`; `just profile VectorFacets` performs the locked build,
flash, capture, marker validation, and artifact attestation.
