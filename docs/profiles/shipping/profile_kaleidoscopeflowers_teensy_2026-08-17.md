# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-08-17, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeFlowers`).
Raw capture: `build/prof/kaleidoscopeflowers_ship.log`.
Replaces `profile_kaleidoscopeflowers_teensy_2026-08-16.md`: re-captured after the
pullback instance-pipeline refactor and the `4d93e449` stage-dispatch fix, on
the same board and with the same capture parameters, to confirm the fix holds
for a look with no surface stage. It does — the per-sample de-inlining that cost
LatticeMelt 16 ms/frame is absent here — but this look shades **1.5% slower**
than the 08-16 baseline rather than faster. See "Comparison to the 08-16
baseline".

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeFlowers 288×144, single-entry playlist, tip `1e717361` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, -D HS_PROFILE_EPOCH_REVS=2400 |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeFlowers profile 260 16 -D HS_PROFILE_EPOCH_REVS=2400` |

Image sections (`arm-none-eabi-size -A`):

```text
.text.code     42,528   .text.progmem  145,672
.text.itcm     25,792   .data            4,800
.bss          310,272   .bss.dma       520,064
```

Exactness cross-check: window frames 1105-1120 root counter cycles ÷ 600 MHz
match the measured wall sum within **0.4 ppm**. Validation reports
258 complete windows, monotonic frame numbers, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py <log> windows` footer): peak frame render
**36.51 ms** at frames 1905-1920, spilled **0/4,128 frames** (0%).

A display window is 62.5 ms. KaleidoscopeFlowers is a strobe effect rendering one
quadrant ≈ 10,368 samples per frame; it sets neither `needs_full_frame()`
nor `persists_pixels()`, so there is no full-canvas multiplier. The worst frame
of the pass keeps **25.99 ms** of the window in hand, so the effect holds
16 fps for the whole capture. `canvas_buffer_wait` is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: each preset holds for 600 frames and then morphs over a
480-frame transition, so one preset owns 1,080 frames. The two regimes below are
captured at the same window indices the 08-16 report used, so the trees compare
line for line; the spread across presets is in the per-preset table.

### Held window (window frames 1345–1360)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  fx_shader_draw       32.40 ms  19.44 Mcyc   51%
  fx_prepare_frame     927.9 us  556.76 kcyc    1%
  fx_advance            1.88 ms   1.13 Mcyc    3%
  fx_timeline_step      42.4 us  25.48 kcyc    0%
  canvas_clear          86.2 us  51.74 kcyc    0%
  canvas_buffer_wait   27.10 ms  16.26 Mcyc   43%
```

Wall min/avg/max = 60.68/62.45/64.23 ms — indistinguishable from the baseline's
60.74/62.43/64.23, so cadence is unchanged and the flywheel is still setting the
frame period. `fx_shader_draw` is 51% of the frame at 32.40 ms against the
baseline's 31.90 ms in this same window. Every scope *outside* the per-sample
path got cheaper: `fx_prepare_frame` 927.9 µs against 944.2 µs and `fx_advance`
1.88 ms against 1.99 ms, both because the per-frame prepare pass now resolves
state once instead of the runtime recomputing it. The regression is confined to
the shade itself.

### Transition window (window frames 1905–1920)

```text
frame                  62.44 ms  37.46 Mcyc  100%
  fx_shader_draw       31.83 ms  19.10 Mcyc   50%
  fx_prepare_frame      1.30 ms  780.86 kcyc    2%
  fx_advance            2.04 ms   1.22 Mcyc    3%
  fx_timeline_step     115.6 us  69.36 kcyc    0%
  canvas_clear          86.5 us  51.91 kcyc    0%
  canvas_buffer_wait   27.05 ms  16.23 Mcyc   43%
```

Wall min/avg/max = 60.35/62.44/64.56 ms. This window holds the pass peak
(36.51 ms). `fx_timeline_step` reads 115.6 µs against the baseline's 53.2 µs:
preset transitions now run through `Animation::Lerp` over the field descriptor
tables rather than a hand-rolled per-field lerp, and the table walk is the extra
62 µs. That is 0.1% of the frame and only during a morph, so it does not move
cadence. `fx_prepare_frame` is 1.30 ms against 1.32 ms.

### Per-preset table

Rows are ranked by clean-hold shader cost. A preset owns the frames it was on
screen for, including the transition that follows it, so the spilled column is
stricter than a hold-only figure. Bucket 0 is the initial hold before the first
`Preset:` marker; open-grid owns two holds because the cycle wraps through it
twice. The markers run 2/3 → 3/3 → 1/3 → 2/3, confirming a full cycle.

| bucket | preset | fx_shader_draw ms | peak render ms | spilled/frames | fps |
|---:|---|--:|--:|--:|--:|
| 2 | open-grid | 32.81 | 36.22 | 0/1,371 | 16 |
| 1 | double-map | 32.25 | 36.14 | 0/1,079 | 16 |
| 0 | double-map (initial hold) | 32.02 | 36.19 | 0/599 | 16 |
| 3 | fine-grid | 31.89 | 36.51 | 0/1,079 | 16 |

Span across presets is 1.03× of shader cost — every preset is green.

## Comparison to the 08-16 baseline

Same board, same capture parameters, same window indices. The 08-16 row was
captured at `69d4751c`, before the pullback refactor.

| preset | 08-16 shader ms | 08-17 shader ms | delta |
|---|--:|--:|--:|
| open-grid | 32.25 | 32.81 | +0.56 |
| double-map | 32.07 | 32.25 | +0.18 |
| fine-grid | 31.39 | 31.89 | +0.50 |
| double-map (initial hold) | 31.32 | 32.02 | +0.70 |
| **peak render** | **35.91** | **36.51** | **+0.60** |

The regression is small (+1.5% mean on the shade, +1.7% on peak render),
consistent across all four buckets, and in the opposite direction from
LatticeMelt, which the same tip shades 11% *faster*. Both directions have the
same cause: the instance pipeline hoists per-frame state out of the per-sample
path. LatticeMelt had expensive hoistable work — resolving the curl loop point
once per frame instead of 10,368 times — so it won big. KaleidoscopeFlowers's stages
(a `Grid` source, a `MirrorTile` warp, an `Equirectangular` projection) have
almost nothing to hoist, so it collects none of the benefit and pays the cost of
threading a wider `Frame` (`{FrameState ctx; PreparedTuple prepared}`) through
an out-of-line `shade`, roughly 29 cycles per sample.

An ELF symbol diff against a baseline build of `69d4751c` confirms no per-sample
helper was de-inlined by the refactor: `Pullback::Source::grid` and
`Pipeline::shade` were **already** out of line at the baseline — that is this
look's normal shape, since it reaches shade through `Scan::Shader::draw_cached`
lambdas and has no `Placed::run`. The only symbols the refactor newly pushed out
of line are `prepare_stages`, `Source::prepare` and `Warp::prepare`, all of which
run once per frame by design. The catastrophic tell — a per-sample helper newly
appearing as its own symbol — is absent.

Both numbers are single captures, so treat ±0.1 ms as noise; the consistency
across four buckets and two window indices is what makes the direction credible.

### Per-pixel figures

The dominant scope spends 1,899 cycles per quadrant sample over the
10,368 samples of a frame, read from the worst clean hold (open-grid); the same
basis gives 1,866 at the 08-16 baseline. The fixed pipeline writes the scan
directly through `Scan::Shader::draw_cached`, so there is no blended-pixel
population and no `filter_blend` subtree in this capture.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.73/2.00/13.66 us  cpu 3.68%
isr_pack         144/frame  min/avg/max 6.33/7.04/9.23 us  cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.65/0.93/5.62 us  cpu 0.21%
```

- Pack plus submit costs 1.15 ms of CPU per frame.
- The LED wire transfer itself stays asynchronous; only marshaling is on the CPU.
- Total ISR CPU share is 5.51%, already absorbed by the render
  counters. With 25.99 ms of margin at the worst frame, no speedup is
  required for cadence.

## Summary ranking

1. `fx_shader_draw` — 51% of the frame, 32.40 ms: the pullback shade over
   10,368 samples. An `Equirectangular` projection through a
   `DodecahedralKaleidoscope` lens, an inner `MirrorTile` warp, and a `Grid`
   source, with a generated-palette colour stage over noise-driven hue.
2. `fx_advance` — 3% of the frame, 1.88 ms: runtime clocks, spatial frames and
   the generated-palette cycler step.
3. `fx_prepare_frame` — 1% of the frame, 928 µs: the frame-context snapshot,
   the hue LUTs, and the pipeline `prepare_stages` pass.

The perf ledger has no WASM or native baseline for this effect, so there is no
cross-target comparison to make yet.

## Caveats

- All scopes absorb live ISR time because CYCCNT free-runs.
- `filter_blend` does not appear in this direct scan; if it were entered it would
  parent under whichever scope reached it first, and its subtree would be hidden
  in windows where that parent had 0 calls.
- Shipping selective O3: `Scan::Shader::draw_cached` is `HS_O3_FN` with cached-
  flash placement, so the whole shade path compiles at `-O3`; the OKLab
  and gamut helpers the generated palette bakes through are `HS_O3_FN` as well.
  The rest of the image keeps the `-Os` base policy.
- Unlike LatticeMelt, this look's `shade` and `Source::grid` are out-of-line at
  every tip measured, so the standalone-symbol regression tell does not apply
  here; the meaningful check is whether a *new* per-sample symbol appears.
- The capture used `-D HS_PROFILE_EPOCH_REVS=2400` to fit the cycle inside one
  epoch. Epoch length changes how long a preset holds, never its per-frame
  cost.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=KaleidoscopeFlowers` and
`HS_PROFILE_WINDOW=16`; `just profile KaleidoscopeFlowers` performs the locked build,
flash, capture, marker validation, and artifact attestation.
