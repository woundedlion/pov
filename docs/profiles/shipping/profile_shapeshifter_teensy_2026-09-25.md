# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-09-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: [raw capture](../review20260925finding77/data/after-ship.log.txt). Replaces `profile_shapeshifter_teensy_2026-08-26.md`. The current roster has nine presets; the old report and skill's four-shape schedule describe earlier versions.

**Finding 77 attribution:** ShapeShifter sets `SAMPLED_RASTER_CONFIG.single_pass = true` in `effects/ShapeShifter.h:548` and passes it to `Plot::rasterize` at line 571. The `if constexpr (SINGLE_PASS)` branch in `core/render/plot/raster.h:648` returns at line 808, before the changed two-pass replay clamp at line 879. Consequently this effect never executes the changed clamp. These measurements describe ShapeShifter's overall image behavior, including run variability and possible compiler/code-layout differences; they do **not** measure the clamp's per-replay-sample execution cost. The hot raster symbols are not all byte-identical (`raster-symbols.json`), so no whole-code identity or statistical-zero claim is made.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel and DMA ISRs live, COM3 |
| Image | `profile`; -Os base with selective-O3 Plot cull/raster and shape-scan regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, clean tip `799ed81d189a2b007cb15dc227704606c72abc4f` |
| Method | `HS_PROFILE`, window 16, 155 s, epoch 1600 revolutions (200 s). Exact runtime frames 2–2458, frame 1 excluded. Complete-window scope/ISR summaries cover 17–2448. |
| Reproduce | From checkout `799ed81d1`: `HS_PROFILE_TREE="$PWD" HS_TEENSY_PORT=COM3 bash tools/profile_one.sh ShapeShifter profile 155 16 '-D HS_PROFILE_EPOCH_REVS=1600'` |

Image size: `FLASH: code:77112, data:152148, headers:8308` / `RAM1: variables:315072, code:49656, padding:15880, free:143680` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: frames 17–32, root 598793948 cycles / 600 = 997989.913 µs versus measured wall sum 997990 µs, **0.09 ppm** difference. The untouched capture passes `parse_profile.py ... validate`, including cycle wrap and monotonic frames.

## Frame cadence

**Runtime aggregate**, exact per-frame telemetry excluding setup: mean render 23.813 ms, peak render 58.395 ms, spilled 0/2457 (0.000%). Complete-window `ss_draw_all` averages 23.465 ms/frame; worst window 2353–2368 is 51.625 ms/frame.

Startup frame 1 renders in 49.821 ms and is excluded. A display interval is 62.5 ms; a quadrant is 10,368 pixels. All nine presets hold 16 fps without a spilled live frame. Peak margin to the interval is 4.105 ms. `ss_buffer_wait` is intentional idle until the next display flip.

## Phase-by-phase readout

The nine presets last 240 frames each, including eight-frame fade-in and fade-out envelopes; parameters snap in the dark. A full cycle is 2160 frames. Clean-hold windows have every frame outside those fade envelopes; boundary windows include fades and may contain adjacent hold frames.

### Clean hold (frames 2353–2368, worst regime window)

```text
frame                      62.41ms 37.44Mcyc 100.0%
  pov_preserve_half        137.1us  82.3kcyc   0.2% x1.0 137.1us/c
  ss_draw_all              51.63ms 30.98Mcyc  82.7%
    ss_plot_dispatch       51.42ms 30.85Mcyc  82.4% x136.9 375.7us/c
  ss_timeline_step          60.1us  36.0kcyc   0.1% x1.0 60.1us/c
  ss_buffer_wait           10.58ms  6.35Mcyc  17.0%
    canvas_clear            84.6us  50.8kcyc   0.1% x1.0 84.6us/c
    canvas_buffer_wait     10.49ms  6.30Mcyc  16.8% x1.0 10493.5us/c
```

Wall min/avg/max = 50.864/62.408/74.523 ms. `ss_draw_all` costs 51.625 ms/frame; render costs 51.829 ms/frame. Remaining wall time is predominantly display synchronization. Parent and child scope costs overlap and must not be added.

### Fade / preset boundary (frames 2161–2176, worst regime window)

```text
frame                      62.60ms 37.56Mcyc 100.0%
  pov_preserve_half        137.6us  82.6kcyc   0.2% x1.0 137.6us/c
  ss_draw_all              51.25ms 30.75Mcyc  81.9%
    ss_plot_dispatch       51.04ms 30.63Mcyc  81.5% x151.8 336.2us/c
  ss_timeline_step          61.7us  37.0kcyc   0.1% x1.0 61.7us/c
  ss_buffer_wait           11.14ms  6.69Mcyc  17.8%
    canvas_clear            84.5us  50.7kcyc   0.1% x1.0 84.5us/c
    canvas_buffer_wait     11.06ms  6.64Mcyc  17.7% x1.0 11059.6us/c
```

Wall min/avg/max = 60.943/62.598/65.003 ms. `ss_draw_all` costs 51.247 ms/frame; render costs 51.453 ms/frame. Remaining wall time is predominantly display synchronization. Parent and child scope costs overlap and must not be added.

### Per-preset table

The capture returns to preset 1 after presets 2–9. The nine rows and README N merge repeated visits after wrap; they count distinct authored presets, not visits. Rows use each preset’s worst clean-hold window, ranked by `ss_draw_all`; exact peak/spill columns include all live frames and following fade transitions. Windows are clean/owned complete windows.

| # | Preset | Windows | Blended px/f | ss_draw_all ms | Render ms | Hold fps | Live peak ms | Spilled |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | Planar star, 208 | 26/29 | unavailable | 51.625 | 51.829 | 16.0 | 58.395 | 0/478 |
| 9 | Flower, 72 | 13/15 | unavailable | 36.910 | 37.118 | 16.0 | 37.705 | 0/240 |
| 4 | Flower, 70 | 13/15 | unavailable | 35.596 | 35.800 | 16.0 | 37.372 | 0/240 |
| 8 | Spherical polygon, 144 / 3.195 sides | 13/15 | unavailable | 24.622 | 24.828 | 16.0 | 27.135 | 0/240 |
| 7 | Spherical polygon, 144 / 4 sides | 13/15 | unavailable | 20.119 | 20.330 | 16.0 | 22.171 | 0/240 |
| 6 | Spherical polygon, 128 | 13/15 | unavailable | 17.586 | 17.793 | 16.0 | 18.395 | 0/240 |
| 2 | Spherical polygon, 74.645 | 15/18 | unavailable | 13.028 | 13.242 | 16.0 | 15.546 | 0/299 |
| 5 | Planar star, 72 | 13/15 | unavailable | 11.251 | 11.470 | 16.0 | 13.804 | 0/240 |
| 3 | Planar star, 43.328 | 13/15 | unavailable | 8.679 | 8.884 | 16.0 | 10.007 | 0/240 |

### Per-pixel figures

This direct-write path emits no `filter_blend` or scan population counter in this timing image. Actual blended pixels and cycles per blended pixel are unavailable; dividing by quadrant size would not measure coverage. No additional per-pixel instrumentation was enabled.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.9/f  0.55/1.65/12.46 us  CPU 3.05%
isr_pack          144.0/f  6.23/6.82/10.05 us  CPU 1.57%
isr_dma_submit    144.0/f  0.59/0.93/2.71 us  CPU 0.21%
```

- Rate/frame and min/mean/max per call are shown; pack plus submit costs 1.117 ms/frame. Submission overhead is smaller than pixel packing.
- LED transmission is asynchronous at 24 MHz: the 72-LED image-plus-black composite is 600 bytes, approximately 230 µs including LPSPI byte framing, within the 434 µs column interval. These counters measure CPU marshaling, not wire occupancy.
- ISR share is 4.84% of elapsed capture windows, leaving approximately 59.48 ms CPU time per display interval. ISR time is already included in render measurements; do not subtract it twice. Neither regime needs a speedup to meet 16 fps.

## Summary ranking

1. `ss_draw_all` — 37.59% of root time, 23.465 ms/frame; inclusive scope.
2. `ss_plot_dispatch` — 37.42% of root time, 23.358 ms/frame; inclusive scope.
3. `pov_preserve_half` — 0.23% of root time, 0.142 ms/frame; inclusive scope.

No matched WASM/native timing capture was used. The separate finding-77 comparison pairs identical frame indices before and after the clamp.

## Caveats

- Counters include ISR time because CYCCNT free-runs.
- No `filter_blend` population is present here; generally it inherits its first active parent and disappears beneath an inactive parent.
- Deep per-pixel profiling is disabled; scope overhead remains in the measurements.
- Shipping uses selective O3 in Plot cull/raster and shape scanning; global O3 is a compiler reference, not the full-roster shipping policy.
- Epoch extension prevents reinitialization; it does not compress dwell or change per-frame work.
- Both source trees were clean; candidate differs from baseline only by finding 77.
- Runtime totals include trailing complete per-frame records beyond the last dumped scope window.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=1600`; `just profile ShapeShifter` routes through the locked harness. Use the explicit command above for the current nine-preset cycle.
