# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-09-25, **-O3**)

Shipping sibling: [selective-O3 report](../shipping/profile_shapeshifter_teensy_2026-09-25.md).

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: [raw capture](../review20260925finding77/data/after-o3.log.txt). Replaces `profile_shapeshifter_teensy_2026-08-26.md`. The current roster has nine presets; the old report and skill's four-shape schedule describe earlier versions.

**Finding 77 attribution:** ShapeShifter sets `SAMPLED_RASTER_CONFIG.single_pass = true` in `effects/ShapeShifter.h:548` and passes it to `Plot::rasterize` at line 571. The `if constexpr (SINGLE_PASS)` branch in `core/render/plot/raster.h:648` returns at line 808, before the changed two-pass replay clamp at line 879. Consequently this effect never executes the changed clamp. These measurements describe ShapeShifter's overall image behavior, including run variability and possible compiler/code-layout differences; they do **not** measure the clamp's per-replay-sample execution cost. The hot raster symbols are not all byte-identical (`raster-symbols.json`), so no whole-code identity or statistical-zero claim is made.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel and DMA ISRs live, COM4 |
| Image | `profile_o3`; global -O3 -ffast-math reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, clean tip `799ed81d189a2b007cb15dc227704606c72abc4f` |
| Method | `HS_PROFILE`, window 16, 155 s, epoch 1600 revolutions (200 s). Exact runtime frames 2–2458, frame 1 excluded. Complete-window scope/ISR summaries cover 17–2448. |
| Reproduce | From checkout `799ed81d1`: `HS_PROFILE_TREE="$PWD" HS_TEENSY_PORT=COM4 bash tools/profile_one.sh ShapeShifter profile_o3 155 16 '-D HS_PROFILE_EPOCH_REVS=1600'` |

Image size: `FLASH: code:106400, data:152028, headers:8836` / `RAM1: variables:315072, code:74152, padding:24152, free:110912` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: frames 17–32, root 598885699 cycles / 600 = 998142.832 µs versus measured wall sum 998144 µs, **1.17 ppm** difference. The untouched capture passes `parse_profile.py ... validate`, including cycle wrap and monotonic frames.

## Frame cadence

**Runtime aggregate**, exact per-frame telemetry excluding setup: mean render 23.806 ms, peak render 60.047 ms, spilled 0/2457 (0.000%). Complete-window `ss_draw_all` averages 23.468 ms/frame; worst window 2289–2304 is 51.255 ms/frame.

Startup frame 1 renders in 48.736 ms and is excluded. A display interval is 62.5 ms; a quadrant is 10,368 pixels. All nine presets hold 16 fps without a spilled live frame. Peak margin to the interval is 2.453 ms. `ss_buffer_wait` is intentional idle until the next display flip.

## Phase-by-phase readout

The nine presets last 240 frames each, including eight-frame fade-in and fade-out envelopes; parameters snap in the dark. A full cycle is 2160 frames. Clean-hold windows have every frame outside those fade envelopes; boundary windows include fades and may contain adjacent hold frames.

### Clean hold (frames 2289–2304, worst regime window)

```text
frame                      62.26ms 37.35Mcyc 100.0%
  pov_preserve_half        138.4us  83.1kcyc   0.2% x1.0 138.4us/c
  ss_draw_all              51.26ms 30.75Mcyc  82.3%
    ss_plot_dispatch       51.04ms 30.62Mcyc  82.0% x137.9 370.0us/c
  ss_timeline_step          57.9us  34.8kcyc   0.1% x1.0 57.9us/c
  ss_buffer_wait           10.80ms  6.48Mcyc  17.3%
    canvas_clear            84.5us  50.7kcyc   0.1% x1.0 84.5us/c
    canvas_buffer_wait     10.71ms  6.43Mcyc  17.2% x1.0 10714.3us/c
```

Wall min/avg/max = 50.429/62.256/76.248 ms. `ss_draw_all` costs 51.255 ms/frame; render costs 51.457 ms/frame. Remaining wall time is predominantly display synchronization. Parent and child scope costs overlap and must not be added.

### Fade / preset boundary (frames 2161–2176, worst regime window)

```text
frame                      62.22ms 37.33Mcyc 100.0%
  pov_preserve_half        138.9us  83.4kcyc   0.2% x1.0 138.9us/c
  ss_draw_all              47.84ms 28.70Mcyc  76.9%
    ss_plot_dispatch       47.61ms 28.57Mcyc  76.5% x165.6 287.5us/c
  ss_timeline_step          48.0us  28.8kcyc   0.1% x1.0 48.0us/c
  ss_buffer_wait           14.20ms  8.52Mcyc  22.8%
    canvas_clear            84.3us  50.6kcyc   0.1% x1.0 84.3us/c
    canvas_buffer_wait     14.11ms  8.47Mcyc  22.7% x1.0 14110.6us/c
```

Wall min/avg/max = 58.948/62.222/65.207 ms. `ss_draw_all` costs 47.835 ms/frame; render costs 48.028 ms/frame. Remaining wall time is predominantly display synchronization. Parent and child scope costs overlap and must not be added.

### Per-preset table

The capture returns to preset 1 after presets 2–9. The nine rows and README N merge repeated visits after wrap; they count distinct authored presets, not visits. Rows use each preset’s worst clean-hold window, ranked by `ss_draw_all`; exact peak/spill columns include all live frames and following fade transitions. Windows are clean/owned complete windows.

| # | Preset | Windows | Blended px/f | ss_draw_all ms | Render ms | Hold fps | Live peak ms | Spilled |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | Planar star, 208 | 26/29 | unavailable | 51.255 | 51.457 | 16.0 | 60.047 | 0/478 |
| 9 | Flower, 72 | 13/15 | unavailable | 35.636 | 35.834 | 16.0 | 36.307 | 0/240 |
| 4 | Flower, 70 | 13/15 | unavailable | 34.891 | 35.085 | 16.0 | 35.961 | 0/240 |
| 8 | Spherical polygon, 144 / 3.195 sides | 13/15 | unavailable | 24.377 | 24.569 | 16.0 | 25.919 | 0/240 |
| 7 | Spherical polygon, 144 / 4 sides | 13/15 | unavailable | 20.035 | 20.239 | 16.0 | 22.917 | 0/240 |
| 6 | Spherical polygon, 128 | 13/15 | unavailable | 17.272 | 17.467 | 16.0 | 18.522 | 0/240 |
| 2 | Spherical polygon, 74.645 | 15/18 | unavailable | 12.953 | 13.154 | 16.0 | 15.508 | 0/299 |
| 5 | Planar star, 72 | 13/15 | unavailable | 9.920 | 10.125 | 16.0 | 11.488 | 0/240 |
| 3 | Planar star, 43.328 | 13/15 | unavailable | 7.684 | 7.884 | 16.0 | 8.928 | 0/240 |

### Per-pixel figures

This direct-write path emits no `filter_blend` or scan population counter in this timing image. Actual blended pixels and cycles per blended pixel are unavailable; dividing by quadrant size would not measure coverage. No additional per-pixel instrumentation was enabled.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.9/f  0.32/1.49/16.74 us  CPU 2.75%
isr_pack          144.0/f  5.99/6.64/14.03 us  CPU 1.53%
isr_dma_submit    144.0/f  0.58/0.93/9.18 us  CPU 0.21%
```

- Rate/frame and min/mean/max per call are shown; pack plus submit costs 1.090 ms/frame. Submission overhead is smaller than pixel packing.
- LED transmission is asynchronous at 24 MHz: the 72-LED image-plus-black composite is 600 bytes, approximately 230 µs including LPSPI byte framing, within the 434 µs column interval. These counters measure CPU marshaling, not wire occupancy.
- ISR share is 4.49% of elapsed capture windows, leaving approximately 59.69 ms CPU time per display interval. ISR time is already included in render measurements; do not subtract it twice. Neither regime needs a speedup to meet 16 fps.

## Summary ranking

1. `ss_draw_all` — 37.59% of root time, 23.468 ms/frame; inclusive scope.
2. `ss_plot_dispatch` — 37.42% of root time, 23.359 ms/frame; inclusive scope.
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

Global -O3 vs selective -O3: mean render 23.813 → 23.806 ms (1.000×); peak 58.395 → 60.047 ms. O3 image minus shipping: FLASH code +29,288 B; ITCM +24,496 B.
