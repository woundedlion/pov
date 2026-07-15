# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw captures:
`build/profile_capture.log` (128-frame windows, 2 epochs),
`build/profile_capture_w32.log` (32-frame windows, 1 epoch + start of next),
`build/profile_capture_isr.log` (with column-ISR/DMA counters).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — the shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `DisplacementField<288, 144>` only (single-entry playlist), current working-tree state (`thickness` default 0.035, no init override) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 128 (resp. 32) frames then reset |
| Reproduce | `just profile DisplacementField` (any roster effect works: `just profile DreamBalls`) |

Image size: FLASH code 57,012 B; ITCM (RAM1 code) 39,288 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
window frames 385–416 root counter = 2,401,307,981 cyc = 4,002,179.97 µs vs
measured `micros()` window sum 4,002,182 µs (Δ ≈ 0.5 ppm).

## Frame cadence (context for every number below)

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering one quadrant
(this segment's 72-row band × the window's 144-column half ≈ 10,370 px —
`filter_blend` confirms ~10.9 k blended pixels/frame). A frame that misses one
window boundary snaps to the next, so wall time quantizes to multiples of
62.5 ms:

- **NOISE dwell / steady BALLS**: render > 62.5 ms → 2-window cadence, **125 ms/frame (8 fps)**
- **early BALLS / flat fade-in**: render < 62.5 ms → 1-window cadence, **62.5 ms/frame (16 fps)**

`df_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately, so the render numbers below are clean.

## Phase-by-phase readout (32-frame windows, per-frame averages)

Phase schedule at this 8 fps cadence: fade-in ≈ frames 0–150, NOISE dwell ≈
150–750, fade-out ≈ 750–900, BALLS from ≈ 900. Epoch (960 revolutions = 120 s)
reconstructs the effect at frame ≈ 1030.

### NOISE dwell — the hot phase (representative window: frames 385–416)

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame, calls/frame,
per-call time (cycles = µs × 600):

```
frame                  125.07 ms  75.04 Mcyc  100%
  df_timeline_step     105.51 ms  63.30 Mcyc   84%
    df_draw_rings      105.41 ms  63.25 Mcyc   84%
      df_ring_scan      86.79 ms  52.07 Mcyc   69%  x35.7  2430 us/ring
      df_lut_bake       12.77 ms   7.66 Mcyc   10%  x35.7   358 us/ring
      df_hue_table_prep  4.02 ms   2.41 Mcyc    3%  x22.1   182 us/call
      df_chunk_cull      1.55 ms   0.93 Mcyc    1%  x41.9    37 us/ring
  df_buffer_wait        19.56 ms  11.74 Mcyc   16%
  df_prepare_fields      0.14 us      86 cyc    0%
```

`df_ring_scan` is the SDF rasterize (`Scan::DistortedRing::draw`);
`df_buffer_wait` is the display-window sync idle.

Wall: min 123.8 / avg 125.1 / max 126.6 ms. Ring accounting per frame: 48
rings → ~42 survive the cap cull, ~36 survive the chunk cull and get baked +
scanned, 0 flat. `df_draw_rings` self time (prefilter, palette, basis) ≈
0.3 ms/frame. Dwell numbers are stable across the whole 600-frame hold
(ring_scan 80–85% of draw_rings in every window of both captures).

### Steady BALLS (representative window: frames 993–1024)

```
frame                  124.95 ms  74.97 Mcyc  100%
  df_timeline_step      79.50 ms  47.70 Mcyc   64%
    df_draw_rings       79.28 ms  47.57 Mcyc   63%
      df_ring_scan      58.42 ms  35.05 Mcyc   47%  x44.8  1303 us/ring
      df_lut_bake       13.24 ms   7.95 Mcyc   11%  x44.8   295 us/ring
      df_hue_table_prep  5.28 ms   3.17 Mcyc    4%  x40.1   132 us/call
      df_chunk_cull      1.58 ms   0.95 Mcyc    1%  x44.9    35 us/ring
      df_flat_scan       0.34 ms   0.21 Mcyc    0%  x 2.7   126 us/ring
  df_buffer_wait        45.44 ms  27.26 Mcyc   36%
  df_prepare_fields      9.3  us    5.6 kcyc    0%
```

(`df_flat_scan` = rings nothing displaces, drawn undistorted.)

Same 125 ms cadence as the dwell but ~26 ms/frame cheaper render: the ball
pool saturates the frame less than the noise field; per-ring scan cost is
about half the dwell's (narrower displaced bands), while slightly more rings
survive the cull (44.8 vs 35.7 — balls only reach a colatitude band, but the
band prefilter keeps far rings, which then bake flat-ish and scan cheap).

### Early BALLS — few balls in flight (window: frames 897–928)

```
frame                   62.82 ms  37.69 Mcyc  100%
  df_timeline_step      34.68 ms  20.81 Mcyc   55%
    df_draw_rings       34.56 ms  20.74 Mcyc   55%
      df_ring_scan      21.05 ms  12.63 Mcyc   34%  x13.4  1570 us/ring
      df_flat_scan       8.61 ms   5.16 Mcyc   14%  x25.0   345 us/ring
        filter_blend                                x10824  181 cyc/blend
      df_lut_bake        2.26 ms   1.36 Mcyc    4%  x13.4   169 us/ring
      df_hue_table_prep  1.96 ms   1.18 Mcyc    3%  x11.9   165 us/call
      df_chunk_cull      0.45 ms   0.27 Mcyc    1%  x14.3    32 us/ring
  df_buffer_wait        28.13 ms  16.88 Mcyc   44%
  df_prepare_fields      1.6  us   1.0  kcyc    0%
```

16 fps cadence (avg 62.8 ms wall, min 38.9). **25.0 rings/frame take the -Os
flat path** — a flat ring (345 µs) is ~4.5× cheaper than a displaced scan
(1.57 ms) in this window. (`filter_blend` counts every blended pixel across
all scan paths, not just its printed parent — see Caveats.)

### Fades

- **Fade-in (epoch start)**: first frames with `master_gain≈0` are the global
  minimum, **17.3–17.6 ms wall** (all 48 rings flat, 48 `df_flat_scan` calls,
  window `frames 1-32` of each epoch). Cost climbs with the amplitude ramp;
  by frames 65–96 the frame is dwell-priced (125 ms).
- **Fade-out (frames 833–896)**: render drops from 67.3 ms/frame to 50.0 ms
  as the noise band narrows; cadence mixes 1- and 2-window frames
  (min 55.7 / max 125.1 ms).

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h): ~10.9 k blended
pixels/frame in every steady window — the 48-ring stack effectively covers
the whole 10,368-px quadrant — at **195.8 cyc (0.33 µs) per blend** (includes
its own scope overhead). In the dwell, `df_ring_scan` ≈ 52.07 Mcyc over
~10.9 k blended px ⇒ ~4,780 cyc per *blended* pixel, i.e. most of the scan
cost is the knot-polyline SDF distance evaluation over the tested band, not
the blend.

## Column-ISR / DMA marshaling cost (ISR-safe counters, `build/profile_capture_isr.log`)

Measured with `HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree
is main-loop-only), steady dwell windows (16 s each). Nested as executed —
the pack and submit run inside the flywheel wake, on the 1-in-8 wakes that
render a column. Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake          18432/s  0.55 / 2.6  / ~180 us  cpu 4.5-5.2%
  isr_pack         2305/s  4.3  / 13.7 / ~180 us  cpu 2.9-3.4%
  isr_dma_submit   2305/s  0.61 / 0.97 / 3-9 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA call itself is ~1 µs of CPU**; the marshaling around it (pixel
  pack) is ~13–15 µs, so the whole per-column CPU cost is **~14–16 µs of the
  434 µs column period (~3.4%)**. The wire transfer (image + strobe-bg
  composite, ~220 µs at 24 MHz SPI) runs asynchronously on eDMA and costs no
  CPU — it never blocks rendering, it only bounds back-to-back submits.
- The pack avg ≈ 170 cyc/pixel is dominated by strided OCRAM reads (one
  display-buffer cache miss per LED, stride = one 288-px row). The ~180 µs
  maxima are preemption artifacts (higher-priority USB/DMA-completion ISRs
  landing inside the scope, inflated by the serial dumps themselves), not
  real pack work.
- Net: the ISR machinery steals **~5% of the chip**. Against the 62.5 ms
  window budget that's ~3.1 ms, leaving **~59 ms of render budget per
  window**. The dwell render (105 ms) needs ~1.8×, steady balls (80 ms)
  ~1.35×, to hit 1-window cadence with that headroom.

## Summary ranking (NOISE dwell, share of the 125 ms frame)

1. `Scan::DistortedRing::draw` (SDF rasterize + shade + blend) — **69%** (86.8 ms)
2. display-window sync (`df_buffer_wait`) — 16% (19.6 ms, idle by design)
3. LUT bake (ball/noise field + hue per column) — 10% (12.8 ms)
4. hue-table prep — 3% (4.0 ms)
5. chunk cull — 1% (1.6 ms)
6. everything else (`prepare_fields`, timeline, prefilter) — <1%

The rasterizer dominates on device even more than in the WASM survey
(perf ledger: NOISE dwell 7.45 ms full-frame on WASM vs ~105 ms/quadrant
here at -Os, ISRs live).

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs
  the flywheel/DMA/USB ISRs that fire inside it. That overhead is the real
  shipping condition (this is the point of segmented-mode profiling), but
  these are not pure-CPU algorithm costs.
- **`filter_blend` tree artifact**: it parents under whichever scope first
  enters it (here `df_flat_scan`), so its printed percentage is nonsense and
  the whole subtree is hidden in windows where no flat ring rendered
  (`log_node` skips zero-count parents). Its calls/cycles are correct.
- Per-pixel `filter_blend` scope overhead inflates the scan counters by
  ~0.3% of frame; the coarser `df_*` scopes are negligible.
- Epoch-boundary windows (`frame` count > window frames, e.g. 250 calls)
  mix the torn-down instance's undumped tail into the new instance's first
  window — per-call averages stay valid, per-frame sums don't.
- `-Os` build: the `__OPTIMIZE_SIZE__` flat-ring path is compiled in
  (it is the shipping Phantasm config); an -O3 profile would differ.
- Run with the uncommitted working-tree `thickness = 0.035` default.

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=<EffectClass>`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake,
  pixel pack, DMA submit) each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile <EffectClass> [seconds]`.
