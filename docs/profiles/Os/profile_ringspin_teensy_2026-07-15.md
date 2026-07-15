# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-07-15)

Point-in-time snapshot at the landed fused trail-scan tip (`56d8c854`);
replaces the 2026-07-14 per-sub-stroke report. Raw captures:
`build/profile_capture_ringspin_fused.log` (128-frame windows, 2 epochs),
`build/profile_capture_ringspin_fused_w32.log` (32-frame windows).

**Headline: the fused scan is device-neutral at `-Os`.** Render, per-blend
cost, and cadence all match the 2026-07-14 baseline within orientation
noise — RingSpin's `-Os` frame is bound by per-blended-pixel work
(~950–960 scan cyc/blend across ~42 k blends of 4× overdraw), not by the
per-pass row/interval/walk overhead the fusion removed. That overhead is
what host builds pay for (host perf_bench: −57% at clang `-Os`, −25% at
`-O3`), so the host win does not transfer here. The fusion does pay on
device at `-O3` — see the twin report
([../O3/profile_ringspin_teensy_2026-07-15.md](../O3/profile_ringspin_teensy_2026-07-15.md)),
where it converts the borderline 16↔8 cadence into a **16 fps lock**.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` (`-D HS_PROFILE_WINDOW=32` for the fine pass) |
| Driver | `POVSegmented<288, 4, 480>` — the shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `RingSpin<288, 144>` at master `56d8c854` — 4 rings × 19-sample trail, each trail frame's ≤4 sub-rings drawn as **one fused `Scan::RingGroup` pass** (76 group passes/frame vs 224 sub-strokes on the legacy path) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 128 (resp. 32) frames then reset |
| Reproduce | `just profile RingSpin` |

Image size: FLASH code 41,908 B; ITCM (RAM1 code) 28,472 B; RAM2 free
4,736 B (+2,080 B FLASH / +2,080 B ITCM vs the pre-fused tip — the
`RingGroup` instantiation).

**Exactness cross-check** — window frames 257–384 (epoch 1): root counter
9,598,585,987 cyc = 15,997,643.3 µs vs measured `micros()` window sum
15,997,662 µs (Δ ≈ 1.2 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × 144-column half ≈ **10,368 px**). Wall time snaps up to a
whole number of 62.5 ms windows.

RingSpin has no preset schedule — cost swings with how much of each ring's
great circle falls in the segment band. Render spans 54–87 ms per 32-frame
window, straddling the 62.5 ms window exactly as the baseline did: locked
2-window stretches (frames 225–384: 125 ms, 8 fps), mixed windows (avg
75–109 ms), effective ≈ **9.4 fps** over the epoch. `rs_buffer_wait` (the
`buffer_free()` spin inside the `Canvas` ctor) is the round-up idle, timed
separately, so render numbers are clean.

## Regime-by-regime readout (per-frame averages)

No phase schedule; the three windows below bracket the orientation-driven
range. Columns: time/frame, cycles/frame, % of frame, calls/frame, per-call.

### Sustained locked 8 fps (128-frame window, frames 257–384)

Identical in both epochs:

```
frame                  124.98 ms  74.99 Mcyc  100%
  rs_draw_rings         67.58 ms  40.55 Mcyc   54%
    rs_ring_scan        66.76 ms  40.05 Mcyc   53%  x76     879 us/group
      filter_blend      11.24 ms   6.74 Mcyc    9%  x41710  162 cyc/blend
  rs_buffer_wait        57.33 ms  34.39 Mcyc   46%
  rs_timeline_step      71.1 us   42.7 kcyc     0%
```

Wall min 123.7 / avg 125.0 / max 126.4 ms. `rs_ring_scan` (now the fused
`Scan::RingGroup::draw`, one call per trail frame) is 99% of render; scan
cost per blended pixel is **960 cyc** — the same as the 2026-07-14
per-sub-stroke baseline (~975), which is the whole story of this report.

### Hottest observed render (32-frame window, frames 33–64)

```
frame                  122.36 ms  73.41 Mcyc  100%
  rs_draw_rings         87.25 ms  52.35 Mcyc   71%
    rs_ring_scan        86.39 ms  51.84 Mcyc   71%  x76    1137 us/group
      filter_blend      16.26 ms   9.75 Mcyc   13%  x61207  159 cyc/blend
  rs_buffer_wait        35.04 ms  21.02 Mcyc   29%
  rs_timeline_step      73.0 us   43.8 kcyc     0%
```

Wall min 92.5 / avg 122.4 / max 129.2 ms. Worst-case orientations push
coverage to 61 k blends/frame (5.9× the quadrant); per-blend scan cost
*drops* to 847 cyc (denser runs amortize the walk), so cost tracks coverage
sub-linearly.

### Mixed cadence (32-frame window, frames 449–480)

```
frame                   75.31 ms  45.19 Mcyc  100%
  rs_draw_rings         54.77 ms  32.86 Mcyc   73%
    rs_ring_scan        54.02 ms  32.41 Mcyc   72%  x76     711 us/group
      filter_blend       9.93 ms   5.96 Mcyc   13%  x39132  152 cyc/blend
  rs_buffer_wait        20.47 ms  12.28 Mcyc   27%
  rs_timeline_step      71.4 us   42.9 kcyc     0%
```

Wall min 48.1 / avg 75.3 / max 123.9 ms — render dips under the window on
favorable orientations, so frames oscillate between 1 and 2 windows.

### Per-pixel figures

Sustained window: 41.7 k blends/frame ≈ **4.0× the 10,368-px quadrant**
(the trail's overdraw) at **162 cyc (0.27 µs) per blend** — down from the
baseline's 183 cyc; the union walk blends a pixel's overlapping sub-rings
back-to-back, which the write path likes. `rs_ring_scan` ≈ 40.05 Mcyc over
41.7 k blended px ⇒ **960 cyc per blended pixel** (847 in the hottest
window) — unchanged from the baseline's ~975, because the per-blend chain
(out-of-line `SDF::Ring::distance<false>` call, `DistanceResult` +
`Fragment` stores, `quintic_kernel`, plot) dominates at GCC `-Os` and the
fusion only removed the per-pass row/interval/walk overhead around it.

## Column-ISR / DMA marshaling cost

Unchanged from the baseline (same driver): sustained 2-window windows show
`isr_wake` 5.7–7.1%, `isr_pack` 3.9–5.4%, `isr_dma_submit` 0.22% CPU;
1-window stretches roughly halve the shares. ISR machinery ⇒ ~3.9 ms per
62.5 ms window ⇒ **~58.6 ms render budget**. The sustained render (67.6 ms)
still needs ~1.15×, and the hottest window (87.3 ms) ~1.49×, to lock
16 fps at `-Os` — the fusion did not close that gap; cheaper per-blend work
(or `-O3`, see the twin) is what does.

## Summary ranking (sustained window, share of the 125 ms frame)

1. `rs_ring_scan` (fused `Scan::RingGroup::draw` × 76) — **53%** (66.8 ms;
   71% / 86.4 ms in the hottest window)
2. display-window sync (`rs_buffer_wait`) — 46% (57.3 ms, idle by design)
3. `filter_blend` (inside the scan) — 9% (11.2 ms)
4. everything else (timeline, trail record, basis, palette) — <0.1%

Host figures for the same change (perf_bench, full 288×144 frame):
`-Os` 11,725 → 4,983 µs/frame (−57%), `-O3` 9,597 → 7,219 µs/frame
(−25%) — clang inlines the per-blend chain, leaving hosts walk-bound where
the M7 at GCC `-Os` is blend-bound. Treat host wins on scan-path effects
as upper bounds until profiled here.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope absorbs the
  flywheel/DMA/USB ISRs firing inside it.
- **`filter_blend` parenting**: nests under `rs_ring_scan` (first entrant);
  calls/cycles are correct, printed percentage is scan-relative.
- **Epoch straddle**: capture-boundary windows report more `frame` calls
  than window frames; skipped for all per-frame figures. The first ~20
  frames of each epoch draw a filling trail.
- The fused path's output diverges from per-sub-stroke draws only by
  AA-tail pixels the per-ring interval clip used to drop (bounded by
  `test_ring_group_matches_sequential`).
- `-Os` shipping config; captured on clean master `56d8c854`.

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=<EffectClass>`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile RingSpin [seconds]`.
