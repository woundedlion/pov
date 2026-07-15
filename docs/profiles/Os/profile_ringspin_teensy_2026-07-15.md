# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-07-15)

Point-in-time snapshot at the slim RingGroup tip (`730b6d2f` — inline
`stroke_alpha` per-blend eval + single covering-ring row intervals);
supersedes the fused-only tip (`56d8c854`) and the 2026-07-14 per-sub-stroke
report. Raw captures: `build/profile_capture_ringspin_slim_full.log`
(128-frame windows, 2 epochs), `build/profile_capture_ringspin_slim_w32.log`
(32-frame windows).

**Headline: the slim rework buys ~29% off per-blend scan cost at `-Os`,
lifting RingSpin from ~9.4 fps effective to ~14 fps** — mean render/frame
102 → 71 ms, and 8 of 27 windows now fully lock 16 fps (none did at the
fused-only tip). The heaviest-coverage orientations still flip to 2-window,
so `-Os` is not a full lock; `-O3` is (see the twin,
[../O3/profile_ringspin_teensy_2026-07-15.md](../O3/profile_ringspin_teensy_2026-07-15.md)).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` (`-D HS_PROFILE_WINDOW=32` for the fine pass) |
| Driver | `POVSegmented<288, 4, 480>` — the shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `RingSpin<288, 144>` at master `730b6d2f` — each trail frame's ≤4 sub-rings drawn as one fused `Scan::RingGroup` pass (76/frame), now with inline `SDF::Ring::stroke_alpha` coverage and one covering-ring row-interval producer |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 128 (resp. 32) frames then reset |
| Reproduce | `just profile RingSpin` |

Image size: FLASH code 42,228 B; ITCM (RAM1 code) 28,792 B; RAM2 free
4,736 B (+320 B FLASH / +320 B ITCM vs the fused-only tip — the inline
`stroke_alpha`).

**Exactness cross-check** — window frames 385–512: root counter
4,800,883,775 cyc = 8,001,472.9 µs vs measured `micros()` window sum
8,001,485 µs (Δ ≈ 1.5 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × 144-column half ≈ **10,368 px**). Wall time snaps up to a
whole number of 62.5 ms windows.

Render swings with the rings' orientation coverage. Over the 260 s capture
(27 128-frame windows): window-mean render 43.5–74.8 ms, wall-mean per
window 62.4–82.5 ms (grand mean **71.2 ms/frame, ≈14 fps effective**). **8
of 27 windows fully lock 16 fps** (wall < 65 ms); the other 19 are mixed
(individual heavy frames snap to 125 ms); **none** stay in the pure 8 fps
regime the fused-only tip showed. `rs_buffer_wait` is the round-up idle,
timed separately.

## Regime-by-regime readout (32-frame windows, per-frame averages)

Columns: time/frame, cycles/frame, % of frame, calls/frame, per-call.

### Locked 16 fps, heavy coverage (frames 97–128)

```
frame                   62.79 ms  37.67 Mcyc  100%
  rs_draw_rings         49.98 ms  29.99 Mcyc   80%
    rs_ring_scan        49.57 ms  29.74 Mcyc   79%  x76     652 us/group
      filter_blend      10.59 ms   6.36 Mcyc   17%  x42522  149 cyc/blend
  rs_buffer_wait        12.27 ms   7.36 Mcyc   20%
  rs_timeline_step      68.8 us   41.3 kcyc     0%
```

Wall min 46.3 / avg 62.8 / max 81.0 ms. The heaviest sustained coverage
(42.5 k blends/frame) renders in 49.6 ms — under the window — at **700 cyc
per blended pixel** (was ~960 at the fused-only tip, **−27%**); leaves
12.3 ms slack, so most frames hold one window.

### Hottest observed render (frames 33–64)

```
frame                  122.50 ms  73.50 Mcyc  100%
  rs_draw_rings         74.78 ms  44.87 Mcyc   61%
    rs_ring_scan        74.78 ms  44.87 Mcyc   61%  x76     984 us/group
      filter_blend      17.60 ms  10.56 Mcyc   14%  x65393  161 cyc/blend
  rs_buffer_wait        46.61 ms  27.97 Mcyc   38%
  rs_timeline_step      72.9 us   43.8 kcyc     0%
```

Wall min 82.4 / avg 122.5 / max 131.9 ms. Peak coverage (65.4 k
blends/frame, 6.3× quadrant, the trail-fill window) at **686 cyc/blend**;
render 74.8 ms still exceeds one window, so this regime stays 2-window even
slimmed — the residual that keeps `-Os` off a full lock.

### Locked 16 fps, light coverage (frames 225–256)

```
frame                   62.84 ms  37.71 Mcyc  100%
  rs_draw_rings         43.87 ms  26.32 Mcyc   70%
    rs_ring_scan        43.49 ms  26.10 Mcyc   69%  x76     572 us/group
      filter_blend       8.96 ms   5.38 Mcyc   14%  x36008  149 cyc/blend
  rs_buffer_wait        18.40 ms  11.04 Mcyc   29%
  rs_timeline_step      68.6 us   41.2 kcyc     0%
```

Wall min 36.4 / avg 62.8 / max 88.8 ms — 725 cyc/blend at the lightest
coverage; render 43.5 ms clears the window comfortably.

### Per-pixel figures

Slim path: **686–725 cyc (1.14–1.21 µs) per blended pixel** across regimes,
down from the fused-only tip's ~960 (−27% to −29%). Two changes contribute:
`SDF::Ring::stroke_alpha` replaces the per-slot out-of-line
`distance<false>()` call + `DistanceResult`/`Fragment` round trip with an
inline coverage eval, and one covering ring computes each row's intervals
once instead of per member (the previously-measured 9.3 ms/frame of
`sqrt`/`acos` collapses ~4×). `filter_blend` itself is ~149 cyc/blend.

## Column-ISR / DMA marshaling cost

Lighter than the fused-only tip in absolute terms because more frames run
1-window: `isr_wake` 4.9–5.1%, `isr_pack` 3.1–3.4%, `isr_dma_submit` 0.21%
CPU in the 2-window stretches; ~half that when locked. ISR ⇒ ~3.5 ms per
62.5 ms window ⇒ ~59 ms render budget. The sustained render (49.6 ms) now
clears it; only the trail-fill peak (74.8 ms) exceeds it — the one regime
still forcing 2-window at `-Os`.

## Summary ranking (heavy locked window, share of the 62.8 ms frame)

1. `rs_ring_scan` (fused `Scan::RingGroup::draw` × 76) — **79%** (49.6 ms)
2. display-window sync (`rs_buffer_wait`) — 20% (12.3 ms)
3. `filter_blend` (inside the scan) — 17% (10.6 ms)
4. everything else — <0.1%

Host figures for the slim change (perf_bench, full 288×144 frame):
`-Os` 11,725 → ~8,050 µs/frame, `-O3` 9,597 → ~6,830 µs/frame (both vs the
2026-07-14 pre-fused baseline; the fused-only tip was 4,983 / 7,219 — the
host `-Os` number *rose* from the fused-only tip because clang had already
inlined the old chain, so the slim eval trades a win it didn't need on host
for the one that matters on the M7). Trust the device counters.

## Caveats

- **ISR time is included**; `filter_blend` parents under `rs_ring_scan`.
- **Epoch straddle**: capture-boundary windows report more `frame` calls
  than window frames; skipped for all per-frame figures.
- The slim covering-ring scan exposes slightly more of the fast-path
  interval underestimate than the per-member scan did — low-alpha
  stroke-edge pixels, bounded by `test_ring_group_matches_sequential`.
- `-Os` shipping config; captured on clean master `730b6d2f`.

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=<EffectClass>`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile RingSpin [seconds]`.
