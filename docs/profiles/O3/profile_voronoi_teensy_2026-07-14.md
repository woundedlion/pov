# Voronoi on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_voronoi_o3.log` (64-frame windows, ~75 s single pass; ISR/DMA
accumulators dumped every window). This is the **-O3** twin of the
shipping-config `-Os` report at `../Os/profile_voronoi_teensy_2026-07-14.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`, N=4) + `-D HS_PROFILE_ENABLE`, but **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — this board strapped as segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `Voronoi<288, 144>` only (single-entry playlist) |
| Method | DWT CYCCNT cycle counters (HS_PROFILE RAII scopes, 600 cyc/µs) + micros() wall clock per draw_frame, dumped every 64 frames |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 54,532 B; ITCM (RAM1 code) 40,264 B; RAM2 free 4,736 B.

**Exactness cross-check** — window frames 65–128 root counter = 2,395,303,416
cyc = 3,992,172.4 µs vs measured `micros()` window sum 3,992,176 µs (Δ ≈ 0.9 ppm).

## Frame cadence

62.5 ms display window; one quadrant ≈ 10,368 px. Voronoi is unusual — no
timeline, no `Scan::Shader`: it open-codes a per-frame prep (site spin +
KD-tree rebuild) and a per-pixel KD-nearest shading loop.

**Cadence changes at -O3.** At -Os the KD-nearest shade (~77 ms) overflowed one
window → 125 ms/frame (8 fps). At -O3 the shade drops to ~57 ms, **under the
62.5 ms window**, so Voronoi runs a steady **16 fps (62.5 ms/frame, 1-window)** —
a full cadence tier gained. `vo_buffer_wait` is the round-up idle.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Steady; representative window frames 65–128. Columns: time/frame, cyc/frame,
% of frame, per-pixel cost on the shade leaf.

```
frame                62.38 ms  37.43 Mcyc  100%
  vo_shade           57.17 ms  34.30 Mcyc   92%  10368 px  3308 cyc/px
  vo_buffer_wait      4.96 ms   2.98 Mcyc    8%
  vo_kdtree           0.18 ms  108 kcyc      0%
  vo_animate          64 us      38 kcyc      0%
```

Wall: min 56.5 / avg 62.4 / max 68.6 ms. `vo_shade` (coarse-grid corner pre-pass
+ per-pixel KD-nearest query and shade) is the entire cost — 57.2 ms shading the
10,368-px quadrant at **3,308 cyc (5.51 µs) per pixel**. The prep is negligible:
the whole per-frame KD-tree rebuild is 0.18 ms (0.3%) and the site spin
`vo_animate` is 64 µs. Now that the shade fits one window, `vo_buffer_wait`
collapses to ~8% (from 38% at -Os) — the effect is barely idle at 16 fps.

### Per-pixel figures

The open-coded shade loop composites in place (no `filter_blend`), so per-pixel
cost is from full quadrant coverage: **10,368 px/frame at 3,308 cyc (5.51 µs)
per pixel** — dominated by the KD-tree nearest-neighbour descent per pixel.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE`, window frames 65–128. rate, per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.51 / 1.50 / 10 us   cpu 2.76%
  isr_pack        2304/s  4.64 / 5.22 / 7.9 us  cpu 1.20%
  isr_dma_submit  2304/s  0.60 / 0.91 / 1.2 us  cpu 0.20%
```

`isr_wake` = whole flywheel ISR; `isr_pack` = 72× packPixel/column;
`isr_dma_submit` = submitFrame (dcache flush + eDMA kick, ~1 µs CPU). Notably,
Voronoi's ISR steal drops to ~2.8% at -O3 (from ~7.9% at -Os): the 16 fps frames
are half the wall time and the packPixel avg falls to its ~5 µs floor (the -Os
long 125 ms frames + heavier memory contention were what inflated it).

## Summary ranking (steady 16 fps frame, share of 62.5 ms)

1. `vo_shade` (per-pixel KD-nearest query + shade) — **92%** (57.2 ms) — the
   entire algorithmic cost.
2. `vo_buffer_wait` (display-window sync idle) — 8% (5.0 ms, idle by design).
3. `vo_kdtree` (per-frame KD-tree rebuild) — 0.3% (0.18 ms).
4. `vo_animate` (site spin) — 0.1% (64 µs).

## -O3 vs -Os

| Metric | -Os | -O3 | Δ |
|---|---|---|---|
| `vo_shade` per pixel | 4,452 cyc | 3,308 cyc | **1.35× faster** |
| `vo_shade` per frame | 76.9 ms | 57.2 ms | 1.35× (−19.7 ms) |
| **Cadence** | **8 fps (125 ms)** | **16 fps (62.5 ms)** | **doubled** |
| ISR steal | ~7.9% | ~2.8% | lower (shorter frames) |
| FLASH code | 33,588 B | 54,532 B | +62.4% |
| ITCM code | 21,832 B | 40,264 B | +84.4% |

Voronoi is one of the standout -O3 wins: the 1.35× shade speedup pulls the
per-pixel KD-nearest render below the 62.5 ms window, **doubling the frame rate
from 8 fps to 16 fps**. The KD-tree rebuild stays negligible at both levels —
the cost is entirely the per-pixel nearest-site queries. Cost: +62% FLASH,
+84% ITCM; -O3 does not ship (the full roster overflows FlexRAM), so this gain
is only reachable as a single-effect profile image.

## Caveats

- **ISR time is included** (CYCCNT free-runs) — the real shipping condition.
- No `filter_blend` subtree: the open-coded pixel loop composites in place; the
  per-pixel figures are derived from the fixed 10,368-px quadrant coverage.
- **No timeline**: Voronoi drives its own prep + shade; `vo_animate`/`vo_kdtree`
  are the prep-stage scopes and `vo_shade` is a single bare scope over all
  pixel work.
- **This is an -O3 build (`profile_o3` env), which does NOT ship**: the shipping
  Phantasm image is -Os because the full roster overflows ITCM at -O3.
- Ran with the uncommitted `HS_PROFILE` instrumentation in `effects/Voronoi.h`.

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Voronoi`, `-D HS_PROFILE_WINDOW=<frames>`).
- `pio run -e profile_o3 -t upload` + `tools/profile_capture.py`.
