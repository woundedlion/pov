# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_voronoi.log` (64-frame windows, ~75 s single pass; ISR/DMA
accumulators dumped every window).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `Voronoi<288, 144>` only (single-entry playlist) |
| Method | DWT CYCCNT cycle counters (HS_PROFILE RAII scopes, 600 cyc/µs) + micros() wall clock per draw_frame, dumped every 64 frames |
| Reproduce | `just profile Voronoi` |

Image size: FLASH code 33,588 B; ITCM (RAM1 code) 21,832 B; RAM2 free 4,736 B.

**Exactness cross-check** — cycle counter vs wall clock agree to sub-ppm:
window frames 65–128 root counter = 4,803,048,776 cyc = 8,005,081.3 µs vs
measured `micros()` window sum 8,005,087 µs (Δ ≈ 0.7 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × 144-col half ≈ **10,368 px**). Voronoi holds steady at
**125 ms/frame (8 fps, 2-window)** — its per-pixel nearest-site shade overflows
one window. `vo_buffer_wait` (the `buffer_free()` spin in the `Canvas` ctor) is
the round-up idle.

Structural note: Voronoi is unusual — it has **no timeline** and does **not**
use `Scan::Shader`. It open-codes its own per-frame prep (site spin + KD-tree
rebuild) and a per-pixel KD-nearest shading loop, so there is no
`vo_timeline_step`; `vo_animate`/`vo_kdtree` stand in for the prep stage.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Steady across the whole capture — one representative window (frames 65–128).
Columns: time/frame, cycles/frame, % of frame, per-pixel cost on leaves.

```
frame                125.08 ms  75.05 Mcyc  100%
  vo_shade            76.93 ms  46.16 Mcyc   61%  10368 px  4452 cyc/px
  vo_buffer_wait      47.56 ms  28.54 Mcyc   38%
  vo_kdtree            0.43 ms  257 kcyc      0%
  vo_animate           0.16 ms   96 kcyc      0%
```

Wall: min 121.5 / avg 125.1 / max 130.2 ms. `vo_shade` (a coarse-grid corner
pre-pass + the per-pixel KD-nearest-site query and shade) is the entire
algorithmic cost — 76.9 ms shading the 10,368-px quadrant at **4,452 cyc
(7.42 µs) per pixel**. The two prep stages are effectively free: the whole
per-frame **KD-tree rebuild is 0.43 ms (0.3%)** and the site spin/renormalize
`vo_animate` is 0.16 ms — i.e. the nearest-site *queries* in the shade loop, not
tree construction, dominate. `vo_buffer_wait` is the 38% round-up idle to the
second display window.

### Per-pixel figures

Voronoi's open-coded shade loop composites in place (no `filter_blend`
counter), so per-pixel cost is derived from full quadrant coverage: it shades
all **10,368 px/frame** at **4,452 cyc (7.42 µs) per pixel** — the same order
as Flyby's expensive regime, dominated by the KD-tree nearest-neighbour descent
per pixel.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` accumulators, window frames 65–128 (8.0 s). Columns: rate,
per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.62 / 4.28 / 156 us  cpu 7.87%
  isr_pack        2304/s  4.25 / 26.78 / 154 us cpu 6.17%
  isr_dma_submit  2304/s  0.64 / 0.98 / 1.3 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s COLUMN_US/8 grid); `isr_pack` = the
72× `packPixel` marshal per column; `isr_dma_submit` = `submitFrame`
(dcache flush + eDMA kick, ~1 µs CPU; wire transfer async on eDMA).
Voronoi shows one of the **highest ISR steals** (~7.9%, pack avg ~27 µs) —
like RingSpin, its memory-bound per-pixel shade contends the OCRAM display
buffer, so the ISR's strided `packPixel` reads miss more; the `min` floor
(pack 4.25 µs) matches every other effect, so the elevated avg is
contention/serial-dump preemption in the long 125 ms frames. Against the
62.5 ms window that ~7.9% is ~4.9 ms, leaving ~57.6 ms render budget; the
77 ms shade needs **~1.34×** to reach 1-window cadence.

## Summary ranking (steady frame, share of 125 ms)

1. `vo_shade` (per-pixel KD-nearest query + shade) — **61%** (76.9 ms) — the
   entire algorithmic cost.
2. `vo_buffer_wait` (display-window sync idle) — 38% (47.6 ms, idle by design).
3. `vo_kdtree` (per-frame KD-tree rebuild) — 0.3% (0.43 ms).
4. `vo_animate` (site spin/renormalize) — 0.1% (0.16 ms).

Voronoi is pure per-pixel-query-bound: the KD-tree rebuild it does every frame
is negligible next to the ~10 k nearest-site lookups. The only lever is the
per-pixel nearest-site cost (site count / KD descent) or the corner pre-pass
coverage.

## Caveats

- **ISR time is included**: CYCCNT free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs — the real shipping condition, and the reason
  `isr_pack` avg is inflated here (memory contention, see above).
- No `filter_blend` subtree: the open-coded pixel loop composites in place; the
  per-pixel figures are derived from the fixed 10,368-px quadrant coverage.
- **No timeline**: Voronoi drives its own prep + shade, so there is no
  `vo_timeline_step`; `vo_animate`/`vo_kdtree` are the prep-stage scopes and
  `vo_shade` is a single bare scope covering all pixel work.
- Epoch-straddle windows (low-min first frame) are skipped for the
  representative.
- `-Os` shipping config; an `-O3` profile would differ.
- Ran with the uncommitted `HS_PROFILE` instrumentation in `effects/Voronoi.h`
  (`vo_buffer_wait`, `vo_animate`, `vo_kdtree`, `vo_shade`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Voronoi`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile Voronoi [seconds]`.
