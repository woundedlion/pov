# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_displacementfield_o3.log` (64-frame windows, single ~80 s pass).
This is the **-O3** twin of the shipping `-Os` report
(`docs/profiles/Os/profile_displacementfield_teensy_2026-07-14.md`); the only
variable between the two is the optimization level.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: same Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but keeps base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forcing `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — the shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `DisplacementField<288, 144>` only (single-entry playlist), current working-tree state (`thickness` default 0.035, no init override) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR counters via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 84,876 B; ITCM (RAM1 code) 65,624 B; RAM2 free 4,736 B.
(At -Os: 57,012 / 39,288 / 4,736 — see the `-O3 vs -Os` section for the cost.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
NOISE-dwell window frames 385–448 root counter = 4,798,131,569 cyc =
7,996,885.9 µs vs measured `micros()` window sum 7,996,891 µs (Δ ≈ 0.6 ppm).

**Capture coverage** — the ~80 s pass (640 frames at the 125 ms cadence below)
spans only the **fade-in ramp (windows 1–128)** and the **NOISE dwell
(windows 129–576)**; the phase schedule's fade-out/BALLS regimes (≈ frame 750+)
were not reached, so there is no steady-BALLS or fade-out window here. The
dwell is the hot phase and is captured cleanly across 7 windows (129–576).

## Frame cadence (context for every number below)

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering one quadrant
(this segment's 72-row band × the window's 144-column half = **10,368 px** —
`filter_blend` measures ~12.1 k blended pixels/frame in the dwell, i.e. the
48-ring stack covers the quadrant with ~1.16× overdraw). A frame that misses a
window boundary snaps to the next, so wall time quantizes to multiples of
62.5 ms:

- **NOISE dwell**: render (`df_timeline_step`) ≈ 78 ms > 62.5 ms →
  2-window cadence, **125 ms/frame (8 fps)**.
- **fade-in ramp**: render climbs from a near-flat opening → mixed 1- and
  2-window cadence (16→8 fps) as the noise band fills in.

`df_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately, so the render numbers below are clean.
At -O3 this idle is large (~46 ms in the dwell) because the faster render
finishes well inside the 125 ms budget — see `-O3 vs -Os`.

## Phase-by-phase readout (64-frame windows, per-frame averages)

### NOISE dwell — the hot phase (representative window: frames 385–448)

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame, calls/frame,
per-call time (cycles = µs × 600):

```
frame                  124.95 ms  74.97 Mcyc  100%
  df_timeline_step      78.52 ms  47.11 Mcyc   63%
    df_draw_rings       78.44 ms  47.06 Mcyc   63%
      df_ring_scan      65.61 ms  39.36 Mcyc   53%  x37.0  1771 us/ring
        filter_blend                                x12077   95 cyc/blend
      df_lut_bake        9.51 ms   5.71 Mcyc    8%  x37.0   257 us/ring
      df_hue_table_prep  1.74 ms   1.04 Mcyc    1%  x22.8    76 us/call
      df_chunk_cull      1.37 ms   0.82 Mcyc    1%  x41.1    33 us/ring
  df_buffer_wait        46.43 ms  27.86 Mcyc   37%
  df_prepare_fields      0.08 us      49 cyc    0%
```

`df_ring_scan` is the SDF rasterize (`Scan::DistortedRing::draw`);
`df_buffer_wait` is the display-window sync idle. **No `df_flat_scan`** —
at -O3 the flat-ring fast path is compiled out (it is `#ifdef __OPTIMIZE_SIZE__`
in `DisplacementField.h`), so undisplaced rings fall through to `df_ring_scan`
with a zero-shift LUT rather than the dedicated flat draw.

Wall: min 123.9 / avg 125.0 / max 126.2 ms. Ring accounting per frame: 48
rings → ~41 survive the cap/chunk cull and get baked + scanned; 0 flat-path.
`df_draw_rings` self time (prefilter, palette, basis) ≈ 0.13 Mcyc ≈ 0.22 ms/frame.
Dwell is stable across all 7 dwell windows (129–576): render 68.9–78.5 ms/frame,
`df_ring_scan` 57.2–65.6 ms (46–53% of frame) — the spread is the animating
noise band, not measurement noise. Note that at -O3 `df_ring_scan` (53%) no
longer dwarfs `df_buffer_wait` (37%) the way it did at -Os (69% vs 16%): the
scan got ~1.3× faster while the 125 ms cadence is fixed, so the reclaimed time
re-emerges as sync idle.

### Fade-in ramp (windows 1–128)

The first two windows are the amplitude ramp from `master_gain≈0`:

- **window 1–64**: wall min 24.0 / avg 61.7 / max 69.0 ms — the opening frames
  are near-flat (all rings displace ~0), so cadence starts at 1-window (16 fps).
  With the -Os flat path gone, these rings still route through `df_ring_scan`
  but at ~712 µs/ring (a cheap zero-shift scan) vs 1771 µs/ring in the dwell,
  and ~42 rings scan/frame (vs ~37 in the dwell — fewer are cap-culled while flat).
- **window 65–128**: wall min 60.4 / avg 88.2 / max 127.6 ms — the noise band
  widens, per-ring scan cost rises and the frame slides to the 2-window
  (8 fps / 125 ms) dwell cadence by the end.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h) parents under
`df_ring_scan` in the dwell (no flat ring renders): ~12,077 blended
pixels/frame at **95.1 cyc (0.16 µs) per blend** (includes its own scope
overhead) — roughly half the -Os per-blend cost. `df_ring_scan` ≈ 39.36 Mcyc
over ~12.1 k blended px ⇒ **~3,259 cyc per blended pixel** (vs ~4,780 at -Os),
i.e. most of the scan cost is still the knot-polyline SDF distance evaluation
over the tested band, not the blend.

## Column-ISR / DMA marshaling cost (ISR-safe counters, same log)

Measured with `HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree
is main-loop-only), representative dwell window 385–448 (8.0 s). Nested as
executed — the pack and submit run inside the flywheel wake, on the 1-in-8
wakes that render a column. Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake         18432/s   0.50 / 2.66 / 139.8 us   cpu 4.3-5.1%
  isr_pack        2304/s   4.65 / 14.5 / 138.1 us   cpu 2.8-3.6%
  isr_dma_submit  2304/s   0.66 / 0.92 / 3.5  us    cpu 0.21%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA call itself is ~0.7–1 µs of CPU**; the marshaling around it (pixel
  pack) is ~4.7–15 µs, so the whole per-column CPU cost is a small fraction of
  the 434 µs column period. The wire transfer runs asynchronously on eDMA and
  costs no CPU — it never blocks rendering, it only bounds back-to-back submits.
- The `min` floors (pack 2789 cyc ≈ 4.65 µs, submit 398 cyc ≈ 0.66 µs) are the
  stable per-call costs; the avg/max inflation (up to ~138 µs) is serial-dump
  preemption landing inside the scope, not real pack work.
- Net: the ISR machinery steals **~4.3–5.1% of the chip** (unchanged from -Os —
  it is display-driver work, not effect code). Against the 62.5 ms window that
  is ~3.1 ms, leaving ~59 ms of render budget per window. The dwell render
  (78.5 ms) needs ~1.32× the window, so it still lands at 2-window / 125 ms
  cadence.

## Summary ranking (NOISE dwell, share of the 125 ms frame)

1. `Scan::DistortedRing::draw` (`df_ring_scan`, SDF rasterize + shade + blend) — **53%** (65.6 ms)
2. display-window sync (`df_buffer_wait`) — 37% (46.4 ms, idle by design)
3. LUT bake (ball/noise field + hue per column) — 8% (9.5 ms)
4. hue-table prep — 1% (1.7 ms)
5. chunk cull — 1% (1.4 ms)
6. everything else (`df_prepare_fields`, timeline, prefilter) — <1%

The rasterizer still dominates the *render* budget on device, but -O3 shrinks
it enough that the fixed display cadence (not the scan) now governs fps — the
21 ms it saves over -Os is entirely absorbed by `df_buffer_wait`.

## -O3 vs -Os

Both builds render the identical scene at the identical **125 ms / 8 fps**
cadence — -O3 does **not** change frame rate here, because even the faster
render (~78 ms) still exceeds one 62.5 ms window, so the saved time only widens
`df_buffer_wait` idle (19.6 → 46.4 ms). -O3 would only lift fps if it dropped
render below 62.5 ms. Per-scope, -O3 is ~1.3–2.4× faster (representative NOISE
dwell, per-frame / per-call):

| metric | -Os | -O3 | speedup |
|---|---|---|---|
| render `df_draw_rings` (ms/frame) | 105.41 | 78.44 | **1.34×** (−27.0 ms) |
| `df_ring_scan` (ms/frame) | 86.79 | 65.61 | **1.32×** (−21.2 ms) |
| `df_ring_scan` per ring (µs) | 2430 | 1771 | 1.37× |
| `df_lut_bake` per ring (µs) | 358 | 257 | 1.39× |
| `df_hue_table_prep` per call (µs) | 182 | 76 | 2.39× |
| `df_chunk_cull` per ring (µs) | 37 | 33 | 1.12× |
| `filter_blend` (cyc/blend) | 195.8 | 95.1 | 2.06× |
| frame wall / cadence | 125 ms / 8 fps | 125 ms / 8 fps | unchanged |
| FLASH code | 57,012 B | 84,876 B | **+27,864 B (+48.9%)** |
| ITCM (RAM1 code) | 39,288 B | 65,624 B | **+26,336 B (+67.0%)** |

The tightest per-pixel/per-column loops (`filter_blend`, `hue_table_prep`)
benefit most from -O3 (~2.1–2.4×); the SDF scan and LUT bake gain ~1.3–1.4×.
One structural difference: the `__OPTIMIZE_SIZE__` flat-ring fast path
(`df_flat_scan`) is compiled out at -O3, so flat rings route through the normal
scan (this shows up in the fade-in, not the noise-saturated dwell).

**-O3 is not the shipping config.** The full 26-effect Phantasm image overflows
FlexRAM/ITCM at -O3 (a single effect already costs +67% ITCM here), which is
exactly why the ship build is `-Os`. -O3 is only viable as this single-effect
profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs
  the flywheel/DMA/USB ISRs that fire inside it. That overhead is the real
  shipping condition (the point of segmented-mode profiling), but these are
  not pure-CPU algorithm costs.
- **`filter_blend` tree artifact**: it parents under whichever scope first
  enters it (here `df_ring_scan`, since no flat ring renders in the dwell), so
  its printed percentage is nonsense and the subtree is hidden in windows where
  its parent scope isn't entered (`log_node` skips zero-count parents). Its
  calls/cycles are correct.
- Per-pixel `filter_blend` scope overhead inflates the scan counters slightly;
  the coarser `df_*` scopes are negligible.
- **No `df_flat_scan` at -O3**: the flat-ring path is `#ifdef __OPTIMIZE_SIZE__`
  only, so it never appears in this build — expected, not a missing scope.
- **Partial phase coverage**: the ~80 s capture reaches only fade-in + NOISE
  dwell; no fade-out/BALLS window was recorded (the schedule's BALLS regime
  starts past the end of the pass). The first window's opening frames
  (`master_gain≈0`) are the global cheapest but are the amplitude ramp, not a
  steady phase.
- **-O3 build**: this is the `profile_o3` env, which does **not** ship. The
  shipping Phantasm image is `-Os` because the full 26-effect roster overflows
  ITCM at -O3.
- Run with the uncommitted working-tree `thickness = 0.035` default and
  uncommitted `HS_PROFILE` instrumentation in `effects/DisplacementField.h`.

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=<EffectClass>`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake,
  pixel pack, DMA submit) each window.
- `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py`
  (the `-Os` twin uses env `profile`, reproducible via `just profile <EffectClass>`).
