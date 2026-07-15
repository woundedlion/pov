# RingSpin on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). This is the **-O3** twin of the
shipping `-Os` report at
`docs/profiles/Os/profile_ringspin_teensy_2026-07-14.md`. Raw capture:
`build/profile_ringspin_o3.log` (64-frame windows, ~64 s single pass, 1,024
frames). The column-ISR/DMA accumulators are dumped every window.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: identical flags (newlib-nano, `USE_DMA_LEDS`, N=4, `tools/phantasm.ld`, `-D HS_PROFILE_ENABLE`) EXCEPT the base is **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forced `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `RingSpin<288, 144>` only (single-entry playlist), current working-tree state (4 rings, `TRAIL_LENGTH = 19`, alpha 0.5, thickness 0.8) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 58,060 B; ITCM (RAM1 code) 42,424 B; RAM2 free 4,736 B.
(The shipping `-Os` image was 39,716 / 26,280.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
steady window frames 961–1024 root `frame` counter = 2,400,246,194 cyc =
4,000,410.3 µs vs measured `micros()` window sum 4,000,413 µs (Δ ≈ −0.67 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × 144-column half ≈ **10,368 px**). Wall time snaps up to a whole
number of 62.5 ms windows.

RingSpin has no preset or spawn schedule — 4 great-circle rings are spawned once
and each random-walks forever. The draw structure is fixed: `deep_tween` over
each ring's 19-sample trail emits a **constant 224 `Scan::Ring::draw`
sub-strokes/frame** (56 per ring), alpha-tapered along the trail. The cost is not
the stroke *count* (constant) but the on-screen *coverage*: it swings with how
much of each ring's great circle currently falls inside this segment's 72-row
colatitude band. So the cadence **jitters** with coverage rather than stepping
through phases — but at **-O3 the whole curve shifts under the 62.5 ms
threshold**, so most of the run is now a single window:

- **Rings drifting to the poles** (low coverage, ~37.5k blends): render dips to
  ~45 ms → comfortably 1-window, **steady 62.5 ms/frame (16 fps)**. This is the
  common state (frames 257–512, 577–832, 961–1024 all run 62.5 ms).
- **Rings in-band** (great circles crossing the band, ~46k blends): render ~55 ms
  — still under 62.5 ms, but sitting right on the threshold, so frame-to-frame
  jitter oscillates between one window (62.5 ms) and two (125 ms); windows
  average ~81 ms.
- **Peak in-band coverage** (~53k blends, frames 833–896): render 64.1 ms just
  *exceeds* the window → falls back to 2-window, **125 ms/frame (8 fps)** — the
  only residual 8 fps state at -O3, and the direct descendant of the -Os steady
  heavy state.

`rs_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor) is the
round-up idle, timed separately.

## Phase-by-phase readout (64-frame windows, per-frame averages)

No discrete phases; the windows below bound the continuous coverage jitter.
Columns: time/frame, cycles/frame, % of frame, calls/frame, per-call cost.

### Rings in-band — borderline 16↔8 fps (representative: frames 513–576)

```
frame                 80.93 ms  48.56 Mcyc  100%
  rs_draw_rings       54.80 ms  32.88 Mcyc   68%
    rs_ring_scan      54.59 ms  32.74 Mcyc   67%  x224  244 us/scan
      filter_blend     6.64 ms   3.99 Mcyc    8%  x45801  87 cyc/blend
  rs_buffer_wait      26.06 ms  15.65 Mcyc   32%
  rs_timeline_step      54 us     32 kcyc     0%
```

Wall: min 51.7 / avg 80.9 / max 131.4 ms. `rs_ring_scan` (the `Scan::Ring::draw`
SDF stroke + AA blend) is still essentially the whole render at 54.6 ms — but
that is now *below* the 62.5 ms window, so the per-frame cadence flickers between
one and two windows (hence the 51.7→131.4 ms spread). `rs_timeline_step` (the
four random-walk steps) is 54 µs. `filter_blend` at 45,801 blends/frame is
**4.42× the 10,368-px quadrant** — the trail massively overdraws itself (56
overlapping strokes per ring stack into the same band).

At **peak in-band coverage** (frames 833–896: 53,405 blends/frame, 5.15×
overdraw) `rs_ring_scan` climbs to 63.9 ms — the one window that still exceeds
62.5 ms and holds the 8 fps 2-window cadence. That window's inflated ISR share
(pack avg 24.7 µs, cpu 5.7 %) is serial-dump preemption of its long frames, not
extra work — its `min` floors match every other window.

### Rings drifting out — steady 16 fps (representative: frames 961–1024)

```
frame                 62.51 ms  37.50 Mcyc  100%
  rs_draw_rings       44.96 ms  26.98 Mcyc   72%
    rs_ring_scan      44.74 ms  26.84 Mcyc   72%  x224  200 us/scan
      filter_blend     5.23 ms   3.14 Mcyc    8%  x37556  84 cyc/blend
  rs_buffer_wait      17.49 ms  10.49 Mcyc   28%
  rs_timeline_step      52 us     31 kcyc     0%
```

Wall: min 56.3 / avg 62.5 / max 69.2 ms. Same 224 strokes, but each covers fewer
band rows (rings near the poles), so per-scan cost drops 244→200 µs and blended
pixels drop 46k→38k. The render (45.0 ms) now sits well under the 62.5 ms window,
so `rs_buffer_wait` (17.5 ms) fills the remainder and every frame holds a single
window — a rock-steady 16 fps (min/max span only ±7 ms). Steady windows 257–320,
321–384, 385–448, 449–512, 577–640, 641–704, 705–768 and 769–832 all reproduce
this 62.5 ms cadence.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in `filter.h`): **45,801
blends/frame** in the in-band window — 4.42× the 10,368-px quadrant from trail
overdraw — at **87 cyc (0.145 µs) per blend** (down from 183 cyc at -Os).
`rs_ring_scan` ≈ 32.74 Mcyc over those 45,801 blended px ⇒ **~715 cyc per blended
pixel** (was ~975 at -Os): most of the ring-scan cost is still the SDF stroke
distance/coverage math, not the blend, but -O3's `-ffast-math` inlining roughly
halves the blend leaf.

## Column-ISR / DMA marshaling cost (`build/profile_ringspin_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative steady window (frames
961–1024, 4.00 s). Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.52 / 1.62 / 201 us  cpu 2.98%
  isr_pack        2304/s  4.66 / 6.22 / 200 us  cpu 1.43%
  isr_dma_submit  2304/s  0.58 / 0.91 / 6.0 us  cpu 0.20%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick is ~0.9 µs of CPU.** The `min` floor (pack 4.66 µs / 2,793 cyc,
  submit 0.58 µs) is the stable per-call cost; the avg/max inflate only in the
  heavy in-band windows (pack avg climbs to 24.7 µs in frames 833–896) from
  serial-dump preemption of the long frames, not extra marshaling work.
- Net: in the steady 16 fps windows the ISR machinery steals **~3 % of the chip**
  (~1.9 ms of the 62.5 ms window), leaving **~60.6 ms render budget**. The
  drifting render (45 ms) clears that with room to spare; the in-band render
  (55 ms) clears it too but with only ~5 ms of margin, which is why per-frame
  jitter tips it over into a second window; the peak in-band render (64 ms)
  exceeds it outright.

## Summary ranking (rings in-band, share of the 80.9 ms frame)

1. `Scan::Ring::draw` (`rs_ring_scan`: SDF stroke + AA blend, 224 sub-strokes) —
   **67 %** (54.6 ms) — the whole render.
2. display-window sync (`rs_buffer_wait`) — 32 % (26.1 ms, idle by design).
3. four random-walk steps (`rs_timeline_step`) — <0.1 % (54 µs).

RingSpin is rasterizer-bound with a self-inflicted **~4.5× overdraw** tax: the
motion-blur trail redraws the band ~4–5 times per frame. -O3 does not change that
structure — it makes each stroke/blend cheaper. The highest-leverage target is
still the trail's `deep_tween` sub-stroke count (224/frame) and the alpha cull
threshold — fewer sub-strokes or an earlier band-cull would cut overdraw
directly. No WASM/native RingSpin figures are recorded in the perf ledger for
comparison.

## -O3 vs -Os

Optimization level is the **only** variable between this run and the `-Os` report
(same env otherwise). -O3 uniformly speeds the rasterizer and lifts the cadence:

| metric | -Os | -O3 | Δ |
|---|---|---|---|
| `rs_ring_scan`, in-band | 76.85 ms | 54.59 ms | **1.41× faster** (−22.3 ms) |
| `rs_ring_scan` per sub-stroke, in-band | 343 µs | 244 µs | 1.41× |
| `rs_ring_scan`, drifting | 59.32 ms | 44.74 ms | 1.33× (−14.6 ms) |
| scan cyc / blended px, in-band | ~975 | ~715 | 1.36× |
| `filter_blend` per blend | 183 cyc | 87 cyc | **2.10× faster** (−96 cyc) |
| in-band cadence | 125 ms (8 fps), steady | 62.5–125 ms, borderline 16↔8 | ↑ toward 16 fps |
| drifting cadence | ~83 ms, jittery | 62.5 ms (16 fps), steady | ↑ to solid 16 fps |
| FLASH code | 39,716 B | 58,060 B | +18,344 B (1.46×) |
| ITCM (RAM1) | 26,280 B | 42,424 B | +16,144 B (1.61×) |
| RAM2 free | 4,736 B | 4,736 B | — |

Headline: **rs_ring_scan is ~1.4× faster** at -O3 (blend leaf alone ~2.1×). More
importantly the in-band render drops from 76.9 ms — inherently 2-window at -Os
(steady 8 fps) — to 54.6 ms, *below* the 62.5 ms window, so **the cadence lifts
from a solid 8 fps toward 16 fps**: the drifting/pole state is now rock-steady
16 fps, the in-band state hovers on the threshold (jittering 16↔8), and only the
peak-coverage window (64.1 ms) still holds 8 fps.

This is **not** the shipping config. -O3 costs +18 KB FLASH and +16 KB ITCM for
this single-effect image; the full 26-effect Phantasm roster overflows FlexRAM /
ITCM at -O3, which is exactly why the ship build is `-Os`. -O3 is viable here
only as a single-effect profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, and the
  reason the in-band `isr_pack` avg inflates in the long-frame windows (see
  above).
- **`filter_blend` parenting**: here it parents cleanly under `rs_ring_scan`
  (which always enters it), so its subtree and counts are valid; its printed % is
  of the frame.
- Per-pixel `filter_blend` scope overhead inflates `rs_ring_scan` by a fraction
  of a percent (~46k calls/frame).
- Epoch/reset windows (a `frame` count with a near-zero first frame, e.g.
  frames 1–64 with min=1,491 µs) mix the reconstructed instance's first frame —
  they are skipped; the clean steady windows (≥8 of them) are used as
  representatives.
- **-O3 build**: this is the `profile_o3` env, which does **not** ship — the
  shipping Phantasm image is `-Os` because the full 26-effect roster overflows
  ITCM at -O3.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/RingSpin.h`
  (`rs_buffer_wait`, `rs_timeline_step`, `rs_draw_rings`, `rs_ring_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=RingSpin`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (the -O3
  env; the shipping `-Os` twin is `just profile RingSpin`).
</content>
</invoke>
