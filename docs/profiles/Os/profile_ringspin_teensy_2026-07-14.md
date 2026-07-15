# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw captures:
`build/profile_ringspin.log` (128-frame windows, ~1.6 epochs),
`build/profile_ringspin_w32.log` (32-frame windows, jitter resolution). The
column-ISR/DMA accumulators are dumped every window in both.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `RingSpin<288, 144>` only (single-entry playlist), current working-tree state (4 rings, `TRAIL_LENGTH = 19`, alpha 0.5, thickness 0.8) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 128 (resp. 32) frames then reset |
| Reproduce | `just profile RingSpin` |

Image size: FLASH code 39,716 B; ITCM (RAM1 code) 26,280 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
coarse window frames 257–384 root counter = 9,599,078,168 cyc = 15,998,463.6 µs
vs measured `micros()` window sum 15,998,483 µs (Δ ≈ 1.2 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × 144-column half ≈ **10,368 px**). Wall time snaps up to a whole
number of 62.5 ms windows.

RingSpin has no preset or spawn schedule — 4 great-circle rings are spawned once
and each random-walks forever. The draw structure is fixed: `deep_tween` over
each ring's 19-sample trail emits a **constant ~224 `Scan::Ring::draw`
sub-strokes/frame** (56 per ring), alpha-tapered along the trail. The cost is
not the stroke *count* (constant) but the on-screen *coverage*: it swings with
how much of each ring's great circle currently falls inside this segment's
72-row colatitude band. So the cadence **jitters** rather than stepping through
phases:

- **Rings in-band** (great circles crossing the segment band): render ~77 ms →
  2-window cadence, **125 ms/frame (8 fps)** — the steady heavy state.
- **Rings drifting to the poles**: render dips to ~59 ms, so individual frames
  fall to 1-window (62.5 ms); windows average ~83 ms with heavy frame-to-frame
  jitter (per-frame min 40–47 ms, max 131 ms).

`rs_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor) is the
round-up idle, timed separately.

## Phase-by-phase readout (32-frame windows, per-frame averages)

No discrete phases; the two windows below bound the continuous jitter. Columns:
time/frame, cycles/frame, % of frame, calls/frame, per-call cost. Epoch (960
revolutions = 120 s) re-seeds the random walks.

### Rings in-band — steady 8 fps (representative: w32 frames 97–128)

```
frame                125.28 ms  75.13 Mcyc  100%
  rs_draw_rings       77.26 ms  46.36 Mcyc   62%
    rs_ring_scan      76.85 ms  46.11 Mcyc   61%  x224   343 us/scan
      filter_blend    14.41 ms   8.65 Mcyc   12%  x47270 183 cyc/blend
  rs_buffer_wait      47.94 ms  28.77 Mcyc   38%
  rs_timeline_step     72 us      43 kcyc     0%
```

Wall: min 123.8 / avg 125.3 / max 126.9 ms. `rs_ring_scan` (the
`Scan::Ring::draw` SDF stroke + AA blend) is essentially the whole render;
`rs_timeline_step` (the four random-walk steps) is 72 µs. `filter_blend` is the
blended-pixel leaf — 47,270 blends/frame is **4.56× the 10,368-px quadrant**, so
the trail massively overdraws itself (56 overlapping strokes per ring stack into
the same band). Only ~12 % of the frame is the blend itself; the other ~50 % of
`rs_ring_scan` is the per-stroke ring SDF + coverage evaluation.

### Rings drifting out — jittery 1–2 window (representative: w32 frames 161–192)

```
frame                 83.28 ms  49.97 Mcyc  100%
  rs_draw_rings       59.72 ms  35.83 Mcyc   72%
    rs_ring_scan      59.32 ms  35.59 Mcyc   71%  x224   265 us/scan
      filter_blend    11.80 ms   7.08 Mcyc   14%  x40104 177 cyc/blend
  rs_buffer_wait      23.50 ms  14.10 Mcyc   28%
  rs_timeline_step     72 us      43 kcyc     0%
```

Wall: min 47.7 / avg 83.3 / max 130.4 ms. Same 224 strokes, but each covers
fewer band rows (rings near the poles), so per-scan cost drops 343→265 µs and
blended pixels drop 47k→40k. The render (59.7 ms) now sits just under the
62.5 ms window, so the cadence oscillates frame-to-frame between one window
(fast frames) and two (slow frames) — the source of the wide min/max spread.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h): **47,270
blends/frame** in the heavy window — 4.56× the 10,368-px quadrant from trail
overdraw — at **183 cyc (0.31 µs) per blend**. `rs_ring_scan` ≈ 46.11 Mcyc over
those 47,270 blended px ⇒ **~975 cyc per blended pixel**: most of the ring-scan
cost is the SDF stroke distance/coverage math, not the blend.

## Column-ISR / DMA marshaling cost (`build/profile_ringspin.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative steady window (coarse
frames 257–384, 16.0 s). Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.60 / 3.40 / 204 us  cpu 6.25%
  isr_pack        2304/s  4.25 / 19.67 / 202 us cpu 4.53%
  isr_dma_submit  2304/s  0.15 / 0.96 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick is ~1 µs of CPU.** RingSpin shows the **highest pack avg** of
  the four effects (19.7 µs vs ~5–14 µs) — its long, memory-bound ring-scan
  frames leave the OCRAM display buffer under heavy cache contention, so the
  ISR's strided `packPixel` reads miss more, and the serial dumps preempt the
  125 ms frames. The `min` floor (pack 4.25 µs) matches the other effects, so
  the elevated avg is contention/preemption, not extra work.
- Net: the ISR machinery steals **~6 % of the chip** here (~3.9 ms of the
  62.5 ms window), leaving **~58.6 ms render budget per window**. The in-band
  render (77 ms) is inherently 2-window; it needs **~1.31×** to hold a single
  window.

## Summary ranking (rings in-band, share of the 125 ms frame)

1. `Scan::Ring::draw` (`rs_ring_scan`: SDF stroke + AA blend, 224 sub-strokes) —
   **61%** (76.9 ms) — the whole render.
2. display-window sync (`rs_buffer_wait`) — 38% (47.9 ms, idle by design).
3. four random-walk steps (`rs_timeline_step`) — <0.1% (72 µs).

RingSpin is rasterizer-bound with a self-inflicted **4.5× overdraw** tax: the
motion-blur trail redraws the band ~4–5 times per frame. The highest-leverage
target is the trail's `deep_tween` sub-stroke count (224/frame) and the alpha
cull threshold — fewer sub-strokes or an earlier band-cull would cut overdraw
directly. No WASM/native RingSpin figures are recorded in the perf ledger for
comparison.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, and the
  reason RingSpin's `isr_pack` avg is inflated (see above).
- **`filter_blend` parenting**: here it parents cleanly under `rs_ring_scan`
  (which always enters it), so its subtree and counts are valid; its printed %
  is of the frame.
- Per-pixel `filter_blend` scope overhead inflates `rs_ring_scan` by a fraction
  of a percent (47k calls/frame).
- Epoch-boundary / reset windows (coarse `frame` count with a near-zero first
  frame, e.g. min=1987 µs) mix the reconstructed instance's first frame — the
  clean steady windows are used as representatives.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/RingSpin.h` (`rs_buffer_wait`, `rs_timeline_step`, `rs_draw_rings`,
  `rs_ring_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=RingSpin`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile RingSpin [seconds]`.
