# DistortedRing on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_distortedring.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `DistortedRing<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile DistortedRing` |

Image size: FLASH code 49,428 B; ITCM (RAM1 code) 33,096 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
window frames 65–128 root counter = 2,398,563,786 cyc = 3,997,606.31 µs vs
measured `micros()` window sum 3,997,614 µs (Δ ≈ 1.9 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). DistortedRing
draws essentially **one** distorted ring per frame (`dr_ring_scan` fires exactly
1×/frame) through the same `Scan::DistortedRing` knot-polyline SDF rasterizer
that DisplacementField stacks ~36 deep — so this effect is its lightweight
sibling: a single ring covering only a slice of the quadrant. Wall time snaps up
to a whole number of 62.5 ms windows.

Render never approaches one window, so the cadence is **locked at 16 fps
(62.5 ms/frame)** for the entire pass — every non-epoch window reads avg
62.4–62.7 ms wall. What *does* move is the render leaf: `dr_ring_scan` breathes
between **6 % and 13 %** of the frame as the ring's displaced band widens and
narrows; the rest of every frame is idle spin.

`dr_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly that round-up idle — **86–
93 % of every frame** here, because the one-ring render leaves the window almost
entirely empty.

## Phase-by-phase readout (64-frame windows, per-frame averages)

The effect has no discrete phases and a single 16 fps cadence; only the ring's
band size oscillates. The two windows below bound that breathing. Columns:
time/frame, cycles/frame, % of frame, calls/frame, per-call cost (cycles =
µs × 600). The epoch-start window (frames 1–64, min ≈ 6,159 µs) is skipped.

### Wide-band frame — render peak (representative: frames 321–384)

```
frame                62.52 ms  37.51 Mcyc  100%
  dr_timeline_step    8.09 ms   4.85 Mcyc   13%
    dr_draw           8.06 ms   4.83 Mcyc   13%
      dr_ring_scan    7.82 ms   4.69 Mcyc   13%  x1     7823 us/ring
        filter_blend  0.41 ms   0.25 Mcyc         x1280  194 cyc/blend
      dr_lut_bake     0.23 ms   0.14 Mcyc    0%  x1      232 us/bake
  dr_buffer_wait     54.43 ms  32.66 Mcyc   87%
  dr_mod_advance      0.39 us    237 cyc     0%
```

Wall: min 54.0 / avg 62.5 / max 71.9 ms. The single `dr_ring_scan`
(`Scan::DistortedRing::draw` — knot-polyline SDF distance → soft stroke →
shade → blend) is the whole render cost; `dr_lut_bake` (the per-column
`shift_lut` sine-sum bake) is a flat ~0.23 ms/frame regardless of band size,
and `dr_mod_advance` (mod-driver + per-wave phase advance) is noise. Even at
its widest the ring covers only ~1,280 of the 10,368 quadrant pixels, so 87 %
of the frame is idle `dr_buffer_wait`.

### Narrow-band frame — render trough (representative: frames 961–1024)

```
frame                62.47 ms  37.48 Mcyc  100%
  dr_timeline_step    3.96 ms   2.38 Mcyc    6%
    dr_draw           3.94 ms   2.36 Mcyc    6%
      dr_ring_scan    3.70 ms   2.22 Mcyc    6%  x1     3700 us/ring
        filter_blend  0.19 ms   0.11 Mcyc         x586   193 cyc/blend
      dr_lut_bake     0.24 ms   0.14 Mcyc    0%  x1      235 us/bake
  dr_buffer_wait     58.51 ms  35.10 Mcyc   94%
  dr_mod_advance      0.39 us    233 cyc     0%
```

Wall: min 59.0 / avg 62.5 / max 66.1 ms. Same one-ring pipeline, but the
displaced band has narrowed to ~586 covered pixels, halving `dr_ring_scan`
(2.22 vs 4.69 Mcyc). The `dr_lut_bake` cost is unchanged (the LUT is baked over
all columns every frame irrespective of band width), so at this trough it is a
larger *relative* slice of the render but still <1 % of the frame. Cadence is
identical — the render never had a chance of crossing 62.5 ms either way.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h) fires once per
blended pixel under `dr_ring_scan`: **1,280 blended px/frame** at the wide-band
peak, **586 px/frame** at the trough — i.e. **6–12 % of the 10,368-px quadrant**,
consistent with a single ring rather than a full field. Cost is a steady
**~193–194 cyc (0.32 µs) per blend** in both regimes. Backing out the scan:
`dr_ring_scan` is ~3,700 cyc per *blended* pixel (4.69 Mcyc ÷ 1,280 wide;
2.22 Mcyc ÷ 586 narrow — both land near 3,700), so the frame-to-frame swing is
purely how many pixels the band covers, not per-pixel cost. Most of that 3,700
cyc is the knot-polyline SDF distance evaluation over the tested band, not the
blend — the same distance kernel that dominates DisplacementField.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative window frames 769–832 (~4.0 s). The pack and
submit run inside the flywheel wake, on the 1-in-8 wakes that render a column.
Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18434/s  0.58 / 1.54 / 43.2 us  cpu 2.83%
  isr_pack        2304/s  4.25 / 4.77 / 41.4 us  cpu 1.09%
  isr_dma_submit  2304/s  0.67 / 0.97 / 1.07 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal is the real
  per-column cost. Whole per-column CPU ≈ **~5 µs of the 434 µs column period
  (~1.2 %)**. The wire transfer runs asynchronously on eDMA and costs no CPU.
- The `min` floor (pack 4.25 µs / 2,551 cyc, submit 0.6 µs) is the stable
  per-call cost; the avg/max inflation (pack up to ~40–90 µs in the noisier
  windows) is serial-dump preemption landing inside the scope, not real pack
  growth.
- Net: the ISR machinery steals **~2.8 % of the chip** = ~1.8 ms of the 62.5 ms
  window, leaving **~60.7 ms of render budget**. The wide-band render (~8 ms)
  uses only ~0.13× of that, so 16 fps is never at risk.

## Summary ranking (share of the 62.5 ms frame)

1. display-window sync (`dr_buffer_wait`) — **87–94 %** (54–59 ms) — idle by
   design; the one-ring render leaves the window nearly empty.
2. `Scan::DistortedRing::draw` (`dr_ring_scan`: knot-polyline SDF + shade +
   blend) — **6–13 %** (3.7–7.8 ms) — the entire real render cost, breathing
   with band width.
3. per-column `shift_lut` bake (`dr_lut_bake`) — <1 % (~0.23 ms, flat).
4. everything else (`dr_mod_advance`, timeline/draw self time) — <1 %.

Unlike its sibling DisplacementField (rasterizer-bound at ~105 ms/quadrant with
~36 stacked rings), DistortedRing is **idle-bound**: a single ring through the
shared SDF is far too cheap to fill a 62.5 ms window, so the frame is dominated
by display-sync spin, not compute.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: it nests under `dr_ring_scan` (the only scope
  that enters it here), so its subtree % is meaningful this time, but its
  per-blend scope overhead inflates the scan counter by a fraction of a percent.
- Epoch-start window (frames 1–64, min ≈ 6,159 µs first frame) mixes the
  instance's cold construction into the window — skipped when picking
  representatives; per-call averages there stay valid, per-frame sums don't.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with uncommitted `HS_PROFILE` instrumentation in
  `effects/DistortedRing.h` (`dr_mod_advance`, `dr_buffer_wait`,
  `dr_timeline_step`, `dr_draw`, `dr_lut_bake`, `dr_ring_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=DistortedRing`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel pack,
  DMA submit) each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile DistortedRing [seconds]`.
