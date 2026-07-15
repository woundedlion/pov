# DistortedRing on-device profile — Teensy 4.0, segmented mode, -O3 (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_distortedring_o3.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm shipping flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but built at **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` — the -O3 twin of the shipping `-Os` profile; only the optimization level differs |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `DistortedRing<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 79,572 B; ITCM (RAM1 code) 62,344 B; RAM2 free 4,736 B.
(vs -Os: FLASH 49,428 B; ITCM 33,096 B; RAM2 free 4,736 B.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
window frames 65–128 root counter = 2,398,848,738 cyc = 3,998,081.23 µs vs
measured `micros()` window sum 3,998,092 µs (Δ ≈ 2.7 ppm).

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
62.4–62.6 ms wall. What *does* move is the render leaf: at -O3 `dr_ring_scan`
breathes between **5 % and 10 %** of the frame as the ring's displaced band
widens and narrows (down from 6–13 % at -Os — the ring got cheaper, the idle
grew); the rest of every frame is idle spin.

`dr_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly that round-up idle — **90–
94 % of every frame** here, because the one-ring render leaves the window almost
entirely empty (a bit emptier than -Os, since -O3 shrank the render).

## Phase-by-phase readout (64-frame windows, per-frame averages)

The effect has no discrete phases and a single 16 fps cadence; only the ring's
band size oscillates. The two windows below bound that breathing (the same
representative frames chosen at -Os, for a like-for-like comparison). Columns:
time/frame, cycles/frame, % of frame, calls/frame, per-call cost (cycles =
µs × 600). The epoch-start window (frames 1–64, min ≈ 4,702 µs) is skipped;
17 clean steady windows follow it.

### Wide-band frame — render peak (representative: frames 321–384)

```
frame                62.51 ms  37.51 Mcyc  100%
  dr_timeline_step    6.30 ms   3.78 Mcyc   10%
    dr_draw           6.27 ms   3.76 Mcyc   10%
      dr_ring_scan    6.05 ms   3.63 Mcyc   10%  x1     6054 us/ring
        filter_blend  0.20 ms   0.12 Mcyc         x1280   94 cyc/blend
      dr_lut_bake     0.21 ms   0.13 Mcyc    0%  x1      215 us/bake
  dr_buffer_wait     56.22 ms  33.73 Mcyc   90%
  dr_mod_advance      0.30 us    179 cyc     0%
```

Wall: min 56.2 / avg 62.5 / max 69.6 ms. The single `dr_ring_scan`
(`Scan::DistortedRing::draw` — knot-polyline SDF distance → soft stroke →
shade → blend) is the whole render cost; `dr_lut_bake` (the per-column
`shift_lut` sine-sum bake) is a flat ~0.21 ms/frame regardless of band size,
and `dr_mod_advance` (mod-driver + per-wave phase advance) is noise. Even at
its widest the ring covers only ~1,280 of the 10,368 quadrant pixels, so 90 %
of the frame is idle `dr_buffer_wait`. This is the same band width (1,280
covered px) as the -Os peak, so it is the cleanest direct A/B: -O3 renders it in
3.63 Mcyc where -Os took 4.69 Mcyc.

### Narrow-band frame — render trough (representative: frames 961–1024)

```
frame                62.46 ms  37.48 Mcyc  100%
  dr_timeline_step    3.52 ms   2.11 Mcyc    6%
    dr_draw           3.49 ms   2.09 Mcyc    6%
      dr_ring_scan    3.28 ms   1.97 Mcyc    5%  x1     3275 us/ring
        filter_blend  0.10 ms   0.06 Mcyc         x662    94 cyc/blend
      dr_lut_bake     0.22 ms   0.13 Mcyc    0%  x1      215 us/bake
  dr_buffer_wait     58.95 ms  35.37 Mcyc   94%
  dr_mod_advance      0.28 us    175 cyc     0%
```

Wall: min 56.8 / avg 62.5 / max 68.1 ms. Same one-ring pipeline, but the
displaced band has narrowed to ~662 covered pixels, roughly halving
`dr_ring_scan` (1.97 vs 3.63 Mcyc). The `dr_lut_bake` cost is unchanged (the LUT
is baked over all columns every frame irrespective of band width), so at this
trough it is a larger *relative* slice of the render but still <1 % of the
frame. Cadence is identical — the render never had a chance of crossing 62.5 ms
either way.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h) fires once per
blended pixel under `dr_ring_scan`: **1,280 blended px/frame** at the wide-band
peak, **662 px/frame** at the trough — i.e. **6–12 % of the 10,368-px quadrant**,
consistent with a single ring rather than a full field. Cost is a steady
**~94 cyc (0.16 µs) per blend** in both regimes — half the -Os cost of ~193 cyc,
the single biggest -O3 win on this path. Backing out the scan: `dr_ring_scan`
is ~2,900 cyc per *blended* pixel (3.63 Mcyc ÷ 1,280 wide → 2,837; 1.97 Mcyc ÷
662 narrow → 2,971), down from ~3,700 cyc at -Os. The frame-to-frame swing is
still purely how many pixels the band covers, not per-pixel cost. Most of that
~2,900 cyc is the knot-polyline SDF distance evaluation over the tested band,
not the blend — the same distance kernel that dominates DisplacementField.

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative window frames 769–832 (~4.0 s). The pack and
submit run inside the flywheel wake, on the 1-in-8 wakes that render a column.
Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18434/s  0.50 / 1.41 / 42.9 us  cpu 2.59%
  isr_pack        2304/s  4.64 / 5.28 / 41.4 us  cpu 1.21%
  isr_dma_submit  2304/s  0.61 / 0.90 / 2.89 us  cpu 0.20%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick itself is ~0.6 µs of CPU**; the pixel-pack marshal is the real
  per-column cost. Whole per-column CPU ≈ **~5.2 µs of the 434 µs column period
  (~1.2 %)**. The wire transfer runs asynchronously on eDMA and costs no CPU.
- The `min` floor (pack 4.64 µs / 2,782 cyc, submit 0.61 µs) is the stable
  per-call cost; the avg/max inflation (pack up to ~40–90 µs in the noisier
  windows) is serial-dump preemption landing inside the scope, not real pack
  growth. The pack floor is marginally higher than -Os's 4.25 µs — the ISR
  marshal barely benefits from -O3 (code layout, not arithmetic-bound), unlike
  the render kernel.
- Net: the ISR machinery steals **~2.6 % of the chip** = ~1.6 ms of the 62.5 ms
  window, leaving **~60.9 ms of render budget**. The wide-band render (~6 ms)
  uses only ~0.10× of that, so 16 fps is never at risk.

## Summary ranking (share of the 62.5 ms frame)

1. display-window sync (`dr_buffer_wait`) — **90–94 %** (56–59 ms) — idle by
   design; the one-ring render leaves the window nearly empty (even emptier at
   -O3 than -Os, because the render shrank).
2. `Scan::DistortedRing::draw` (`dr_ring_scan`: knot-polyline SDF + shade +
   blend) — **5–10 %** (3.3–6.1 ms) — the entire real render cost, breathing
   with band width.
3. per-column `shift_lut` bake (`dr_lut_bake`) — <1 % (~0.21 ms, flat).
4. everything else (`dr_mod_advance`, timeline/draw self time) — <1 %.

Unlike its sibling DisplacementField (rasterizer-bound with ~36 stacked rings),
DistortedRing is **idle-bound**: a single ring through the shared SDF is far too
cheap to fill a 62.5 ms window, so the frame is dominated by display-sync spin,
not compute. -O3 makes the ring cheaper still — the speedup lands entirely in
the small render term and is completely absorbed by `dr_buffer_wait`; the frame
cadence does not move.

## -O3 vs -Os

Both builds are locked at the same **16 fps / 62.5 ms** cadence — this effect is
idle-bound, so -O3 buys no frames, only a smaller render term (and a larger
idle-spin term to match). The gains are real but hidden below the display-sync
floor:

| Metric | -Os | -O3 | -O3 win |
|---|---|---|---|
| `dr_ring_scan`, wide-band peak (1,280 px, same band) | 4.69 Mcyc / 7.82 ms | 3.63 Mcyc / 6.05 ms | **1.29×**, −1.77 ms |
| scan cost per blended pixel | ~3,700 cyc | ~2,900 cyc | ~1.28× |
| `filter_blend` per-blend | ~193 cyc | ~94 cyc | **2.06×** |
| frame wall / cadence | 62.5 ms / 16 fps | 62.5 ms / 16 fps | no change (idle-bound) |
| FLASH code | 49,428 B | 79,572 B | +30,144 B (+61 %) |
| ITCM (RAM1 code) | 33,096 B | 62,344 B | +29,248 B (+88 %) |

Headline: the render hot spot (`dr_ring_scan`) is **~1.29× faster at -O3**
(−1.77 ms at the peak band), and the per-pixel blend is a full **2× faster** —
but because the ring already fits ~10× inside the 62.5 ms window, none of that
reaches the frame rate; it just widens `dr_buffer_wait` from 87–94 % to 90–94 %.

**-O3 is not the shipping config.** The single-effect profile image fits, but
the full 26-effect Phantasm roster overflows FlexRAM/ITCM at -O3 (the ITCM code
alone nearly doubles here, +29 KB for one effect) — which is exactly why the
ship build is -Os. -O3 is only viable as this single-effect profiling image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: it nests under `dr_ring_scan` (the only scope
  that enters it here), so its subtree % is meaningful this time, but its
  per-blend scope overhead inflates the scan counter by a fraction of a percent.
- Epoch-start window (frames 1–64, min ≈ 4,702 µs first frame) mixes the
  instance's cold construction into the window — skipped when picking
  representatives; per-call averages there stay valid, per-frame sums don't.
- **`-O3` build** (the `profile_o3` env): this does **not** ship. The shipping
  Phantasm image is `-Os` because the full 26-effect roster overflows ITCM at
  -O3; this single-effect image is the only place -O3 fits.
- Run with uncommitted `HS_PROFILE` instrumentation in
  `effects/DistortedRing.h` (`dr_mod_advance`, `dr_buffer_wait`,
  `dr_timeline_step`, `dr_draw`, `dr_lut_bake`, `dr_ring_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=DistortedRing`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel pack,
  DMA submit) each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py`. (The
  `just profile <Effect>` recipe drives the shipping `-Os` `profile` env; there
  is no `just` recipe for the -O3 twin.)
