# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw captures:
`build/profile_islamicstars.log` (128-frame windows, ~2 epochs),
`build/profile_islamicstars_w32.log` (32-frame windows, phase resolution). The
column-ISR/DMA accumulators are dumped every window in both.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `IslamicStars<288, 144>` only (single-entry playlist), current working-tree state (`TerminatorSweep` segue, ripple amp 0.15, 4-ripple burst) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 128 (resp. 32) frames then reset |
| Reproduce | `just profile IslamicStars` |

Image size: FLASH code 69,156 B; ITCM (RAM1 code) 40,360 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
w32 window frames 449–480 root counter = 2,398,954,818 cyc = 3,998,258.0 µs vs
measured `micros()` window sum 3,998,260 µs (Δ ≈ 0.5 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × 144-column half ≈ **10,368 px**). Wall time snaps up to a whole
number of 62.5 ms windows.

The whole render is `Scan::Mesh::draw` (`is_mesh_scan`) — the exact per-face SDF
rasterize (deliberately *not* the congruence-class LUT, which ripple/segue
deformation would mis-shade). Cost tracks two things: the **shape's on-screen
face pixels** (a new solid every ~304 frames, F = 362…542 faces of 2,016…4,320
face indices) and the **ripple deformation** mid-life. So cadence steps per
shape and per choreography stage:

- **Fade-in / still hold, small shape** (cube, F=434): render ~51 ms →
  1-window, **62.5 ms/frame (16 fps)**.
- **Full hold / ripple, mid shapes** (truncatedIcosahedron, F=452): render
  ~68–72 ms → 2-window, **125 ms/frame (8 fps)**.
- **Heavy shapes** (dodecahedron-ambo F=362 / truncatedIcosahedron-ambo F=542,
  big faces): mesh-scan ~79 ms → sustained 8 fps.
- **Segue crossings**: the per-face terminator sweep drops most faces, so those
  frames are the **cheapest** (per-frame min 17–20 ms).

`is_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor) is the
round-up idle, timed separately. `is_ripple_prepare` (the per-frame ripple pool
advance) is ~27 ns — free.

## Phase-by-phase readout (32-frame windows, per-frame averages)

Per-shape choreography (`spawn_shape`): fade-in (72 f) → still (16 f) → ripple
burst (4 ripples, ~128 f) → still (16 f) → fade-out (72 f), duration 304 f. The
`is_mesh_scan` subtree exposes pre-existing render scopes `scan_mesh_raster`
(the per-face rasterize) and `scan_face_setup` (per-face SDF/basis setup).
Columns: time/frame, cycles/frame, % of frame, calls/frame, per-call cost.
Epoch (960 revolutions = 120 s) reconstructs the effect.

### Small shape, still hold — 16 fps (representative: w32 frames 33–64, cube F=434)

```
frame                  62.76 ms  37.66 Mcyc  100%
  is_timeline_step     51.74 ms  31.04 Mcyc   82%
    is_draw_shape      51.70 ms  31.02 Mcyc   82%
      is_mesh_scan     51.31 ms  30.79 Mcyc   82%
        scan_mesh_raster 45.80 ms 27.48 Mcyc  73%  x434   105 us/face
          filter_blend  2.32 ms   1.39 Mcyc    4%  x7790  179 cyc/blend
        scan_face_setup  5.42 ms   3.25 Mcyc    9%  x434   12 us/face
      is_face_offsets    0.25 ms  151 kcyc     0%  x1
      is_mesh_transform  0.13 ms   78 kcyc     0%  x1
  is_buffer_wait       11.02 ms   6.61 Mcyc   18%
  is_ripple_prepare      27 ns     27 cyc      0%
```

Wall: min 59.3 / avg 62.8 / max 67.0 ms. `scan_mesh_raster` (per-face SDF
distance + shade) is 73 % of the frame; `scan_face_setup` adds 9 %; the actual
`filter_blend` is only 4 %. `is_mesh_transform` (orient + ripple) is 0.13 ms
here — no ripple active in the hold. All 434 faces are rasterized every frame
(x434 = the cube's full face count; its 2,016 face-vertex indices are the
per-face SDF workload, not extra faces).

### Mid shape ripple burst — 8 fps (representative: w32 frames 449–480, truncatedIcosahedron F=452)

```
frame                 124.95 ms  74.97 Mcyc  100%
  is_timeline_step     71.95 ms  43.17 Mcyc   58%
    is_draw_shape      71.90 ms  43.14 Mcyc   58%
      is_mesh_scan     68.18 ms  40.91 Mcyc   55%
        scan_mesh_raster 60.30 ms 36.18 Mcyc  48%  x452   133 us/face
          filter_blend  5.22 ms   3.13 Mcyc    4%  x16336 192 cyc/blend
        scan_face_setup  7.79 ms   4.67 Mcyc    6%  x452   17 us/face
      is_mesh_transform  3.39 ms   2.03 Mcyc    3%  x1  (ripple deform)
      is_face_offsets    0.34 ms  202 kcyc     0%  x1
  is_buffer_wait       52.99 ms  31.80 Mcyc   42%
  is_ripple_prepare      27 ns     27 cyc      0%
```

Wall: min 123.9 / avg 124.9 / max 126.1 ms. During the ripple burst,
`is_mesh_transform` (the `RippleTransformer` displacing every vertex + orient)
rises to 3.4 ms (2.7 %), and the rippled mesh presents steeper/larger faces so
`scan_mesh_raster` per-face cost climbs 105→133 µs and blended pixels rise
7.8k→16.3k. The render (72 ms) crosses the window boundary into 8 fps. The
heavier dodecahedron-ambo (F=362, big faces) sustains this 8 fps for its whole
hold with mesh-scan ~79 ms/frame (coarse frames 641–768).

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h): **7,790
blends/frame** in the small-shape hold, **16,336** in the mid-shape ripple —
0.75× and 1.58× the 10,368-px quadrant (the mesh fills a fraction of the band
in the hold, more when rippled) — at **179–192 cyc (0.30–0.32 µs) per blend**.
`scan_mesh_raster` ≈ 27–36 Mcyc over those blends ⇒ **~2,200–3,530 cyc per
blended pixel**: the cost is overwhelmingly the exact per-face SDF distance
evaluation, not the blend (only ~4 % of frame).

## Column-ISR / DMA marshaling cost (`build/profile_islamicstars_w32.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative 16 fps window (w32
frames 33–64, 2.0 s). Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.58 / 1.63 / 148 us  cpu 3.00%
  isr_pack        2304/s  4.25 / 5.45 / 146 us  cpu 1.25%
  isr_dma_submit  2304/s  0.11 / 0.97 / 1.2 us  cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick is ~1 µs of CPU.** IslamicStars has the **lowest base ISR
  cost** of the four (3.00 % vs Flyby 3.3 %, DreamBalls 2.8 %, RingSpin 6.3 %) —
  clean mesh-scan memory access. In the 8 fps ripple windows the avg/max inflate
  (isr_wake cpu 7.5 %, pack avg ~15 µs) from serial-dump preemption of the
  125 ms frames; the `min` floor (pack 4.25 µs) is stable.
- Net: the ISR machinery steals **~3 % of the chip** (~1.9 ms of the 62.5 ms
  window), leaving **~60.6 ms of render budget per window**. The mid-shape
  ripple render (72 ms) needs **~1.19×** to hold 1-window cadence; the
  dodecahedron-ambo holds (79 ms) need ~1.3×; the small-shape hold already fits.

## Summary ranking (mid-shape ripple, share of the 125 ms frame)

1. `scan_mesh_raster` (per-face SDF rasterize + shade) — **48%** (60.3 ms).
2. display-window sync (`is_buffer_wait`) — 42% (53.0 ms, idle by design).
3. `scan_face_setup` (per-face SDF/basis setup) — 6% (7.8 ms).
4. `is_mesh_transform` (ripple vertex deform + orient) — 3% (3.4 ms).
5. `is_face_offsets` (per-face segue sweep ordering) — <1% (0.34 ms).
6. `is_ripple_prepare` — negligible (27 ns).

IslamicStars is rasterizer-bound on the exact per-face SDF path: raster + setup
together are ~54 % of the frame; the ripple deformation it explicitly renders on
the exact path (not the LUT) is only ~3 %. This matches the perf ledger's
"IslamicStars is rasterizer-bound" finding — the WASM survey put it at ~9.3 ms
full-frame after the LUT-drop + convex-path win; on device at -Os with ISRs live
it is ~68–79 ms per **quadrant** in the 8 fps holds. The lever remains the
per-face raster (face count × on-screen face size), as the ledger's dead-end
list already established for spatial schemes at this face radius.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition.
- **`filter_blend` parenting**: here it parents cleanly under
  `scan_mesh_raster`, so its subtree and counts are valid; its printed % is of
  the frame.
- `scan_mesh_raster` / `scan_face_setup` are **pre-existing** render-side
  `HS_PROFILE` scopes inside `Scan::Mesh::draw` (only IslamicStars, of these
  four, uses that path), not added by this instrumentation.
- Segue-crossing frames drop most faces (terminator sweep) and are the cheapest;
  per-shape holds and ripple bursts are used as representatives.
- Epoch-boundary / segue windows (coarse `frame` count with a low first frame,
  e.g. min=17 ms) mix a torn-down shape's tail — skipped when picking
  representatives.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/IslamicStars.h` (`is_buffer_wait`, `is_ripple_prepare`,
  `is_timeline_step`, `is_draw_shape`, `is_mesh_transform`, `is_face_offsets`,
  `is_mesh_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=IslamicStars`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile IslamicStars [seconds]`.
