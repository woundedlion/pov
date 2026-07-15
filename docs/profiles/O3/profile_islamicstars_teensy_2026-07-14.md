# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-07-14, **-O3**)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_islamicstars_o3.log` (64-frame windows, single ~71 s pass). This
is the **-O3** twin of the shipping `-Os` report; the only variable
between the two is the optimization level.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: same Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but keeps base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forcing `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — the shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `IslamicStars<288, 144>` only (single-entry playlist), current working-tree state (`TerminatorSweep` segue, ripple amp 0.15, 4-ripple burst) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset; column-ISR counters via `HS_ISR_PROFILE` |
| Reproduce | `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 101,372 B; ITCM (RAM1 code) 56,696 B; RAM2 free 4,736 B.
(At -Os: 69,156 / 40,360 / 4,736 — see the `-O3 vs -Os` section for the cost.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
8 fps window frames 705–768 root counter = 4,802,035,626 cyc = 8,003,392.7 µs vs
measured `micros()` window sum 8,003,398 µs (Δ ≈ 0.66 ppm).

**Capture coverage** — the ~71 s pass (896 frames across 14 windows) spans three
full shapes plus a 4th spawned at the very end: **cube** `…hk675_expand5` (F=434,
windows 1–4), **truncatedIcosahedron** `…chamfer63` (F=452, windows 5–8),
**dodecahedron-ambo** `…hk66` (F=362 with big faces, windows 9–13), and
`truncatedIcosahedron-ambo` (F=542) spawned at the end with no window recorded.
Each shape's fade-in → still → 4-ripple burst → still → fade-out choreography is
seen; spawn/segue windows (low `min`, high `max`) are skipped when picking
representatives.

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × 144-column half ≈ **10,368 px**). Wall time snaps up to a whole
number of 62.5 ms windows.

The whole render is `Scan::Mesh::draw` (`is_mesh_scan`) — the exact per-face SDF
rasterize (deliberately *not* the congruence-class LUT, which ripple/segue
deformation would mis-shade). Cost tracks the **shape's on-screen face pixels**
(face count × on-screen face size) plus the **ripple deformation** mid-life. At
-O3 the render is fast enough that the mid-size shapes now land inside a single
window, so cadence steps *up* versus -Os:

- **cube (F=434)** — render (`is_timeline_step`) ≈ 42.6 ms < 62.5 ms →
  1-window, **62.5 ms/frame (16 fps)** through its whole life.
- **truncatedIcosahedron (F=452)** — render ≈ 50.8 ms in the ripple burst,
  still < 62.5 ms → **16 fps** (this shape was **8 fps at -Os**; see `-O3 vs -Os`).
- **dodecahedron-ambo (F=362, big faces)** — the fewer-but-larger faces push
  render to ≈ 69.4 ms > 62.5 ms → 2-window, **125 ms/frame (8 fps)**. This is the
  one regime -O3 does not lift: the render still crosses the window boundary.
- **Segue / spawn crossings** — the per-face terminator sweep drops most faces,
  so those frames are the **cheapest** (per-frame min 24–33 ms).

`is_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor) is the
round-up idle, timed separately. `is_ripple_prepare` (the per-frame ripple pool
advance) is ~12 ns — free.

## Phase-by-phase readout (64-frame windows, per-frame averages)

The `is_mesh_scan` subtree exposes pre-existing render scopes `scan_mesh_raster`
(the per-face rasterize) and `scan_face_setup` (per-face SDF/basis setup).
Columns: time/frame, cycles/frame, % of frame, calls/frame, per-call cost. Each
parent includes its children. Three regimes are captured across ≥ 2 steady
windows each (cube 2–4, truncatedIcosahedron 6–8, dodecahedron-ambo 11–12).

### cube, 16 fps hold+ripple (representative: frames 65–128, F=434)

```
frame                  62.51 ms  37.51 Mcyc  100%
  is_timeline_step     42.64 ms  25.59 Mcyc   68%
    is_draw_shape      42.62 ms  25.57 Mcyc   68%
      is_mesh_scan     42.11 ms  25.27 Mcyc   67%
        scan_mesh_raster 37.84 ms 22.70 Mcyc  61%  x434    87.2 us/face
          filter_blend  2.40 ms   1.44 Mcyc    4%  x16593  86.8 cyc/blend
        scan_face_setup  4.22 ms   2.53 Mcyc    7%  x434     9.7 us/face
      is_face_offsets    0.12 ms   75 kcyc     0%  x1
      is_mesh_transform  0.23 ms  227 kcyc     1%  x1
  is_buffer_wait       19.86 ms  11.92 Mcyc   32%
  is_ripple_prepare      12 ns      7 cyc      0%
```

Wall: min 60.5 / avg 62.5 / max 64.4 ms. The cube renders in ~42.6 ms — well
inside one 62.5 ms window — so it holds 16 fps the whole time, spending 32% of
the frame in `is_buffer_wait` idle. `scan_mesh_raster` (per-face SDF distance +
shade) is 61% of the frame; `scan_face_setup` adds 7%; `filter_blend` is 4%. All
434 faces are rasterized every frame (x434 = the cube's full face count; its
2,016 face-vertex indices are the per-face SDF workload, not extra faces).
`is_mesh_transform` is 0.23 ms with the mild mid-life ripple active.

### truncatedIcosahedron ripple burst, 16 fps (representative: frames 449–512, F=452)

```
frame                  62.46 ms  37.47 Mcyc  100%
  is_timeline_step     50.76 ms  30.45 Mcyc   81%
    is_draw_shape      50.72 ms  30.43 Mcyc   81%
      is_mesh_scan     49.05 ms  29.43 Mcyc   79%
        scan_mesh_raster 43.35 ms 26.01 Mcyc  69%  x452    95.9 us/face
          filter_blend  2.33 ms   1.40 Mcyc    4%  x16278  86.0 cyc/blend
        scan_face_setup  5.64 ms   3.38 Mcyc    9%  x452    12.5 us/face
      is_mesh_transform  1.53 ms   0.92 Mcyc    2%  x1  (ripple deform)
      is_face_offsets    0.14 ms   86 kcyc     0%  x1
  is_buffer_wait       11.70 ms   7.02 Mcyc   19%
  is_ripple_prepare      12 ns      7 cyc      0%
```

Wall: min 60.4 / avg 62.5 / max 64.8 ms. This is the direct counterpart of the
-Os "mid-shape ripple burst" representative — which ran at **8 fps / 125 ms**. At
-O3 the same rippled mesh renders in 50.8 ms, inside one window, so it holds
**16 fps**. `is_mesh_transform` (the `RippleTransformer` displacing every vertex +
orient) is 1.53 ms during the burst; the rippled mesh presents steeper faces so
`scan_mesh_raster` per-face cost is 95.9 µs (vs 87.2 µs for the cube hold) and
blended pixels rise to 16.3k.

### dodecahedron-ambo, 8 fps (representative: frames 705–768, F=362 big faces)

```
frame                 125.05 ms  75.03 Mcyc  100%
  is_timeline_step     69.44 ms  41.67 Mcyc   55%
    is_draw_shape      69.41 ms  41.64 Mcyc   55%
      is_mesh_scan     67.27 ms  40.36 Mcyc   54%
        scan_mesh_raster 61.32 ms 36.79 Mcyc  49%  x362   169.4 us/face
          filter_blend  2.51 ms   1.51 Mcyc    2%  x15870  94.9 cyc/blend
        scan_face_setup  5.91 ms   3.54 Mcyc    5%  x362    16.3 us/face
      is_mesh_transform  1.99 ms   1.19 Mcyc    2%  x1  (ripple deform)
      is_face_offsets    0.14 ms   85 kcyc     0%  x1
  is_buffer_wait       55.61 ms  33.37 Mcyc   44%
  is_ripple_prepare      12 ns      7 cyc      0%
```

Wall: min 123.5 / avg 125.1 / max 126.5 ms. Fewer faces than the
truncatedIcosahedron (362 vs 452) but each is **much larger** (169.4 µs/face vs
95.9 µs). That pushes render to 69.4 ms > 62.5 ms, so the frame crosses into
2-window / **8 fps** and 44% of the frame becomes `is_buffer_wait` idle. This is
the one shape whose per-face size -O3 cannot shrink under the window; the F=542
`truncatedIcosahedron-ambo` spawned at the end would behave the same.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h, parenting under
`scan_mesh_raster`): **16,278 blends/frame** in the truncatedIcosahedron ripple
(1.57× the 10,368-px quadrant — the rippled mesh overdraws the band), **16,593**
for the cube hold, **15,870** for the dodecahedron-ambo 8 fps frame — at **86–95
cyc (0.14–0.16 µs) per blend** (roughly half the -Os per-blend cost).
`scan_mesh_raster` ≈ 22.7–36.8 Mcyc over those blends ⇒ **~1,370–2,320 cyc per
blended pixel** (vs ~2,200–3,530 at -Os): the cost is still overwhelmingly the
exact per-face SDF distance evaluation, not the blend (only ~2–4% of frame).

## Column-ISR / DMA marshaling cost (same log)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative 16 fps window (frames
449–512, 4.0 s). Columns: rate, per-call min/avg/max (µs), CPU share:

```
isr_wake         18432/s  0.55 / 1.73 / 174.7 us  cpu 3.19%
  isr_pack        2304/s  4.76 / 6.85 / 172.9 us  cpu 1.57%
  isr_dma_submit  2304/s  0.61 / 0.93 / 1.27 us   cpu 0.21%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick is ~0.9 µs of CPU.** IslamicStars keeps its clean mesh-scan
  memory access at -O3: base ISR steal is ~3.2% (essentially unchanged from -Os's
  3.00% — this is display-driver work, not effect code). The `min` floors (pack
  2856 cyc ≈ 4.76 µs, submit 364 cyc ≈ 0.61 µs) are the stable per-call costs;
  avg/max inflate from serial-dump preemption of long frames.
- In the 8 fps big-face windows the avg/max balloon (frames 705–768: isr_wake cpu
  7.21%, isr_pack avg ~14.6 µs) — that is serial-dump preemption landing inside
  the 125 ms frames, not real pack work; the `min` floor is identical.
- Net: the ISR machinery steals **~3.2% of the chip** (~2.0 ms of the 62.5 ms
  window), leaving **~60.5 ms of render budget per window**. The mid-shape ripple
  render (50.8 ms) now fits 1-window (16 fps); only the big-face
  dodecahedron-ambo render (69.4 ms) needs ~1.15× and stays 8 fps.

## Summary ranking (dodecahedron-ambo 8 fps, share of the 125 ms frame)

1. `scan_mesh_raster` (per-face SDF rasterize + shade) — **49%** (61.3 ms).
2. display-window sync (`is_buffer_wait`) — 44% (55.6 ms, idle by design).
3. `scan_face_setup` (per-face SDF/basis setup) — 5% (5.9 ms).
4. `is_mesh_transform` (ripple vertex deform + orient) — 2% (2.0 ms).
5. `is_face_offsets` (per-face segue sweep ordering) — <1% (0.14 ms).
6. `is_ripple_prepare` — negligible (12 ns).

IslamicStars is rasterizer-bound on the exact per-face SDF path at -O3 too: raster
+ setup together are ~54% of the frame; the ripple deform it explicitly renders on
the exact path (not the LUT) is only ~2%. This matches the perf ledger's
"IslamicStars is rasterizer-bound" finding — the lever remains the per-face raster
(face count × on-screen face size), exactly what makes the big-face
dodecahedron-ambo the sole shape that stays 8 fps here.

## -O3 vs -Os

The headline: **-O3 roughly halves render time on the per-face SDF path, which
lifts the mid-size shapes from 8 fps to 16 fps** — because their render drops
below the 62.5 ms window boundary. Only the big-face dodecahedron-ambo still
exceeds one window and stays 8 fps at both levels. Comparing the shared
"mid-shape ripple burst" representative (truncatedIcosahedron, F=452, per-frame /
per-call):

| metric | -Os | -O3 | speedup |
|---|---|---|---|
| render `is_timeline_step` (ms/frame) | 71.95 | 50.76 | **1.42×** (−21.2 ms) |
| `is_mesh_scan` (ms/frame) | 68.18 | 49.05 | 1.39× (−19.1 ms) |
| `scan_mesh_raster` (ms/frame) | 60.30 | 43.35 | **1.39×** (−17.0 ms) |
| `scan_mesh_raster` per face (µs) | 133 | 95.9 | **1.39×** |
| `scan_face_setup` per face (µs) | 17 | 12.5 | 1.36× |
| `is_mesh_transform` (ripple, ms) | 3.39 | 1.53 | 2.22× |
| `filter_blend` (cyc/blend) | 192 | 86.0 | **2.23×** |
| scan cyc per blended px | ~3,530 | ~1,598 | 2.21× |
| frame wall / cadence | 125 ms / 8 fps | 62.5 ms / **16 fps** | **2× fps** |
| FLASH code | 69,156 B | 101,372 B | **+32,216 B (+46.6%)** |
| ITCM (RAM1 code) | 40,360 B | 56,696 B | **+16,336 B (+40.5%)** |

The tightest per-pixel/per-vertex loops (`filter_blend`, `is_mesh_transform`)
gain most from -O3 (~2.2×); the SDF scan and per-face setup gain ~1.4×. Because
the mid-shape render clears the 62.5 ms window at -O3, the reclaimed time does not
disappear into idle the way it does for a shape stuck above the boundary — it
converts directly into a doubled frame rate. The dodecahedron-ambo (F=362 big
faces) is the exception: its 69.4 ms render still exceeds one window, so there the
saved time re-emerges as `is_buffer_wait` idle (44% of the 125 ms frame) and
cadence stays 8 fps.

**-O3 is not the shipping config.** The full 26-effect Phantasm image overflows
FlexRAM/ITCM at -O3 (a single effect already costs +40.5% ITCM here), which is
exactly why the ship build is `-Os`. -O3 is only viable as this single-effect
profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, but not
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: it parents cleanly under `scan_mesh_raster` in
  every steady window here, so its subtree and counts are valid; the printed % is
  of the frame (recomputed, not the log's parent-relative %).
- `scan_mesh_raster` / `scan_face_setup` are **pre-existing** render-side
  `HS_PROFILE` scopes inside `Scan::Mesh::draw` (only IslamicStars uses that
  path), not added by this instrumentation.
- Segue-crossing / spawn frames drop most faces (terminator sweep) and are the
  cheapest; per-shape holds and ripple bursts are used as representatives.
- Spawn/epoch windows (coarse `frame` count with a low first frame, e.g. min=24.7
  or 33.1 ms) mix a torn-down shape's tail — skipped when picking representatives.
- **-O3 build**: this is the `profile_o3` env, which does **not** ship. The
  shipping Phantasm image is `-Os` because the full 26-effect roster overflows
  ITCM at -O3.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/IslamicStars.h` (`is_buffer_wait`, `is_ripple_prepare`,
  `is_timeline_step`, `is_draw_shape`, `is_mesh_transform`, `is_face_offsets`,
  `is_mesh_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=IslamicStars`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py`
  (the `-Os` twin uses env `profile`, reproducible via `just profile IslamicStars`).
