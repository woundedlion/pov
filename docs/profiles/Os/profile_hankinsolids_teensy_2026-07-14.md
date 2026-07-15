# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_hankinsolids.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `HankinSolids<288, 144>` only (single-entry playlist), current working-tree state (default carousel of Platonic/Archimedean solids) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile HankinSolids` |

Image size: FLASH code 65,540 B; ITCM (RAM1 code) 41,368 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
steady window frames 641–704 root counter = 2,405,350,795 cyc =
4,008,917.99 µs vs measured `micros()` window sum 4,008,920 µs
(Δ ≈ 0.5 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). HankinSolids is
the mesh/carousel sibling of IslamicStars: a solid is projected and each face
is SDF-rasterized (`Scan::Mesh::draw`), so the frame cost is essentially
*covered pixels × per-pixel SDF+blend*. Wall time snaps up to a whole number
of 62.5 ms windows.

The effect is **phased**: a new solid cycles in periodically, and each
transition is a **crossfade** where two meshes are drawn into the same frame.
That shows up two ways in the log:

- **Held single shape** — one `hk_draw_mesh` per frame (64 calls/window),
  render 31–48 ms → 1-window cadence, **62.5 ms/frame (16 fps)**.
- **Crossfade** — `hk_draw_mesh` call count rises above 64/window (e.g. **77**)
  as two shapes overlap; the frames that draw *both* meshes cost ≈ 2× and spill
  into a second window → **125 ms/frame (8 fps)**, and the rasterized face
  count (`scan_mesh_raster` calls) jumps. Large-solid crossfades push a few
  frames as high as **137–192 ms**.

Wall breathes from **~61 ms** (16 fps, small/held shape) to **~75–137 ms**
(8 fps, larger shape or crossfade). `hk_buffer_wait` (the `buffer_free()` spin
inside the `Canvas` ctor waiting for the display flip) is timed separately, so
the render numbers below are the clean draw cost; it *is* the round-up idle by
design.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame, calls/frame,
per-call cost (cycles = µs × 600). `scan_mesh_raster`, `scan_face_setup` and
`filter_blend` are **pre-existing render-side scopes** inside `Scan::Mesh::draw`
(not this effect's own instrumentation).

### Steady 16 fps hold — held solid (representative: frames 641–704)

```
frame                 62.64 ms  37.58 Mcyc  100%
  hk_timeline_step    42.14 ms  25.28 Mcyc   67%
    hk_draw_mesh      41.74 ms  25.05 Mcyc   67%  x1
      hk_mesh_scan    41.71 ms  25.02 Mcyc   67%  x1
        scan_mesh_raster 40.50 ms 24.30 Mcyc 65%  x62    653 us/face
          filter_blend   3.58 ms  2.15 Mcyc  6%  x12480  172 cyc/blend
        scan_face_setup  1.20 ms  0.72 Mcyc  2%  x62      19 us/face
      hk_mesh_transform   36 us     22 kcyc  0%  x1
    hk_update_hankin     174 us    105 kcyc  0%  x1
  hk_buffer_wait      20.50 ms  12.30 Mcyc   33%
```

Wall: min 58.1 / avg 62.6 / max 73.1 ms. One mesh per frame (a 62-face solid),
render 41.7 ms clears the 62.5 ms window → 16 fps, with 20.5 ms of
`hk_buffer_wait` idle rounding out the window. `scan_mesh_raster` (the per-face
SDF rasterize) is essentially all of it; `hk_mesh_transform` (camera transform)
and `hk_update_hankin` (per-frame mesh re-eval) together are **0.3 %** of the
frame.

### Crossfade — two solids overlapping (representative: frames 65–128)

```
frame                 75.25 ms  45.15 Mcyc  100%
  hk_timeline_step    40.62 ms  24.37 Mcyc   54%
    hk_draw_mesh      40.56 ms  24.33 Mcyc   54%  x1.20
      hk_mesh_scan    40.55 ms  24.33 Mcyc   54%  x1.20
        scan_mesh_raster 40.27 ms 24.16 Mcyc 54%  x15.5  2593 us/face
          filter_blend   3.99 ms  2.40 Mcyc  5%  x13529  177 cyc/blend
        scan_face_setup  0.28 ms  0.17 Mcyc  0%  x15.5    18 us/face
      hk_mesh_transform  9.6 us   5.7 kcyc   0%  x1.20
    hk_update_hankin    25.9 us    16 kcyc   0%  x0.75
  hk_buffer_wait      34.63 ms  20.78 Mcyc   46%
```

Wall: min 55.7 / avg 75.3 / **max 137.0 ms**. This is a small-shape crossfade
(octahedron → next solid). `hk_draw_mesh` fires **77 times over 64 frames** —
the 13 overlap frames each draw *two* meshes (~80 ms → 8 fps / 125 ms), the
rest draw one (16 fps); the per-frame counter averages the two cadences, so the
tree understates the peak (wall min 55.7 = single-mesh frames, max 137.0 =
double-mesh frames spilling two windows). Note the inversion vs the steady hold:
only ~15.5 faces/frame but **2.59 ms each** — a few huge faces (octahedron
covers the quadrant in 8) cost as much per pixel as the 62 tiny faces of a held
Archimedean solid. The largest-solid crossfades (frames 577–640, 705–768) reach
~190 faces/frame and avg 112–116 ms.

Across every window the story is identical: `scan_mesh_raster` is **94–99 %** of
`hk_draw_mesh`; `hk_mesh_transform` ≤ 0.2 % and `hk_update_hankin` ≤ 0.4 % of
the frame in all 14 windows.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h, under
`scan_mesh_raster`): **12,480 blended px/frame** in the steady hold and
**13,529** in the crossfade — both *exceed* the 10,368-px quadrant because
projected mesh faces overlap and each covered pixel is blended once per covering
face (1.20–1.30× coverage). Cost per blend is **172–177 cyc (~0.29 µs)**,
stable. But `scan_mesh_raster` spends **~1,790–1,950 cyc per *blended* pixel**
(40.3 Mcyc over ~12.5–13.5 k blends), i.e. ~10× the blend itself — most of the
scan cost is the per-face SDF distance evaluation over the tested band, not the
composite. The blend is ~9 % of the scan; the SDF is the other ~91 %.

## Column-ISR / DMA marshaling cost (`build/profile_hankinsolids.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative steady window (frames 641–704, 4.0 s). The pack
and submit run inside the flywheel wake, on the 1-in-8 wakes that render a
column. Columns: rate, per-call min/avg/max µs, CPU share:

```
isr_wake         18416/s  0.59 / 1.72 / ~171 us  cpu 3.2%
  isr_pack        2304/s  4.37 / 6.42 / ~169 us  cpu 1.5%
  isr_dma_submit  2304/s  0.64 / 0.96 / 1.2 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, nests
pack+submit); `isr_pack` = the 72× `packPixel` marshal, once per column;
`isr_dma_submit` = `submitFrame` (dcache flush + eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **7–15 µs of the 434 µs
  column period (~3.2 %)**. The wire transfer runs asynchronously on eDMA and
  costs no CPU.
- The `min` floor (pack **4.37 µs / 2,621 cyc**, submit ~0.6 µs) is the stable
  per-call cost; avg/max inflate in the long crossfade windows — isr_wake cpu
  climbs to **5.3–7.6 %** and pack avg to ~15 µs in frames 385–640 — because the
  serial dumps and higher-priority USB/DMA-completion ISRs preempt inside the
  scope during the slow 8 fps frames, not because pack work grew.
- Net: the ISR machinery steals **~3.2 %** of the chip in steady state (up to
  ~7.6 % transiently). Against the 62.5 ms window that is ~2.0 ms, leaving
  **~60 ms of render budget per window**. The held-solid render (41.7 ms) fits
  one window (16 fps); a crossfade double-mesh frame (~80 ms) needs ~1.3× → 8 fps.

## Summary ranking (crossfade regime, share of the 75.3 ms frame)

1. `scan_mesh_raster` (per-face SDF rasterize + shade + blend) — **54 %**
   (40.3 ms) — **99 % of active render**; the whole algorithmic cost.
2. display-window sync (`hk_buffer_wait`) — 46 % (34.6 ms, idle by design;
   inflated here by the mixed 16/8 fps cadence rounding up).
3. `scan_face_setup` (per-face setup) — 0.4 % (0.28 ms).
4. everything else — `hk_mesh_transform` (camera) 0.01 ms, `hk_update_hankin`
   (mesh re-eval) 0.03 ms, `hk_timeline_step` self — **< 0.1 %**.

In the steady hold the same leaf is **65 %** (40.5 ms) of the 62.6 ms frame,
buffer_wait 33 %. Either way `scan_mesh_raster` **is** the effect: mesh
transform and per-frame Hankin re-evaluation are ~0 %, exactly as expected for a
rasterizer-bound mesh carousel. Like its sibling IslamicStars, the only levers
are on the per-face SDF path (distance eval, face-band cull, supersample), not
geometry or camera. No committed WASM/native HankinSolids figures exist for
comparison.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`scan_mesh_raster` / `scan_face_setup` / `filter_blend` are render-side**
  scopes pre-existing in `Scan::Mesh::draw`, not added by this effect; the
  effect's own scopes are the `hk_*` nodes.
- **`filter_blend` counts every covering face's blend**, so blended px/frame
  exceeds the 10,368-px quadrant (face overlap); its per-blend cost is correct,
  its printed % is a fraction of its immediate parent.
- Per-pixel `filter_blend` scope overhead inflates the scan counters slightly;
  the coarse `hk_*` scopes are negligible.
- The **crossfade windows** mix 16 fps single-mesh and 8 fps double-mesh frames
  (77 draw calls / 64 frames); per-frame averages are valid means but understate
  the double-mesh peak — read the wall min/max for the spread.
- Epoch-start window (frames 1–64) shows a low `min` (32.1 ms) as the first
  solid ramps in — **skipped** when picking representatives.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/HankinSolids.h` (`hk_buffer_wait`, `hk_timeline_step`,
  `hk_update_hankin`, `hk_draw_mesh`, `hk_mesh_transform`, `hk_mesh_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=HankinSolids`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel pack,
  DMA submit) each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile HankinSolids [seconds]`.
