# HankinSolids on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_hankinsolids_o3.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the **-O3** twin of the shipping `-Os` profile: identical Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`, `-D HS_PROFILE_ENABLE`) but keeps the base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forcing `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `HankinSolids<288, 144>` only (single-entry playlist), current working-tree state (default carousel of Platonic/Archimedean solids) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code **94,060 B**; ITCM (RAM1 code) **59,256 B**; RAM2 free 4,736 B.
(-Os shipping build: 65,540 / 41,368 / 4,736 — see the [-O3 vs -Os](#-o3-vs--os) section.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm.
Verified in **three** clean windows:

- steady hold, frames 641–704: root `frame` = 2,404,573,135 cyc = 4,007,621.89 µs
  vs `micros()` window sum 4,007,628 µs (Δ ≈ 1.5 ppm);
- crossfade, frames 65–128: 2,403,753,790 cyc = 4,006,256.3 µs vs sum
  4,006,256 µs (Δ ≈ 0.08 ppm);
- held solid, frames 321–384: 2,402,076,762 cyc = 4,003,461.3 µs vs sum
  4,003,475 µs (Δ ≈ 3.4 ppm).

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
  render ~24–33 ms → 1-window cadence, **62.5 ms/frame (16 fps)**.
- **Crossfade** — `hk_draw_mesh` call count rises to **77/window** as two shapes
  overlap; the 13 overlap frames draw *two* meshes. At -O3 a **small**-solid
  double-mesh frame (~49 ms render) still fits one 62.5 ms window → stays
  **16 fps**, but a **large**-solid double-mesh frame (~99 ms) spills a second
  window → **125 ms/frame (8 fps)**. The rasterized face count
  (`scan_mesh_raster` calls) jumps into the thousands during large crossfades.

Wall breathes from **~44 ms** (16 fps, small/held shape) to **~130–150 ms**
(8 fps, largest-solid crossfade). `hk_buffer_wait` (the `buffer_free()` spin
inside the `Canvas` ctor waiting for the display flip) is timed separately, so
the render numbers below are the clean draw cost; it *is* the round-up idle by
design — and at -O3 it grows, because the faster render leaves more slack inside
the fixed window.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Per-frame averages, nested as in the counter tree (each parent includes its
children). Columns: time/frame, cycles/frame, % of frame, calls/frame,
per-call cost (cycles = µs × 600). `scan_mesh_raster`, `scan_face_setup` and
`filter_blend` are **pre-existing render-side scopes** inside `Scan::Mesh::draw`
(not this effect's own instrumentation).

### Steady 16 fps hold — held solid (representative: frames 641–704)

```
frame                 62.62 ms  37.57 Mcyc  100%
  hk_timeline_step    32.71 ms  19.62 Mcyc   52%
    hk_draw_mesh      32.46 ms  19.48 Mcyc   52%  x1
      hk_mesh_scan    32.45 ms  19.47 Mcyc   52%  x1
        scan_mesh_raster 31.51 ms 18.91 Mcyc 50%  x62    508 us/face
          filter_blend   1.82 ms  1.09 Mcyc  3%  x12480   87 cyc/blend
        scan_face_setup  0.93 ms  0.56 Mcyc  1%  x62      15 us/face
      hk_mesh_transform   11 us   6.7 kcyc   0%  x1
    hk_update_hankin     102 us    61 kcyc   0%  x1
  hk_buffer_wait      29.91 ms  17.95 Mcyc   48%
```

Wall: min 53.8 / avg 62.6 / max 74.9 ms. One mesh per frame (the 62-face
truncated icosidodecahedron), render 32.5 ms clears the 62.5 ms window → 16 fps,
with **29.9 ms of `hk_buffer_wait` idle** rounding out the window (vs 20.5 ms at
-Os — the win shows up as idle inside a window-locked cadence). `scan_mesh_raster`
(the per-face SDF rasterize) is essentially all of it; `hk_mesh_transform`
(camera transform) and `hk_update_hankin` (per-frame mesh re-eval) together are
**0.3 %** of the frame.

### Crossfade — two solids overlapping (representative: frames 65–128)

```
frame                 62.60 ms  37.56 Mcyc  100%
  hk_timeline_step    29.41 ms  17.65 Mcyc   47%
    hk_draw_mesh      29.36 ms  17.62 Mcyc   47%  x1.20
      hk_mesh_scan    29.36 ms  17.62 Mcyc   47%  x1.20
        scan_mesh_raster 29.14 ms 17.49 Mcyc 47%  x15.5  1876 us/face
          filter_blend   1.96 ms  1.17 Mcyc  3%  x13536   87 cyc/blend
        scan_face_setup  0.21 ms  0.13 Mcyc  0%  x15.5    14 us/face
      hk_mesh_transform  2.8 us   1.7 kcyc   0%  x1.20
    hk_update_hankin    16.3 us   9.8 kcyc   0%  x0.75
  hk_buffer_wait      33.19 ms  19.91 Mcyc   53%
```

Wall: min 44.3 / avg 62.6 / **max 89.9 ms**. This is a small-shape crossfade
(octahedron → truncatedOctahedron). `hk_draw_mesh` fires **77 times over 64
frames** — 13 overlap frames each draw *two* meshes. Per-mesh render is
**24.4 ms**, so a single-mesh frame ≈ 24 ms and a double-mesh frame ≈ 49 ms —
**both under the 62.5 ms window**, which is the key -O3 change: the whole
crossfade holds **16 fps** (wall avg 62.6 = one window/frame), whereas at -Os
these ~80 ms double-mesh frames spilled a second window and dropped to 8 fps.
The max 89.9 ms is per-frame flip jitter, not a doubled cadence. Note the
per-face inversion vs the steady hold: only ~15.5 faces/frame but **1.88 ms
each** — a few huge faces (octahedron covers the quadrant in 8) cost as much per
pixel as the 62 tiny faces of a held solid. The **largest**-solid crossfades
(frames 577–640, 705–768: ~190 faces/frame) still spill to 8 fps, avg 80.8–83.2 ms
(vs 112–116 ms at -Os).

Across every window the story is identical: `scan_mesh_raster` is **94–99 %** of
`hk_draw_mesh`; `hk_mesh_transform` ≤ 0.1 % and `hk_update_hankin` ≤ 0.4 % of
the frame in all 15 windows.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h, under
`scan_mesh_raster`): **12,480 blended px/frame** in the steady hold and
**13,536** in the crossfade — both *exceed* the 10,368-px quadrant because
projected mesh faces overlap and each covered pixel is blended once per covering
face (1.20–1.30× coverage). Cost per blend is **87 cyc (~0.15 µs)**, stable —
**about half the 172–177 cyc of the -Os build** (the blend is the biggest single
-O3 win; see below). `scan_mesh_raster` spends **~1,290–1,515 cyc per *blended*
pixel** (vs ~1,790–1,950 at -Os), i.e. still ~15× the blend itself — most of the
scan cost is the per-face SDF distance evaluation over the tested band, not the
composite. At -O3 the blend is ~6 % of the scan; the SDF is the other ~94 %.

## Column-ISR / DMA marshaling cost (`build/profile_hankinsolids_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative steady window (frames 641–704, 4.0 s). The pack
and submit run inside the flywheel wake, on the 1-in-8 wakes that render a
column. Columns: rate, per-call min/avg/max µs, CPU share:

```
isr_wake         18416/s  0.46 / 1.69 / ~169 us  cpu 3.1%
  isr_pack        2304/s  4.65 / 6.68 / ~167 us  cpu 1.5%
  isr_dma_submit  2304/s  0.61 / 0.91 / 1.3 us   cpu 0.21%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, nests
pack+submit); `isr_pack` = the 72× `packPixel` marshal, once per column;
`isr_dma_submit` = `submitFrame` (dcache flush + eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **7–15 µs of the 434 µs
  column period (~3.1 %)**. The wire transfer runs asynchronously on eDMA and
  costs no CPU.
- The `min` floor (pack **4.65 µs / 2,789 cyc**, submit ~0.6 µs) is the stable
  per-call cost; avg/max inflate in the long crossfade windows — isr_wake cpu
  climbs to **5.4–5.9 %** and pack avg to ~10–11 µs in frames 449–640/833–960 —
  because the serial dumps and higher-priority USB/DMA-completion ISRs preempt
  inside the scope during the slow 8 fps frames, not because pack work grew.
- The ISR path is **essentially unchanged from -Os** (pack floor 4.65 vs 4.37 µs,
  steady cpu 3.1 vs 3.2 %): it is fixed-cost pack/DMA marshal that the -O3
  render speedup does not touch. Net: the ISR machinery steals **~3.1 %** of the
  chip in steady state (up to ~5.9 % transiently), ~1.9 ms against the 62.5 ms
  window, leaving **~60 ms of render budget per window**.

## Summary ranking (crossfade regime, share of the 62.6 ms frame)

1. `scan_mesh_raster` (per-face SDF rasterize + shade + blend) — **47 %**
   (29.1 ms) — **99 % of active render**; the whole algorithmic cost.
2. display-window sync (`hk_buffer_wait`) — 53 % (33.2 ms, idle by design; at
   -O3 it *exceeds* render because the faster draw leaves more slack inside the
   fixed window).
3. `scan_face_setup` (per-face setup) — 0.3 % (0.21 ms).
4. everything else — `hk_mesh_transform` (camera) 0.003 ms, `hk_update_hankin`
   (mesh re-eval) 0.02 ms, `hk_timeline_step` self — **< 0.1 %**.

In the steady hold the same leaf is **50 %** (31.5 ms) of the 62.6 ms frame,
buffer_wait 48 %. Either way `scan_mesh_raster` **is** the effect: mesh
transform and per-frame Hankin re-evaluation are ~0 %, exactly as expected for a
rasterizer-bound mesh carousel. Like its sibling IslamicStars, the only levers
are on the per-face SDF path (distance eval, face-band cull, supersample), not
geometry or camera.

## -O3 vs -Os

The only variable between this run and the shipping-config report
(`docs/profiles/Os/profile_hankinsolids_teensy_2026-07-14.md`) is the
optimization level. Both windows are the **same shapes** (identical face and
blend counts — 62 faces / 12,480 blends steady, 15.5 faces / ~13.5 k blends in
the octahedron crossfade), so the deltas are a clean A/B on the SDF+blend path:

| Metric | -Os | -O3 | speedup |
|---|---|---|---|
| `scan_mesh_raster` per-face SDF (steady hold) | 653 µs | 508 µs | **1.29×** |
| scan cyc per *blended* pixel (steady) | 1,947 | 1,515 | **1.29×** |
| `filter_blend` per-blend | 172 cyc | 87 cyc | **1.97×** |
| `scan_mesh_raster` per-face (crossfade) | 2,593 µs | 1,876 µs | **1.38×** |
| render `hk_draw_mesh`/frame (steady) | 41.7 ms | 32.5 ms | 1.29× (−9.2 ms) |
| FLASH code | 65,540 B | 94,060 B | +28,520 (+43.5 %) |
| ITCM code | 41,368 B | 59,256 B | +17,888 (+43.2 %) |

Headline: **-O3 speeds the per-face SDF rasterize by ~1.29×** (508 vs 653 µs/face
in the steady hold; up to 1.38× on the big-face crossfade) and **nearly halves
the per-blend composite (172 → 87 cyc, 1.97×)** — the blend, being a tight
fixed-shape kernel, benefits most from `-ffast-math`/vectorization while the
branch-heavy SDF distance loop gains a more modest 1.3×.

**Regime crossing a window boundary:** the small-shape (octahedron) crossfade.
At -Os its ~80 ms double-mesh frames spilled a second 62.5 ms window and dropped
to a mixed 8/16 fps (wall avg 75.3 ms). At -O3 the double-mesh render falls to
~49 ms — under the 62.5 ms window — so the whole crossfade **holds 16 fps** (wall
avg 62.6 ms). Large-solid crossfades (~99 ms double-mesh) still spill to 8 fps in
both builds, but their wall drops from ~112–116 ms to ~81–83 ms. In the
window-locked **steady hold** the cadence is unchanged (16 fps both builds); the
-O3 win there materializes as **+9.4 ms more `hk_buffer_wait` idle** (29.9 vs
20.5 ms), not a higher frame rate.

**Cost:** -O3 is **not** the shipping config. FLASH grows +43.5 % and ITCM
+43.2 % for this single-effect image; the full 26-effect Phantasm roster
overflows FlexRAM/ITCM at -O3, which is exactly why the ship build is -Os. -O3
is only viable here as a single-effect profiling image.

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
- The **crossfade windows** mix 16 fps single-mesh and (for large solids) 8 fps
  double-mesh frames (77 draw calls / 64 frames); per-frame averages are valid
  means but can understate the double-mesh peak — read the wall min/max for the
  spread.
- Epoch-start window (frames 1–64) shows a low `min` (23.1 ms) as the first
  solid ramps in — **skipped** when picking representatives.
- **-O3 build** (the `profile_o3` env) — this **does not ship**. The shipping
  Phantasm image is `-Os` because the full 26-effect roster overflows ITCM at
  -O3; this single-effect -O3 image is a profiling artifact only.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/HankinSolids.h` (`hk_buffer_wait`, `hk_timeline_step`,
  `hk_update_hankin`, `hk_draw_mesh`, `hk_mesh_transform`, `hk_mesh_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=HankinSolids`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel pack,
  DMA submit) each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py`
  (the -O3 twin; the -Os `just profile HankinSolids` recipe builds the shipping
  `profile` env).
