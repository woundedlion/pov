# IslamicStars per-shape on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw captures:
`build/profile_islamicstars_pershape.log` (-Os) and
`build/profile_islamicstars_pershape_o3.log` (-O3), 16-frame windows, ~210 s
single pass each. This report times **every one of the 24 Islamic star-pattern
solids** individually, with per-shape geometry metadata, at both optimization
levels. It complements the aggregate IslamicStars reports under `Os/` and `O3/`
(those pick one representative shape; this one measures the whole roster).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` (-Os, shipping) and `profile_o3` (-O3) PlatformIO envs; Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `IslamicStars<288, 144>`, cycling all 24 solids from `Solids::Collections::get_islamic_solids()` |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE`, 600 cyc/µs), 16-frame windows; each `Spawning Shape:` serial line tags the windows that follow it |
| Carousel | `Trans Speed = 4` (new slider) set from the harness (`-D HS_PROFILE_TRANS_SPEED=4`) so all 24 shapes walk past inside one epoch |
| Reproduce | `pio run -e profile -t upload` (or `profile_o3`) with `PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=IslamicStars -D HS_PROFILE_WINDOW=16 -D HS_PROFILE_TRANS_SPEED=4'`, then `tools/profile_capture.py --seconds 210` |

Image size (single-effect profile image): **-Os** FLASH code 69,428 B, ITCM
40,504 B; **-O3** FLASH 101,596 B, ITCM 56,872 B (+46 % / +40 %). RAM2 free
4,736 B both.

**Exactness cross-check** — window frames 1-16 (-Os, shape 0) root counter
581,351,051 cyc = 968,918.4 µs vs measured `micros()` window sum 968,923 µs
(Δ ≈ 4.7 ppm; the 16-frame window's smaller sum widens the ppm vs the 64-frame
runs but the counters still track the wall clock exactly).

## Method — how each shape is isolated

IslamicStars shows one solid at a time; the `TerminatorSweep` segue is a
**within-shape** reveal/conceal (no two-shape crossfade). The per-shape cost is
read from the shape's **fullest clean hold**: among its 16-frame windows, those
whose `scan_mesh_raster` call count is exactly **16·F** (every face drawn once
per frame), taking the one with maximum `is_mesh_scan`. Sweep-in/out frames draw
only the revealed part of the mesh (calls < 16·F), and during the spawn
crossover the outgoing sprite is still drawing (calls > 16·F); both kinds of
window mis-state the shape's own cost and are excluded. Every shape has 3+
clean windows in each capture.

The **`Trans Speed` slider** (added to IslamicStars: divides every per-shape
stage length — fade, still holds, ripple span — by its value, default 1 =
shipping cadence) only shortens how long each shape *holds*. It does **not**
change the per-frame render cost, so the `is_mesh_scan` ms and the resulting fps
below are the real shipping figures; `Trans Speed = 4` is used purely to fit all
24 shapes into one 210 s capture. The harness sets it after `init()` via
`updateParameter("Trans Speed", 4)` (a no-op for effects without that param).

The counter tree per frame:

```
frame
  is_timeline_step
    is_draw_shape
      is_mesh_scan          <- per-shape rasterize cost (the headline column)
        scan_mesh_raster     x = faces drawn/frame (one call per drawn face)
          filter_blend       x = AA-blended pixels/frame
        scan_face_setup
      is_face_offsets        per-face terminator ordering (small)
      is_mesh_transform      ripple + orient (small)
  is_ripple_prepare
  is_buffer_wait
```

## Per-shape table (ranked by -Os `is_mesh_scan`, the mesh rasterize cost)

`V`/`F` = solid vertex/face count (matching the solids tool and the `Spawning
Shape` log). `I` = face-vertex indices (Σ face side counts = 2·E) — the flat
index workload the rasterizer walks. `blend/f` = `filter_blend` calls
per frame = AA-blended pixels (coverage, -Os window). `scan` = `is_mesh_scan`
ms/frame. `O3×` = -Os ÷ -O3 mesh-scan. `render` = frame − buffer_wait (-Os).
`fps` = observed cadence (16/8/5 fps = 1/2/3 display windows) at each level.
Every value is read from the shape's fullest clean hold (see Method), where
`scan_mesh_raster` runs exactly F calls per frame — one per face, so a separate
faces-per-frame column would duplicate `F`.

| # | Shape | V | F | I | blend/f | scan Os | scan O3 | O3× | render Os | fps Os | fps O3 |
|--:|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` | 3240 | 1082 | 8640 | 18058 | 133.3 | 110.1 | 1.21× | 139.5 | 5 | 8 |
| 2 | `truncatedIcosahedron_ambo_relax_truncate001_hankin73` | 1620 | 542 | 4320 | 17436 | 110.2 | 92.6 | 1.19× | 111.8 | 8 | 8 |
| 3 | `truncatedIcosahedron_ambo_relax_truncate001_hankin59` | 1620 | 542 | 4320 | 15809 | 110.1 | 92.1 | 1.20× | 114.8 | 8 | 8 |
| 4 | `truncatedIcosidodecahedron_bevel5_relax_hk77` | 2160 | 722 | 5760 | 18728 | 99.6 | 80.7 | 1.23× | 102.8 | 8 | 8 |
| 5 | `truncatedIcosahedron_ambo_relax_truncate33_hk64` | 1620 | 542 | 4320 | 16867 | 99.0 | 83.2 | 1.19× | 104.1 | 8 | 8 |
| 6 | `truncatedIcosahedron_hk54_ambo_hk72` | 1620 | 542 | 4320 | 17219 | 98.8 | 84.1 | 1.17× | 103.9 | 8 | 8 |
| 7 | `truncatedOctahedron_gyro_kis_hk17` | 1620 | 542 | 4320 | 18094 | 95.3 | 79.6 | 1.20× | 99.0 | 8 | 8 |
| 8 | `icosahedron_snub_relax_truncate033_hankin62` | 1350 | 452 | 3600 | 16564 | 85.7 | 70.9 | 1.21× | 89.9 | 8 | 8 |
| 9 | `dodecahedron_ambo_bevel33_relax_hk66` | 1080 | 362 | 2880 | 15714 | 79.7 | 66.9 | 1.19× | 82.7 | 8 | 8 |
| 10 | `rhombicuboctahedron_hk63_ambo_hk63` | 864 | 290 | 2304 | 14797 | 74.6 | 64.3 | 1.16× | 76.7 | 8 | 8 |
| 11 | `dodecahedron_hk72_ambo_dual_hk20` | 540 | 182 | 1440 | 14631 | 73.3 | 53.6 | 1.37× | 73.9 | 8 | **16** |
| 12 | `truncatedIcosidodecahedron_truncate50d_ambo_dual` | 542 | 540 | 2160 | 19334 | 72.8 | 48.3 | 1.51× | 74.4 | 8 | **16** |
| 13 | `dodecahedron_hk54_ambo_hk72` | 540 | 182 | 1440 | 14793 | 70.6 | 53.5 | 1.32× | 72.5 | 8 | **16** |
| 14 | `icosahedron_bevel033_hk59` | 540 | 182 | 1440 | 13564 | 68.1 | 51.5 | 1.32× | 70.0 | 8 | **16** |
| 15 | `icosahedron_ambo_truncate033_hankin59` | 540 | 182 | 1440 | 14014 | 68.0 | 51.2 | 1.33× | 69.9 | 8 | **16** |
| 16 | `truncatedIcosahedron_hk58_chamfer63` | 990 | 452 | 2880 | 16144 | 66.5 | 49.1 | 1.35× | 70.0 | 8 | **16** |
| 17 | `cube_relax_bevel33_relax_hk675_expand5` | 576 | 434 | 2016 | 16645 | 62.1 | 44.0 | 1.41× | 64.3 | 8 | **16** |
| 18 | `octahedron_hk17_ambo_hk73` | 216 | 74 | 576 | 13108 | 53.8 | 42.8 | 1.26× | 54.0 | 16 | 16 |
| 19 | `octahedron_hk34_ambo_hk72` | 216 | 74 | 576 | 13089 | 51.3 | 41.6 | 1.23× | 51.5 | 16 | 16 |
| 20 | `snubDodecahedron_truncate5d_ambo_dual` | 452 | 450 | 1800 | 16658 | 47.8 | 35.9 | 1.33× | 49.3 | 16 | 16 |
| 21 | `truncatedIcosahedron_truncate50d_ambo_dual` | 272 | 270 | 1080 | 15813 | 45.6 | 35.1 | 1.30× | 46.2 | 16 | 16 |
| 22 | `dodecahedron_bevel2_relax_gyro` | 542 | 360 | 1800 | 15670 | 45.0 | 34.0 | 1.32× | 46.7 | 16 | 16 |
| 23 | `icosahedron_kis_gyro` | 272 | 180 | 900 | 14488 | 40.1 | 29.9 | 1.34× | 40.9 | 16 | 16 |
| 24 | `icosidodecahedron_truncate5d_ambo_dual` | 182 | 180 | 720 | 14092 | 38.7 | 28.6 | 1.35× | 38.9 | 16 | 16 |

**Bold fps O3** = a shape that crosses from 8 fps to 16 fps at -O3.

## Analysis

- **Range**: mesh-scan spans **38.7 → 133.3 ms** at -Os (a 3.4× spread across
  the roster) and 28.6 → 110.1 ms at -O3. Seven shapes (rows 18–24: the
  octahedron pair, the truncate-dual family, `dodecahedron_bevel2_relax_gyro`,
  `icosahedron_kis_gyro`) hold 16 fps at -Os; the other 17 are 8 fps (2-window)
  shapes except the heaviest (`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`,
  1,082 faces), a **5 fps** (3-window) shape.
- **No single geometry column predicts cost.** AA coverage is nearly flat
  (`blend/f` 13.1–19.3 k for every shape — the patterns are all space-filling),
  so blended-pixel count does not differentiate the roster. At equal F the cost
  splits ~1.8× by face complexity: the 182-face hankin dodecahedra (I=1,440,
  ~8 sides/face) cost 68–73 ms while the 180-face `icosahedron_kis_gyro` /
  `icosidodecahedron_truncate5d_ambo_dual` (I=720–900, 4–5 sides/face) cost
  38.7–40.1 ms; likewise row 9 (F=362, I=2,880) at 79.7 ms vs row 22 (F=360,
  I=1,800) at 45.0 ms. The octahedron pair pays 51–54 ms with only 74 faces —
  few huge faces probe large pixel areas. Per-face setup is the F-proportional
  part: `scan_face_setup` is 22.6 ms/frame on shape 1 (1,082 faces, -Os),
  ~17 % of its mesh-scan.
- **-O3 buys a mean 1.28× (range 1.16–1.51×)** on mesh-scan. The ten heaviest
  shapes (rows 1–10, all ≥74.6 ms) gain the least (1.16–1.23×); the lighter
  half gains 1.23–1.51×, topped by 12 `..truncate50d_ambo_dual` (1.51×) and 17
  `cube..hk675` (1.41×). **Seven shapes cross 8 → 16 fps at -O3** (rows 11–17)
  and the heaviest drops 5 → 8 fps — cutting the sub-16-fps roster from 17
  shapes to 10.

## Caveats

- **`is_mesh_scan` absorbs live ISR time** (`CYCCNT` free-runs) — the real
  shipping condition, not pure-CPU cost.
- **Representative window = fullest clean hold.** Per shape the
  max-`is_mesh_scan` window among those with exactly 16·F `scan_mesh_raster`
  calls is used; the ripple burst may or may not overlap it, but the
  ripple deforms the mesh without changing the visible-face count, so the scan
  cost is stable (ripple lands in `is_mesh_transform`, separately tracked and
  <1 % of the frame). `fps` is the worst-case (fully-visible) cadence for the
  shape.
- **`Trans Speed = 4`** compresses the dwell to fit the roster in one capture; it
  does not change per-frame cost (verified: the same shapes read the same
  mesh-scan ms at Trans Speed 1 in the aggregate reports). Ripple/fade timing is
  4× faster here than shipping.
- The `Trans Speed` slider + the `HS_PROFILE_TRANS_SPEED` harness knob are new
  working-tree code (`effects/IslamicStars.h`, `targets/Profile/Profile.ino`);
  the slider defaults to 1 (shipping cadence unchanged) and is a real control,
  not profiling-only scaffolding.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/IslamicStars.h` (`is_*` scopes) and `core/render/scan.h`
  (`scan_mesh_raster`/`scan_face_setup`).

## Harness

- `pio run -e profile -t upload` / `profile_o3` with
  `PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=IslamicStars -D HS_PROFILE_WINDOW=16 -D HS_PROFILE_TRANS_SPEED=4'`,
  then `python tools/profile_capture.py --seconds 210 --out build/profile_islamicstars_pershape[_o3].log`.
- Per-shape attribution: `scratchpad/parse_islamic.py` correlates each
  `Spawning Shape:` line with the windows that follow and picks the
  max-`is_mesh_scan` window among the shape's clean holds (`scan_mesh_raster`
  calls == 16·F).
