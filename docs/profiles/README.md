# On-device effect profiles — Teensy 4.0, segmented mode

On-device timing for the **38 effects in the Phantasm image**, captured on
bench-attached Teensy 4.0 boards with the real segmented POV driver
(`POVSegmented<288, 4, 480>`), DMA LEDs, and live flywheel/DMA ISRs. Each
effect renders one 288×144 image quadrant (about 10,368 pixels); the 62.5 ms
display window makes cadence quantize to 16, 8, 5.3 fps, and below.

## Capture configurations

Timing sets contain one canonical report per effect, with optional supplemental
reports identified by underscore-separated variant suffixes before `_teensy_`.
Each effect/variant pair has one dated report per set. Supplemental reports use
the same title, date, section, roster, and index checks as canonical reports;
they do not satisfy the required canonical shipping coverage.

### Shipping selective-O3

The [`profile` report set](shipping/README.md) uses the shipping `-Os` image
with landed `HS_O3` hot-loop regions active.

### Global-O3 reference

The [`profile_o3` report set](O3/README.md) replaces `-Os` globally with
`-O3 -ffast-math`. It is a single-effect optimization ceiling, not a shippable
full-roster image.

## Paired shipping/O3 captures

Rows rank by shipping spill fraction, then shipping peak render. Both peaks
are worst-frame render, never wall time; spilled is the number of frames whose
render exceeded one 62.5 ms window. Colours are strict per config: 🟢 zero spill, 🟡 under 25% spill, 🔴 at least 25% spill. MeshFeedback and KaleidoscopeStainedGlass use adjacent-source pairs:
`63268c376` O3 versus `20ca3cb48` shipping; their image deltas include source
changes. Other image deltas use each O3 capture's own same-source shipping pair, including
retired shipping captures. They match the currently linked shipping report only
where that report belongs to the same capture pair. GSReactionDiffusion’s O3
column and image deltas predate the October 5 shipping optimization.

| Effect | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---|--:|--:|--:|--:|--:|--:|---|
| [Raymarch](shipping/profile_raymarch_teensy_2026-09-28.md) / [O3](O3/profile_raymarch_teensy_2026-09-20.md) | `rm_shader_draw` | 🟢 54.87 | 🟢 56.04 | 🟢 0/1087 (0.00%) | 🟢 0/1736 (0.00%) | +12,568 B | +7,904 B | ship 2026-09-28 19:09<br>O3 2026-09-20 23:07 |
| [MindSplatter](shipping/profile_mindsplatter_architecture_teensy_2026-10-01.md) / [O3](O3/profile_mindsplatter_architecture_teensy_2026-10-01.md) § ● | `msp_draw_particles` | 🟢 54.457 (8) | 🟢 54.539 (8) | 🟢 0/1736 (0.00%) | 🟢 0/1736 (0.00%) | +21,976 B | +17,200 B | ship 2026-10-01 09:15<br>O3 2026-10-01 09:22 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-09-29.md) / [O3](O3/profile_shapeshifter_teensy_2026-09-25.md) § | `ss_draw_all` | 🟢 53.22 (9) | 🟢 60.05 (9) | 🟢 0/2448 (0.00%) | 🟢 0/2457 (0.00%) | +29,288 B | +24,496 B | ship 2026-09-29 22:38<br>O3 2026-09-25 07:29 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-09-28.md) / [O3](O3/profile_hopffibration_teensy_2026-08-26.md) | `hf_render_trails` | 🟢 51.52 | 🟢 46.52 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +19,616 B | +18,272 B | ship 2026-09-28 19:02<br>O3 2026-08-26 01:27 |
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-09-28.md) / [O3](O3/profile_islamicstars_teensy_2026-09-28.md) § | `is_timeline_step` | 🟢 50.51 (23) | 🟢 49.112 (23) | 🟢 0/3327 (0.00%) | 🟢 0/3336 (0.00%) | +23,904 B | +8,240 B | ship 2026-09-28 18:40<br>O3 2026-09-28 16:50 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-09-28.md) / [O3](O3/profile_bzreactiondiffusion_teensy_2026-08-26.md) | `bz_render` | 🟢 48.84 | 🟢 48.65 | 🟢 0/2047 (0.00%) | 🟢 0/2048 (0%) | +12,760 B | +10,224 B | ship 2026-09-28 19:12<br>O3 2026-08-26 01:16 |
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-09-29.md) / [O3](O3/profile_meshfeedback_teensy_2026-08-26.md) § | `mf_feedback_flush` | 🟢 48.61 (12) | 🟢 58.32 (13) | 🟢 0/6688 (0.00%) | 🟢 0/6688 (0%) | +34,144 B | +21,776 B | ship 2026-09-29 10:03<br>O3 2026-08-26 02:10 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-09-28.md) / [O3](O3/profile_ringspin_teensy_2026-09-24.md) | `rs_draw_rings` | 🟢 47.04 | 🟢 50.750 | 🟢 0/1087 (0.00%) | 🟢 0/1087 (0.0%) | +15,400 B | +12,880 B | ship 2026-09-28 19:13<br>O3 2026-09-24 20:18 |
| [AshCloud](shipping/profile_ashcloud_teensy_2026-09-28.md) / [O3](O3/profile_ashcloud_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 44.57 | 🔴 79.81 | 🟢 0/1087 (0.00%) | 🔴 544/544 (100%) | +16,608 B | +12,112 B | ship 2026-09-28 17:33<br>O3 2026-08-26 02:45 |
| [KaleidoscopeStainedGlass](shipping/profile_kaleidoscopestainedglass_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopestainedglass_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 43.50 | 🟢 46.99 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +13,784 B | +11,632 B | ship 2026-09-28 17:35<br>O3 2026-08-26 02:51 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-09-28.md) / [O3](O3/profile_dreamballs_teensy_2026-08-26.md) § | `db_timeline_step` | 🟢 42.66 (10) | 🟢 34.32 (11) | 🟢 0/3647 (0.00%) | 🟢 0/3648 (0%) | +26,896 B | +12,352 B | ship 2026-09-28 18:45<br>O3 2026-08-26 02:19 |
| [MermaidSkin](shipping/profile_mermaidskin_teensy_2026-09-28.md) / [O3](O3/profile_mermaidskin_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 39.57 | 🟢 54.55 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +16,544 B | +11,952 B | ship 2026-09-28 18:49<br>O3 2026-08-26 02:43 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-10-05.md) / [O3](O3/profile_gsreactiondiffusion_teensy_2026-09-30.md) ● | `grd_shader_draw` | 🟢 39.271 | 🔴 279.686 | 🟢 0/2032 (0.00%) | 🔴 538/538 (100.00%) | +16,024 B | +10,832 B | ship 2026-10-05 00:38<br>O3 2026-09-30 19:15 |
| [KaleidoscopeHexOil](shipping/profile_kaleidoscopehexoil_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopehexoil_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 38.51 (2) | 🟢 38.94 (3) | 🟢 0/2207 (0.00%) | 🟢 0/2208 (0%) | +13,456 B | +10,608 B | ship 2026-09-28 18:54<br>O3 2026-08-26 02:49 |
| [HyperLattice](shipping/profile_hyperlattice_architecture_teensy_2026-10-01.md) / [O3](O3/profile_hyperlattice_architecture_teensy_2026-10-01.md) § ● | `hl_shader_draw` | 🟢 36.839 (3) | 🟢 36.975 (3) | 🟢 0/2696 (0.00%) | 🟢 0/2696 (0.00%) | +9,616 B | +4,912 B | ship 2026-10-01 09:19<br>O3 2026-10-01 09:26 |
| [LatticeMelt](shipping/profile_latticemelt_architecture_teensy_2026-10-01.md) / [O3](O3/profile_latticemelt_architecture_teensy_2026-10-01.md) § ● | `fx_shader_draw` | 🟢 36.583 (2) | 🟢 38.705 (2) | 🟢 0/1736 (0.00%) | 🟢 0/1736 (0.00%) | +25,752 B | +19,472 B | ship 2026-10-01 09:38<br>O3 2026-10-01 09:49 |
| [ChromaticLichen](shipping/profile_chromaticlichen_architecture_teensy_2026-10-01.md) / [O3](O3/profile_chromaticlichen_architecture_teensy_2026-10-01.md) ● | `fx_shader_draw` | 🟢 35.773 | 🟢 35.218 | 🟢 0/1096 (0.00%) | 🟢 0/1096 (0.00%) | +15,264 B | +8,672 B | ship 2026-10-01 09:40<br>O3 2026-10-01 09:51 |
| [KaleidoscopeMandala](shipping/profile_kaleidoscopemandala_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopemandala_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 35.40 (2) | 🟢 36.86 (3) | 🟢 0/2367 (0.00%) | 🟢 0/2368 (0%) | +13,656 B | +11,440 B | ship 2026-09-28 18:39<br>O3 2026-08-26 03:17 |
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-09-29.md) / [O3](O3/profile_displacementfield_teensy_2026-09-19.md) | `df_timeline_step` | 🟢 35.27 | 🟢 57.45 | 🟢 0/2368 (0.00%) | 🟢 0/1088 (0%) | +25,856 B | +22,144 B | ship 2026-09-29 09:09<br>O3 2026-09-19 22:19 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-09-28.md) / [O3](O3/profile_hankinsolids_teensy_2026-09-28.md) § | `hk_timeline_step` | 🟢 34.81 (19) | 🟢 32.865 (18) | 🟢 0/4447 (0.00%) | 🟢 0/4137 (0.00%) | +17,088 B | +2,144 B | ship 2026-09-28 19:27<br>O3 2026-09-28 16:45 |
| [KaleidoscopeFlowers](shipping/profile_kaleidoscopeflowers_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopeflowers_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 32.52 (3) | 🟢 33.38 (4) | 🟢 0/4127 (0.00%) | 🟢 0/4128 (0%) | +14,576 B | +11,712 B | ship 2026-09-28 19:07<br>O3 2026-08-26 03:03 |
| [KaleidoscopeSmooth](shipping/profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md) / [O3](O3/profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md) § ● | `fx_shader_draw` | 🟢 32.230 (4) | 🟢 30.314 (4) | 🟢 0/4136 (0.00%) | 🟢 0/4137 (0.00%) | +15,048 B | +8,624 B | ship 2026-10-01 11:05<br>O3 2026-10-01 11:19 |
| [KaleidoscopeHexBright](shipping/profile_kaleidoscopehexbright_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopehexbright_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 31.64 (2) | 🟢 32.41 (3) | 🟢 0/2367 (0.00%) | 🟢 0/2368 (0%) | +14,520 B | +11,712 B | ship 2026-09-28 19:02<br>O3 2026-08-26 03:14 |
| [Comets](shipping/profile_comets_teensy_2026-09-28.md) / [O3](O3/profile_comets_teensy_2026-08-26.md) § | `cm_draw_trail` | 🟢 30.66 (12) | 🟢 28.71 (13) | 🟢 0/4127 (0.00%) | 🟢 0/4128 (0%) | +15,432 B | +12,944 B | ship 2026-09-28 18:28<br>O3 2026-08-26 02:02 |
| [KaleidoscopePentBright](shipping/profile_kaleidoscopepentbright_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopepentbright_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 27.24 | 🟢 29.53 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +14,464 B | +11,712 B | ship 2026-09-28 18:51<br>O3 2026-08-26 02:46 |
| [AlienBrain](shipping/profile_alienbrain_teensy_2026-09-28.md) / [O3](O3/profile_alienbrain_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 25.89 (4) | 🟢 27.90 (5) | 🟢 0/4767 (0.00%) | 🟢 0/4768 (0%) | +14,560 B | +11,712 B | ship 2026-09-28 18:30<br>O3 2026-08-26 02:24 |
| [GridSpace](shipping/profile_gridspace_teensy_2026-09-28.md) / [O3](O3/profile_gridspace_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 25.31 | 🟢 25.39 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +14,424 B | +11,712 B | ship 2026-09-28 17:38<br>O3 2026-08-26 02:33 |
| [KaleidoscopeHexSoft](shipping/profile_kaleidoscopehexsoft_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopehexsoft_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 24.52 | 🟢 24.58 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +14,496 B | +11,712 B | ship 2026-09-28 18:32<br>O3 2026-08-26 02:26 |
| [CosmicEyeball](shipping/profile_cosmiceyeball_teensy_2026-09-28.md) / [O3](O3/profile_cosmiceyeball_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 24.11 | 🟢 22.97 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +12,552 B | +10,320 B | ship 2026-09-28 19:09<br>O3 2026-08-26 03:05 |
| [Fishbowl](shipping/profile_fishbowl_teensy_2026-09-28.md) / [O3](O3/profile_fishbowl_teensy_2026-08-26.md) | `fish_build_vertices` | 🟢 23.98 | 🟢 21.29 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +27,856 B | +26,400 B | ship 2026-09-28 19:14<br>O3 2026-08-26 01:18 |
| [AlienOcean](shipping/profile_alienocean_teensy_2026-09-28.md) / [O3](O3/profile_alienocean_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 22.56 | 🟢 24.98 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-09-28 18:34<br>O3 2026-08-26 02:28 |
| [MobiusGrid](shipping/profile_mobiusgrid_teensy_2026-09-28.md) / [O3](O3/profile_mobiusgrid_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 22.05 (2) | 🟢 21.44 (3) | 🟢 0/2687 (0.00%) | 🟢 0/2688 (0%) | +13,448 B | +10,560 B | ship 2026-09-28 19:05<br>O3 2026-08-26 01:30 |
| [GnomonicStars](shipping/profile_gnomonicstars_teensy_2026-09-28.md) / [O3](O3/profile_gnomonicstars_teensy_2026-08-26.md) | `gn_draw_stars` | 🟢 21.86 | 🟢 26.29 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +12,184 B | +11,104 B | ship 2026-09-28 19:18<br>O3 2026-08-26 01:22 |
| [AlienCore](shipping/profile_aliencore_teensy_2026-09-28.md) / [O3](O3/profile_aliencore_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 18.03 | 🟢 20.10 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-09-28 18:36<br>O3 2026-08-26 02:30 |
| [PetalFlow](shipping/profile_petalflow_teensy_2026-09-28.md) / [O3](O3/profile_petalflow_teensy_2026-08-26.md) | `pf_draw_rings` | 🟢 13.61 | 🟢 10.77 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +22,272 B | +21,008 B | ship 2026-09-28 19:07<br>O3 2026-08-26 01:32 |
| [SphericalHarmonics](shipping/profile_sphericalharmonics_teensy_2026-09-28.md) / [O3](O3/profile_sphericalharmonics_teensy_2026-08-26.md) § | `sh_rasterize` | 🟢 12.89 (24) | 🟢 11.70 (24) | 🟢 0/3487 (0.00%) | 🟢 0/3488 (0%) | +8,016 B | +6,400 B | ship 2026-09-28 19:00<br>O3 2026-08-26 01:57 |
| [Voronoi](shipping/profile_voronoi_teensy_2026-09-28.md) / [O3](O3/profile_voronoi_teensy_2026-08-26.md) | `vo_shade` | 🟢 8.51 | 🟢 7.71 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +15,568 B | +12,688 B | ship 2026-09-28 19:15<br>O3 2026-08-26 01:39 |
| [RingShower](shipping/profile_ringshower_teensy_2026-09-28.md) / [O3](O3/profile_ringshower_teensy_2026-08-26.md) | `rsh_draw_rings` | 🟢 4.33 | 🟢 3.86 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +16,336 B | +15,136 B | ship 2026-09-28 19:11<br>O3 2026-08-26 01:36 |

● GSReactionDiffusion shipping refreshed 2026-10-05. Architecture captures refreshed 2026-10-01: MindSplatter, HyperLattice, LatticeMelt, ChromaticLichen and KaleidoscopeSmooth. Setup frame 1 is excluded. Cycler buckets include every live owner-attributed transition frame; each report lists per-preset peaks, fractions and the worst preset. Captured times are local raw-log mtimes.

§ Cyclers carry one aligned line per parser-owned colour bucket, worst first;
(N) counts parser-owned preset, shape, or mode entries in that colour bucket.
Bucket frames include the following transition, so they are stricter than clean holds.
Captured timestamps are local raw-log mtimes.

ShapeShifter’s O3 column predates entry 1’s increase from 208 to 288 contours
and the chord walk; its workload differs from the linked shipping capture.

## Memory captures

[Arena high-water measurements](memory/arena_high_water.md) come from a host
probe and are independent of the on-device timing tables.

## Supplemental experimental presets

Historical fixed opt-in captures; these do not change the normal firmware
roster ranking. Triangular has since been removed.

| Preset | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| HyperLattice Triangular (supporting artifact removed) | `hl_shader_draw` | 🔴 134.620 | 🔴 141.587 | 🔴 441/441 (100%) | 🔴 366/366 (100%) | +10,440 B | +8,208 B | ship 2026-09-27 21:59<br>O3 2026-09-27 22:01 |
| HyperLattice Octet 3D | `hl_shader_draw` | 🔴 96.869 | 🔴 95.710 | 🔴 548/548 (100%) | 🔴 548/548 (100%) | +11,688 B | +8,128 B | ship 2026-09-27 22:50<br>O3 2026-09-27 22:53 |

Octet 3D refreshed on 2026-09-27 after the single-owner strut correction.
These fixed-preset captures retain the original oscillating camera path.
The earlier Octet measurements (supporting artifact removed)
record 3D shipping peak 103.733 ms (548/548 spills) before that correction.
The 4D capture remains historical: 549.828 ms peak, 120/120 spills,
captured 2026-09-27 at 22:18; no 4D global-O3 capture.

Cubic Wide Flight (supporting artifact removed)
has a separate fixed-preset shipping check: peak 50.804 ms, 0/616 spills,
captured 2026-09-27 22:25.

## HyperLattice experimental cycle

Opt-in image (`HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1`), full nine-preset cycle
with family segues (lerp within a pattern and view, fade through black across), 345 s on COM3; setup frame 1 excluded. A bucket opens when
its preset's parameters are adopted: the morph into it, or the dark midpoint
of a fade. Supersedes the octet and shell supplements
below. [Report](shipping/profile_hyperlattice_teensy_2026-09-29.md#per-preset-table).

| Preset | Ship peak ms | Ship spilled | Captured |
|---|--:|--:|---|
| 5: Octet 4D Flight | 🟢 55.87 | 🟢 0/559 (0.00%) | 2026-09-29 16:11 |
| 8: Shell 4D Flight | 🟢 36.76 | 🟢 0/559 (0.00%) | 2026-09-29 16:11 |
| 2: Hypercube Flight | 🟢 36.70 | 🟢 0/559 (0.00%) | 2026-09-29 16:11 |
| 3: Octet Flight | 🟢 32.23 | 🟢 0/439 (0.00%) | 2026-09-29 16:11 |
| 4: Octet Wide Flight | 🟢 31.75 | 🟢 0/679 (0.00%) | 2026-09-29 16:11 |
| 7: Shell Close Flight | 🟢 31.22 | 🟢 0/679 (0.00%) | 2026-09-29 16:11 |
| 1: Cubic Wide Flight | 🟢 29.64 | 🟢 0/817 (0.00%) | 2026-09-29 16:11 |
| 6: Shell Flight | 🟢 28.03 | 🟢 0/439 (0.00%) | 2026-09-29 16:11 |
| 0: Cubic Flight | 🟢 26.41 | 🟢 0/758 (0.00%) | 2026-09-29 16:11 |

## Octet optimization supplement

Historical captures: [2026-09-27 shipping](shipping/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-27 global-O3](O3/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-28 global-O3 preset 3](O3/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md),
and [2026-09-28 global-O3 preset 5](O3/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md).

| Preset | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---:|---:|---:|---:|---:|---:|---|
| [Octet 4D (index 4 at capture time)](shipping/profile_hyperlattice_octet4d_teensy_2026-09-28.md) | 🔴 143.318 | — | 🔴 233/233 (100.00%) | — | — | — | ship 2026-09-28 09:53 |
| [3: Octet 3D](shipping/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) | 🟢 45.688 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:56 |
| [Octet wide (index 5 at capture time)](shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) | 🟢 45.242 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:53 |

Updated 2026-09-28 09:53–09:56 after the Octet 58 ms Tiers 1–2 (`1bb522211`). Opt-in held presets, startup excluded; 4D on COM3, 3D on COM4. Earlier optimization and 4D evidence (supporting artifact removed). No global-O3 twins or size deltas were captured for this code.

## Shell flight supplement

Sphere-only renderer; two authored fixed presets on COM4, startup excluded. [Shipping report](shipping/profile_hyperlattice_shell_flight_teensy_2026-09-29.md) and [global-O3 report](O3/profile_hyperlattice_shell_flight_teensy_2026-09-29.md).

| Preset | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured local |
|---|--:|--:|--:|--:|--:|--:|---|
| Shell Flight (index 5 at capture time) | 🔴 123.553 | 🔴 117.363 | 🔴 228/228 (100.00%) | 🔴 228/228 (100.00%) | +13,352 B | +4,736 B | ship 2026-09-29 08:03<br>O3 2026-09-29 08:04 |
| Shell Close Flight (index 6 at capture time) | 🔴 125.186 | 🔴 121.337 | 🔴 228/228 (100.00%) | 🔴 228/228 (100.00%) | +13,352 B | +4,736 B | ship 2026-09-29 08:06<br>O3 2026-09-29 08:07 |

## Earlier captures

Preserved canonical snapshots; current architecture rankings use the dated variants above. Each variant records its exact checkpoint source and paired baseline.

- MindSplatter: [shipping 2026-09-29](shipping/profile_mindsplatter_teensy_2026-09-29.md), [O3 2026-08-26](O3/profile_mindsplatter_teensy_2026-08-26.md).
- HyperLattice: [shipping 2026-09-29](shipping/profile_hyperlattice_teensy_2026-09-29.md), [O3 2026-09-27](O3/profile_hyperlattice_teensy_2026-09-27.md).
- LatticeMelt: [shipping 2026-09-28](shipping/profile_latticemelt_teensy_2026-09-28.md), [O3 2026-08-26](O3/profile_latticemelt_teensy_2026-08-26.md).
- ChromaticLichen: [shipping 2026-09-28](shipping/profile_chromaticlichen_teensy_2026-09-28.md), [O3 2026-08-26](O3/profile_chromaticlichen_teensy_2026-08-26.md).
- KaleidoscopeSmooth: [shipping 2026-09-28](shipping/profile_kaleidoscopesmooth_teensy_2026-09-28.md), [O3 2026-08-26](O3/profile_kaleidoscopesmooth_teensy_2026-08-26.md).

## Artifact retention and source reachability

Only memory, shipping and global-O3 reports and their indexes are checked in.
Raw captures, source patches, build logs, ELF/map files and experimental campaign
artifacts are local outputs; they are no longer retained in this archive.
The profile wrapper writes these outputs under ignored `build/prof/`.
Historical hashes identify captures but do not guarantee reproducible source.

Capture `d25dd85dee17` maps to landed `43858e5fb`; only test headers differ,
with no firmware-relevant tree difference. Capture `0156d0d7490355` maps to
landed `3fc710350`. Capture `7baf3cc4307` maps to landed `daf63d812` on a
different base (319 differing files); it cannot be treated as the same tree.
Capture `0df961b818ae` is unrecoverable. These historical orphan hashes are
not reachable from the published branch. The KaleidoscopeSmooth shipping
capture additionally used an unretained source patch, so its recorded base
and patch hash cannot reconstruct the captured source from a fresh clone.

For the 2026-08-26 O3 cycler captures, the parser counted the unmarked
startup run as an extra entry 0. Their parenthesized entry counts are the
authored preset count plus one; raw logs are no longer available to rederive
these buckets.
