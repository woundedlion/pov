# On-device effect profiles — Teensy 4.0, segmented mode

On-device timing for the **38 effects in the Phantasm image**, captured on
bench-attached Teensy 4.0 boards with the real segmented POV driver
(`POVSegmented<288, 4, 480>`), DMA LEDs, and live flywheel/DMA ISRs. Each
effect renders one 288×144 image quadrant (about 10,368 pixels); the 62.5 ms
display window makes cadence quantize to 16, 8, 5.3 fps, and below.

## Capture configurations

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
render exceeded one 62.5 ms window. Colours are strict per config: 🟢 zero spill, 🟡 under 25% spill, 🔴 at least 25% spill. Image deltas are raw global-O3
minus shipping bytes from each pair's own image-size reports.

| Effect | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---|--:|--:|--:|--:|--:|--:|---|
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-09-28.md) / [O3](O3/profile_meshfeedback_teensy_2026-08-26.md) § | `mf_feedback_flush` | 🟢 61.87 (12) | 🟢 58.32 (13) | 🟢 0/6687 (0.00%) | 🟢 0/6688 (0%) | +34,144 B | +21,776 B | ship 2026-09-28 18:36<br>O3 2026-08-26 02:10 |
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-09-28.md) / [O3](O3/profile_displacementfield_teensy_2026-09-19.md) | `df_timeline_step` | 🟢 59.42 | 🟢 57.45 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +25,856 B | +22,144 B | ship 2026-09-28 19:16<br>O3 2026-09-19 22:19 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-09-28.md) / [O3](O3/profile_shapeshifter_teensy_2026-09-25.md) § | `ss_draw_all` | 🟢 58.24 (9) | 🟢 60.05 (9) | 🟢 0/2447 (0.00%) | 🟢 0/2457 (0.00%) | +29,288 B | +24,496 B | ship 2026-09-28 18:48<br>O3 2026-09-25 07:29 |
| [HyperLattice](shipping/profile_hyperlattice_teensy_2026-09-28.md) / [O3](O3/profile_hyperlattice_teensy_2026-09-27.md) § | `hl_shader_draw` | 🟢 56.12 (3) | 🟢 54.817 (2) | 🟢 0/2687 (0.00%) | 🟢 0/1576 (0.00%) | +9,736 B | +8,240 B | ship 2026-09-28 18:43<br>O3 2026-09-27 01:05 |
| [MindSplatter](shipping/profile_mindsplatter_teensy_2026-09-28.md) / [O3](O3/profile_mindsplatter_teensy_2026-08-26.md) § | `msp_draw_particles` | 🟢 56.11 (8) | 🟢 52.79 (9) | 🟢 0/1727 (0.00%) | 🟢 0/1728 (0%) | +22,800 B | +20,400 B | ship 2026-09-28 18:51<br>O3 2026-08-26 07:45 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-09-28.md) / [O3](O3/profile_gsreactiondiffusion_teensy_2026-08-26.md) | `grd_render` | 🟢 55.53 | 🟢 56.16 | 🟢 0/2047 (0.00%) | 🟢 0/2048 (0%) | +13,224 B | +11,280 B | ship 2026-09-28 19:21<br>O3 2026-08-26 01:25 |
| [Raymarch](shipping/profile_raymarch_teensy_2026-09-28.md) / [O3](O3/profile_raymarch_teensy_2026-09-20.md) | `rm_shader_draw` | 🟢 54.87 | 🟢 56.04 | 🟢 0/1087 (0.00%) | 🟢 0/1736 (0.00%) | +12,568 B | +7,904 B | ship 2026-09-28 19:09<br>O3 2026-09-20 23:07 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-09-28.md) / [O3](O3/profile_hopffibration_teensy_2026-08-26.md) | `hf_render_trails` | 🟢 51.52 | 🟢 46.52 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +19,616 B | +18,272 B | ship 2026-09-28 19:02<br>O3 2026-08-26 01:27 |
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-09-28.md) / [O3](O3/profile_islamicstars_teensy_2026-09-28.md) § | `is_timeline_step` | 🟢 50.51 (23) | 🟢 49.112 (23) | 🟢 0/3327 (0.00%) | 🟢 0/3336 (0.00%) | +23,904 B | +8,240 B | ship 2026-09-28 18:40<br>O3 2026-09-28 16:50 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-09-28.md) / [O3](O3/profile_bzreactiondiffusion_teensy_2026-08-26.md) | `bz_render` | 🟢 48.84 | 🟢 48.65 | 🟢 0/2047 (0.00%) | 🟢 0/2048 (0%) | +12,760 B | +10,224 B | ship 2026-09-28 19:12<br>O3 2026-08-26 01:16 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-09-28.md) / [O3](O3/profile_ringspin_teensy_2026-09-24.md) | `rs_draw_rings` | 🟢 47.04 | 🟢 50.750 | 🟢 0/1087 (0.00%) | 🟢 0/1087 (0.0%) | +15,400 B | +12,880 B | ship 2026-09-28 19:13<br>O3 2026-09-24 20:18 |
| [AshCloud](shipping/profile_ashcloud_teensy_2026-09-28.md) / [O3](O3/profile_ashcloud_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 44.57 | 🔴 79.81 | 🟢 0/1087 (0.00%) | 🔴 544/544 (100%) | +16,608 B | +12,112 B | ship 2026-09-28 17:33<br>O3 2026-08-26 02:45 |
| [KaleidoscopeStainedGlass](shipping/profile_kaleidoscopestainedglass_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopestainedglass_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 43.50 | 🟢 46.99 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +13,784 B | +11,632 B | ship 2026-09-28 17:35<br>O3 2026-08-26 02:51 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-09-28.md) / [O3](O3/profile_dreamballs_teensy_2026-08-26.md) § | `db_timeline_step` | 🟢 42.66 (10) | 🟢 34.32 (11) | 🟢 0/3647 (0.00%) | 🟢 0/3648 (0%) | +26,896 B | +12,352 B | ship 2026-09-28 18:45<br>O3 2026-08-26 02:19 |
| [MermaidSkin](shipping/profile_mermaidskin_teensy_2026-09-28.md) / [O3](O3/profile_mermaidskin_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 39.57 | 🟢 54.55 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +16,544 B | +11,952 B | ship 2026-09-28 18:49<br>O3 2026-08-26 02:43 |
| [KaleidoscopeHexOil](shipping/profile_kaleidoscopehexoil_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopehexoil_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 38.51 (2) | 🟢 38.94 (3) | 🟢 0/2207 (0.00%) | 🟢 0/2208 (0%) | +13,456 B | +10,608 B | ship 2026-09-28 18:54<br>O3 2026-08-26 02:49 |
| [LatticeMelt](shipping/profile_latticemelt_teensy_2026-09-28.md) / [O3](O3/profile_latticemelt_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 37.18 (2) | 🔴 104.75 (3) | 🟢 0/1727 (0.00%) | 🔴 1824/1824 (100%) | +16,592 B | +11,952 B | ship 2026-09-28 18:45<br>O3 2026-08-26 03:22 |
| [ChromaticLichen](shipping/profile_chromaticlichen_teensy_2026-09-28.md) / [O3](O3/profile_chromaticlichen_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 36.63 | 🟢 61.87 | 🟢 0/1087 (0.00%) | 🟢 0/1088 (0%) | +16,576 B | +11,952 B | ship 2026-09-28 18:47<br>O3 2026-08-26 02:41 |
| [KaleidoscopeMandala](shipping/profile_kaleidoscopemandala_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopemandala_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 35.40 (2) | 🟢 36.86 (3) | 🟢 0/2367 (0.00%) | 🟢 0/2368 (0%) | +13,656 B | +11,440 B | ship 2026-09-28 18:39<br>O3 2026-08-26 03:17 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-09-28.md) / [O3](O3/profile_hankinsolids_teensy_2026-09-28.md) § | `hk_timeline_step` | 🟢 34.81 (19) | 🟢 32.865 (18) | 🟢 0/4447 (0.00%) | 🟢 0/4137 (0.00%) | +17,088 B | +2,144 B | ship 2026-09-28 19:27<br>O3 2026-09-28 16:45 |
| [KaleidoscopeSmooth](shipping/profile_kaleidoscopesmooth_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopesmooth_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 32.93 (4) | 🟢 32.58 (5) | 🟢 0/4127 (0.00%) | 🟢 0/4128 (0%) | +14,560 B | +11,712 B | ship 2026-09-28 18:59<br>O3 2026-08-26 02:56 |
| [KaleidoscopeFlowers](shipping/profile_kaleidoscopeflowers_teensy_2026-09-28.md) / [O3](O3/profile_kaleidoscopeflowers_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 32.52 (3) | 🟢 33.38 (4) | 🟢 0/4127 (0.00%) | 🟢 0/4128 (0%) | +14,576 B | +11,712 B | ship 2026-09-28 19:07<br>O3 2026-08-26 03:03 |
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

Shipping columns refreshed on 2026-09-28 after the finding-10 placement change `97eb0bf78`; setup frame 1 is excluded. The O3 columns and the FLASH/ITCM deltas are each pair's original global-O3 capture, so they predate this shipping image.

§ Cyclers carry one aligned line per parser-owned colour bucket, worst first;
(N) counts parser-owned preset, shape, or mode entries in that colour bucket.
Bucket frames include the following transition, so they are stricter than clean holds.
Captured timestamps are local raw-log mtimes.

## Memory captures

[Arena high-water measurements](memory/arena_high_water.md) come from a host
probe and are independent of the on-device timing tables.


● HyperLattice refreshed on 2026-09-27: both presets and transitions; setup frame 1 excluded.

## Supplemental experimental presets

Historical fixed opt-in captures; these do not change the normal firmware
roster ranking. Triangular was subsequently selected for removal.

| Preset | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| [HyperLattice Triangular](hyperlattice_triangular_2026-09-27.md) | `hl_shader_draw` | 🔴 134.620 | 🔴 141.587 | 🔴 441/441 (100%) | 🔴 366/366 (100%) | +10,440 B | +8,208 B | ship 2026-09-27 21:59<br>O3 2026-09-27 22:01 |
| [HyperLattice Octet 3D ●](shipping/profile_hyperlattice_teensy_2026-09-28.md#supplemental-octet-3d-single-owner-correction) | `hl_shader_draw` | 🔴 96.869 | 🔴 95.710 | 🔴 548/548 (100%) | 🔴 548/548 (100%) | +11,688 B | +8,128 B | ship 2026-09-27 22:50<br>O3 2026-09-27 22:53 |

● Octet 3D refreshed on 2026-09-27 after the single-owner strut correction.
These fixed-preset captures retain the original oscillating camera path.
The [earlier Octet measurements](hyperlattice_experimental_presets_2026-09-27.md#octet-validation-and-device-measurements)
record 3D shipping peak 103.733 ms (548/548 spills) before that correction.
The 4D capture remains historical: 549.828 ms peak, 120/120 spills,
captured 2026-09-27 at 22:18; no 4D global-O3 capture.

[Cubic Wide Flight](hyperlattice_experimental_presets_2026-09-27.md#cubic-wide-flight)
has a separate fixed-preset shipping check: peak 50.804 ms, 0/616 spills,
captured 2026-09-27 22:25. The earlier HyperLattice cycle rows above cover
its two original presets; they predate this third default preset.

## Octet optimization supplement

Historical captures: [2026-09-27 shipping](shipping/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-27 global-O3](O3/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-28 global-O3 preset 3](O3/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md),
and [2026-09-28 global-O3 preset 5](O3/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md).

Timing sets contain one canonical report per effect, with optional supplemental
reports identified by underscore-separated variant suffixes before `_teensy_`.
Each effect/variant pair has one dated report per set. Supplemental reports use
the same title, date, section, roster, and index checks as canonical reports;
they do not satisfy the required canonical shipping coverage.

| Preset | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---:|---:|---:|---:|---:|---:|---|
| [4: Octet 4D](shipping/profile_hyperlattice_octet4d_teensy_2026-09-28.md) ● | 🔴 143.318 | — | 🔴 233/233 (100.00%) | — | — | — | ship 2026-09-28 09:53 |
| [3: Octet 3D](shipping/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) ● | 🟢 45.688 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:56 |
| [5: Octet wide](shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) ● | 🟢 45.242 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:53 |

● Updated 2026-09-28 09:53–09:56 after the Octet 58 ms Tiers 1–2 (`1bb522211`). Opt-in held presets, startup excluded; 4D on COM3, 3D on COM4. [Earlier optimization and 4D evidence](hyperlattice_octet_optimization_2026-09-27.md). No global-O3 twins or size deltas were captured for this code.
