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
render exceeded one 62.5 ms window. Colours are strict per config: 🟢 zero
spill, 🟡 under 25%, 🔴 25% or more. Image deltas are raw global-O3
minus shipping bytes from each pair's own image-size reports.

| Effect | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured |
|---|---|--:|--:|--:|--:|--:|--:|---|
| [MindSplatter](shipping/profile_mindsplatter_teensy_2026-08-26.md) / [O3](O3/profile_mindsplatter_teensy_2026-08-26.md) § ● | `msp_draw_particles` | 🟡 67.05 (2)<br>🟢 60.54 (7) | 🟡 66.86 (1)<br>🟢 58.78 (8) | 🟡 39/424 (9.2%)<br>🟢 0/1272 (0.0%) | 🟡 10/318 (3.1%)<br>🟢 0/1394 (0%) | +22,760 B | +20,368 B | ship 2026-08-26 01:52<br>O3 2026-08-26 01:49 |
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-08-26.md) / [O3](O3/profile_islamicstars_teensy_2026-08-26.md) § ● | `is_timeline_step` | 🟡 64.18 (1)<br>🟢 61.87 (22) | 🟡 63.74 (1)<br>🟢 61.36 (22) | 🟡 2/156 (1.3%)<br>🟢 0/3172 (0.0%) | 🟡 1/156 (0.6%)<br>🟢 0/3172 (0%) | +23,240 B | +12,368 B | ship 2026-08-26 02:17<br>O3 2026-08-26 02:14 |
| [Raymarch](shipping/profile_raymarch_teensy_2026-08-26.md) / [O3](O3/profile_raymarch_teensy_2026-08-26.md) ● | `rm_shader_draw` | 🟢 62.22 | 🟢 58.80 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +10,216 B | +6,688 B | ship 2026-08-26 01:37<br>O3 2026-08-26 01:34 |
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-08-26.md) / [O3](O3/profile_meshfeedback_teensy_2026-08-26.md) § ● | `mf_feedback_flush` | 🟢 58.30 (13) | 🟢 58.32 (13) | 🟢 0/6688 (0.0%) | 🟢 0/6688 (0%) | +34,144 B | +21,776 B | ship 2026-08-26 03:31<br>O3 2026-08-26 02:10 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-08-26.md) / [O3](O3/profile_shapeshifter_teensy_2026-08-26.md) § ● | `ss_draw_all` | 🟢 58.13 (10) | 🟢 59.80 (10) | 🟢 0/2448 (0.0%) | 🟢 0/2448 (0%) | +27,272 B | +23,440 B | ship 2026-08-26 03:28<br>O3 2026-08-26 03:11 |
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-08-26.md) / [O3](O3/profile_displacementfield_teensy_2026-08-26.md) ● | `df_timeline_step` | 🟢 57.40 | 🟢 57.48 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +25,512 B | +21,920 B | ship 2026-08-26 01:24<br>O3 2026-08-26 01:20 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-08-26.md) / [O3](O3/profile_gsreactiondiffusion_teensy_2026-08-26.md) ● | `grd_render` | 🟢 55.26 | 🟢 56.16 | 🟢 0/2048 (0.0%) | 🟢 0/2048 (0%) | +13,224 B | +11,280 B | ship 2026-08-26 01:28<br>O3 2026-08-26 01:25 |
| [HyperLattice](shipping/profile_hyperlattice_teensy_2026-08-26.md) / [O3](O3/profile_hyperlattice_teensy_2026-08-26.md) § ● | `hl_shader_draw` | 🟢 51.19 (4) | 🟢 49.92 (4) | 🟢 0/2688 (0.0%) | 🟢 0/2688 (0%) | +10,152 B | +7,728 B | ship 2026-08-26 03:32<br>O3 2026-08-26 02:37 |
| [AshCloud](shipping/profile_ashcloud_teensy_2026-08-26.md) / [O3](O3/profile_ashcloud_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 50.09 | 🔴 79.81 | 🟢 0/1088 (0.0%) | 🔴 544/544 (100%) | +16,608 B | +12,112 B | ship 2026-08-26 02:47<br>O3 2026-08-26 02:45 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-08-26.md) / [O3](O3/profile_ringspin_teensy_2026-08-26.md) ● | `rs_draw_rings` | 🟢 49.42 | 🟢 50.81 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +15,080 B | +13,088 B | ship 2026-08-26 01:41<br>O3 2026-08-26 01:38 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-08-26.md) / [O3](O3/profile_bzreactiondiffusion_teensy_2026-08-26.md) ● | `bz_render` | 🟢 48.91 | 🟢 48.65 | 🟢 0/2048 (0.0%) | 🟢 0/2048 (0%) | +12,760 B | +10,224 B | ship 2026-08-26 01:20<br>O3 2026-08-26 01:16 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-08-26.md) / [O3](O3/profile_hopffibration_teensy_2026-08-26.md) ● | `hf_render_trails` | 🟢 48.50 | 🟢 46.52 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +19,616 B | +18,272 B | ship 2026-08-26 01:30<br>O3 2026-08-26 01:27 |
| [KaleidoscopeStainedGlass](shipping/profile_kaleidoscopestainedglass_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopestainedglass_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 47.91 | 🟢 46.99 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +13,784 B | +11,632 B | ship 2026-08-26 03:23<br>O3 2026-08-26 02:51 |
| [MermaidSkin](shipping/profile_mermaidskin_teensy_2026-08-26.md) / [O3](O3/profile_mermaidskin_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 45.60 | 🟢 54.55 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +16,544 B | +11,952 B | ship 2026-08-26 02:46<br>O3 2026-08-26 02:43 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-08-26.md) / [O3](O3/profile_hankinsolids_teensy_2026-08-26.md) § ● | `hk_timeline_step` | 🟢 45.01 (19) | 🟢 42.04 (19) | 🟢 0/3328 (0.0%) | 🟢 0/3328 (0%) | +17,752 B | +8,752 B | ship 2026-08-26 01:56<br>O3 2026-08-26 01:53 |
| [LatticeMelt](shipping/profile_latticemelt_teensy_2026-08-26.md) / [O3](O3/profile_latticemelt_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 43.73 (3) | 🔴 104.75 (3) | 🟢 0/1728 (0.0%) | 🔴 1824/1824 (100%) | +16,592 B | +11,952 B | ship 2026-08-26 02:42<br>O3 2026-08-26 03:22 |
| [ChromaticLichen](shipping/profile_chromaticlichen_teensy_2026-08-26.md) / [O3](O3/profile_chromaticlichen_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 42.87 | 🟢 61.87 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +16,576 B | +11,952 B | ship 2026-08-26 02:44<br>O3 2026-08-26 02:41 |
| [KaleidoscopeMandala](shipping/profile_kaleidoscopemandala_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopemandala_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 40.29 (3) | 🟢 36.86 (3) | 🟢 0/2368 (0.0%) | 🟢 0/2368 (0%) | +13,656 B | +11,440 B | ship 2026-08-26 03:22<br>O3 2026-08-26 03:17 |
| [KaleidoscopeHexOil](shipping/profile_kaleidoscopehexoil_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopehexoil_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 39.04 (3) | 🟢 38.94 (3) | 🟢 0/2208 (0.0%) | 🟢 0/2208 (0%) | +13,456 B | +10,608 B | ship 2026-08-26 02:52<br>O3 2026-08-26 02:49 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-08-26.md) / [O3](O3/profile_dreamballs_teensy_2026-08-26.md) § ● | `db_timeline_step` | 🟢 38.73 (11) | 🟢 34.32 (11) | 🟢 0/3648 (0.0%) | 🟢 0/3648 (0%) | +26,896 B | +12,352 B | ship 2026-08-26 02:21<br>O3 2026-08-26 02:19 |
| [KaleidoscopeFlowers](shipping/profile_kaleidoscopeflowers_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopeflowers_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 36.69 (4) | 🟢 33.38 (4) | 🟢 0/4128 (0.0%) | 🟢 0/4128 (0%) | +14,576 B | +11,712 B | ship 2026-08-26 03:06<br>O3 2026-08-26 03:03 |
| [KaleidoscopeSmooth](shipping/profile_kaleidoscopesmooth_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopesmooth_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 35.85 (5) | 🟢 32.58 (5) | 🟢 0/4128 (0.0%) | 🟢 0/4128 (0%) | +14,560 B | +11,712 B | ship 2026-08-26 02:59<br>O3 2026-08-26 02:56 |
| [KaleidoscopeHexBright](shipping/profile_kaleidoscopehexbright_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopehexbright_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 35.55 (3) | 🟢 32.41 (3) | 🟢 0/2368 (0.0%) | 🟢 0/2368 (0%) | +14,520 B | +11,712 B | ship 2026-08-26 03:19<br>O3 2026-08-26 03:14 |
| [Comets](shipping/profile_comets_teensy_2026-08-26.md) / [O3](O3/profile_comets_teensy_2026-08-26.md) § ● | `cm_draw_trail` | 🟢 33.91 (13) | 🟢 28.71 (13) | 🟢 0/4128 (0.0%) | 🟢 0/4128 (0%) | +15,432 B | +12,944 B | ship 2026-08-26 02:05<br>O3 2026-08-26 02:02 |
| [KaleidoscopePentBright](shipping/profile_kaleidoscopepentbright_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopepentbright_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 33.37 | 🟢 29.53 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +14,464 B | +11,712 B | ship 2026-08-26 02:49<br>O3 2026-08-26 02:46 |
| [AlienBrain](shipping/profile_alienbrain_teensy_2026-08-26.md) / [O3](O3/profile_alienbrain_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 31.61 (5) | 🟢 27.90 (5) | 🟢 0/4768 (0.0%) | 🟢 0/4768 (0%) | +14,560 B | +11,712 B | ship 2026-08-26 02:27<br>O3 2026-08-26 02:24 |
| [GnomonicStars](shipping/profile_gnomonicstars_teensy_2026-08-26.md) / [O3](O3/profile_gnomonicstars_teensy_2026-08-26.md) ● | `gn_draw_stars` | 🟢 29.64 | 🟢 26.29 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +12,184 B | +11,104 B | ship 2026-08-26 01:25<br>O3 2026-08-26 01:22 |
| [GridSpace](shipping/profile_gridspace_teensy_2026-08-26.md) / [O3](O3/profile_gridspace_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 28.38 | 🟢 25.39 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +14,424 B | +11,712 B | ship 2026-08-26 02:36<br>O3 2026-08-26 02:33 |
| [AlienOcean](shipping/profile_alienocean_teensy_2026-08-26.md) / [O3](O3/profile_alienocean_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 27.79 | 🟢 24.98 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-08-26 02:31<br>O3 2026-08-26 02:28 |
| [KaleidoscopeHexSoft](shipping/profile_kaleidoscopehexsoft_teensy_2026-08-26.md) / [O3](O3/profile_kaleidoscopehexsoft_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 27.30 | 🟢 24.58 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +14,496 B | +11,712 B | ship 2026-08-26 02:29<br>O3 2026-08-26 02:26 |
| [CosmicEyeball](shipping/profile_cosmiceyeball_teensy_2026-08-26.md) / [O3](O3/profile_cosmiceyeball_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 24.01 | 🟢 22.97 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +12,552 B | +10,320 B | ship 2026-08-26 03:07<br>O3 2026-08-26 03:05 |
| [Fishbowl](shipping/profile_fishbowl_teensy_2026-08-26.md) / [O3](O3/profile_fishbowl_teensy_2026-08-26.md) ● | `fish_build_vertices` | 🟢 23.22 | 🟢 21.29 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +27,856 B | +26,400 B | ship 2026-08-26 01:22<br>O3 2026-08-26 01:18 |
| [MobiusGrid](shipping/profile_mobiusgrid_teensy_2026-08-26.md) / [O3](O3/profile_mobiusgrid_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 22.37 (3) | 🟢 21.44 (3) | 🟢 0/2688 (0.0%) | 🟢 0/2688 (0%) | +13,448 B | +10,560 B | ship 2026-08-26 01:34<br>O3 2026-08-26 01:30 |
| [AlienCore](shipping/profile_aliencore_teensy_2026-08-26.md) / [O3](O3/profile_aliencore_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 21.09 | 🟢 20.10 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-08-26 02:32<br>O3 2026-08-26 02:30 |
| [SphericalHarmonics](shipping/profile_sphericalharmonics_teensy_2026-08-26.md) / [O3](O3/profile_sphericalharmonics_teensy_2026-08-26.md) § ● | `sh_rasterize` | 🟢 12.64 (24) | 🟢 11.70 (24) | 🟢 0/3488 (0.0%) | 🟢 0/3488 (0%) | +8,016 B | +6,400 B | ship 2026-08-26 02:00<br>O3 2026-08-26 01:57 |
| [PetalFlow](shipping/profile_petalflow_teensy_2026-08-26.md) / [O3](O3/profile_petalflow_teensy_2026-08-26.md) ● | `pf_draw_rings` | 🟢 11.85 | 🟢 10.77 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +22,272 B | +21,008 B | ship 2026-08-26 01:35<br>O3 2026-08-26 01:32 |
| [Voronoi](shipping/profile_voronoi_teensy_2026-08-26.md) / [O3](O3/profile_voronoi_teensy_2026-08-26.md) ● | `vo_shade` | 🟢 8.96 | 🟢 7.71 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +15,568 B | +12,688 B | ship 2026-08-26 01:43<br>O3 2026-08-26 01:39 |
| [RingShower](shipping/profile_ringshower_teensy_2026-08-26.md) / [O3](O3/profile_ringshower_teensy_2026-08-26.md) ● | `rsh_draw_rings` | 🟢 3.98 | 🟢 3.86 | 🟢 0/1088 (0.0%) | 🟢 0/1088 (0%) | +16,336 B | +15,136 B | ship 2026-08-26 01:39<br>O3 2026-08-26 01:36 |

§ Cyclers carry one aligned line per parser-owned colour bucket, worst first;
the count after a bucket peak is the number of presets in that bucket. Bucket
frames include the following transition, so they are stricter than clean holds.
● Capture refreshed on 2026-08-26. Timestamps are local raw-log mtimes.

## Retired captures

The last device captures for retired effects remain archived as
[Flyby](retired/profile_flyby_teensy_2026-07-27.md) and
[Liquid2D](retired/profile_liquid2d_teensy_2026-07-25.md).

## Memory captures

[Arena high-water measurements](memory/arena_high_water.md) come from a host
probe and are independent of the on-device timing tables.
