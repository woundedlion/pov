# Global-O3 reference profiles

Ranked captures covering all **38 effects** in the 288×144 Phantasm
image. These are single-effect `profile_o3` compiler-reference images built
with global `-O3 -ffast-math`; they are not a shippable full-roster image.

Peak is worst-frame render time (never wall). Spilled is the fraction of frames
whose render exceeded one 62.5 ms display interval. Rows rank by spill fraction, then peak render.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [LatticeMelt](profile_latticemelt_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🔴 104.75 (3) | 🔴 1824/1824 (100%) | 2026-08-26 03:22 |
| [AshCloud](profile_ashcloud_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🔴 79.81 | 🔴 544/544 (100%) | 2026-08-26 02:45 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 61.87 | 🟢 0/1088 (0%) | 2026-08-26 02:41 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-08-26.md)§ ● | `ss_draw_all` | 🟢 59.80 (10) | 🟢 0/2448 (0%) | 2026-08-26 03:11 |
| [Raymarch](profile_raymarch_teensy_2026-08-26.md) ● | `rm_shader_draw` | 🟢 58.80 | 🟢 0/1088 (0%) | 2026-08-26 01:34 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-08-26.md)§ ● | `mf_feedback_flush` | 🟢 58.32 (13) | 🟢 0/6688 (0%) | 2026-08-26 02:10 |
| [DisplacementField](profile_displacementfield_teensy_2026-08-26.md) ● | `df_timeline_step` | 🟢 57.48 | 🟢 0/1088 (0%) | 2026-08-26 01:20 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-08-26.md) ● | `grd_render` | 🟢 56.16 | 🟢 0/2048 (0%) | 2026-08-26 01:25 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 54.55 | 🟢 0/1088 (0%) | 2026-08-26 02:43 |
| [MindSplatter](profile_mindsplatter_teensy_2026-08-26.md)§ ● | `msp_draw_particles` | 🟢 52.79 (9) | 🟢 0/1728 (0%) | 2026-08-26 07:45 |
| [RingSpin](profile_ringspin_teensy_2026-08-26.md) ● | `rs_draw_rings` | 🟢 50.81 | 🟢 0/1088 (0%) | 2026-08-26 01:38 |
| [IslamicStars](profile_islamicstars_teensy_2026-08-26.md)§ ● | `is_timeline_step` | 🟢 50.80 (23) | 🟢 0/3328 (0%) | 2026-08-26 07:53 |
| [HyperLattice](profile_hyperlattice_teensy_2026-08-26.md)§ ● | `hl_shader_draw` | 🟢 49.92 (4) | 🟢 0/2688 (0%) | 2026-08-26 02:37 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-08-26.md) ● | `bz_render` | 🟢 48.65 | 🟢 0/2048 (0%) | 2026-08-26 01:16 |
| [KaleidoscopeStainedGlass](profile_kaleidoscopestainedglass_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 46.99 | 🟢 0/1088 (0%) | 2026-08-26 02:51 |
| [HopfFibration](profile_hopffibration_teensy_2026-08-26.md) ● | `hf_render_trails` | 🟢 46.52 | 🟢 0/1088 (0%) | 2026-08-26 01:27 |
| [HankinSolids](profile_hankinsolids_teensy_2026-08-26.md)§ ● | `hk_timeline_step` | 🟢 42.04 (19) | 🟢 0/3328 (0%) | 2026-08-26 01:53 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🟢 38.94 (3) | 🟢 0/2208 (0%) | 2026-08-26 02:49 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🟢 36.86 (3) | 🟢 0/2368 (0%) | 2026-08-26 03:17 |
| [DreamBalls](profile_dreamballs_teensy_2026-08-26.md)§ ● | `db_timeline_step` | 🟢 34.32 (11) | 🟢 0/3648 (0%) | 2026-08-26 02:19 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🟢 33.38 (4) | 🟢 0/4128 (0%) | 2026-08-26 03:03 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🟢 32.58 (5) | 🟢 0/4128 (0%) | 2026-08-26 02:56 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🟢 32.41 (3) | 🟢 0/2368 (0%) | 2026-08-26 03:14 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 29.53 | 🟢 0/1088 (0%) | 2026-08-26 02:46 |
| [Comets](profile_comets_teensy_2026-08-26.md)§ ● | `cm_draw_trail` | 🟢 28.71 (13) | 🟢 0/4128 (0%) | 2026-08-26 02:02 |
| [AlienBrain](profile_alienbrain_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🟢 27.90 (5) | 🟢 0/4768 (0%) | 2026-08-26 02:24 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-08-26.md) ● | `gn_draw_stars` | 🟢 26.29 | 🟢 0/1088 (0%) | 2026-08-26 01:22 |
| [GridSpace](profile_gridspace_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 25.39 | 🟢 0/1088 (0%) | 2026-08-26 02:33 |
| [AlienOcean](profile_alienocean_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 24.98 | 🟢 0/1088 (0%) | 2026-08-26 02:28 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 24.58 | 🟢 0/1088 (0%) | 2026-08-26 02:26 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 22.97 | 🟢 0/1088 (0%) | 2026-08-26 03:05 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-08-26.md)§ ● | `fx_shader_draw` | 🟢 21.44 (3) | 🟢 0/2688 (0%) | 2026-08-26 01:30 |
| [Fishbowl](profile_fishbowl_teensy_2026-08-26.md) ● | `fish_build_vertices` | 🟢 21.29 | 🟢 0/1088 (0%) | 2026-08-26 01:18 |
| [AlienCore](profile_aliencore_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 20.10 | 🟢 0/1088 (0%) | 2026-08-26 02:30 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-08-26.md)§ ● | `sh_rasterize` | 🟢 11.70 (24) | 🟢 0/3488 (0%) | 2026-08-26 01:57 |
| [PetalFlow](profile_petalflow_teensy_2026-08-26.md) ● | `pf_draw_rings` | 🟢 10.77 | 🟢 0/1088 (0%) | 2026-08-26 01:32 |
| [Voronoi](profile_voronoi_teensy_2026-08-26.md) ● | `vo_shade` | 🟢 7.71 | 🟢 0/1088 (0%) | 2026-08-26 01:39 |
| [RingShower](profile_ringshower_teensy_2026-08-26.md) ● | `rsh_draw_rings` | 🟢 3.86 | 🟢 0/1088 (0%) | 2026-08-26 01:36 |

§ Cyclers carry one aligned line per colour bucket, worst first. The count in
parentheses is presets in that bucket; bucket frames include transitions and are
therefore stricter than the clean-hold table in each report.

**Colour:** 🟢 no spills; 🔴 one or more spilled frames.
**●** — capture refreshed 2026-08-26. Captured timestamps are raw-log local mtimes.
