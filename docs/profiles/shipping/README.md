# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
38 effects in `HS_PHANTASM_EFFECT_LIST`. Peak is worst-frame render
time; spilled counts frames whose render exceeded the 62.5 ms display window.
Rows rank by worst spill fraction, then peak render. 🟢 is zero spill, 🟡 is
under 25%, and 🔴 is 25% or more. Cyclers (§) use parser-owned cadence buckets.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [MindSplatter](profile_mindsplatter_teensy_2026-08-26.md) § ● | `msp_draw_particles` | 🟡 67.05 (2)<br>🟢 60.54 (7) | 🟡 39/424 (9.2%)<br>🟢 0/1272 (0.0%) | 2026-08-26 01:52 |
| [IslamicStars](profile_islamicstars_teensy_2026-08-26.md) § ● | `is_timeline_step` | 🟡 64.18 (1)<br>🟢 61.87 (22) | 🟡 2/156 (1.3%)<br>🟢 0/3172 (0.0%) | 2026-08-26 02:17 |
| [Raymarch](profile_raymarch_teensy_2026-08-26.md) ● | `rm_shader_draw` | 🟢 62.22 | 🟢 0/1088 (0.0%) | 2026-08-26 01:37 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-08-26.md) § ● | `mf_feedback_flush` | 🟢 58.30 (13) | 🟢 0/6688 (0.0%) | 2026-08-26 03:31 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-08-26.md) § ● | `ss_draw_all` | 🟢 58.13 (10) | 🟢 0/2448 (0.0%) | 2026-08-26 03:28 |
| [DisplacementField](profile_displacementfield_teensy_2026-08-26.md) ● | `df_timeline_step` | 🟢 57.40 | 🟢 0/1088 (0.0%) | 2026-08-26 01:24 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-08-26.md) ● | `grd_render` | 🟢 55.26 | 🟢 0/2048 (0.0%) | 2026-08-26 01:28 |
| [HyperLattice](profile_hyperlattice_teensy_2026-08-26.md) § ● | `hl_shader_draw` | 🟢 51.19 (4) | 🟢 0/2688 (0.0%) | 2026-08-26 03:32 |
| [AshCloud](profile_ashcloud_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 50.09 | 🟢 0/1088 (0.0%) | 2026-08-26 02:47 |
| [RingSpin](profile_ringspin_teensy_2026-08-26.md) ● | `rs_draw_rings` | 🟢 49.42 | 🟢 0/1088 (0.0%) | 2026-08-26 01:41 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-08-26.md) ● | `bz_render` | 🟢 48.91 | 🟢 0/2048 (0.0%) | 2026-08-26 01:20 |
| [HopfFibration](profile_hopffibration_teensy_2026-08-26.md) ● | `hf_render_trails` | 🟢 48.50 | 🟢 0/1088 (0.0%) | 2026-08-26 01:30 |
| [KaleidoscopeStainedGlass](profile_kaleidoscopestainedglass_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 47.91 | 🟢 0/1088 (0.0%) | 2026-08-26 03:23 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 45.60 | 🟢 0/1088 (0.0%) | 2026-08-26 02:46 |
| [HankinSolids](profile_hankinsolids_teensy_2026-08-26.md) § ● | `hk_timeline_step` | 🟢 45.01 (19) | 🟢 0/3328 (0.0%) | 2026-08-26 01:56 |
| [LatticeMelt](profile_latticemelt_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 43.73 (3) | 🟢 0/1728 (0.0%) | 2026-08-26 02:42 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 42.87 | 🟢 0/1088 (0.0%) | 2026-08-26 02:44 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 40.29 (3) | 🟢 0/2368 (0.0%) | 2026-08-26 03:22 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 39.04 (3) | 🟢 0/2208 (0.0%) | 2026-08-26 02:52 |
| [DreamBalls](profile_dreamballs_teensy_2026-08-26.md) § ● | `db_timeline_step` | 🟢 38.73 (11) | 🟢 0/3648 (0.0%) | 2026-08-26 02:21 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 36.69 (4) | 🟢 0/4128 (0.0%) | 2026-08-26 03:06 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 35.85 (5) | 🟢 0/4128 (0.0%) | 2026-08-26 02:59 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 35.55 (3) | 🟢 0/2368 (0.0%) | 2026-08-26 03:19 |
| [Comets](profile_comets_teensy_2026-08-26.md) § ● | `cm_draw_trail` | 🟢 33.91 (13) | 🟢 0/4128 (0.0%) | 2026-08-26 02:05 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 33.37 | 🟢 0/1088 (0.0%) | 2026-08-26 02:49 |
| [AlienBrain](profile_alienbrain_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 31.61 (5) | 🟢 0/4768 (0.0%) | 2026-08-26 02:27 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-08-26.md) ● | `gn_draw_stars` | 🟢 29.64 | 🟢 0/1088 (0.0%) | 2026-08-26 01:25 |
| [GridSpace](profile_gridspace_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 28.38 | 🟢 0/1088 (0.0%) | 2026-08-26 02:36 |
| [AlienOcean](profile_alienocean_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 27.79 | 🟢 0/1088 (0.0%) | 2026-08-26 02:31 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 27.30 | 🟢 0/1088 (0.0%) | 2026-08-26 02:29 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 24.01 | 🟢 0/1088 (0.0%) | 2026-08-26 03:07 |
| [Fishbowl](profile_fishbowl_teensy_2026-08-26.md) ● | `fish_build_vertices` | 🟢 23.22 | 🟢 0/1088 (0.0%) | 2026-08-26 01:22 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 22.37 (3) | 🟢 0/2688 (0.0%) | 2026-08-26 01:34 |
| [AlienCore](profile_aliencore_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 21.09 | 🟢 0/1088 (0.0%) | 2026-08-26 02:32 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-08-26.md) § ● | `sh_rasterize` | 🟢 12.64 (24) | 🟢 0/3488 (0.0%) | 2026-08-26 02:00 |
| [PetalFlow](profile_petalflow_teensy_2026-08-26.md) ● | `pf_draw_rings` | 🟢 11.85 | 🟢 0/1088 (0.0%) | 2026-08-26 01:35 |
| [Voronoi](profile_voronoi_teensy_2026-08-26.md) ● | `vo_shade` | 🟢 8.96 | 🟢 0/1088 (0.0%) | 2026-08-26 01:43 |
| [RingShower](profile_ringshower_teensy_2026-08-26.md) ● | `rsh_draw_rings` | 🟢 3.98 | 🟢 0/1088 (0.0%) | 2026-08-26 01:39 |

**● captured 2026-08-26.** Captured timestamps are local raw-log mtimes. Every row
links to the report generated from that exact log and its provenance sidecar.

For cyclers, each `<br>`-joined line is one colour bucket, worst first. Counts
in parentheses are parser ownership buckets; spill fractions include the
transition following a preset and are therefore stricter than clean holds.

- **MindSplatter**: 9 parser ownership buckets spanning 23.11–67.05 ms; the sequence closes back to its first entry.
- **IslamicStars**: 23 parser ownership buckets spanning 21.66–64.18 ms; the sequence closes back to its first entry.
- **MeshFeedback**: 13 parser ownership buckets spanning 47.02–58.30 ms; the sequence closes back to its first entry.
- **ShapeShifter**: 10 parser ownership buckets spanning 9.86–58.13 ms; the sequence closes back to its first entry.
- **HyperLattice**: 4 parser ownership buckets spanning 41.77–51.19 ms; the sequence closes back to its first entry.
- **HankinSolids**: 19 parser ownership buckets spanning 19.66–45.01 ms; the sequence closes back to its first entry.
- **LatticeMelt**: 3 parser ownership buckets spanning 42.00–43.73 ms; the sequence closes back to its first entry.
- **KaleidoscopeMandala**: 3 parser ownership buckets spanning 35.91–40.29 ms; the sequence closes back to its first entry.
- **KaleidoscopeHexOil**: 3 parser ownership buckets spanning 38.63–39.04 ms; the sequence closes back to its first entry.
- **DreamBalls**: 11 parser ownership buckets spanning 13.21–38.73 ms; the sequence closes back to its first entry.
- **KaleidoscopeFlowers**: 4 parser ownership buckets spanning 35.92–36.69 ms; the sequence closes back to its first entry.
- **KaleidoscopeSmooth**: 5 parser ownership buckets spanning 31.33–35.85 ms; the sequence closes back to its first entry.
- **KaleidoscopeHexBright**: 3 parser ownership buckets spanning 34.66–35.55 ms; the sequence closes back to its first entry.
- **Comets**: 13 parser ownership buckets spanning 16.67–33.91 ms; the sequence closes back to its first entry.
- **AlienBrain**: 5 parser ownership buckets spanning 29.12–31.61 ms; the sequence closes back to its first entry.
- **MobiusGrid**: 3 parser ownership buckets spanning 21.73–22.37 ms; the sequence closes back to its first entry.
- **SphericalHarmonics**: 24 parser ownership buckets spanning 8.05–12.64 ms; the sequence closes back to its first entry.

Retired-effect captures live in [`../retired/`](../retired/). Shipping reports
in this directory correspond exactly to `HS_PHANTASM_EFFECT_LIST`.
