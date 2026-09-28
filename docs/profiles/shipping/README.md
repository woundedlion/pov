# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
38 effects in `HS_PHANTASM_EFFECT_LIST`. Peak is worst-frame render
time; spilled counts frames whose render exceeded the 62.5 ms display window.
Rows rank by any spill, then peak render. 🟢 is zero spill, 🟡 is under 25% spill and 🔴 is at least 25% spill. Cyclers (§) use parser-owned cadence buckets.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [IslamicStars](profile_islamicstars_teensy_2026-09-24.md) § | `is_timeline_step` | 🟢 60.695 (23) | 🟢 0/3327 (0.0%) | 2026-09-24 11:01 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-09-25.md) § ● | `ss_draw_all` | 🟢 58.40 (9) | 🟢 0/2457 (0.00%) | 2026-09-25 07:33 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-08-26.md) § | `mf_feedback_flush` | 🟢 58.30 (13) | 🟢 0/6688 (0.0%) | 2026-08-26 03:31 |
| [DisplacementField](profile_displacementfield_teensy_2026-09-19.md) | `df_timeline_step` | 🟢 58.18 | 🟢 0/1088 (0.0%) | 2026-09-19 22:17 |
| [Raymarch](profile_raymarch_teensy_2026-09-20.md) | `rm_shader_draw` | 🟢 56.07 | 🟢 0/1736 (0.00%) | 2026-09-20 22:51 |
| [HyperLattice](profile_hyperlattice_teensy_2026-09-27.md) § ● | `hl_shader_draw` | 🟢 55.509 (2) | 🟢 0/1576 (0.00%) | 2026-09-27 01:02 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-08-26.md) | `grd_render` | 🟢 55.26 | 🟢 0/2048 (0.0%) | 2026-08-26 01:28 |
| [MindSplatter](profile_mindsplatter_teensy_2026-08-26.md) § | `msp_draw_particles` | 🟢 52.77 (9) | 🟢 0/1728 (0.0%) | 2026-08-26 07:40 |
| [AshCloud](profile_ashcloud_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 50.09 | 🟢 0/1088 (0.0%) | 2026-08-26 02:47 |
| [RingSpin](profile_ringspin_teensy_2026-09-24.md) | `rs_draw_rings` | 🟢 49.920 | 🟢 0/1087 (0.0%) | 2026-09-24 20:15 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-08-26.md) | `bz_render` | 🟢 48.91 | 🟢 0/2048 (0.0%) | 2026-08-26 01:20 |
| [HopfFibration](profile_hopffibration_teensy_2026-08-26.md) | `hf_render_trails` | 🟢 48.50 | 🟢 0/1088 (0.0%) | 2026-08-26 01:30 |
| [KaleidoscopeStainedGlass](profile_kaleidoscopestainedglass_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 47.91 | 🟢 0/1088 (0.0%) | 2026-08-26 03:23 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 45.60 | 🟢 0/1088 (0.0%) | 2026-08-26 02:46 |
| [HankinSolids](profile_hankinsolids_teensy_2026-08-26.md) § | `hk_timeline_step` | 🟢 45.01 (19) | 🟢 0/3328 (0.0%) | 2026-08-26 01:56 |
| [LatticeMelt](profile_latticemelt_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 43.73 (3) | 🟢 0/1728 (0.0%) | 2026-08-26 02:42 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 42.87 | 🟢 0/1088 (0.0%) | 2026-08-26 02:44 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 40.29 (3) | 🟢 0/2368 (0.0%) | 2026-08-26 03:22 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 39.04 (3) | 🟢 0/2208 (0.0%) | 2026-08-26 02:52 |
| [DreamBalls](profile_dreamballs_teensy_2026-08-26.md) § | `db_timeline_step` | 🟢 38.73 (11) | 🟢 0/3648 (0.0%) | 2026-08-26 02:21 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 36.69 (4) | 🟢 0/4128 (0.0%) | 2026-08-26 03:06 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 35.85 (5) | 🟢 0/4128 (0.0%) | 2026-08-26 02:59 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 35.55 (3) | 🟢 0/2368 (0.0%) | 2026-08-26 03:19 |
| [Comets](profile_comets_teensy_2026-08-26.md) § | `cm_draw_trail` | 🟢 33.91 (13) | 🟢 0/4128 (0.0%) | 2026-08-26 02:05 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 33.37 | 🟢 0/1088 (0.0%) | 2026-08-26 02:49 |
| [AlienBrain](profile_alienbrain_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 31.61 (5) | 🟢 0/4768 (0.0%) | 2026-08-26 02:27 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-08-26.md) | `gn_draw_stars` | 🟢 29.64 | 🟢 0/1088 (0.0%) | 2026-08-26 01:25 |
| [GridSpace](profile_gridspace_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 28.38 | 🟢 0/1088 (0.0%) | 2026-08-26 02:36 |
| [AlienOcean](profile_alienocean_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 27.79 | 🟢 0/1088 (0.0%) | 2026-08-26 02:31 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 27.30 | 🟢 0/1088 (0.0%) | 2026-08-26 02:29 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 24.01 | 🟢 0/1088 (0.0%) | 2026-08-26 03:07 |
| [Fishbowl](profile_fishbowl_teensy_2026-08-26.md) | `fish_build_vertices` | 🟢 23.22 | 🟢 0/1088 (0.0%) | 2026-08-26 01:22 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 22.37 (3) | 🟢 0/2688 (0.0%) | 2026-08-26 01:34 |
| [AlienCore](profile_aliencore_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 21.09 | 🟢 0/1088 (0.0%) | 2026-08-26 02:32 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-08-26.md) § | `sh_rasterize` | 🟢 12.64 (24) | 🟢 0/3488 (0.0%) | 2026-08-26 02:00 |
| [PetalFlow](profile_petalflow_teensy_2026-08-26.md) | `pf_draw_rings` | 🟢 11.85 | 🟢 0/1088 (0.0%) | 2026-08-26 01:35 |
| [Voronoi](profile_voronoi_teensy_2026-08-26.md) | `vo_shade` | 🟢 8.96 | 🟢 0/1088 (0.0%) | 2026-08-26 01:43 |
| [RingShower](profile_ringshower_teensy_2026-08-26.md) | `rsh_draw_rings` | 🟢 3.98 | 🟢 0/1088 (0.0%) | 2026-08-26 01:39 |

● ShapeShifter refreshed on 2026-09-25: N counts nine distinct presets, merging repeated visits after wrap. Its runtime buckets exclude setup frame 1 and include transitions.

IslamicStars refreshed on 2026-09-24 (unlanded finding 2 candidate); Raymarch on 2026-09-20; DisplacementField on 2026-09-19; RingSpin on 2026-09-24; other rows on 2026-08-26.
Captured timestamps are local raw-log mtimes.
Each row links to the report generated from its capture log.

For cyclers, each row summarizes the parser-owned preset, shape, or mode entries.
(N) gives the entry count; the linked report contains the individual buckets.
Spill fractions include the transition following an entry and are stricter than clean holds.

- **MindSplatter**: 9 parser ownership buckets spanning 21.61–52.77 ms; the sequence closes back to its first entry.
- **IslamicStars**: 23 shape ownership buckets including transitions; setup frame 1 excluded. The unlanded finding 2 candidate completes the full cycle. The controlled baseline comparison is no longer retained.
- **MeshFeedback**: 13 parser ownership buckets spanning 47.02–58.30 ms; the sequence closes back to its first entry.
- **ShapeShifter**: 9 distinct presets with live peaks spanning 10.007–58.395 ms; repeated visits after wrap are merged.
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

Shipping reports in this directory correspond exactly to
`HS_PHANTASM_EFFECT_LIST`.


● HyperLattice refreshed on 2026-09-27: both presets and transitions; setup frame 1 excluded.

## Supplemental experimental presets

Historical fixed opt-in captures; these do not change the normal firmware
roster ranking. Triangular was subsequently selected for removal.

| Preset | Dominant scope | Peak render ms | Spilled | Captured |
| --- | --- | ---: | ---: | --- |
| [HyperLattice Triangular](../hyperlattice_triangular_2026-09-27.md) | `hl_shader_draw` | 🔴 134.620 | 🔴 441/441 (100%) | 2026-09-27 21:59 |
| [HyperLattice Octet 3D ●](profile_hyperlattice_teensy_2026-09-27.md#supplemental-octet-3d-single-owner-correction) | `hl_shader_draw` | 🔴 96.869 | 🔴 548/548 (100%) | 2026-09-27 22:50 |

● Octet 3D refreshed on 2026-09-27 after the single-owner strut correction.
These fixed-preset captures retain the original oscillating camera path.
The [earlier Octet measurements](../hyperlattice_experimental_presets_2026-09-27.md#octet-validation-and-device-measurements)
record 3D shipping peak 103.733 ms (548/548 spills) before that correction.
The 4D capture remains historical: 549.828 ms peak, 120/120 spills,
captured 2026-09-27 at 22:18; no 4D global-O3 capture.

[Cubic Wide Flight](../hyperlattice_experimental_presets_2026-09-27.md#cubic-wide-flight)
has a separate fixed-preset shipping check: peak 50.804 ms, 0/616 spills,
captured 2026-09-27 22:25. The earlier HyperLattice cycle rows above cover
its two original presets; they predate this third default preset.

## Octet optimization supplement

| Preset | Peak ms | Spilled | Captured |
|---|---:|---:|---|
| [3: Octet 3D](profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) ● | 🔴 76.566 | 🔴 549/549 (100.00%) | 2026-09-28 00:17 |
| [5: Octet wide](profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) ● | 🔴 76.425 | 🔴 549/549 (100.00%) | 2026-09-28 00:22 |

● Updated 2026-09-28. Opt-in held presets, with startup excluded. [Earlier optimization and 4D evidence](../hyperlattice_octet_optimization_2026-09-27.md).
