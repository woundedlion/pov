# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
38 effects in `HS_PHANTASM_EFFECT_LIST`. Peak is worst-frame render
time; spilled counts frames whose render exceeded the 62.5 ms display window.
Rows rank by any spill, then peak render. 🟢 is zero spill, 🟡 is under 25% spill and 🔴 is at least 25% spill. Cyclers (§) use parser-owned cadence buckets.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [HyperLattice](profile_hyperlattice_teensy_2026-09-29.md) § | `hl_shader_draw` | 🟢 36.65 (3) | 🟢 0/1887 (0.00%) | 2026-09-29 14:54 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-09-28.md) | `grd_render` | 🟢 55.53 | 🟢 0/2047 (0.00%) | 2026-09-28 19:21 |
| [MindSplatter](profile_mindsplatter_teensy_2026-09-29.md) § | `msp_draw_particles` | 🟢 54.94 (8) | 🟢 0/1727 (0.00%) | 2026-09-29 14:48 |
| [Raymarch](profile_raymarch_teensy_2026-09-28.md) | `rm_shader_draw` | 🟢 54.87 | 🟢 0/1087 (0.00%) | 2026-09-28 19:09 |
| [HopfFibration](profile_hopffibration_teensy_2026-09-28.md) | `hf_render_trails` | 🟢 51.52 | 🟢 0/1087 (0.00%) | 2026-09-28 19:02 |
| [IslamicStars](profile_islamicstars_teensy_2026-09-28.md) § | `is_timeline_step` | 🟢 50.51 (23) | 🟢 0/3327 (0.00%) | 2026-09-28 18:40 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-09-28.md) | `bz_render` | 🟢 48.84 | 🟢 0/2047 (0.00%) | 2026-09-28 19:12 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-09-29.md) § | `mf_feedback_flush` | 🟢 48.61 (12) | 🟢 0/6688 (0.00%) | 2026-09-29 10:03 |
| [RingSpin](profile_ringspin_teensy_2026-09-28.md) | `rs_draw_rings` | 🟢 47.04 | 🟢 0/1087 (0.00%) | 2026-09-28 19:13 |
| [AshCloud](profile_ashcloud_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 44.57 | 🟢 0/1087 (0.00%) | 2026-09-28 17:33 |
| [KaleidoscopeStainedGlass](profile_kaleidoscopestainedglass_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 43.50 | 🟢 0/1087 (0.00%) | 2026-09-28 17:35 |
| [DreamBalls](profile_dreamballs_teensy_2026-09-28.md) § | `db_timeline_step` | 🟢 42.66 (10) | 🟢 0/3647 (0.00%) | 2026-09-28 18:45 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-09-29.md) ● § | `ss_draw_all` | 🟢 40.45 (9) | 🟢 0/2448 (0.00%) | 2026-09-29 19:01 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 39.57 | 🟢 0/1087 (0.00%) | 2026-09-28 18:49 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 38.51 (2) | 🟢 0/2207 (0.00%) | 2026-09-28 18:54 |
| [LatticeMelt](profile_latticemelt_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 37.18 (2) | 🟢 0/1727 (0.00%) | 2026-09-28 18:45 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 36.63 | 🟢 0/1087 (0.00%) | 2026-09-28 18:47 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 35.40 (2) | 🟢 0/2367 (0.00%) | 2026-09-28 18:39 |
| [DisplacementField](profile_displacementfield_teensy_2026-09-29.md) | `df_timeline_step` | 🟢 35.27 | 🟢 0/2368 (0.00%) | 2026-09-29 09:09 |
| [HankinSolids](profile_hankinsolids_teensy_2026-09-28.md) § | `hk_timeline_step` | 🟢 34.81 (19) | 🟢 0/4447 (0.00%) | 2026-09-28 19:27 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 32.93 (4) | 🟢 0/4127 (0.00%) | 2026-09-28 18:59 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 32.52 (3) | 🟢 0/4127 (0.00%) | 2026-09-28 19:07 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 31.64 (2) | 🟢 0/2367 (0.00%) | 2026-09-28 19:02 |
| [Comets](profile_comets_teensy_2026-09-28.md) § | `cm_draw_trail` | 🟢 30.66 (12) | 🟢 0/4127 (0.00%) | 2026-09-28 18:28 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 27.24 | 🟢 0/1087 (0.00%) | 2026-09-28 18:51 |
| [AlienBrain](profile_alienbrain_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 25.89 (4) | 🟢 0/4767 (0.00%) | 2026-09-28 18:30 |
| [GridSpace](profile_gridspace_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 25.31 | 🟢 0/1087 (0.00%) | 2026-09-28 17:38 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 24.52 | 🟢 0/1087 (0.00%) | 2026-09-28 18:32 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 24.11 | 🟢 0/1087 (0.00%) | 2026-09-28 19:09 |
| [Fishbowl](profile_fishbowl_teensy_2026-09-28.md) | `fish_build_vertices` | 🟢 23.98 | 🟢 0/1087 (0.00%) | 2026-09-28 19:14 |
| [AlienOcean](profile_alienocean_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 22.56 | 🟢 0/1087 (0.00%) | 2026-09-28 18:34 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 22.05 (2) | 🟢 0/2687 (0.00%) | 2026-09-28 19:05 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-09-28.md) | `gn_draw_stars` | 🟢 21.86 | 🟢 0/1087 (0.00%) | 2026-09-28 19:18 |
| [AlienCore](profile_aliencore_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 18.03 | 🟢 0/1087 (0.00%) | 2026-09-28 18:36 |
| [PetalFlow](profile_petalflow_teensy_2026-09-28.md) | `pf_draw_rings` | 🟢 13.61 | 🟢 0/1087 (0.00%) | 2026-09-28 19:07 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-09-28.md) § | `sh_rasterize` | 🟢 12.89 (24) | 🟢 0/3487 (0.00%) | 2026-09-28 19:00 |
| [Voronoi](profile_voronoi_teensy_2026-09-28.md) | `vo_shade` | 🟢 8.51 | 🟢 0/1087 (0.00%) | 2026-09-28 19:15 |
| [RingShower](profile_ringshower_teensy_2026-09-28.md) | `rsh_draw_rings` | 🟢 4.33 | 🟢 0/1087 (0.00%) | 2026-09-28 19:11 |

Every row was refreshed on 2026-09-28 after the finding-10 placement change `97eb0bf78`. AshCloud, KaleidoscopeStainedGlass and GridSpace come from the placement A/B's variant-6 arm, which carries the same kernel attributes. Setup frame 1 is excluded from every peak and spill figure; each report states its frame-1 render. Boards are named per report.

For cyclers, each row summarizes the parser-owned preset, shape, or mode entries.
(N) gives the entry count; the linked report contains the individual buckets.
Spill fractions include the transition following an entry and are stricter than clean holds.

- **MeshFeedback**: 12 parser ownership buckets spanning 40.88–48.61 ms.
- **ShapeShifter**: 9 parser ownership buckets spanning 6.87–40.45 ms.
- **HyperLattice**: 3 parser ownership buckets spanning 29.60–36.65 ms (setup frame excluded).
- **MindSplatter**: 8 parser ownership buckets spanning 22.36–54.94 ms.
- **IslamicStars**: 23 parser ownership buckets spanning 19.92–50.51 ms.
- **DreamBalls**: 10 parser ownership buckets spanning 14.41–42.66 ms.
- **KaleidoscopeHexOil**: 2 parser ownership buckets spanning 38.03–38.51 ms.
- **LatticeMelt**: 2 parser ownership buckets spanning 37.18–37.18 ms.
- **KaleidoscopeMandala**: 2 parser ownership buckets spanning 34.88–35.40 ms.
- **HankinSolids**: 19 parser ownership buckets spanning 14.52–34.81 ms.
- **KaleidoscopeSmooth**: 4 parser ownership buckets spanning 28.62–32.93 ms.
- **KaleidoscopeFlowers**: 3 parser ownership buckets spanning 32.31–32.52 ms.
- **KaleidoscopeHexBright**: 2 parser ownership buckets spanning 31.44–31.64 ms.
- **Comets**: 12 parser ownership buckets spanning 16.48–30.66 ms.
- **AlienBrain**: 4 parser ownership buckets spanning 25.80–25.89 ms.
- **MobiusGrid**: 2 parser ownership buckets spanning 22.05–22.05 ms.
- **SphericalHarmonics**: 24 parser ownership buckets spanning 8.30–12.89 ms.

Shipping reports in this directory correspond exactly to
`HS_PHANTASM_EFFECT_LIST`.


● HyperLattice refreshed on 2026-09-29 after its optimization campaign: all three presets and their segues; setup frame 1 excluded.

## Supplemental experimental presets

Historical fixed opt-in captures; these do not change the normal firmware
roster ranking. Triangular was subsequently selected for removal.

| Preset | Dominant scope | Peak render ms | Spilled | Captured |
| --- | --- | ---: | ---: | --- |
| [HyperLattice Triangular](../hyperlattice_triangular_2026-09-27.md) | `hl_shader_draw` | 🔴 134.620 | 🔴 441/441 (100%) | 2026-09-27 21:59 |
| [HyperLattice Octet 3D ●](profile_hyperlattice_teensy_2026-09-29.md) | `hl_shader_draw` | 🔴 96.869 | 🔴 548/548 (100%) | 2026-09-27 22:50 |

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

[2026-09-27 Octet snapshot](profile_hyperlattice_octet_teensy_2026-09-27.md).

| Preset | Peak ms | Spilled | Captured |
|---|---:|---:|---|
| [4: Octet 4D](profile_hyperlattice_octet4d_teensy_2026-09-28.md) ● | 🔴 143.318 | 🔴 233/233 (100.00%) | 2026-09-28 09:53 |
| [3: Octet 3D](profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) ● | 🟢 45.688 | 🟢 0/1096 (0.00%) | 2026-09-28 09:56 |
| [5: Octet wide](profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) ● | 🟢 45.242 | 🟢 0/1096 (0.00%) | 2026-09-28 09:53 |

● Updated 2026-09-28 09:53–09:56 after the Octet 58 ms Tiers 1–2 (`1bb522211`). Opt-in held presets, startup excluded; 4D on COM3, 3D on COM4. [Earlier optimization and 4D evidence](../hyperlattice_octet_optimization_2026-09-27.md).

## Shell flight supplement

Sphere-only renderer; two authored fixed presets on COM4, startup excluded. [Paired surface report](profile_hyperlattice_shell_flight_teensy_2026-09-29.md).

| Preset | Peak render ms | Spilled/live frames | Captured local |
|---|--:|--:|---|
| 5: Shell Flight | 🔴 123.553 | 🔴 228/228 (100.00%) | 2026-09-29 08:03 |
| 6: Shell Close Flight | 🔴 125.186 | 🔴 228/228 (100.00%) | 2026-09-29 08:06 |
