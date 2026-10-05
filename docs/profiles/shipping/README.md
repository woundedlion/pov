# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
38 effects in `HS_PHANTASM_EFFECT_LIST`. Peak is worst-frame render
time; spilled counts frames whose render exceeded the 62.5 ms display window.
Rows rank by any spill, then peak render. 🟢 is zero spill, 🟡 is under 25% spill and 🔴 is at least 25% spill. Cyclers (§) use parser-owned cadence buckets.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [Raymarch](profile_raymarch_teensy_2026-09-28.md) | `rm_shader_draw` | 🟢 54.87 | 🟢 0/1087 (0.00%) | 2026-09-28 19:09 |
| [MindSplatter](profile_mindsplatter_architecture_teensy_2026-10-01.md) § ● | `msp_draw_particles` | 🟢 54.457 (8) | 🟢 0/1736 (0.00%) | 2026-10-01 09:15 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-09-29.md) § | `ss_draw_all` | 🟢 53.22 (9) | 🟢 0/2448 (0.00%) | 2026-09-29 22:38 |
| [HopfFibration](profile_hopffibration_teensy_2026-09-28.md) | `hf_render_trails` | 🟢 51.52 | 🟢 0/1087 (0.00%) | 2026-09-28 19:02 |
| [IslamicStars](profile_islamicstars_teensy_2026-09-28.md) § | `is_timeline_step` | 🟢 50.51 (23) | 🟢 0/3327 (0.00%) | 2026-09-28 18:40 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-09-28.md) | `bz_render` | 🟢 48.84 | 🟢 0/2047 (0.00%) | 2026-09-28 19:12 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-09-29.md) § | `mf_feedback_flush` | 🟢 48.61 (12) | 🟢 0/6688 (0.00%) | 2026-09-29 10:03 |
| [RingSpin](profile_ringspin_teensy_2026-09-28.md) | `rs_draw_rings` | 🟢 47.04 | 🟢 0/1087 (0.00%) | 2026-09-28 19:13 |
| [AshCloud](profile_ashcloud_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 44.57 | 🟢 0/1087 (0.00%) | 2026-09-28 17:33 |
| [KaleidoscopeStainedGlass](profile_kaleidoscopestainedglass_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 43.50 | 🟢 0/1087 (0.00%) | 2026-09-28 17:35 |
| [DreamBalls](profile_dreamballs_teensy_2026-09-28.md) § | `db_timeline_step` | 🟢 42.66 (10) | 🟢 0/3647 (0.00%) | 2026-09-28 18:45 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-09-28.md) | `fx_shader_draw` | 🟢 39.57 | 🟢 0/1087 (0.00%) | 2026-09-28 18:49 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-10-05.md) ● | `grd_shader_draw` | 🟢 39.371 | 🟢 0/2032 (0.00%) | 2026-10-05 12:20 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 38.51 (2) | 🟢 0/2207 (0.00%) | 2026-09-28 18:54 |
| [HyperLattice](profile_hyperlattice_architecture_teensy_2026-10-01.md) § ● | `hl_shader_draw` | 🟢 36.839 (3) | 🟢 0/2696 (0.00%) | 2026-10-01 09:19 |
| [LatticeMelt](profile_latticemelt_architecture_teensy_2026-10-01.md) § ● | `fx_shader_draw` | 🟢 36.583 (2) | 🟢 0/1736 (0.00%) | 2026-10-01 09:38 |
| [ChromaticLichen](profile_chromaticlichen_architecture_teensy_2026-10-01.md) ● | `fx_shader_draw` | 🟢 35.773 | 🟢 0/1096 (0.00%) | 2026-10-01 09:40 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 35.40 (2) | 🟢 0/2367 (0.00%) | 2026-09-28 18:39 |
| [DisplacementField](profile_displacementfield_teensy_2026-09-29.md) | `df_timeline_step` | 🟢 35.27 | 🟢 0/2368 (0.00%) | 2026-09-29 09:09 |
| [HankinSolids](profile_hankinsolids_teensy_2026-09-28.md) § | `hk_timeline_step` | 🟢 34.81 (19) | 🟢 0/4447 (0.00%) | 2026-09-28 19:27 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-09-28.md) § | `fx_shader_draw` | 🟢 32.52 (3) | 🟢 0/4127 (0.00%) | 2026-09-28 19:07 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md) § ● | `fx_shader_draw` | 🟢 32.230 (4) | 🟢 0/4136 (0.00%) | 2026-10-01 11:05 |
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

● GSReactionDiffusion shipping refreshed 2026-10-05. Architecture captures refreshed 2026-10-01: MindSplatter, HyperLattice, LatticeMelt, ChromaticLichen and KaleidoscopeSmooth. Setup frame 1 is excluded. Cycler buckets include every live owner-attributed transition frame; each report lists per-preset peaks, fractions and the worst preset. Captured times are local raw-log mtimes.

The shipping set was refreshed on 2026-09-28 after the kernel-placement change `97eb0bf78`. AshCloud, KaleidoscopeStainedGlass and GridSpace come from the placement A/B's variant-6 arm, which carries the same kernel attributes. Setup frame 1 is excluded from every peak and spill figure; each report states its frame-1 render. Boards are named per report.

For cyclers, each row summarizes the parser-owned preset, shape, or mode entries.
(N) gives the entry count; the linked report contains the individual buckets.
Spill fractions include the transition following an entry and are stricter than clean holds.

- **MeshFeedback**: 12 parser ownership buckets spanning 40.88–48.61 ms.
- **ShapeShifter**: 9 parser ownership buckets spanning 6.30–53.22 ms (entry 1 at 288 contours).
- **HyperLattice**: 3 parser ownership buckets spanning 26.686–36.839 ms (setup frame excluded).
- **MindSplatter**: 8 parser ownership buckets spanning 22.224–54.457 ms (setup frame excluded).
- **IslamicStars**: 23 parser ownership buckets spanning 19.92–50.51 ms.
- **DreamBalls**: 10 parser ownership buckets spanning 14.41–42.66 ms.
- **KaleidoscopeHexOil**: 2 parser ownership buckets spanning 38.03–38.51 ms.
- **LatticeMelt**: 2 parser ownership buckets spanning 36.535–36.583 ms (setup frame excluded).
- **KaleidoscopeMandala**: 2 parser ownership buckets spanning 34.88–35.40 ms.
- **HankinSolids**: 19 parser ownership buckets spanning 14.52–34.81 ms.
- **KaleidoscopeSmooth**: 4 parser ownership buckets spanning 28.432–32.230 ms (setup frame excluded).
- **KaleidoscopeFlowers**: 3 parser ownership buckets spanning 32.31–32.52 ms.
- **KaleidoscopeHexBright**: 2 parser ownership buckets spanning 31.44–31.64 ms.
- **Comets**: 12 parser ownership buckets spanning 16.48–30.66 ms.
- **AlienBrain**: 4 parser ownership buckets spanning 25.80–25.89 ms.
- **MobiusGrid**: 2 parser ownership buckets spanning 22.05–22.05 ms.
- **SphericalHarmonics**: 24 parser ownership buckets spanning 8.30–12.89 ms.

Each effect in `HS_PHANTASM_EFFECT_LIST` has exactly one un-suffixed report
here; suffixed reports are supplements or preserved captures.

## Supplemental experimental presets

Historical fixed opt-in captures; these do not change the normal firmware
roster ranking. Triangular has since been removed.

| Preset | Dominant scope | Peak render ms | Spilled | Captured |
| --- | --- | ---: | ---: | --- |
| HyperLattice Triangular (supporting artifact removed) | `hl_shader_draw` | 🔴 134.620 | 🔴 441/441 (100%) | 2026-09-27 21:59 |
| HyperLattice Octet 3D (supporting artifact removed) | `hl_shader_draw` | 🔴 96.869 | 🔴 548/548 (100%) | 2026-09-27 22:50 |

Octet 3D refreshed on 2026-09-27 after the single-owner strut correction.
These fixed-preset captures retain the original oscillating camera path.
The earlier Octet measurements (supporting artifact removed)
record 3D shipping peak 103.733 ms (548/548 spills) before that correction.
The 4D capture remains historical: 549.828 ms peak, 120/120 spills,
captured 2026-09-27 at 22:18; no 4D global-O3 capture.

Cubic Wide Flight (supporting artifact removed)
has a separate fixed-preset shipping check: peak 50.804 ms, 0/616 spills,
captured 2026-09-27 22:25.

## Octet optimization supplement

[2026-09-27 Octet snapshot](profile_hyperlattice_octet_teensy_2026-09-27.md).

| Preset | Peak ms | Spilled | Captured |
|---|---:|---:|---|
| [Octet 4D (index 4 at capture time)](profile_hyperlattice_octet4d_teensy_2026-09-28.md) ● | 🔴 143.318 | 🔴 233/233 (100.00%) | 2026-09-28 09:53 |
| [3: Octet 3D](profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) ● | 🟢 45.688 | 🟢 0/1096 (0.00%) | 2026-09-28 09:56 |
| [Octet wide (index 5 at capture time)](profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) ● | 🟢 45.242 | 🟢 0/1096 (0.00%) | 2026-09-28 09:53 |

Updated 2026-09-28 09:53–09:56 after the Octet 58 ms Tiers 1–2 (`1bb522211`). Opt-in held presets, startup excluded; 4D on COM3, 3D on COM4. Earlier optimization and 4D evidence (supporting artifact removed).

## Shell flight supplement

Sphere-only renderer; two authored fixed presets on COM4, startup excluded. [Paired surface report](profile_hyperlattice_shell_flight_teensy_2026-09-29.md).

| Preset | Peak render ms | Spilled/live frames | Captured local |
|---|--:|--:|---|
| Shell Flight (index 5 at capture time) | 🔴 123.553 | 🔴 228/228 (100.00%) | 2026-09-29 08:03 |
| Shell Close Flight (index 6 at capture time) | 🔴 125.186 | 🔴 228/228 (100.00%) | 2026-09-29 08:06 |

## Earlier captures

Preserved canonical snapshots; current architecture rankings use the dated variants above. Each variant records its exact checkpoint source and paired baseline.

- [MindSplatter 2026-09-29](profile_mindsplatter_teensy_2026-09-29.md).
- [HyperLattice 2026-09-29](profile_hyperlattice_teensy_2026-09-29.md).
- [LatticeMelt 2026-09-28](profile_latticemelt_teensy_2026-09-28.md).
- [ChromaticLichen 2026-09-28](profile_chromaticlichen_teensy_2026-09-28.md).
- [KaleidoscopeSmooth 2026-09-28](profile_kaleidoscopesmooth_teensy_2026-09-28.md).
