# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
38 effects in `HS_PHANTASM_EFFECT_LIST`. Peak is worst-frame render
time; spilled counts frames whose render exceeded the 62.5 ms display window.
Rows rank by any spill, then peak render. 🟢 is zero spill; 🔴 is any spill. Cyclers (§) use parser-owned cadence buckets.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [MindSplatter](profile_mindsplatter_teensy_2026-10-06.md) § | `msp_draw_particles` | 🟢 56.60 (8) | 🟢 0/1727 (0.00%) | 2026-10-06 18:32 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-10-06.md) § | `ss_draw_all` | 🟢 54.27 (9) | 🟢 0/2447 (0.00%) | 2026-10-06 18:29 |
| [Raymarch](profile_raymarch_teensy_2026-10-06.md) | `rm_shader_draw` | 🟢 51.95 | 🟢 0/1087 (0.00%) | 2026-10-06 18:17 |
| [HyperLattice](profile_hyperlattice_teensy_2026-10-06.md) § | `hl_shader_draw` | 🟢 51.77 (9) | 🟢 0/5487 (0.00%) | 2026-10-06 20:43 |
| [HopfFibration](profile_hopffibration_teensy_2026-10-06.md) | `hf_render_trails` | 🟢 48.61 | 🟢 0/1087 (0.00%) | 2026-10-06 18:03 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-10-06.md) § | `mf_feedback_flush` | 🟢 48.53 (12) | 🟢 0/6687 (0.00%) | 2026-10-06 18:18 |
| [IslamicStars](profile_islamicstars_teensy_2026-10-07.md) § ● | `is_timeline_step` | 🟢 48.510 (23) | 🟢 0/3327 (0.00%) | 2026-10-07 15:16 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-10-06.md) | `bz_render` | 🟢 48.34 | 🟢 0/2047 (0.00%) | 2026-10-06 19:05 |
| [RingSpin](profile_ringspin_teensy_2026-10-06.md) | `rs_draw_rings` | 🟢 47.16 | 🟢 0/1087 (0.00%) | 2026-10-06 18:22 |
| [KaleidoscopeStainedGlass](profile_kaleidoscopestainedglass_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 43.00 | 🟢 0/1087 (0.00%) | 2026-10-06 18:50 |
| [AshCloud](profile_ashcloud_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 42.67 | 🟢 0/1087 (0.00%) | 2026-10-06 19:01 |
| [DreamBalls](profile_dreamballs_teensy_2026-10-06.md) § | `db_timeline_step` | 🟢 41.44 (10) | 🟢 0/3647 (0.00%) | 2026-10-06 18:28 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-10-06.md) | `grd_shader_draw` | 🟢 39.43 | 🟢 0/2047 (0.00%) | 2026-10-06 19:15 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 39.18 | 🟢 0/1087 (0.00%) | 2026-10-06 18:59 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 37.08 (2) | 🟢 0/2207 (0.00%) | 2026-10-06 18:47 |
| [LatticeMelt](profile_latticemelt_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 35.83 (2) | 🟢 0/1727 (0.00%) | 2026-10-06 18:55 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 35.53 | 🟢 0/1087 (0.00%) | 2026-10-06 18:57 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 34.97 (2) | 🟢 0/2367 (0.00%) | 2026-10-06 18:45 |
| [HankinSolids](profile_hankinsolids_teensy_2026-10-07.md) § ● | `hk_timeline_step` | 🟢 34.681 (19) | 🟢 0/4447 (0.00%) | 2026-10-07 14:18 |
| [DisplacementField](profile_displacementfield_teensy_2026-10-06.md) | `df_timeline_step` | 🟢 33.66 | 🟢 0/2527 (0.00%) | 2026-10-06 19:24 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 32.02 (4) | 🟢 0/4127 (0.00%) | 2026-10-06 18:55 |
| [Comets](profile_comets_teensy_2026-10-06.md) § | `cm_draw_trail` | 🟢 31.92 (12) | 🟢 0/4127 (0.00%) | 2026-10-06 18:07 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 31.42 (3) | 🟢 0/4127 (0.00%) | 2026-10-06 19:04 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 30.03 (2) | 🟢 0/2367 (0.00%) | 2026-10-06 18:58 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 27.04 | 🟢 0/1087 (0.00%) | 2026-10-06 18:44 |
| [AlienBrain](profile_alienbrain_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 25.11 (4) | 🟢 0/4767 (0.00%) | 2026-10-06 18:35 |
| [GridSpace](profile_gridspace_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 24.57 | 🟢 0/1087 (0.00%) | 2026-10-06 18:48 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 23.88 | 🟢 0/1087 (0.00%) | 2026-10-06 18:37 |
| [Fishbowl](profile_fishbowl_teensy_2026-10-06.md) | `fish_build_vertices` | 🟢 23.78 | 🟢 0/1087 (0.00%) | 2026-10-06 19:07 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 23.63 | 🟢 0/1087 (0.00%) | 2026-10-06 19:06 |
| [AlienOcean](profile_alienocean_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 22.37 | 🟢 0/1087 (0.00%) | 2026-10-06 18:39 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-10-07.md) | `gn_draw_stars` | 🟢 21.84 | 🟢 0/1087 (0.00%) | 2026-10-07 13:42 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-10-06.md) § | `fx_shader_draw` | 🟢 21.77 (2) | 🟢 0/2687 (0.00%) | 2026-10-06 18:08 |
| [AlienCore](profile_aliencore_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 17.86 | 🟢 0/1087 (0.00%) | 2026-10-06 18:42 |
| [PetalFlow](profile_petalflow_teensy_2026-10-06.md) | `pf_draw_rings` | 🟢 13.29 | 🟢 0/1087 (0.00%) | 2026-10-06 18:12 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-10-06.md) § | `sh_rasterize` | 🟢 12.92 (24) | 🟢 0/3487 (0.00%) | 2026-10-06 18:42 |
| [Voronoi](profile_voronoi_teensy_2026-10-06.md) | `vo_shade` | 🟢 8.00 | 🟢 0/1087 (0.00%) | 2026-10-06 18:25 |
| [RingShower](profile_ringshower_teensy_2026-10-06.md) | `rsh_draw_rings` | 🟢 4.32 | 🟢 0/1087 (0.00%) | 2026-10-06 18:20 |

The full-roster sweep was captured on 2026-10-06 at `e2f5b0a3d` across COM3
and COM4, except HyperLattice, re-captured after its 4D optimization at
`a13c149d4`; its linked report records the 51.774 ms peak. GnomonicStars was
re-captured on October 7; its report records the source and clock variant.
● marks the October 7 IslamicStars and HankinSolids optimized face-distance
measurements. Each report names its board and source snapshot; setup frame 1
is excluded from peak and spill figures.

HankinSolids, HyperLattice and DisplacementField were re-captured with longer
budgets than `tools/profile_sweep.sh` assigns them, so that each capture covers its
full cycle; their reports give the commands.

For cyclers, each row summarizes the parser-owned preset, shape, or mode entries.
(N) gives the entry count; the linked report contains the individual buckets.
Spill fractions include the transition following an entry and are stricter than clean holds.

- **HyperLattice**: 9 parser ownership buckets spanning 26.22–51.77 ms.
- **MindSplatter**: 8 parser ownership buckets spanning 22.99–56.60 ms.
- **ShapeShifter**: 9 parser ownership buckets spanning 6.28–54.27 ms.
- **IslamicStars**: 23 parser ownership buckets spanning 20.32–48.51 ms.
- **MeshFeedback**: 12 parser ownership buckets spanning 40.71–48.53 ms.
- **DreamBalls**: 10 parser ownership buckets spanning 14.27–41.44 ms.
- **KaleidoscopeHexOil**: 2 parser ownership buckets spanning 36.83–37.08 ms.
- **LatticeMelt**: 2 parser ownership buckets spanning 35.77–35.83 ms.
- **KaleidoscopeMandala**: 2 parser ownership buckets spanning 34.40–34.97 ms.
- **HankinSolids**: 19 parser ownership buckets spanning 14.569–34.681 ms (18 named shapes plus startup).
- **KaleidoscopeSmooth**: 4 parser ownership buckets spanning 27.36–32.02 ms.
- **Comets**: 12 parser ownership buckets spanning 17.05–31.92 ms.
- **KaleidoscopeFlowers**: 3 parser ownership buckets spanning 30.78–31.42 ms.
- **KaleidoscopeHexBright**: 2 parser ownership buckets spanning 29.85–30.03 ms.
- **AlienBrain**: 4 parser ownership buckets spanning 24.86–25.11 ms.
- **MobiusGrid**: 2 parser ownership buckets spanning 21.68–21.77 ms.
- **SphericalHarmonics**: 24 parser ownership buckets spanning 8.36–12.92 ms.

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

Preserved architecture snapshots, superseded in the ranking by the 2026-10-06 sweep.

- [MindSplatter 2026-10-01](profile_mindsplatter_architecture_teensy_2026-10-01.md).
- [HyperLattice 2026-10-01](profile_hyperlattice_architecture_teensy_2026-10-01.md).
- [LatticeMelt 2026-10-01](profile_latticemelt_architecture_teensy_2026-10-01.md).
- [ChromaticLichen 2026-10-01](profile_chromaticlichen_architecture_teensy_2026-10-01.md).
- [KaleidoscopeSmooth 2026-10-01](profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md).

Experimental clock characterization: [GnomonicStars bounded-channel comparison](profile_gnomonicstars_clock_experiment_teensy_2026-10-07.md). The candidate remains outside the shipping implementation; its numbers are not substituted into the ranked row.
