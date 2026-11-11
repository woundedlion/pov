# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
39 effects in `HS_PHANTASM_EFFECT_LIST`. All 39 have canonical reports: 39 completed captures. Steady peak is worst-frame render excluding setup frame 1; setup render is listed separately. Steady spilled counts frames whose render exceeded
the 62.5 ms display window. Rows rank by spill fraction, then steady peak render.
🟢 zero spill; 🔴 any spill.

| Effect | Dominant scope | Peak steady ms | Setup frame 1 ms | Spilled | Captured |
|---|---|--:|--:|--:|---|
| [IslamicStars](profile_islamicstars_teensy_2026-11-05.md) § ● | `is_timeline_step` | 🟢 59.311 (16) | 16.664 | 🟢 0/4863 (0.00%) | 2026-11-04 23:39 |
| [AshCloud](profile_ashcloud_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 58.057 | 80.859 | 🟢 0/927 (0.00%) | 2026-11-05 08:38 |
| [MindSplatter](profile_mindsplatter_teensy_2026-11-05.md) § ● | `msp_draw_particles` | 🟢 57.896 (8) | 0.389 | 🟢 0/1247 (0.00%) | 2026-11-04 23:16 |
| [HopfFibration](profile_hopffibration_teensy_2026-11-05.md) ● | `hf_render_trails` | 🟢 57.199 | 0.646 | 🟢 0/927 (0.00%) | 2026-11-04 23:01 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-11-05.md) § ● | `ss_draw_all` | 🟢 57.031 (9) | 78.200 | 🟢 0/4271 (0.00%) | 2026-11-04 23:13 |
| [RingSpin](profile_ringspin_teensy_2026-11-05.md) ● | `rs_draw_rings` | 🟢 53.875 | 3.399 | 🟢 0/927 (0.00%) | 2026-11-04 23:07 |
| [HyperLattice](profile_hyperlattice_teensy_2026-11-05.md) § ● | `hl_shader_draw` | 🟢 53.649 (9) | 43.112 | 🟢 0/4767 (0.00%) | 2026-11-05 00:03 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-11-05.md) § ● | `mf_feedback_flush` | 🟢 53.505 (12) | 6.712 | 🟢 0/2863 (0.00%) | 2026-11-04 23:32 |
| [Raymarch](profile_raymarch_teensy_2026-11-05.md) ● | `rm_shader_draw` | 🟢 52.693 | 93.551 | 🟢 0/447 (0.00%) | 2026-11-04 23:06 |
| [DreamBalls](profile_dreamballs_teensy_2026-11-05.md) § ● | `db_timeline_step` | 🟢 52.024 (10) | 52.869 | 🟢 0/3167 (0.00%) | 2026-11-04 23:43 |
| [Millefiori](profile_millefiori_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 51.758 | 78.463 | 🟢 0/927 (0.00%) | 2026-11-04 23:58 |
| [GridWarp](profile_gridwarp_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 49.281 | 74.267 | 🟢 0/927 (0.00%) | 2026-11-05 00:19 |
| [Iridescence](profile_iridescence_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 49.094 | 91.161 | 🟢 0/927 (0.00%) | 2026-11-05 00:22 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-11-05.md) ● | `grd_render` | 🟢 48.568 | 24.927 | 🟢 0/447 (0.00%) | 2026-11-04 22:59 |
| [Involution](profile_involution_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 47.347 | 83.449 | 🟢 0/927 (0.00%) | 2026-11-05 00:15 |
| [Galaxies](profile_galaxies_teensy_2026-11-05.md) ● | `gx_particle_step` | 🟢 45.126 | 1.482 | 🟢 0/927 (0.00%) | 2026-11-04 22:44 |
| [DisplacementField](profile_displacementfield_teensy_2026-11-07.md) § ● | `df_timeline_step` | 🟢 45.107 (2) | 18.485 | 🟢 0/2367 (0.00%) | 2026-11-07 23:11 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-11-05.md) ● | `bz_render` | 🟢 44.574 | 76.467 | 🟢 0/927 (0.00%) | 2026-11-04 22:48 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-11-05.md) § ● | `fx_shader_draw` | 🟢 42.739 (2) | 64.707 | 🟢 0/1647 (0.00%) | 2026-11-05 00:26 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-11-05.md) § ● | `fx_shader_draw` | 🟢 42.737 (2) | 62.398 | 🟢 0/1631 (0.00%) | 2026-11-04 23:53 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 41.164 | 64.523 | 🟢 0/447 (0.00%) | 2026-11-05 00:07 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 40.138 | 63.665 | 🟢 0/447 (0.00%) | 2026-11-05 00:05 |
| [GridWorms](profile_gridworms_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 37.926 | 56.868 | 🟢 0/927 (0.00%) | 2026-11-04 23:56 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-11-05.md) § ● | `fx_shader_draw` | 🟢 36.868 (3) | 57.037 | 🟢 0/2735 (0.00%) | 2026-11-05 00:39 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-11-05.md) § ● | `fx_shader_draw` | 🟢 36.811 (2) | 53.101 | 🟢 0/1631 (0.00%) | 2026-11-05 00:35 |
| [LiquidWindows](profile_liquidwindows_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 36.585 | 57.965 | 🟢 0/927 (0.00%) | 2026-11-05 00:17 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-11-05.md) § ● | `fx_shader_draw` | 🟢 36.376 (4) | 56.382 | 🟢 0/3807 (0.00%) | 2026-11-05 00:32 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-11-05.md) § ● | `fx_shader_draw` | 🟢 32.874 (2) | 53.328 | 🟢 0/1647 (0.00%) | 2026-11-04 23:03 |
| [AlienBrain](profile_alienbrain_teensy_2026-11-05.md) § ● | `fx_shader_draw` | 🟢 32.584 (4) | 53.448 | 🟢 0/3807 (0.00%) | 2026-11-04 23:48 |
| [Comets](profile_comets_teensy_2026-11-05.md) § ● | `cm_draw_trail` | 🟢 30.637 (12) | 0.860 | 🟢 0/1887 (0.00%) | 2026-11-04 23:28 |
| [HankinSolids](profile_hankinsolids_teensy_2026-11-05.md) § ● | `hk_timeline_step` | 🟢 30.408 (19) | 19.120 | 🟢 0/2687 (0.00%) | 2026-11-04 23:22 |
| [Fishbowl](profile_fishbowl_teensy_2026-11-05.md) ● | `fish_build_vertices` | 🟢 28.932 | 1.755 | 🟢 0/927 (0.00%) | 2026-11-04 22:52 |
| [GridSpace](profile_gridspace_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 28.401 | 41.152 | 🟢 0/927 (0.00%) | 2026-11-04 23:55 |
| [AlienOcean](profile_alienocean_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 26.435 | 48.819 | 🟢 0/447 (0.00%) | 2026-11-04 23:49 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-11-05.md) ● | `gn_draw_stars` | 🟢 26.237 | 15.236 | 🟢 0/447 (0.00%) | 2026-11-04 22:58 |
| [PetalFlow](profile_petalflow_teensy_2026-11-05.md) ● | `pf_draw_rings` | 🟢 25.961 | 20.617 | 🟢 0/447 (0.00%) | 2026-11-04 23:04 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 25.370 | 43.803 | 🟢 0/927 (0.00%) | 2026-11-05 00:40 |
| [AlienCore](profile_aliencore_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 25.172 | 40.847 | 🟢 0/447 (0.00%) | 2026-11-04 23:50 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-11-05.md) § ● | `sh_rasterize` | 🟢 12.951 (24) | 21.277 | 🟢 0/1503 (0.00%) | 2026-11-04 23:24 |

Per-effect reports retain preset coverage and strict validator diagnostics. Per-frame render telemetry supports cadence figures without a claim of an unobserved wrap or fully validated cycle-counter exactness.

● marks shipping captures updated in the 2026-11-05 sweep or the DisplacementField re-capture on 2026-11-07.

§ Cyclers carry one aligned line per parser ownership colour bucket, red first;
(N) counts observed ownership entries. Buckets include transitions and exclude setup frame 1.
Show duration covers the authored show segment; a segment can end before its
last departure wraps. The per-effect report states observed coverage.

## Supplemental experimental presets

Historical fixed opt-in captures; these do not change the normal firmware
roster ranking. Triangular has since been removed.

| Preset | Dominant scope | Peak render ms | Spilled | Captured |
| --- | --- | ---: | ---: | --- |
| HyperLattice Triangular (supporting artifact removed) | `hl_shader_draw` | 🔴 134.620 | 🔴 441/441 (100%) | 2026-09-27 21:59 |
| HyperLattice Octet 3D (supporting artifact removed) | `hl_shader_draw` | 🔴 96.869 | 🔴 548/548 (100%) | 2026-09-27 22:50 |

See the main index's [historical fixed-preset notes](../README.md#historical-fixed-preset-captures), including the separate Cubic Wide Flight shipping check.

## Octet optimization supplement

[2026-09-27 Octet snapshot](profile_hyperlattice_octet_teensy_2026-09-27.md).

| Preset | Peak ms | Spilled | Captured |
|---|---:|---:|---|
| [Octet 4D (index 3 at capture time)](profile_hyperlattice_octet4d_teensy_2026-09-28.md) | 🔴 143.318 | 🔴 233/233 (100.00%) | 2026-09-28 09:53 |
| [Octet 3D (index 2 at capture time)](profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) | 🟢 45.688 | 🟢 0/1096 (0.00%) | 2026-09-28 09:56 |
| [Octet wide (index 4 at capture time)](profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) | 🟢 45.242 | 🟢 0/1096 (0.00%) | 2026-09-28 09:53 |

Updated 2026-09-28 09:53–09:56 after the Octet 58 ms Tiers 1–2 (`1bb522211`). Opt-in held presets, startup excluded; 4D on COM3, 3D on COM4. Earlier optimization and 4D evidence (supporting artifact removed).

## Shell flight supplement

Sphere-only renderer; two authored fixed presets on COM4, startup excluded. [Paired surface report](profile_hyperlattice_shell_flight_teensy_2026-09-29.md).

| Preset | Peak render ms | Spilled/live frames | Captured local |
|---|--:|--:|---|
| Shell Flight (index 5 at capture time) | 🔴 123.553 | 🔴 228/228 (100.00%) | 2026-09-29 08:03 |
| Shell Close Flight (index 6 at capture time) | 🔴 125.186 | 🔴 228/228 (100.00%) | 2026-09-29 08:06 |

## Earlier captures

Preserved architecture snapshots, superseded in the ranking by the current canonical captures.

- [MindSplatter 2026-10-01](profile_mindsplatter_architecture_teensy_2026-10-01.md).
- [HyperLattice 2026-10-01](profile_hyperlattice_architecture_teensy_2026-10-01.md).
- [ChromaticLichen 2026-10-01](profile_chromaticlichen_architecture_teensy_2026-10-01.md).
- [KaleidoscopeSmooth 2026-10-01](profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md).

Experimental clock characterization: [GnomonicStars bounded-channel comparison](profile_gnomonicstars_clock_experiment_teensy_2026-10-07.md). Historical paired measurements; the current canonical report measures the shipping bounded channel phases.

Integer direct-AA validation: [ShapeShifter paired capture](profile_shapeshifter_integer_aa_teensy_2026-10-10.md). Historical paired measurements at their recorded source snapshots.


## Temporal sampling comparison

[GridSpace matched A/B](profile_gridspace_teensy_2026-10-17.md) preserves the historical same-source temporal-on/off measurements. Current shipping timing is in the canonical report above.

## Palette-pattern recipe qualification

[GridSpace recipe overrides, 2026-10-29](profile_gridspace_palette_patterns_teensy_2026-10-29.md) qualify selected source, mechanic and evaluation combinations separately from the normal GridSpace shipping row.
