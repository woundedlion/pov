# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
36 effects in `HS_PHANTASM_EFFECT_LIST`. All 36 have canonical reports, including
Galaxies. Peak is worst-frame render; spilled counts frames whose render exceeded
the 62.5 ms display window. Rows rank by spill fraction, then peak render.
🟢 zero spill; 🔴 any spill.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [IslamicStars](profile_islamicstars_teensy_2026-10-22.md) § | `is_timeline_step` | 🔴 67.393 (1)<br>🟢 57.235 (15) | 🔴 33/352 (9.38%)<br>🟢 0/4447 (0.00%) | 2026-10-22 19:39 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 61.132 | 🟢 0/447 (0.00%) | 2026-10-22 20:03 |
| [MindSplatter](profile_mindsplatter_teensy_2026-10-22.md) § | `msp_draw_particles` | 🟢 58.709 (8) | 🟢 0/1247 (0.00%) | 2026-10-22 19:20 |
| [Raymarch](profile_raymarch_teensy_2026-10-22.md) | `rm_shader_draw` | 🟢 57.674 | 🟢 0/447 (0.00%) | 2026-10-22 19:10 |
| [HopfFibration](profile_hopffibration_teensy_2026-10-22.md) | `hf_render_trails` | 🟢 57.358 | 🟢 0/927 (0.00%) | 2026-10-22 19:05 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-10-22.md) § | `ss_draw_all` | 🟢 55.606 (9) | 🟢 0/4271 (0.00%) | 2026-10-22 19:18 |
| [HyperLattice](profile_hyperlattice_teensy_2026-10-22.md) § | `hl_shader_draw` | 🟢 55.080 (9) | 🟢 0/4767 (0.00%) | 2026-10-22 20:01 |
| [DreamBalls](profile_dreamballs_teensy_2026-10-22.md) § | `db_timeline_step` | 🟢 54.443 (10) | 🟢 0/3167 (0.00%) | 2026-10-22 19:43 |
| [RingSpin](profile_ringspin_teensy_2026-10-22.md) | `rs_draw_rings` | 🟢 53.562 | 🟢 0/927 (0.00%) | 2026-10-22 19:12 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-10-22.md) § | `mf_feedback_flush` | 🟢 53.043 (12) | 🟢 0/2863 (0.00%) | 2026-10-22 19:34 |
| [DisplacementField](profile_displacementfield_teensy_2026-10-22.md) § | `df_timeline_step` | 🟢 51.350 (2) | 🟢 0/1887 (0.00%) | 2026-10-22 19:02 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 49.679 (2) | 🟢 0/1647 (0.00%) | 2026-10-22 20:09 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-10-22.md) | `bz_render` | 🟢 48.271 | 🟢 0/927 (0.00%) | 2026-10-22 18:57 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-10-22.md) | `grd_render` | 🟢 48.119 | 🟢 0/447 (0.00%) | 2026-10-22 19:04 |
| [AshCloud](profile_ashcloud_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 46.831 | 🟢 0/927 (0.00%) | 2026-10-22 20:05 |
| [Galaxies](profile_galaxies_teensy_2026-10-22.md) | `gx_particle_step` | 🟢 44.932 | 🟢 0/927 (0.00%) | 2026-10-22 18:55 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 40.095 | 🟢 0/447 (0.00%) | 2026-10-22 20:02 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 39.963 (2) | 🟢 0/1631 (0.00%) | 2026-10-22 19:54 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 33.982 (3) | 🟢 0/2735 (0.00%) | 2026-10-22 20:19 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 33.729 (2) | 🟢 0/1631 (0.00%) | 2026-10-22 20:16 |
| [HankinSolids](profile_hankinsolids_teensy_2026-10-22.md) § | `hk_timeline_step` | 🟢 33.483 (18) | 🟢 0/2687 (0.00%) | 2026-10-22 19:23 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 33.333 (4) | 🟢 0/3807 (0.00%) | 2026-10-22 20:14 |
| [Comets](profile_comets_teensy_2026-10-22.md) § | `cm_draw_trail` | 🟢 30.550 (12) | 🟢 0/1887 (0.00%) | 2026-10-22 19:30 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 30.419 (2) | 🟢 0/1647 (0.00%) | 2026-10-22 19:08 |
| [PetalFlow](profile_petalflow_teensy_2026-10-22.md) | `pf_draw_rings` | 🟢 30.379 | 🟢 0/447 (0.00%) | 2026-10-22 19:09 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 28.921 | 🟢 0/927 (0.00%) | 2026-10-22 20:21 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 28.723 | 🟢 0/927 (0.00%) | 2026-10-22 20:07 |
| [Fishbowl](profile_fishbowl_teensy_2026-10-22.md) | `fish_build_vertices` | 🟢 28.480 | 🟢 0/927 (0.00%) | 2026-10-22 18:59 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 27.376 | 🟢 0/927 (0.00%) | 2026-10-22 19:50 |
| [AlienBrain](profile_alienbrain_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 26.382 (4) | 🟢 0/3807 (0.00%) | 2026-10-22 19:48 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-10-22.md) | `gn_draw_stars` | 🟢 26.284 | 🟢 0/447 (0.00%) | 2026-10-22 19:03 |
| [GridSpace](profile_gridspace_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 25.957 | 🟢 0/927 (0.00%) | 2026-10-22 19:55 |
| [AlienOcean](profile_alienocean_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 25.517 | 🟢 0/447 (0.00%) | 2026-10-22 19:51 |
| [AlienCore](profile_aliencore_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 23.648 | 🟢 0/447 (0.00%) | 2026-10-22 19:52 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-10-22.md) § | `sh_rasterize` | 🟢 13.292 (24) | 🟢 0/1503 (0.00%) | 2026-10-22 19:28 |
| [Voronoi](profile_voronoi_teensy_2026-10-22.md) | `vo_shade` | 🟢 9.243 | 🟢 0/447 (0.00%) | 2026-10-22 19:13 |

All 36 canonical shipping captures were refreshed on COM4 at source
`aa61e787859bf31d69057b436b08ad8a751f1bb4`. Captures use the current Phantasm
show-segment durations and shipping parameter defaults. Setup frame 1 is excluded
from peak and spill figures. Captured times are local raw-log mtimes.

The root-cycle/wall-sum cross-check exceeds the parser's 5 ppm limit in all 36
captures; each report retains its measured discrepancy. Effect identity, frame
monotonicity and complete-window telemetry checks pass. Rankings use the logged
per-frame render values.

§ Cyclers carry one aligned line per parser-owned colour bucket, red first;
(N) counts ownership entries. Bucket frames include the following transition.
Each report gives per-entry timing, complete-window coverage and parser validation
diagnostics. A segment can end before its final transition wraps, and spills can
limit the entries observed; the tables do not imply full-cycle validation.

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

Preserved architecture snapshots, superseded in the ranking by the current canonical captures.

- [MindSplatter 2026-10-01](profile_mindsplatter_architecture_teensy_2026-10-01.md).
- [HyperLattice 2026-10-01](profile_hyperlattice_architecture_teensy_2026-10-01.md).
- [ChromaticLichen 2026-10-01](profile_chromaticlichen_architecture_teensy_2026-10-01.md).
- [KaleidoscopeSmooth 2026-10-01](profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md).

Experimental clock characterization: [GnomonicStars bounded-channel comparison](profile_gnomonicstars_clock_experiment_teensy_2026-10-07.md). Historical paired measurements; the current canonical report measures the shipping bounded channel phases.

Integer direct-AA validation: [ShapeShifter paired capture](profile_shapeshifter_integer_aa_teensy_2026-10-10.md). Historical paired measurements at their recorded source snapshots.


## Temporal sampling comparison

[GridSpace matched A/B](profile_gridspace_teensy_2026-10-17.md) preserves the historical same-source temporal-on/off measurements cited by the rendering specification. Current shipping timing is in the canonical report above.
