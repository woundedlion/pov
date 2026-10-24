# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
36 effects in `HS_PHANTASM_EFFECT_LIST`. All 36 have canonical reports: 36 completed captures. Steady peak is worst-frame render excluding setup frame 1; setup render is listed separately. Steady spilled counts frames whose render exceeded
the 62.5 ms display window. Rows rank by spill fraction, then steady peak render.
🟢 zero spill; 🔴 any spill.

| Effect | Dominant scope | Peak steady ms | Setup frame 1 ms | Spilled | Captured |
|---|---|--:|--:|--:|---|
| [MindSplatter](profile_mindsplatter_teensy_2026-10-23.md) § | `msp_draw_particles` | 🟢 59.437 (8) | 0.379 | 🟢 0/1247 (0.00%) | 2026-10-23 13:40 |
| [IslamicStars](profile_islamicstars_teensy_2026-10-23.md) § | `is_timeline_step` | 🟢 58.966 (16) | 16.644 | 🟢 0/4863 (0.00%) | 2026-10-23 13:57 |
| [Raymarch](profile_raymarch_teensy_2026-10-23.md) | `rm_shader_draw` | 🟢 57.751 | 103.890 | 🟢 0/447 (0.00%) | 2026-10-23 13:30 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-10-23.md) § | `ss_draw_all` | 🟢 57.537 (9) | 47.865 | 🟢 0/4271 (0.00%) | 2026-10-23 17:21 |
| [HopfFibration](profile_hopffibration_teensy_2026-10-23.md) | `hf_render_trails` | 🟢 56.256 | 0.696 | 🟢 0/927 (0.00%) | 2026-10-23 13:26 |
| [HyperLattice](profile_hyperlattice_teensy_2026-10-23.md) § | `hl_shader_draw` | 🟢 55.035 (9) | 43.131 | 🟢 0/4767 (0.00%) | 2026-10-23 14:18 |
| [DreamBalls](profile_dreamballs_teensy_2026-10-23.md) § | `db_timeline_step` | 🟢 54.311 (10) | 49.756 | 🟢 0/3167 (0.00%) | 2026-10-23 14:00 |
| [RingSpin](profile_ringspin_teensy_2026-10-23.md) | `rs_draw_rings` | 🟢 53.743 | 3.396 | 🟢 0/927 (0.00%) | 2026-10-23 13:32 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-10-23.md) § | `mf_feedback_flush` | 🟢 52.652 (12) | 5.659 | 🟢 0/2863 (0.00%) | 2026-10-23 13:51 |
| [DisplacementField](profile_displacementfield_teensy_2026-10-23.md) § | `df_timeline_step` | 🟢 51.257 (2) | 18.953 | 🟢 0/1887 (0.00%) | 2026-10-23 13:23 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-10-23.md) | `bz_render` | 🟢 48.217 | 78.649 | 🟢 0/927 (0.00%) | 2026-10-23 13:19 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-10-23.md) | `grd_render` | 🟢 48.112 | 24.883 | 🟢 0/447 (0.00%) | 2026-10-23 13:25 |
| [AshCloud](profile_ashcloud_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 47.180 | 73.108 | 🟢 0/927 (0.00%) | 2026-10-23 14:21 |
| [Galaxies](profile_galaxies_teensy_2026-10-23.md) | `gx_particle_step` | 🟢 45.223 | 1.494 | 🟢 0/927 (0.00%) | 2026-10-23 13:17 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 43.705 (2) | 69.694 | 🟢 0/1647 (0.00%) | 2026-10-23 14:25 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 42.073 | 70.074 | 🟢 0/447 (0.00%) | 2026-10-23 14:20 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 40.428 (2) | 55.678 | 🟢 0/1631 (0.00%) | 2026-10-23 14:11 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 40.003 | 65.897 | 🟢 0/447 (0.00%) | 2026-10-23 14:19 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 34.712 (3) | 53.438 | 🟢 0/2735 (0.00%) | 2026-10-23 14:35 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 34.565 (4) | 53.018 | 🟢 0/3807 (0.00%) | 2026-10-23 14:29 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 33.137 (2) | 49.091 | 🟢 0/1631 (0.00%) | 2026-10-23 14:31 |
| [HankinSolids](profile_hankinsolids_teensy_2026-10-23.md) § | `hk_timeline_step` | 🟢 31.644 (19) | 19.276 | 🟢 0/2687 (0.00%) | 2026-10-23 17:16 |
| [Fishbowl](profile_fishbowl_teensy_2026-10-23.md) | `fish_build_vertices` | 🟢 31.232 | 1.811 | 🟢 0/927 (0.00%) | 2026-10-23 13:20 |
| [Comets](profile_comets_teensy_2026-10-23.md) § | `cm_draw_trail` | 🟢 30.701 (12) | 0.897 | 🟢 0/1887 (0.00%) | 2026-10-23 13:47 |
| [PetalFlow](profile_petalflow_teensy_2026-10-23.md) | `pf_draw_rings` | 🟢 30.286 | 21.618 | 🟢 0/447 (0.00%) | 2026-10-23 13:29 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 27.973 | 48.788 | 🟢 0/927 (0.00%) | 2026-10-23 14:23 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 27.761 (2) | 42.938 | 🟢 0/1647 (0.00%) | 2026-10-23 13:28 |
| [AlienBrain](profile_alienbrain_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 26.745 (4) | 48.180 | 🟢 0/3807 (0.00%) | 2026-10-23 14:05 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 26.624 | 42.046 | 🟢 0/927 (0.00%) | 2026-10-23 14:06 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 26.390 | 48.097 | 🟢 0/927 (0.00%) | 2026-10-23 14:36 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-10-23.md) | `gn_draw_stars` | 🟢 26.228 | 15.374 | 🟢 0/447 (0.00%) | 2026-10-23 13:24 |
| [GridSpace](profile_gridspace_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 25.495 | 38.937 | 🟢 0/927 (0.00%) | 2026-10-23 14:12 |
| [AlienOcean](profile_alienocean_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 24.092 | 44.614 | 🟢 0/447 (0.00%) | 2026-10-23 14:07 |
| [AlienCore](profile_aliencore_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 21.483 | 38.170 | 🟢 0/447 (0.00%) | 2026-10-23 14:08 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-10-23.md) § | `sh_rasterize` | 🟢 13.296 (24) | 17.983 | 🟢 0/1503 (0.00%) | 2026-10-23 13:45 |
| [Voronoi](profile_voronoi_teensy_2026-10-23.md) | `vo_shade` | 🟢 9.242 | 13.222 | 🟢 0/447 (0.00%) | 2026-10-23 13:33 |

Shipping rows use source `49bb199f40a655c23983e73ef6739d47ceb7908c` on COM4, except HankinSolids and ShapeShifter, re-captured on COM4 at `bda58b65d3fc317f5d6ccbe29888014738a0ee25`. Each capture uses its Phantasm show-segment budget. Capture times are local raw-log mtimes.
The kaleidoscope effects run with column strobe off.

Shipping steady figures use frames 2 through the last complete counter window. Setup frame 1 is separate. The measured scope is per-frame render telemetry (`frame` minus buffer wait), under the shipping segment clip.
Each report separately gives trailing raw-frame telemetry, exclusions, preset
coverage and strict validator diagnostics. Root-cycle/wall mismatches remain
explicit; logged per-frame render supports cadence measurements without a claim
of fully validated cycle-counter exactness.

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
