# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
36 effects in `HS_PHANTASM_EFFECT_LIST`. All 36 have canonical reports: 35 completed captures and the failed
HankinSolids capture. Peak is worst-frame render; spilled counts frames whose render exceeded
the 62.5 ms display window. The failed capture appears first; measured rows rank by spill fraction, then peak render.
🟢 zero spill; 🔴 any spill.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [HankinSolids](profile_hankinsolids_teensy_2026-10-23.md) § | unavailable (arena trap) | 🔴 unavailable (trap) | 🔴 unavailable (trap) | 2026-10-23 (aborted) |
| [ShapeShifter](profile_shapeshifter_teensy_2026-10-23.md) § | `ss_draw_all` | 🔴 70.224 (2)<br>🟢 57.503 (7) | 🔴 134/958 (13.99%)<br>🟢 0/3042 (0.00%) | 2026-10-23 13:38 |
| [Raymarch](profile_raymarch_teensy_2026-10-23.md) | `rm_shader_draw` | 🔴 103.890 | 🔴 1/448 (0.22%) | 2026-10-23 13:30 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-10-23.md) | `fx_shader_draw` | 🔴 70.074 | 🔴 1/448 (0.22%) | 2026-10-23 14:20 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-10-23.md) | `fx_shader_draw` | 🔴 65.897 | 🔴 1/448 (0.22%) | 2026-10-23 14:19 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-10-23.md) | `bz_render` | 🔴 78.649 | 🔴 1/928 (0.11%) | 2026-10-23 13:19 |
| [AshCloud](profile_ashcloud_teensy_2026-10-23.md) | `fx_shader_draw` | 🔴 73.108 | 🔴 1/928 (0.11%) | 2026-10-23 14:21 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-10-23.md) § | `fx_shader_draw` | 🔴 69.694 (1)<br>🟢 43.705 (1) | 🔴 1/599 (0.17%)<br>🟢 0/1049 (0.00%) | 2026-10-23 14:25 |
| [MindSplatter](profile_mindsplatter_teensy_2026-10-23.md) § | `msp_draw_particles` | 🟢 59.437 (8) | 🟢 0/1248 (0.00%) | 2026-10-23 13:40 |
| [IslamicStars](profile_islamicstars_teensy_2026-10-23.md) § | `is_timeline_step` | 🟢 58.966 (16) | 🟢 0/4864 (0.00%) | 2026-10-23 13:57 |
| [HopfFibration](profile_hopffibration_teensy_2026-10-23.md) | `hf_render_trails` | 🟢 56.256 | 🟢 0/928 (0.00%) | 2026-10-23 13:26 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 55.678 (2) | 🟢 0/1632 (0.00%) | 2026-10-23 14:11 |
| [HyperLattice](profile_hyperlattice_teensy_2026-10-23.md) § | `hl_shader_draw` | 🟢 55.035 (9) | 🟢 0/4768 (0.00%) | 2026-10-23 14:18 |
| [DreamBalls](profile_dreamballs_teensy_2026-10-23.md) § | `db_timeline_step` | 🟢 54.311 (10) | 🟢 0/3168 (0.00%) | 2026-10-23 14:00 |
| [RingSpin](profile_ringspin_teensy_2026-10-23.md) | `rs_draw_rings` | 🟢 53.743 | 🟢 0/928 (0.00%) | 2026-10-23 13:32 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 53.438 (3) | 🟢 0/2736 (0.00%) | 2026-10-23 14:35 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 53.018 (4) | 🟢 0/3808 (0.00%) | 2026-10-23 14:29 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-10-23.md) § | `mf_feedback_flush` | 🟢 52.652 (12) | 🟢 0/2864 (0.00%) | 2026-10-23 13:51 |
| [DisplacementField](profile_displacementfield_teensy_2026-10-23.md) § | `df_timeline_step` | 🟢 51.257 (2) | 🟢 0/1888 (0.00%) | 2026-10-23 13:23 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 49.091 (2) | 🟢 0/1632 (0.00%) | 2026-10-23 14:31 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 48.788 | 🟢 0/928 (0.00%) | 2026-10-23 14:23 |
| [AlienBrain](profile_alienbrain_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 48.180 (4) | 🟢 0/3808 (0.00%) | 2026-10-23 14:05 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-10-23.md) | `grd_render` | 🟢 48.112 | 🟢 0/448 (0.00%) | 2026-10-23 13:25 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 48.097 | 🟢 0/928 (0.00%) | 2026-10-23 14:36 |
| [Galaxies](profile_galaxies_teensy_2026-10-23.md) | `gx_particle_step` | 🟢 45.223 | 🟢 0/928 (0.00%) | 2026-10-23 13:17 |
| [AlienOcean](profile_alienocean_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 44.614 | 🟢 0/448 (0.00%) | 2026-10-23 14:07 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 42.938 (2) | 🟢 0/1648 (0.00%) | 2026-10-23 13:28 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 42.046 | 🟢 0/928 (0.00%) | 2026-10-23 14:06 |
| [GridSpace](profile_gridspace_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 38.937 | 🟢 0/928 (0.00%) | 2026-10-23 14:12 |
| [AlienCore](profile_aliencore_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 38.170 | 🟢 0/448 (0.00%) | 2026-10-23 14:08 |
| [Fishbowl](profile_fishbowl_teensy_2026-10-23.md) | `fish_build_vertices` | 🟢 31.232 | 🟢 0/928 (0.00%) | 2026-10-23 13:20 |
| [Comets](profile_comets_teensy_2026-10-23.md) § | `cm_draw_trail` | 🟢 30.701 (12) | 🟢 0/1888 (0.00%) | 2026-10-23 13:47 |
| [PetalFlow](profile_petalflow_teensy_2026-10-23.md) | `pf_draw_rings` | 🟢 30.286 | 🟢 0/448 (0.00%) | 2026-10-23 13:29 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-10-23.md) | `gn_draw_stars` | 🟢 26.228 | 🟢 0/448 (0.00%) | 2026-10-23 13:24 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-10-23.md) § | `sh_rasterize` | 🟢 17.983 (24) | 🟢 0/1504 (0.00%) | 2026-10-23 13:45 |
| [Voronoi](profile_voronoi_teensy_2026-10-23.md) | `vo_shade` | 🟢 13.222 | 🟢 0/448 (0.00%) | 2026-10-23 13:33 |

All shipping rows use source `49bb199f40a655c23983e73ef6739d47ceb7908c` on COM4, with each capture
sized to its Phantasm show segment. Capture times are local raw-log mtimes.
The kaleidoscope effects run with column strobe off.
The failed capture's date identifies the attempt; its raw-log mtime is unavailable.

Shipping figures include setup and all frames in complete counter windows.
Each report separately gives trailing raw-frame telemetry, exclusions, preset
coverage and strict validator diagnostics. Root-cycle/wall mismatches remain
explicit; logged per-frame render supports cadence measurements without a claim
of fully validated cycle-counter exactness.

§ Cyclers carry one aligned line per parser ownership colour bucket, red first;
(N) counts observed ownership entries. Buckets include setup and transitions.
Failed captures appear first with unavailable timing. Show duration covers the authored show segment; a segment can end before its
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
