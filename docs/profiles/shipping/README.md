# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
36 effects in `HS_PHANTASM_EFFECT_LIST`. All 36 have canonical reports, including
Galaxies. Peak is worst-frame render; spilled counts frames whose render exceeded
the 62.5 ms display window. Rows rank by spill fraction, then peak render.
🟢 zero spill; 🔴 any spill.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [IslamicStars](profile_islamicstars_teensy_2026-10-23.md) § | `is_timeline_step` | 🟢 61.337 (16 + stress) | 🟢 0/3862 (0.00%) | 2026-10-23 08:07 |
| [MindSplatter](profile_mindsplatter_teensy_2026-10-23.md) § | `msp_draw_particles` | 🟢 59.311 (8) | 🟢 0/1295 (0.00%) | 2026-10-23 08:36 |
| [Raymarch](profile_raymarch_teensy_2026-10-22.md) | `rm_shader_draw` | 🟢 57.674 | 🟢 0/447 (0.00%) | 2026-10-22 19:10 |
| [HopfFibration](profile_hopffibration_teensy_2026-10-23.md) | `hf_render_trails` | 🟢 56.986 | 🟢 0/927 (0.00%) | 2026-10-23 08:44 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-10-23.md) § | `ss_draw_all` | 🟢 56.941 (9 observed; no wrap) | 🟢 0/4287 (0.00%) | 2026-10-23 08:42 |
| [HyperLattice](profile_hyperlattice_teensy_2026-10-22.md) § | `hl_shader_draw` | 🟢 55.080 (9) | 🟢 0/4767 (0.00%) | 2026-10-22 20:01 |
| [DreamBalls](profile_dreamballs_teensy_2026-10-22.md) § | `db_timeline_step` | 🟢 54.443 (10) | 🟢 0/3167 (0.00%) | 2026-10-22 19:43 |
| [RingSpin](profile_ringspin_teensy_2026-10-22.md) | `rs_draw_rings` | 🟢 53.562 | 🟢 0/927 (0.00%) | 2026-10-22 19:12 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-10-22.md) § | `mf_feedback_flush` | 🟢 53.043 (12) | 🟢 0/2863 (0.00%) | 2026-10-22 19:34 |
| [DisplacementField](profile_displacementfield_teensy_2026-10-22.md) § | `df_timeline_step` | 🟢 51.350 (2) | 🟢 0/1887 (0.00%) | 2026-10-22 19:02 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-10-22.md) | `bz_render` | 🟢 48.271 | 🟢 0/927 (0.00%) | 2026-10-22 18:57 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-10-22.md) | `grd_render` | 🟢 48.119 | 🟢 0/447 (0.00%) | 2026-10-22 19:04 |
| [AshCloud](profile_ashcloud_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 46.831 | 🟢 0/927 (0.00%) | 2026-10-22 20:05 |
| [Galaxies](profile_galaxies_teensy_2026-10-23.md) | `gx_particle_step` | 🟢 45.496 | 🟢 0/927 (0.00%) | 2026-10-23 08:41 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-10-23.md) § | `fx_shader_draw` | 🟢 43.790 (scoped) | 🟢 0/416 (0.00%) | 2026-10-22 23:19 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-10-23.md) | `fx_shader_draw` | 🟢 42.116 (scoped) | 🟢 0/416 (0.00%) | 2026-10-22 23:18 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 40.095 | 🟢 0/447 (0.00%) | 2026-10-22 20:02 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 39.963 (2) | 🟢 0/1631 (0.00%) | 2026-10-22 19:54 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 33.982 (3) | 🟢 0/2735 (0.00%) | 2026-10-22 20:19 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 33.729 (2) | 🟢 0/1631 (0.00%) | 2026-10-22 20:16 |
| [HankinSolids](profile_hankinsolids_teensy_2026-10-22.md) § | `hk_timeline_step` | 🟢 33.483 (18) | 🟢 0/2687 (0.00%) | 2026-10-22 19:23 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 33.333 (4) | 🟢 0/3807 (0.00%) | 2026-10-22 20:14 |
| [Fishbowl](profile_fishbowl_teensy_2026-10-23.md) | `fish_build_vertices` | 🟢 31.100 | 🟢 0/927 (0.00%) | 2026-10-23 08:39 |
| [Comets](profile_comets_teensy_2026-10-22.md) § | `cm_draw_trail` | 🟢 30.550 (12) | 🟢 0/1887 (0.00%) | 2026-10-22 19:30 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 30.419 (2) | 🟢 0/1647 (0.00%) | 2026-10-22 19:08 |
| [PetalFlow](profile_petalflow_teensy_2026-10-23.md) | `pf_draw_rings` | 🟢 30.159 | 🟢 0/447 (0.00%) | 2026-10-23 08:46 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 28.921 | 🟢 0/927 (0.00%) | 2026-10-22 20:21 |
| [KaleidoscopePentBright](profile_kaleidoscopepentbright_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 28.723 | 🟢 0/927 (0.00%) | 2026-10-22 20:07 |
| [KaleidoscopeHexSoft](profile_kaleidoscopehexsoft_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 27.376 | 🟢 0/927 (0.00%) | 2026-10-22 19:50 |
| [AlienBrain](profile_alienbrain_teensy_2026-10-22.md) § | `fx_shader_draw` | 🟢 26.382 (4) | 🟢 0/3807 (0.00%) | 2026-10-22 19:48 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-10-22.md) | `gn_draw_stars` | 🟢 26.284 | 🟢 0/447 (0.00%) | 2026-10-22 19:03 |
| [GridSpace](profile_gridspace_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 25.957 | 🟢 0/927 (0.00%) | 2026-10-22 19:55 |
| [AlienOcean](profile_alienocean_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 25.517 | 🟢 0/447 (0.00%) | 2026-10-22 19:51 |
| [AlienCore](profile_aliencore_teensy_2026-10-22.md) | `fx_shader_draw` | 🟢 23.648 | 🟢 0/447 (0.00%) | 2026-10-22 19:52 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-10-22.md) § | `sh_rasterize` | 🟢 13.292 (24) | 🟢 0/1503 (0.00%) | 2026-10-22 19:28 |
| [Voronoi](profile_voronoi_teensy_2026-10-22.md) | `vo_shade` | 🟢 9.243 | 🟢 0/447 (0.00%) | 2026-10-22 19:13 |

The October 22 full-roster captures used COM4 at source
`aa61e787859bf31d69057b436b08ad8a751f1bb4`. Refreshed reports identify their
newer sources, capture scenarios and coverage individually. Scoped checks are
marked in the table; they do not inherit a preceding report's full-cycle coverage.
IslamicStars retains the higher fixed-shape stress peak alongside its complete
carousel. Captured times are local raw-log mtimes.

The six temporal-effect reports use accepted `a6222f23`/`70fe9492` captures
that predate the PointRoute equator-cell speed-bound correction in engine
revision `10df63d5432710f2ab22f1f726b7127ad0d6242e`. Their timing and size
figures remain captured-source measurements; the correction has no device timing
or performance claim.

The independently verified image at `f9c5c5efa6765489888265cc62524c68b3fe4e97`
has RAM1 code 171496 B, variables 312736 B, FLASH code 993776 B and FLASH data
851824 B. Its `.text.itcm` is 171488 B across six code banks, with 25112 B code
padding and no added pins. Against `4e75923c`, RAM1 code rises 2160 B
(2.109 KiB), variables are unchanged, FLASH data rises 160 B (0.156 KiB), and
FLASH code rises 77136 B (75.328 KiB). This image's footprint is separate from
the capture-specific sizes retained in each report.

This revision stores prepared face metadata and tables together in draw-scoped
arena storage; `SDF::Face` is an 8-byte generation-checked handle on 32-bit
targets. The device captures predate this storage change and retain their
captured-source timing measurements; this revision has no new device timing
claim.

Each report states its setup/warmup exclusion, complete-window denominator and
strict validator diagnostics. Root-cycle/wall mismatches remain explicit;
rankings use logged per-frame render telemetry and do not imply fully validated
cycle-counter exactness.

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
