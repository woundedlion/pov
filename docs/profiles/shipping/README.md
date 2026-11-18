# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
42 effects in `targets/Phantasm/phantasm_playlist.def`. All 42 have canonical reports. Steady peak is worst-frame render within the scheduled show slot, excluding setup frame 1; setup render is listed separately. Steady spilled counts frames whose render exceeded
the 62.5 ms display window. Rows rank by spill fraction, then steady peak render.
🟢 zero spill; 🔴 any spill.

| Effect | Dominant scope | Peak steady ms | Setup frame 1 ms | Spilled | Board | Captured |
|---|---|--:|--:|--:|---|---|
| [RiverStones](profile_riverstones_teensy_2026-11-18.md) ● | `fx_shader_draw` | 🟢 58.932 | 70.477 | 🟢 0/472 (0.00%) | COM4 | 2026-11-18 14:55 |
| [IslamicStars](profile_islamicstars_teensy_2026-11-17.md) § ● | `is_timeline_step`<br>worst: dodecahedron_hk35_ambo_hk62_ambo_relax_hk42, 58.059 ms / 16 fps | 🟢 58.059 (16) | 16.537 | 🟢 0/4895 (0.00%) | COM3 | 2026-11-17 09:09 |
| [GridWarp](profile_gridwarp_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 56.343 | 80.625 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 09:41 |
| [AshCloud](profile_ashcloud_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 56.216 | 85.379 | 🟢 0/959 (0.00%) | COM4 | 2026-11-17 09:39 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-11-17.md) § ● | `ss_draw_all`<br>worst: 8, 56.021 ms / 16 fps | 🟢 56.021 (9) | 76.786 | 🟢 0/4303 (0.00%) | COM3 | 2026-11-17 09:00 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-11-17.md) § ● | `mf_feedback_flush`<br>worst: WigglingWormhole, 55.672 ms / 16 fps | 🟢 55.672 (12) | 6.470 | 🟢 0/2895 (0.00%) | COM3 | 2026-11-17 09:23 |
| [HopfFibration](profile_hopffibration_teensy_2026-11-17.md) ● | `hf_render_trails` | 🟢 55.635 | 0.672 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 08:45 |
| [HyperLattice](profile_hyperlattice_teensy_2026-11-17.md) § ● | `hl_shader_draw`<br>worst: octet-4d-flight, 53.311 ms / 16 fps | 🟢 53.311 (9) | 45.039 | 🟢 0/4799 (0.00%) | COM4 | 2026-11-17 09:08 |
| [Millefiori](profile_millefiori_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 52.832 | 79.856 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 09:35 |
| [Raymarch](profile_raymarch_teensy_2026-11-17.md) ● | `rm_shader_draw` | 🟢 52.759 | 94.040 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 08:52 |
| [JewelMelt](profile_jewelmelt_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 52.547 | 79.906 | 🟢 0/959 (0.00%) | COM4 | 2026-11-17 09:45 |
| [DreamBalls](profile_dreamballs_teensy_2026-11-17.md) § ● | `db_timeline_step`<br>worst: 9, 51.363 ms / 16 fps | 🟢 51.363 (10) | 49.416 | 🟢 0/3199 (0.00%) | COM3 | 2026-11-17 09:19 |
| [RingSpin](profile_ringspin_teensy_2026-11-17.md) ● | `rs_draw_rings` | 🟢 51.117 | 3.279 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 08:54 |
| [MindSplatter](profile_mindsplatter_teensy_2026-11-17.md) § ● | `msp_draw_particles`<br>worst: 7, 50.810 ms / 16 fps | 🟢 50.810 (8) | 0.372 | 🟢 0/1279 (0.00%) | COM3 | 2026-11-17 09:03 |
| [KaleidoscopeMelt](profile_kaleidoscopemelt_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 50.074 | 78.699 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 09:45 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-11-17.md) § ● | `fx_shader_draw`<br>worst: Kaleidoscope Hex Oil 2, 48.354 ms / 16 fps | 🟢 48.354 (2) | 68.013 | 🟢 0/1679 (0.00%) | COM3 | 2026-11-17 09:26 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-11-17.md) ● | `grd_render` | 🟢 47.331 | 24.651 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 08:43 |
| [Iridescence](profile_iridescence_teensy_2026-11-18.md) ● | `fx_shader_draw` | 🟢 47.242 | 88.278 | 🟢 0/479 (0.00%) | COM3 | 2026-11-18 11:18 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 45.898 | 68.306 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 09:37 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-11-17.md) § ● | `fx_shader_draw`<br>worst: Cup Hue, 45.881 ms / 16 fps | 🟢 45.881 (2) | 65.217 | 🟢 0/1679 (0.00%) | COM4 | 2026-11-17 09:25 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 45.730 | 67.214 | 🟢 0/479 (0.00%) | COM4 | 2026-11-17 09:37 |
| [Involution](profile_involution_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 44.738 | 82.497 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 09:39 |
| [DisplacementField](profile_displacementfield_teensy_2026-11-17.md) § ● | `df_timeline_step`<br>worst: noise, 44.616 ms / 16 fps | 🟢 44.616 (2) | 18.450 | 🟢 0/1935 (0.00%) | COM3 | 2026-11-17 08:38 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-11-17.md) ● | `bz_render` | 🟢 44.329 | 76.068 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 08:32 |
| [Galaxies](profile_galaxies_teensy_2026-11-17.md) ● | `gx_particle_step` | 🟢 43.803 | 1.477 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 08:28 |
| [LiquidWindows](profile_liquidwindows_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 42.732 | 63.633 | 🟢 0/479 (0.00%) | COM4 | 2026-11-17 09:41 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-11-17.md) § ● | `fx_shader_draw`<br>worst: Fine Grid, 41.644 ms / 16 fps | 🟢 41.644 (3) | 61.934 | 🟢 0/2767 (0.00%) | COM4 | 2026-11-17 09:18 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-11-17.md) § ● | `fx_shader_draw`<br>worst: Hex Twin Wave, 41.633 ms / 16 fps | 🟢 41.633 (2) | 57.993 | 🟢 0/1679 (0.00%) | COM4 | 2026-11-17 09:29 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-11-17.md) § ● | `fx_shader_draw`<br>worst: Direct Grid, 40.288 ms / 16 fps | 🟢 40.288 (4) | 63.892 | 🟢 0/3839 (0.00%) | COM3 | 2026-11-17 09:14 |
| [GridWorms](profile_gridworms_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 39.075 | 54.282 | 🟢 0/479 (0.00%) | COM4 | 2026-11-17 09:35 |
| [AlienBrain](profile_alienbrain_teensy_2026-11-17.md) § ● | `fx_shader_draw`<br>worst: Alien Brain 3, 37.447 ms / 16 fps | 🟢 37.447 (4) | 57.420 | 🟢 0/3839 (0.00%) | COM4 | 2026-11-17 09:13 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-11-17.md) § ● | `fx_shader_draw`<br>worst: Mobius Grid 2, 33.522 ms / 16 fps | 🟢 33.522 (2) | 55.942 | 🟢 0/1679 (0.00%) | COM3 | 2026-11-17 08:48 |
| [HankinSolids](profile_hankinsolids_teensy_2026-11-17.md) § ● | `hk_timeline_step`<br>worst: truncatedCube, 32.980 ms / 16 fps | 🟢 32.980 (18) | 18.945 | 🟢 0/2719 (0.00%) | COM4 | 2026-11-17 09:22 |
| [GridSpace](profile_gridspace_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 32.794 | 46.539 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 09:33 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 30.734 | 48.727 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 09:43 |
| [Comets](profile_comets_teensy_2026-11-17.md) § ● | `cm_draw_trail`<br>worst: 5, 29.869 ms / 16 fps | 🟢 29.869 (12) | 0.869 | 🟢 0/1919 (0.00%) | COM3 | 2026-11-17 09:29 |
| [Fishbowl](profile_fishbowl_teensy_2026-11-17.md) ● | `fish_build_vertices` | 🟢 28.256 | 1.681 | 🟢 0/959 (0.00%) | COM3 | 2026-11-17 08:34 |
| [AlienCore](profile_aliencore_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 27.530 | 42.245 | 🟢 0/479 (0.00%) | COM4 | 2026-11-17 09:33 |
| [AlienOcean](profile_alienocean_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 27.306 | 49.224 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 09:31 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-11-17.md) ● | `gn_draw_stars` | 🟢 26.078 | 15.179 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 08:40 |
| [PetalFlow](profile_petalflow_teensy_2026-11-17.md) ● | `pf_draw_rings` | 🟢 25.420 | 19.712 | 🟢 0/479 (0.00%) | COM3 | 2026-11-17 08:50 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-11-17.md) § ● | `sh_rasterize`<br>worst: l=4, m=4, 12.841 ms / 16 fps | 🟢 12.841 (24) | 21.305 | 🟢 0/1535 (0.00%) | COM4 | 2026-11-17 09:31 |

● captures updated 2026-11-17, 2026-11-18.

Shipping results use only the authored show budget, including in-show transitions.
Setup frame 1 remains separate. Post-show capture padding is excluded from every
peak, spill denominator, phase readout and per-preset table.
A full cycle means the scheduled show slot; an additional wrap is not required.

§ Shipping cyclers carry one aligned line per preset colour bucket, red first;
(N) counts presets. Buckets include transitions and exclude setup frame 1.
Per-preset cost spans and hold-window selection are in the linked reports.

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
