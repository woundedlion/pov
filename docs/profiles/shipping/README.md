# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the
43 effects in `targets/Phantasm/phantasm_playlist.def`. All 43 have canonical reports. Steady peak is worst-frame render within the scheduled show slot, excluding setup frame 1; setup render is listed separately. Steady spilled counts frames whose render exceeded
the 62.5 ms display window. Rows rank by spill fraction, then steady peak render.
🟢 zero spill; 🔴 any spill.

The earlier [bidirectional Filigree capture](profile_filigree_bidirectional_teensy_2026-11-20.md) includes both transition directions.

| Effect | Dominant scope | Peak steady ms | Setup frame 1 ms | Spilled | Board | Show s | Captured |
|---|---|--:|--:|--:|---|--:|---|
| [IslamicStars](profile_islamicstars_teensy_2026-11-20.md) § ● | `is_timeline_step`<br>worst: dodecahedron_hk35_ambo_hk62_ambo_relax_hk42, 58.986 ms / 16 fps | 🟢 58.986 (16) | 16.611 | 🟢 0/4880 (0.00%) | COM3 | 306 | 2026-11-20 15:00 |
| [RiverStones](profile_riverstones_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: river-stones → river-stones-2, 58.139 ms / 16 fps | 🟢 58.139 (2) | 69.048 | 🟢 0/1670 (0.00%) | COM3 | 105 | 2026-11-20 15:27 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-11-20.md) § ● | `ss_draw_all`<br>worst: 8, 56.029 ms / 16 fps | 🟢 56.029 (9) | 76.265 | 🟢 0/4289 (0.00%) | COM4 | 269 | 2026-11-20 15:33 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-11-20.md) § ● | `mf_feedback_flush`<br>worst: WigglingWormhole, 55.879 ms / 16 fps | 🟢 55.879 (12) | 6.551 | 🟢 0/2886 (0.00%) | COM4 | 181 | 2026-11-20 15:40 |
| [AshCloud](profile_ashcloud_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 55.513 | 83.423 | 🟢 0/950 (0.00%) | COM3 | 60 | 2026-11-20 14:53 |
| [HopfFibration](profile_hopffibration_teensy_2026-11-20.md) ● | `hf_render_trails` | 🟢 55.490 | 0.672 | 🟢 0/945 (0.00%) | COM4 | 60 | 2026-11-20 16:15 |
| [Filigree](profile_filigree_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: jade-weave → ember-lattice, 55.155 ms / 16 fps | 🟢 55.155 (2) | 78.329 | 🟢 0/1666 (0.00%) | COM4 | 105 | 2026-11-20 15:45 |
| [Millefiori](profile_millefiori_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 53.861 | 82.441 | 🟢 0/950 (0.00%) | COM4 | 60 | 2026-11-20 16:04 |
| [JewelMelt](profile_jewelmelt_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 53.778 | 83.605 | 🟢 0/949 (0.00%) | COM4 | 60 | 2026-11-20 16:01 |
| [HyperLattice](profile_hyperlattice_teensy_2026-11-20.md) § ● | `hl_shader_draw`<br>worst: octet-4d-flight, 53.279 ms / 16 fps | 🟢 53.279 (9) | 45.645 | 🟢 0/4786 (0.00%) | COM4 | 300 | 2026-11-20 15:27 |
| [Raymarch](profile_raymarch_teensy_2026-11-20.md) ● | `rm_shader_draw` | 🟢 52.092 | 93.629 | 🟢 0/469 (0.00%) | COM4 | 30 | 2026-11-20 16:05 |
| [KaleidoscopeMelt](profile_kaleidoscopemelt_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 51.489 | 81.763 | 🟢 0/950 (0.00%) | COM4 | 60 | 2026-11-20 16:18 |
| [DreamBalls](profile_dreamballs_teensy_2026-11-20.md) § ● | `db_timeline_step`<br>worst: 9, 51.441 ms / 16 fps | 🟢 51.441 (10) | 51.219 | 🟢 0/3187 (0.00%) | COM4 | 200 | 2026-11-20 15:36 |
| [RingSpin](profile_ringspin_teensy_2026-11-20.md) ● | `rs_draw_rings` | 🟢 50.813 | 3.257 | 🟢 0/948 (0.00%) | COM4 | 60 | 2026-11-20 16:20 |
| [GridWarp](profile_gridwarp_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 50.732 | 78.390 | 🟢 0/949 (0.00%) | COM4 | 60 | 2026-11-20 16:14 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-11-20.md) ● | `bz_render` | 🟢 50.177 | 87.527 | 🟢 0/950 (0.00%) | COM4 | 60 | 2026-11-20 15:54 |
| [KaleidoscopeHexOil](profile_kaleidoscopehexoil_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: kaleidoscope-hex-oil → kaleidoscope-hex-oil-2, 49.700 ms / 16 fps | 🟢 49.700 (2) | 71.154 | 🟢 0/1668 (0.00%) | COM4 | 105 | 2026-11-20 15:47 |
| [KaleidoscopeMandala](profile_kaleidoscopemandala_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: wave-mirror → cup-hue, 47.556 ms / 16 fps | 🟢 47.556 (2) | 67.029 | 🟢 0/1668 (0.00%) | COM4 | 105 | 2026-11-20 16:07 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-11-20.md) ● | `grd_render` | 🟢 47.186 | 24.668 | 🟢 0/945 (0.00%) | COM4 | 60 | 2026-11-20 15:56 |
| [MermaidSkin](profile_mermaidskin_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 45.788 | 69.784 | 🟢 0/948 (0.00%) | COM4 | 60 | 2026-11-20 16:19 |
| [Iridescence](profile_iridescence_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 45.752 | 84.370 | 🟢 0/949 (0.00%) | COM4 | 60 | 2026-11-20 16:17 |
| [ChromaticLichen](profile_chromaticlichen_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 45.706 | 68.909 | 🟢 0/948 (0.00%) | COM4 | 60 | 2026-11-20 16:10 |
| [MindSplatter](profile_mindsplatter_teensy_2026-11-20.md) § ● | `msp_draw_particles`<br>worst: 6 → 7, 45.451 ms / 16 fps | 🟢 45.451 (8) | 0.423 | 🟢 0/1263 (0.00%) | COM4 | 80 | 2026-11-20 16:09 |
| [Involution](profile_involution_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 45.232 | 81.961 | 🟢 0/949 (0.00%) | COM4 | 60 | 2026-11-20 16:00 |
| [DisplacementField](profile_displacementfield_teensy_2026-11-20.md) § ● | `df_timeline_step`<br>worst: noise, 44.765 ms / 16 fps | 🟢 44.765 (2) | 18.584 | 🟢 0/1921 (0.00%) | COM3 | 121 | 2026-11-20 15:17 |
| [KaleidoscopeFlowers](profile_kaleidoscopeflowers_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: open-grid → fine-grid, 44.353 ms / 16 fps | 🟢 44.353 (3) | 63.854 | 🟢 0/2755 (0.00%) | COM3 | 173 | 2026-11-20 15:14 |
| [Galaxies](profile_galaxies_teensy_2026-11-20.md) ● | `gx_particle_step` | 🟢 43.659 | 1.470 | 🟢 0/949 (0.00%) | COM4 | 60 | 2026-11-20 16:13 |
| [KaleidoscopeHexBright](profile_kaleidoscopehexbright_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: hex-twin-wave, 42.290 ms / 16 fps | 🟢 42.290 (2) | 59.706 | 🟢 0/1667 (0.00%) | COM3 | 105 | 2026-11-20 15:22 |
| [LiquidWindows](profile_liquidwindows_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 42.228 | 64.841 | 🟢 0/947 (0.00%) | COM4 | 60 | 2026-11-20 17:20 |
| [KaleidoscopeSmooth](profile_kaleidoscopesmooth_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: coupled-grid, 41.435 ms / 16 fps | 🟢 41.435 (4) | 65.479 | 🟢 0/3825 (0.00%) | COM3 | 240 | 2026-11-20 15:10 |
| [GridWorms](profile_gridworms_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 38.612 | 52.447 | 🟢 0/945 (0.00%) | COM4 | 60 | 2026-11-20 15:59 |
| [AlienBrain](profile_alienbrain_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: alien-brain-2 → alien-brain-3, 37.554 ms / 16 fps | 🟢 37.554 (4) | 58.378 | 🟢 0/3825 (0.00%) | COM3 | 240 | 2026-11-20 15:05 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: mobius-grid → mobius-grid-2, 34.612 ms / 16 fps | 🟢 34.612 (2) | 59.131 | 🟢 0/1666 (0.00%) | COM4 | 105 | 2026-11-20 15:49 |
| [HankinSolids](profile_hankinsolids_teensy_2026-11-20.md) § ● | `hk_timeline_step`<br>worst: truncatedCuboctahedron, 32.517 ms / 16 fps | 🟢 32.517 (18) | 19.007 | 🟢 0/2703 (0.00%) | COM4 | 170 | 2026-11-20 15:43 |
| [GridSpace](profile_gridspace_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 32.296 | 45.993 | 🟢 0/943 (0.00%) | COM4 | 60 | 2026-11-20 15:21 |
| [Comets](profile_comets_teensy_2026-11-20.md) § ● | `cm_draw_trail`<br>worst: 5, 29.874 ms / 16 fps | 🟢 29.874 (12) | 0.879 | 🟢 0/1903 (0.00%) | COM3 | 120 | 2026-11-20 15:19 |
| [AlienOcean](profile_alienocean_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 28.276 | 51.217 | 🟢 0/463 (0.00%) | COM4 | 30 | 2026-11-20 16:05 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-11-20.md) ● | `gn_draw_stars` | 🟢 28.079 | 15.235 | 🟢 0/928 (0.00%) | COM4 | 60 | 2026-11-20 15:57 |
| [Fishbowl](profile_fishbowl_teensy_2026-11-20.md) ● | `fish_build_vertices` | 🟢 27.892 | 1.653 | 🟢 0/941 (0.00%) | COM4 | 60 | 2026-11-20 16:12 |
| [AlienCore](profile_aliencore_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 27.527 | 41.722 | 🟢 0/941 (0.00%) | COM4 | 60 | 2026-11-20 15:53 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-11-20.md) ● | `fx_shader_draw` | 🟢 27.513 | 47.003 | 🟢 0/943 (0.00%) | COM4 | 60 | 2026-11-20 15:55 |
| [PetalFlow](profile_petalflow_teensy_2026-11-20.md) ● | `pf_draw_rings` | 🟢 25.346 | 19.772 | 🟢 0/459 (0.00%) | COM4 | 30 | 2026-11-20 16:21 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-11-20.md) § ● | `sh_rasterize`<br>worst: morph toward l=4, m=4, 12.368 ms / 16 fps | 🟢 12.368 (24) | 20.591 | 🟢 0/1519 (0.00%) | COM4 | 96 | 2026-11-20 15:51 |

● The 2026-11-20 sweep measured 43 configured slots totaling 4,476 seconds at `6cca5fb7f`. LiquidWindows was remeasured 2026-11-20 17:20 at `e76722518`.

Shipping results use each configured show duration from the shipping firmware. The denominator counts frames confirmed complete by that elapsed-time boundary. Each report gives the possible additional boundary frames and their measured costs.
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
