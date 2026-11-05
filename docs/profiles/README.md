# On-device effect profiles — Teensy 4.0, segmented mode

On-device timing for the **39 effects in the Phantasm image**, captured on
bench-attached Teensy 4.0 boards with the real segmented POV driver
(`POVSegmented<288, 4, 480>`), DMA LEDs, and live flywheel/DMA ISRs. Each
effect renders one 288×144 image quadrant (about 10,368 pixels); the 62.5 ms
display window makes cadence quantize to 16, 8, 5.3 fps, and below.

Shipping reports cover all 39 effects: 39 completed captures.

Recorded capture revisions retain their original IDs. The
[source revision crosswalk](source_revisions.json) preserves each one-to-one
metadata rewrite in `current_commit`. For a published-branch source tree, use
`published_commit` when present, otherwise `current_commit`; these commits have
identical trees. Targets `013bf2db1b`, `449e8f9c5f` and `75bd9d5f4b`
resolve only from retained local branches. None has an identical-tree twin on
the published branch.

The current Profile harness captures ISR counters and their timestamp together
with interrupts masked. ISR CPU shares use the separately logged `ISR window`
interval between successive snapshots, including report output between them.
The frame header's elapsed window excludes the preceding report output.

## Capture configurations

The default profiling epoch is one hour. Images built before 2026-10-11
(`e70ec46de`) ignored `HS_PROFILE_EPOCH_REVS`, so older reports listing an
override ran on the default epoch.

Timing sets contain one canonical report per effect, with optional supplemental
reports identified by underscore-separated variant suffixes before `_teensy_`.
Each effect/variant pair has one dated report per set. Supplemental reports use
the same title, date, section, roster, and index checks as canonical reports;
they do not satisfy the required canonical shipping coverage. The linked
October 17 GridSpace comparison retains its original filename without a variant
suffix and is treated as the `temporal_comparison` supplement.

### Shipping selective-O3

The [`profile` report set](shipping/README.md) uses the shipping `-Os` image
with landed `HS_O3` hot-loop regions active.

### Global-O3 reference

The [`profile_o3` report set](O3/README.md) replaces `-Os` globally with
`-O3 -ffast-math`. It is a single-effect optimization ceiling, not a shippable
full-roster image.

## Shipping/O3 captures

● marks shipping captures updated in the 2026-11-05 sweep.

Rows rank by shipping spill fraction, then shipping steady peak render. Shipping peak and spilled exclude setup frame 1, listed separately; subsequent transitions remain included. Both peaks are render, never wall time; spilled counts frames exceeding 62.5 ms.
Shipping colours are binary: 🟢 zero spill; 🔴 any spill. O3 columns retain their
original values and colours and predate the per-column temporal work. O3 code-size
deltas and timestamps belong to each O3 capture's original source pair, including
MeshFeedback's adjacent-source pair and GSReactionDiffusion's pre-optimization pair.

| Effect | Dominant scope | Ship peak steady ms | Ship setup frame 1 ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---|--:|--:|--:|--:|--:|--:|--:|---|
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-11-05.md) / [O3](O3/profile_islamicstars_teensy_2026-09-28.md) § ● | `is_timeline_step` | 🟢 59.311 (16) | 16.664 | 🟢 49.112 (23) | 🟢 0/4863 (0.00%) | 🟢 0/3336 (0.00%) | +23,904 B | +8,240 B | ship 2026-11-04 23:39<br>O3 2026-09-28 16:50 |
| [MindSplatter](shipping/profile_mindsplatter_teensy_2026-11-05.md) / [O3](O3/profile_mindsplatter_architecture_teensy_2026-10-01.md) § ● | `msp_draw_particles` | 🟢 57.896 (8) | 0.389 | 🟢 54.539 (8) | 🟢 0/1247 (0.00%) | 🟢 0/1736 (0.00%) | +21,976 B | +17,200 B | ship 2026-11-04 23:16<br>O3 2026-10-01 09:22 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-11-05.md) / [O3](O3/profile_hopffibration_teensy_2026-08-26.md) ● | `hf_render_trails` | 🟢 57.199 | 0.646 | 🟢 46.52 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +19,616 B | +18,272 B | ship 2026-11-04 23:01<br>O3 2026-08-26 01:27 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-11-05.md) / [O3](O3/profile_shapeshifter_teensy_2026-09-25.md) § ● | `ss_draw_all` | 🟢 57.031 (9) | 78.200 | 🟢 60.05 (9) | 🟢 0/4271 (0.00%) | 🟢 0/2457 (0.00%) | +29,288 B | +24,496 B | ship 2026-11-04 23:13<br>O3 2026-09-25 07:29 |
| [AshCloud](shipping/profile_ashcloud_teensy_2026-11-05.md) / [O3](O3/profile_ashcloud_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 56.969 | 81.521 | 🔴 79.81 | 🟢 0/927 (0.00%) | 🔴 544/544 (100%) | — | — | ship 2026-11-05 00:10<br>O3 2026-08-26 02:45 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-11-05.md) / [O3](O3/profile_ringspin_teensy_2026-09-24.md) ● | `rs_draw_rings` | 🟢 53.875 | 3.399 | 🟢 50.750 | 🟢 0/927 (0.00%) | 🟢 0/1087 (0.0%) | +15,400 B | +12,880 B | ship 2026-11-04 23:07<br>O3 2026-09-24 20:18 |
| [HyperLattice](shipping/profile_hyperlattice_teensy_2026-11-05.md) / [O3](O3/profile_hyperlattice_architecture_teensy_2026-10-01.md) § ● | `hl_shader_draw` | 🟢 53.649 (9) | 43.112 | 🟢 36.975 (3) | 🟢 0/4767 (0.00%) | 🟢 0/2696 (0.00%) | +9,616 B | +4,912 B | ship 2026-11-05 00:03<br>O3 2026-10-01 09:26 |
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-11-05.md) / [O3](O3/profile_meshfeedback_teensy_2026-08-26.md) § ● | `mf_feedback_flush` | 🟢 53.505 (12) | 6.712 | 🟢 58.32 (13) | 🟢 0/2863 (0.00%) | 🟢 0/6688 (0%) | +34,144 B | +21,776 B | ship 2026-11-04 23:32<br>O3 2026-08-26 02:10 |
| [Raymarch](shipping/profile_raymarch_teensy_2026-11-05.md) / [O3](O3/profile_raymarch_teensy_2026-09-20.md) ● | `rm_shader_draw` | 🟢 52.693 | 93.551 | 🟢 56.04 | 🟢 0/447 (0.00%) | 🟢 0/1736 (0.00%) | +12,568 B | +7,904 B | ship 2026-11-04 23:06<br>O3 2026-09-20 23:07 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-11-05.md) / [O3](O3/profile_dreamballs_teensy_2026-08-26.md) § ● | `db_timeline_step` | 🟢 52.024 (10) | 52.869 | 🟢 34.32 (11) | 🟢 0/3167 (0.00%) | 🟢 0/3648 (0%) | — | — | ship 2026-11-04 23:43<br>O3 2026-08-26 02:19 |
| [Millefiori](shipping/profile_millefiori_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 51.758 | 78.463 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-11-04 23:58 |
| [GridWarp](shipping/profile_gridwarp_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 49.281 | 74.267 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-11-05 00:19 |
| [Iridescence](shipping/profile_iridescence_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 49.094 | 91.161 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-11-05 00:22 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-11-05.md) / [O3](O3/profile_gsreactiondiffusion_teensy_2026-09-30.md) ● | `grd_render` | 🟢 48.568 | 24.927 | 🔴 279.686 | 🟢 0/447 (0.00%) | 🔴 538/538 (100.00%) | +16,024 B | +10,832 B | ship 2026-11-04 22:59<br>O3 2026-09-30 19:15 |
| [Involution](shipping/profile_involution_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 47.347 | 83.449 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-11-05 00:15 |
| [Galaxies](shipping/profile_galaxies_teensy_2026-11-05.md) ● | `gx_particle_step` | 🟢 45.126 | 1.482 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-11-04 22:44 |
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-11-05.md) / [O3](O3/profile_displacementfield_teensy_2026-09-19.md) § ● | `df_timeline_step` | 🟢 44.911 (2) | 18.817 | 🟢 57.45 | 🟢 0/1887 (0.00%) | 🟢 0/1088 (0%) | +25,856 B | +22,144 B | ship 2026-11-04 22:56<br>O3 2026-09-19 22:19 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-11-05.md) / [O3](O3/profile_bzreactiondiffusion_teensy_2026-08-26.md) ● | `bz_render` | 🟢 44.574 | 76.467 | 🟢 48.65 | 🟢 0/927 (0.00%) | 🟢 0/2048 (0%) | +12,760 B | +10,224 B | ship 2026-11-04 22:48<br>O3 2026-08-26 01:16 |
| [KaleidoscopeHexOil](shipping/profile_kaleidoscopehexoil_teensy_2026-11-05.md) / [O3](O3/profile_kaleidoscopehexoil_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 42.739 (2) | 64.707 | 🟢 38.94 (3) | 🟢 0/1647 (0.00%) | 🟢 0/2208 (0%) | +13,456 B | +10,608 B | ship 2026-11-05 00:26<br>O3 2026-08-26 02:49 |
| [KaleidoscopeMandala](shipping/profile_kaleidoscopemandala_teensy_2026-11-05.md) / [O3](O3/profile_kaleidoscopemandala_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 42.737 (2) | 62.398 | 🟢 36.86 (3) | 🟢 0/1631 (0.00%) | 🟢 0/2368 (0%) | +13,656 B | +11,440 B | ship 2026-11-04 23:53<br>O3 2026-08-26 03:17 |
| [MermaidSkin](shipping/profile_mermaidskin_teensy_2026-11-05.md) / [O3](O3/profile_mermaidskin_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 41.164 | 64.523 | 🟢 54.55 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +16,544 B | +11,952 B | ship 2026-11-05 00:07<br>O3 2026-08-26 02:43 |
| [ChromaticLichen](shipping/profile_chromaticlichen_teensy_2026-11-05.md) / [O3](O3/profile_chromaticlichen_architecture_teensy_2026-10-01.md) ● | `fx_shader_draw` | 🟢 40.138 | 63.665 | 🟢 35.218 | 🟢 0/447 (0.00%) | 🟢 0/1096 (0.00%) | +15,264 B | +8,672 B | ship 2026-11-05 00:05<br>O3 2026-10-01 09:51 |
| [GridWorms](shipping/profile_gridworms_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 37.926 | 56.868 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-11-04 23:56 |
| [KaleidoscopeFlowers](shipping/profile_kaleidoscopeflowers_teensy_2026-11-05.md) / [O3](O3/profile_kaleidoscopeflowers_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 36.868 (3) | 57.037 | 🟢 33.38 (4) | 🟢 0/2735 (0.00%) | 🟢 0/4128 (0%) | +14,576 B | +11,712 B | ship 2026-11-05 00:39<br>O3 2026-08-26 03:03 |
| [KaleidoscopeHexBright](shipping/profile_kaleidoscopehexbright_teensy_2026-11-05.md) / [O3](O3/profile_kaleidoscopehexbright_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 36.811 (2) | 53.101 | 🟢 32.41 (3) | 🟢 0/1631 (0.00%) | 🟢 0/2368 (0%) | +14,520 B | +11,712 B | ship 2026-11-05 00:35<br>O3 2026-08-26 03:14 |
| [LiquidWindows](shipping/profile_liquidwindows_teensy_2026-11-05.md) ● | `fx_shader_draw` | 🟢 36.585 | 57.965 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-11-05 00:17 |
| [KaleidoscopeSmooth](shipping/profile_kaleidoscopesmooth_teensy_2026-11-05.md) / [O3](O3/profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md) § ● | `fx_shader_draw` | 🟢 36.376 (4) | 56.382 | 🟢 30.314 (4) | 🟢 0/3807 (0.00%) | 🟢 0/4137 (0.00%) | +15,048 B | +8,624 B | ship 2026-11-05 00:32<br>O3 2026-10-01 11:19 |
| [MobiusGrid](shipping/profile_mobiusgrid_teensy_2026-11-05.md) / [O3](O3/profile_mobiusgrid_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 32.874 (2) | 53.328 | 🟢 21.44 (3) | 🟢 0/1647 (0.00%) | 🟢 0/2688 (0%) | +13,448 B | +10,560 B | ship 2026-11-04 23:03<br>O3 2026-08-26 01:30 |
| [AlienBrain](shipping/profile_alienbrain_teensy_2026-11-05.md) / [O3](O3/profile_alienbrain_teensy_2026-08-26.md) § ● | `fx_shader_draw` | 🟢 32.584 (4) | 53.448 | 🟢 27.90 (5) | 🟢 0/3807 (0.00%) | 🟢 0/4768 (0%) | +14,560 B | +11,712 B | ship 2026-11-04 23:48<br>O3 2026-08-26 02:24 |
| [Comets](shipping/profile_comets_teensy_2026-11-05.md) / [O3](O3/profile_comets_teensy_2026-08-26.md) § ● | `cm_draw_trail` | 🟢 30.637 (12) | 0.860 | 🟢 28.71 (13) | 🟢 0/1887 (0.00%) | 🟢 0/4128 (0%) | +15,432 B | +12,944 B | ship 2026-11-04 23:28<br>O3 2026-08-26 02:02 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-11-05.md) / [O3](O3/profile_hankinsolids_teensy_2026-09-28.md) § ● | `hk_timeline_step` | 🟢 30.408 (19) | 19.120 | 🟢 32.865 (18) | 🟢 0/2687 (0.00%) | 🟢 0/4137 (0.00%) | +17,088 B | +2,144 B | ship 2026-11-04 23:22<br>O3 2026-09-28 16:45 |
| [Fishbowl](shipping/profile_fishbowl_teensy_2026-11-05.md) / [O3](O3/profile_fishbowl_teensy_2026-08-26.md) ● | `fish_build_vertices` | 🟢 28.932 | 1.755 | 🟢 21.29 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +27,856 B | +26,400 B | ship 2026-11-04 22:52<br>O3 2026-08-26 01:18 |
| [GridSpace](shipping/profile_gridspace_teensy_2026-11-05.md) / [O3](O3/profile_gridspace_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 28.401 | 41.152 | 🟢 25.39 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +14,424 B | +11,712 B | ship 2026-11-04 23:55<br>O3 2026-08-26 02:33 |
| [AlienOcean](shipping/profile_alienocean_teensy_2026-11-05.md) / [O3](O3/profile_alienocean_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 26.435 | 48.819 | 🟢 24.98 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-11-04 23:49<br>O3 2026-08-26 02:28 |
| [GnomonicStars](shipping/profile_gnomonicstars_teensy_2026-11-05.md) / [O3](O3/profile_gnomonicstars_teensy_2026-08-26.md) ● | `gn_draw_stars` | 🟢 26.237 | 15.236 | 🟢 26.29 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +12,184 B | +11,104 B | ship 2026-11-04 22:58<br>O3 2026-08-26 01:22 |
| [PetalFlow](shipping/profile_petalflow_teensy_2026-11-05.md) / [O3](O3/profile_petalflow_teensy_2026-08-26.md) ● | `pf_draw_rings` | 🟢 25.961 | 20.617 | 🟢 10.77 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +22,272 B | +21,008 B | ship 2026-11-04 23:04<br>O3 2026-08-26 01:32 |
| [CosmicEyeball](shipping/profile_cosmiceyeball_teensy_2026-11-05.md) / [O3](O3/profile_cosmiceyeball_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 25.370 | 43.803 | 🟢 22.97 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +12,552 B | +10,320 B | ship 2026-11-05 00:40<br>O3 2026-08-26 03:05 |
| [AlienCore](shipping/profile_aliencore_teensy_2026-11-05.md) / [O3](O3/profile_aliencore_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 25.172 | 40.847 | 🟢 20.10 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-11-04 23:50<br>O3 2026-08-26 02:30 |
| [SphericalHarmonics](shipping/profile_sphericalharmonics_teensy_2026-11-05.md) / [O3](O3/profile_sphericalharmonics_teensy_2026-08-26.md) § ● | `sh_rasterize` | 🟢 12.951 (24) | 21.277 | 🟢 11.70 (24) | 🟢 0/1503 (0.00%) | 🟢 0/3488 (0%) | +8,016 B | +6,400 B | ship 2026-11-04 23:24<br>O3 2026-08-26 01:57 |

Per-effect reports retain preset coverage and strict validator diagnostics. Per-frame render telemetry supports cadence figures without a claim of an unobserved wrap or fully validated cycle-counter exactness.

§ Cyclers carry one aligned line per parser ownership colour bucket, red first;
(N) counts observed ownership entries. Buckets include transitions and exclude setup frame 1.
Show duration covers the authored show segment; a segment can end before its
last departure wraps. The per-effect report states observed coverage.


## Memory captures

[Arena high-water snapshots (July 2026, host probe)](memory/arena_high_water.md) come from a host
probe and are independent of the on-device timing tables.

## Historical fixed-preset captures

Historical fixed opt-in captures; these do not change the normal firmware
roster ranking. Triangular has since been removed.

| Preset | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| HyperLattice Triangular (supporting artifact removed) | `hl_shader_draw` | 🔴 134.620 | 🔴 141.587 | 🔴 441/441 (100%) | 🔴 366/366 (100%) | +10,440 B | +8,208 B | ship 2026-09-27 21:59<br>O3 2026-09-27 22:01 |
| HyperLattice Octet 3D | `hl_shader_draw` | 🔴 96.869 | 🔴 95.710 | 🔴 548/548 (100%) | 🔴 548/548 (100%) | +11,688 B | +8,128 B | ship 2026-09-27 22:50<br>O3 2026-09-27 22:53 |

Octet 3D refreshed on 2026-09-27 after the single-owner strut correction.
These fixed-preset captures retain the original oscillating camera path.
The earlier Octet measurements (supporting artifact removed)
record 3D shipping peak 103.733 ms (548/548 spills) before that correction.
The 4D capture remains historical: 549.828 ms peak, 120/120 spills,
captured 2026-09-27 at 22:18, with no global-O3 twin. Later 4D global-O3
results appear in the [Octet optimization supplement](#octet-optimization-supplement).

Cubic Wide Flight (supporting artifact removed)
has a separate fixed-preset shipping check: peak 50.804 ms, 0/616 spills,
captured 2026-09-27 22:25.

## HyperLattice shipping segment

Current nine-preset shipping segment on COM4. The canonical
[report](shipping/profile_hyperlattice_teensy_2026-11-05.md#per-preset-table) gives the observed entry buckets
and validation diagnostics. The fixed-preset supplements below retain their
historical source snapshots and do not describe current shipping performance.

## Octet optimization supplement

Historical captures: [2026-09-27 shipping](shipping/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-27 global-O3](O3/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-28 global-O3 (index 2 at capture time)](O3/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md),
and [2026-09-28 global-O3 (index 4 at capture time)](O3/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md).

| Preset | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---:|---:|---:|---:|---:|---:|---|
| [Octet 4D (index 3 at capture time)](shipping/profile_hyperlattice_octet4d_teensy_2026-09-28.md) | 🔴 143.318 | — | 🔴 233/233 (100.00%) | — | — | — | ship 2026-09-28 09:53 |
| [Octet 3D (index 2 at capture time)](shipping/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) | 🟢 45.688 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:56 |
| [Octet wide (index 4 at capture time)](shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) | 🟢 45.242 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:53 |

Updated 2026-09-28 09:53–09:56 after the Octet 58 ms Tiers 1–2 (`1bb522211`). Opt-in held presets, startup excluded; 4D on COM3, 3D on COM4. Earlier optimization and 4D evidence (supporting artifact removed). No global-O3 twins or size deltas were captured for this code.

## Shell flight supplement

Sphere-only renderer; two authored fixed presets on COM4, startup excluded. [Shipping report](shipping/profile_hyperlattice_shell_flight_teensy_2026-09-29.md) and [global-O3 report](O3/profile_hyperlattice_shell_flight_teensy_2026-09-29.md).

| Preset | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured local |
|---|--:|--:|--:|--:|--:|--:|---|
| Shell Flight (index 5 at capture time) | 🔴 123.553 | 🔴 117.363 | 🔴 228/228 (100.00%) | 🔴 228/228 (100.00%) | +13,352 B | +4,736 B | ship 2026-09-29 08:03<br>O3 2026-09-29 08:04 |
| Shell Close Flight (index 6 at capture time) | 🔴 125.186 | 🔴 121.337 | 🔴 228/228 (100.00%) | 🔴 228/228 (100.00%) | +13,352 B | +4,736 B | ship 2026-09-29 08:06<br>O3 2026-09-29 08:07 |

## Earlier captures

Preserved snapshots: the 2026-10-01 shipping architecture captures (superseded by the
current canonical captures) and the canonical O3 captures that predate the paired architecture O3 rows.

- MindSplatter: [shipping architecture 2026-10-01](shipping/profile_mindsplatter_architecture_teensy_2026-10-01.md), [O3 2026-08-26](O3/profile_mindsplatter_teensy_2026-08-26.md).
- HyperLattice: [shipping architecture 2026-10-01](shipping/profile_hyperlattice_architecture_teensy_2026-10-01.md), [O3 2026-09-27](O3/profile_hyperlattice_teensy_2026-09-27.md).
- ChromaticLichen: [shipping architecture 2026-10-01](shipping/profile_chromaticlichen_architecture_teensy_2026-10-01.md), [O3 2026-08-26](O3/profile_chromaticlichen_teensy_2026-08-26.md).
- KaleidoscopeSmooth: [shipping architecture 2026-10-01](shipping/profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md), [O3 2026-08-26](O3/profile_kaleidoscopesmooth_teensy_2026-08-26.md).

## Artifact retention and source reachability

The shipping, global-O3 and excluded report sets and their indexes, the memory
snapshot, and the source revision crosswalk are checked in.
Raw captures, source patches, build logs, ELF/map files and experimental campaign
artifacts are local outputs; they are no longer retained in this archive.
The profile wrapper writes these outputs under ignored `build/prof/`.
Historical hashes identify captures but do not guarantee reproducible source.

The historical landed IDs below resolve to published-branch source trees through
[source_revisions.json](source_revisions.json), using `published_commit` when
present and otherwise `current_commit`.

Octet capture tip `1bb522211` is unpublished; landed `35f0348b7` has an
identical tree. The three Tier commits `29305dfe2`, `f130ee632` and
`1bb522211` landed patch-identical as `dac3ddde3`, `c88880d53` and
`35f0348b7`, respectively. Per-report capture tips retain their original
source identity.

Capture `d25dd85dee17` maps to landed `064cc88a5`; only test headers differ,
with no firmware-relevant tree difference. Capture `0156d0d7490355` maps to
landed `2b9b2a15c`. Capture `7baf3cc4307` maps to landed `35420268b` on a
different base (319 differing files); it cannot be treated as the same tree.
Capture `17390a7a1` is unpublished: it is `d8bfc08f7` plus one profiling scope.
The canonical shipping source `49bb199f4` is unpublished. Its landed
equivalent `b49317e21` differs in firmware: `walk_uniform_edge_columns()`
has 32-byte alignment in `core/render/plot/column_walk.h`, and
`HankinSolids::SCRATCH_A_BYTES` is 26 KiB instead of 24 KiB. It is not a
firmware-identical substitute for the captured source.
Capture `0df961b818ae` is unrecoverable. The unpublished capture hashes
are not reachable from the published branch. The historical KaleidoscopeSmooth shipping
capture additionally used an unretained source patch, so its recorded base
and patch hash cannot reconstruct the captured source from a fresh clone.

For the 2026-08-26 O3 cycler captures, the parser counted the unmarked
startup run as an extra entry 0. Their parenthesized entry counts are the
authored preset count plus one, except SphericalHarmonics, whose mode markers
start at 1 and whose entry count equals its authored preset count. Raw logs
are no longer available to rederive these buckets.

Experimental clock characterization: [GnomonicStars bounded-channel comparison](shipping/profile_gnomonicstars_clock_experiment_teensy_2026-10-07.md). Historical matched clock measurements at the recorded sources; the current canonical report measures shipping bounded channel phases.

Integer direct-AA validation: [ShapeShifter paired capture](shipping/profile_shapeshifter_integer_aa_teensy_2026-10-10.md). Both images preserve 16 fps throughout the cycle.

## Effects excluded from the Phantasm image

JewelMelt, KaleidoscopeMelt, KaleidoscopeStainedGlass and RingShower remain in the full roster. KaleidoscopeMelt’s current capture covers all four segment roles; the other rows retain their archived measurements. ● marks captures taken after the 2026-10-06 re-capture.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [JewelMelt selective o3](excluded/profile_jewelmelt_selective_o3_teensy_2026-11-03.md) ● | `fx_shader_draw` | 🔴 77.642 | 🔴 351/351 (100.00%) | 2026-11-03 22:26 |
| [KaleidoscopeMelt selective o3](excluded/profile_kaleidoscopemelt_selective_o3_teensy_2026-10-31.md) ● | `fx_shader_draw` | 🟢 55.668 | 🟢 0/8188 (0.00%) | 2026-10-31 03:55 |
| [KaleidoscopeStainedGlass selective o3](excluded/profile_kaleidoscopestainedglass_selective_o3_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 43.00 | 🟢 0/1087 (0.00%) | 2026-10-06 18:50 |
| [KaleidoscopeStainedGlass global o3](excluded/profile_kaleidoscopestainedglass_global_o3_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 46.99 | 🟢 0/1088 (0%) | 2026-08-26 02:51 |
| [RingShower selective o3](excluded/profile_ringshower_selective_o3_teensy_2026-10-06.md) | `rsh_draw_rings` | 🟢 4.32 | 🟢 0/1087 (0.00%) | 2026-10-06 18:20 |
| [RingShower global o3](excluded/profile_ringshower_global_o3_teensy_2026-08-26.md) | `rsh_draw_rings` | 🟢 3.86 | 🟢 0/1088 (0%) | 2026-08-26 01:36 |

KaleidoscopeMelt aggregates ongoing frames from roles 0–3; its report lists startup frame 1 and full-pass spills for every role.

Low-resolution-only (96×20) effects:

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [MobiusRings Hole A/B](excluded/profile_mobiusrings_teensy_2026-10-13.md) | `mr_draw_grid` | 🟢 12.30 | 🟢 0/1088 (0.00%) | 2026-10-13 |

## Temporal sampling comparison

[GridSpace matched A/B](shipping/profile_gridspace_teensy_2026-10-17.md) preserves the historical same-source temporal-on/off measurements cited by the rendering specification. Current shipping timing is in the canonical report above.

## Palette-pattern recipe qualification

[GridSpace recipe overrides, 2026-10-29](shipping/profile_gridspace_palette_patterns_teensy_2026-10-29.md) under the segmented driver. Complete counter windows include setup; red startup-only rows can meet the separate 250 ms handoff budget. Ordinary GridSpace shipping numbers above are unchanged.

| Recipe | Peak ms | Spilled | Captured local |
|---|--:|--:|---|
| Single noise / hue / baked | 🟢 58.221 | 🟢 0/1088 (0.00%) | 2026-10-29 11:42 |
| Directional sine / opacity / direct | 🟢 46.303 | 🟢 0/1088 (0.00%) | 2026-10-29 11:45 |
| Smooth grid / lightness lift / direct | 🟢 59.895 | 🟢 0/1088 (0.00%) | 2026-10-29 11:44 |
| Interference / hue / direct | 🟢 55.702 | 🟢 0/1088 (0.00%) | 2026-10-29 11:49 |
| Octave noise / lightness lift / direct | 🔴 109.018 | 🔴 352/352 (100.00%) | 2026-10-29 11:47 |
| Single noise / lightness lift / direct | 🔴 70.284 | 🔴 1/1088 (0.09%) | 2026-10-29 11:51 |
| Octave noise / lightness lift / baked | 🔴 76.126 | 🔴 271/544 (49.82%) | 2026-10-29 11:55 |
| Single noise / hue / direct | 🔴 66.274 | 🔴 1/1088 (0.09%) | 2026-10-29 11:58 |
| Octave noise / hue / baked | 🔴 64.235 | 🔴 1/1088 (0.09%) | 2026-10-29 12:01 |
| Octave noise / opacity / baked | 🟢 58.521 | 🟢 0/1088 (0.00%) | 2026-10-29 12:03 |
| Single noise / lightness lift / baked | 🔴 63.194 | 🔴 1/1088 (0.09%) | 2026-10-29 11:59 |
