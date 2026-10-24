# On-device effect profiles — Teensy 4.0, segmented mode

On-device timing for the **36 effects in the Phantasm image**, captured on
bench-attached Teensy 4.0 boards with the real segmented POV driver
(`POVSegmented<288, 4, 480>`), DMA LEDs, and live flywheel/DMA ISRs. Each
effect renders one 288×144 image quadrant (about 10,368 pixels); the 62.5 ms
display window makes cadence quantize to 16, 8, 5.3 fps, and below.

Shipping reports cover all 36 effects: 36 completed captures.

The current Profile harness captures ISR counters and their timestamp together
with interrupts masked. ISR CPU shares use the separately logged `ISR window`
interval between successive snapshots, including report output between them.
The frame header's elapsed window excludes the preceding report output.

## Capture configurations

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

Rows rank by shipping spill fraction, then shipping steady peak render. Shipping peak and spilled exclude setup frame 1, listed separately; subsequent transitions remain included. Both peaks are render, never wall time; spilled counts frames exceeding 62.5 ms.
Shipping colours are binary: 🟢 zero spill; 🔴 any spill. O3 columns retain their
original values and colours and predate the per-column temporal work. O3 code-size
deltas and timestamps belong to each O3 capture's original source pair, including
MeshFeedback's adjacent-source pair and GSReactionDiffusion's pre-optimization pair.

| Effect | Dominant scope | Ship peak steady ms | Ship setup frame 1 ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---|--:|--:|--:|--:|--:|--:|--:|---|
| [MindSplatter](shipping/profile_mindsplatter_teensy_2026-10-23.md) / [O3](O3/profile_mindsplatter_architecture_teensy_2026-10-01.md) § | `msp_draw_particles` | 🟢 59.437 (8) | 0.379 | 🟢 54.539 (8) | 🟢 0/1247 (0.00%) | 🟢 0/1736 (0.00%) | +21,976 B | +17,200 B | ship 2026-10-23 13:40<br>O3 2026-10-01 09:22 |
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-10-23.md) / [O3](O3/profile_islamicstars_teensy_2026-09-28.md) § | `is_timeline_step` | 🟢 58.966 (16) | 16.644 | 🟢 49.112 (23) | 🟢 0/4863 (0.00%) | 🟢 0/3336 (0.00%) | +23,904 B | +8,240 B | ship 2026-10-23 13:57<br>O3 2026-09-28 16:50 |
| [Raymarch](shipping/profile_raymarch_teensy_2026-10-23.md) / [O3](O3/profile_raymarch_teensy_2026-09-20.md) | `rm_shader_draw` | 🟢 57.751 | 103.890 | 🟢 56.04 | 🟢 0/447 (0.00%) | 🟢 0/1736 (0.00%) | +12,568 B | +7,904 B | ship 2026-10-23 13:30<br>O3 2026-09-20 23:07 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-10-23.md) / [O3](O3/profile_shapeshifter_teensy_2026-09-25.md) § | `ss_draw_all` | 🟢 57.537 (9) | 47.865 | 🟢 60.05 (9) | 🟢 0/4271 (0.00%) | 🟢 0/2457 (0.00%) | +29,288 B | +24,496 B | ship 2026-10-23 17:21<br>O3 2026-09-25 07:29 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-10-23.md) / [O3](O3/profile_hopffibration_teensy_2026-08-26.md) | `hf_render_trails` | 🟢 56.256 | 0.696 | 🟢 46.52 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +19,616 B | +18,272 B | ship 2026-10-23 13:26<br>O3 2026-08-26 01:27 |
| [HyperLattice](shipping/profile_hyperlattice_teensy_2026-10-23.md) / [O3](O3/profile_hyperlattice_architecture_teensy_2026-10-01.md) § | `hl_shader_draw` | 🟢 55.035 (9) | 43.131 | 🟢 36.975 (3) | 🟢 0/4767 (0.00%) | 🟢 0/2696 (0.00%) | +9,616 B | +4,912 B | ship 2026-10-23 14:18<br>O3 2026-10-01 09:26 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-10-23.md) / [O3](O3/profile_dreamballs_teensy_2026-08-26.md) § | `db_timeline_step` | 🟢 54.311 (10) | 49.756 | 🟢 34.32 (11) | 🟢 0/3167 (0.00%) | 🟢 0/3648 (0%) | +26,896 B | +12,352 B | ship 2026-10-23 14:00<br>O3 2026-08-26 02:19 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-10-23.md) / [O3](O3/profile_ringspin_teensy_2026-09-24.md) | `rs_draw_rings` | 🟢 53.743 | 3.396 | 🟢 50.750 | 🟢 0/927 (0.00%) | 🟢 0/1087 (0.0%) | +15,400 B | +12,880 B | ship 2026-10-23 13:32<br>O3 2026-09-24 20:18 |
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-10-23.md) / [O3](O3/profile_meshfeedback_teensy_2026-08-26.md) § | `mf_feedback_flush` | 🟢 52.652 (12) | 5.659 | 🟢 58.32 (13) | 🟢 0/2863 (0.00%) | 🟢 0/6688 (0%) | +34,144 B | +21,776 B | ship 2026-10-23 13:51<br>O3 2026-08-26 02:10 |
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-10-23.md) / [O3](O3/profile_displacementfield_teensy_2026-09-19.md) § | `df_timeline_step` | 🟢 51.257 (2) | 18.953 | 🟢 57.45 | 🟢 0/1887 (0.00%) | 🟢 0/1088 (0%) | +25,856 B | +22,144 B | ship 2026-10-23 13:23<br>O3 2026-09-19 22:19 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-10-23.md) / [O3](O3/profile_bzreactiondiffusion_teensy_2026-08-26.md) | `bz_render` | 🟢 48.217 | 78.649 | 🟢 48.65 | 🟢 0/927 (0.00%) | 🟢 0/2048 (0%) | +12,760 B | +10,224 B | ship 2026-10-23 13:19<br>O3 2026-08-26 01:16 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-10-23.md) / [O3](O3/profile_gsreactiondiffusion_teensy_2026-09-30.md) | `grd_render` | 🟢 48.112 | 24.883 | 🔴 279.686 | 🟢 0/447 (0.00%) | 🔴 538/538 (100.00%) | +16,024 B | +10,832 B | ship 2026-10-23 13:25<br>O3 2026-09-30 19:15 |
| [AshCloud](shipping/profile_ashcloud_teensy_2026-10-23.md) / [O3](O3/profile_ashcloud_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 47.180 | 73.108 | 🔴 79.81 | 🟢 0/927 (0.00%) | 🔴 544/544 (100%) | +16,608 B | +12,112 B | ship 2026-10-23 14:21<br>O3 2026-08-26 02:45 |
| [Galaxies](shipping/profile_galaxies_teensy_2026-10-23.md) | `gx_particle_step` | 🟢 45.223 | 1.494 | — | 🟢 0/927 (0.00%) | — | — | — | ship 2026-10-23 13:17 |
| [KaleidoscopeHexOil](shipping/profile_kaleidoscopehexoil_teensy_2026-10-23.md) / [O3](O3/profile_kaleidoscopehexoil_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 43.705 (2) | 69.694 | 🟢 38.94 (3) | 🟢 0/1647 (0.00%) | 🟢 0/2208 (0%) | +13,456 B | +10,608 B | ship 2026-10-23 14:25<br>O3 2026-08-26 02:49 |
| [MermaidSkin](shipping/profile_mermaidskin_teensy_2026-10-23.md) / [O3](O3/profile_mermaidskin_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 42.073 | 70.074 | 🟢 54.55 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +16,544 B | +11,952 B | ship 2026-10-23 14:20<br>O3 2026-08-26 02:43 |
| [KaleidoscopeMandala](shipping/profile_kaleidoscopemandala_teensy_2026-10-23.md) / [O3](O3/profile_kaleidoscopemandala_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 40.428 (2) | 55.678 | 🟢 36.86 (3) | 🟢 0/1631 (0.00%) | 🟢 0/2368 (0%) | +13,656 B | +11,440 B | ship 2026-10-23 14:11<br>O3 2026-08-26 03:17 |
| [ChromaticLichen](shipping/profile_chromaticlichen_teensy_2026-10-23.md) / [O3](O3/profile_chromaticlichen_architecture_teensy_2026-10-01.md) | `fx_shader_draw` | 🟢 40.003 | 65.897 | 🟢 35.218 | 🟢 0/447 (0.00%) | 🟢 0/1096 (0.00%) | +15,264 B | +8,672 B | ship 2026-10-23 14:19<br>O3 2026-10-01 09:51 |
| [KaleidoscopeFlowers](shipping/profile_kaleidoscopeflowers_teensy_2026-10-23.md) / [O3](O3/profile_kaleidoscopeflowers_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 34.712 (3) | 53.438 | 🟢 33.38 (4) | 🟢 0/2735 (0.00%) | 🟢 0/4128 (0%) | +14,576 B | +11,712 B | ship 2026-10-23 14:35<br>O3 2026-08-26 03:03 |
| [KaleidoscopeSmooth](shipping/profile_kaleidoscopesmooth_teensy_2026-10-23.md) / [O3](O3/profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md) § | `fx_shader_draw` | 🟢 34.565 (4) | 53.018 | 🟢 30.314 (4) | 🟢 0/3807 (0.00%) | 🟢 0/4137 (0.00%) | +15,048 B | +8,624 B | ship 2026-10-23 14:29<br>O3 2026-10-01 11:19 |
| [KaleidoscopeHexBright](shipping/profile_kaleidoscopehexbright_teensy_2026-10-23.md) / [O3](O3/profile_kaleidoscopehexbright_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 33.137 (2) | 49.091 | 🟢 32.41 (3) | 🟢 0/1631 (0.00%) | 🟢 0/2368 (0%) | +14,520 B | +11,712 B | ship 2026-10-23 14:31<br>O3 2026-08-26 03:14 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-10-23.md) / [O3](O3/profile_hankinsolids_teensy_2026-09-28.md) § | `hk_timeline_step` | 🟢 31.644 (19) | 19.276 | 🟢 32.865 (18) | 🟢 0/2687 (0.00%) | 🟢 0/4137 (0.00%) | +17,088 B | +2,144 B | ship 2026-10-23 17:16<br>O3 2026-09-28 16:45 |
| [Fishbowl](shipping/profile_fishbowl_teensy_2026-10-23.md) / [O3](O3/profile_fishbowl_teensy_2026-08-26.md) | `fish_build_vertices` | 🟢 31.232 | 1.811 | 🟢 21.29 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +27,856 B | +26,400 B | ship 2026-10-23 13:20<br>O3 2026-08-26 01:18 |
| [Comets](shipping/profile_comets_teensy_2026-10-23.md) / [O3](O3/profile_comets_teensy_2026-08-26.md) § | `cm_draw_trail` | 🟢 30.701 (12) | 0.897 | 🟢 28.71 (13) | 🟢 0/1887 (0.00%) | 🟢 0/4128 (0%) | +15,432 B | +12,944 B | ship 2026-10-23 13:47<br>O3 2026-08-26 02:02 |
| [PetalFlow](shipping/profile_petalflow_teensy_2026-10-23.md) / [O3](O3/profile_petalflow_teensy_2026-08-26.md) | `pf_draw_rings` | 🟢 30.286 | 21.618 | 🟢 10.77 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +22,272 B | +21,008 B | ship 2026-10-23 13:29<br>O3 2026-08-26 01:32 |
| [KaleidoscopePentBright](shipping/profile_kaleidoscopepentbright_teensy_2026-10-23.md) / [O3](O3/profile_kaleidoscopepentbright_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 27.973 | 48.788 | 🟢 29.53 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +14,464 B | +11,712 B | ship 2026-10-23 14:23<br>O3 2026-08-26 02:46 |
| [MobiusGrid](shipping/profile_mobiusgrid_teensy_2026-10-23.md) / [O3](O3/profile_mobiusgrid_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 27.761 (2) | 42.938 | 🟢 21.44 (3) | 🟢 0/1647 (0.00%) | 🟢 0/2688 (0%) | +13,448 B | +10,560 B | ship 2026-10-23 13:28<br>O3 2026-08-26 01:30 |
| [AlienBrain](shipping/profile_alienbrain_teensy_2026-10-23.md) / [O3](O3/profile_alienbrain_teensy_2026-08-26.md) § | `fx_shader_draw` | 🟢 26.745 (4) | 48.180 | 🟢 27.90 (5) | 🟢 0/3807 (0.00%) | 🟢 0/4768 (0%) | +14,560 B | +11,712 B | ship 2026-10-23 14:05<br>O3 2026-08-26 02:24 |
| [KaleidoscopeHexSoft](shipping/profile_kaleidoscopehexsoft_teensy_2026-10-23.md) / [O3](O3/profile_kaleidoscopehexsoft_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 26.624 | 42.046 | 🟢 24.58 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +14,496 B | +11,712 B | ship 2026-10-23 14:06<br>O3 2026-08-26 02:26 |
| [CosmicEyeball](shipping/profile_cosmiceyeball_teensy_2026-10-23.md) / [O3](O3/profile_cosmiceyeball_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 26.390 | 48.097 | 🟢 22.97 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +12,552 B | +10,320 B | ship 2026-10-23 14:36<br>O3 2026-08-26 03:05 |
| [GnomonicStars](shipping/profile_gnomonicstars_teensy_2026-10-23.md) / [O3](O3/profile_gnomonicstars_teensy_2026-08-26.md) | `gn_draw_stars` | 🟢 26.228 | 15.374 | 🟢 26.29 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +12,184 B | +11,104 B | ship 2026-10-23 13:24<br>O3 2026-08-26 01:22 |
| [GridSpace](shipping/profile_gridspace_teensy_2026-10-23.md) / [O3](O3/profile_gridspace_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 25.495 | 38.937 | 🟢 25.39 | 🟢 0/927 (0.00%) | 🟢 0/1088 (0%) | +14,424 B | +11,712 B | ship 2026-10-23 14:12<br>O3 2026-08-26 02:33 |
| [AlienOcean](shipping/profile_alienocean_teensy_2026-10-23.md) / [O3](O3/profile_alienocean_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 24.092 | 44.614 | 🟢 24.98 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-10-23 14:07<br>O3 2026-08-26 02:28 |
| [AlienCore](shipping/profile_aliencore_teensy_2026-10-23.md) / [O3](O3/profile_aliencore_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 21.483 | 38.170 | 🟢 20.10 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | ship 2026-10-23 14:08<br>O3 2026-08-26 02:30 |
| [SphericalHarmonics](shipping/profile_sphericalharmonics_teensy_2026-10-23.md) / [O3](O3/profile_sphericalharmonics_teensy_2026-08-26.md) § | `sh_rasterize` | 🟢 13.296 (24) | 17.983 | 🟢 11.70 (24) | 🟢 0/1503 (0.00%) | 🟢 0/3488 (0%) | +8,016 B | +6,400 B | ship 2026-10-23 13:45<br>O3 2026-08-26 01:57 |
| [Voronoi](shipping/profile_voronoi_teensy_2026-10-23.md) / [O3](O3/profile_voronoi_teensy_2026-08-26.md) | `vo_shade` | 🟢 9.242 | 13.222 | 🟢 7.71 | 🟢 0/447 (0.00%) | 🟢 0/1088 (0%) | +15,568 B | +12,688 B | ship 2026-10-23 13:33<br>O3 2026-08-26 01:39 |

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
captured 2026-09-27 at 22:18; no 4D global-O3 capture.

Cubic Wide Flight (supporting artifact removed)
has a separate fixed-preset shipping check: peak 50.804 ms, 0/616 spills,
captured 2026-09-27 22:25.

## HyperLattice shipping segment

Current nine-preset shipping segment on COM4. The canonical
[report](shipping/profile_hyperlattice_teensy_2026-10-23.md#per-preset-table) gives the observed entry buckets
and validation diagnostics. The fixed-preset supplements below retain their
historical source snapshots and do not describe current shipping performance.

## Octet optimization supplement

Historical captures: [2026-09-27 shipping](shipping/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-27 global-O3](O3/profile_hyperlattice_octet_teensy_2026-09-27.md),
[2026-09-28 global-O3 preset 3](O3/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md),
and [2026-09-28 global-O3 preset 5](O3/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md).

| Preset | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Captured |
|---|---:|---:|---:|---:|---:|---:|---|
| [Octet 4D (index 4 at capture time)](shipping/profile_hyperlattice_octet4d_teensy_2026-09-28.md) | 🔴 143.318 | — | 🔴 233/233 (100.00%) | — | — | — | ship 2026-09-28 09:53 |
| [3: Octet 3D](shipping/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) | 🟢 45.688 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:56 |
| [Octet wide (index 5 at capture time)](shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) | 🟢 45.242 | — | 🟢 0/1096 (0.00%) | — | — | — | ship 2026-09-28 09:53 |

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

Only memory, shipping and global-O3 reports and their indexes are checked in.
Raw captures, source patches, build logs, ELF/map files and experimental campaign
artifacts are local outputs; they are no longer retained in this archive.
The profile wrapper writes these outputs under ignored `build/prof/`.
Historical hashes identify captures but do not guarantee reproducible source.

Capture `d25dd85dee17` maps to landed `43858e5fb`; only test headers differ,
with no firmware-relevant tree difference. Capture `0156d0d7490355` maps to
landed `3fc710350`. Capture `7baf3cc4307` maps to landed `daf63d812` on a
different base (319 differing files); it cannot be treated as the same tree.
Capture `17390a7a1` is unpublished: it is `e7f19b0d1` plus one profiling scope.
Capture `d5ca81403` landed as `c40c1def8`.
Capture `0df961b818ae` is unrecoverable. These historical orphan hashes are
not reachable from the published branch. The historical KaleidoscopeSmooth shipping
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

JewelMelt, KaleidoscopeMelt, KaleidoscopeStainedGlass and RingShower remain in the full roster. Their archived captures precede removal from the Phantasm image. ● marks captures taken after the 2026-10-06 re-capture.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [JewelMelt selective o3](excluded/profile_jewelmelt_selective_o3_teensy_2026-10-14.md) ● | `fx_shader_draw` | 🟢 46.35 | 🟢 0/1087 (0.00%) | 2026-10-14 20:16 |
| [KaleidoscopeMelt selective o3](excluded/profile_kaleidoscopemelt_selective_o3_teensy_2026-10-11.md) ● | `fx_shader_draw` | 🟢 46.26 | 🟢 0/1087 (0.00%) | 2026-10-11 11:38 |
| [KaleidoscopeStainedGlass selective o3](excluded/profile_kaleidoscopestainedglass_selective_o3_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 43.00 | 🟢 0/1087 (0.00%) | 2026-10-06 18:50 |
| [KaleidoscopeStainedGlass global o3](excluded/profile_kaleidoscopestainedglass_global_o3_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 46.99 | 🟢 0/1088 (0%) | 2026-08-26 02:51 |
| [RingShower selective o3](excluded/profile_ringshower_selective_o3_teensy_2026-10-06.md) | `rsh_draw_rings` | 🟢 4.32 | 🟢 0/1087 (0.00%) | 2026-10-06 18:20 |
| [RingShower global o3](excluded/profile_ringshower_global_o3_teensy_2026-08-26.md) | `rsh_draw_rings` | 🟢 3.86 | 🟢 0/1088 (0%) | 2026-08-26 01:36 |

Low-resolution-only (96×20) effects:

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [MobiusRings Hole A/B](excluded/profile_mobiusrings_teensy_2026-10-13.md) | `mr_draw_grid` | 🟢 12.30 | 🟢 0/1088 (0.00%) | 2026-10-13 |

## Temporal sampling comparison

[GridSpace matched A/B](shipping/profile_gridspace_teensy_2026-10-17.md) preserves the historical same-source temporal-on/off measurements cited by the rendering specification. Current shipping timing is in the canonical report above.
