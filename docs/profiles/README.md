# On-device effect profiles — Teensy 4.0, segmented mode

On-device timing for the **43 effects in the Phantasm image**, captured on
bench-attached Teensy 4.0 boards with the real segmented POV driver
(`POVSegmented<288, 4, 480>`), DMA LEDs, and live flywheel/DMA ISRs. Each display window uses one 288×144 image quadrant (about 10,368 pixels).
MeshFeedback renders the full 41,472-pixel canvas for feedback. The 62.5 ms
display window makes cadence quantize to 16, 8, 5.3 fps, and below.

Shipping reports cover all 43 effects with 43 completed captures.

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

The archive checker validates setup-frame values against each report's
`README setup frame 1` value and capture timestamps against retained capture
metadata. Six historical reports lack timestamp metadata:
[ShapeShifter O3](O3/profile_shapeshifter_teensy_2026-09-25.md),
[Raymarch O3](O3/profile_raymarch_teensy_2026-09-20.md),
[Octet 4D](shipping/profile_hyperlattice_octet4d_teensy_2026-09-28.md),
[Octet 3D](shipping/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md),
[Octet wide](shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md),
and [MobiusRings](excluded/profile_mobiusrings_teensy_2026-10-13.md).
For these fixed historical report paths only, capture cells receive
calendar/time format checks. Their timestamp cells remain historical index
metadata rather than values verified against the report or raw capture.

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

The [`profile_o3` report set](O3/README.md) replaces `-Os` with global `-O3`;
both images share the base `-ffast-math -fno-finite-math-only` flags. It is a
single-effect optimization ceiling, not a shippable full-roster image.

## Shipping/O3 captures

Rows rank by in-show shipping spill fraction, then in-show live peak render. Shipping peak and spilled exclude setup frame 1, listed separately; subsequent transitions remain included. Both peaks are render, never wall time; spilled counts frames exceeding 62.5 ms.
Shipping colours are binary: 🟢 zero spill; 🔴 any spill. IslamicStars O3 is a fresh matched capture of the same 306-second show slot.
Other O3 columns retain their original historical acquisition scopes, values
and colours; they have not been reanalyzed to the current show scope. O3 code-size
deltas and timestamps belong to each O3 capture's original source pair, including
MeshFeedback's adjacent-source pair and GSReactionDiffusion's pre-optimization pair.

| Effect | Dominant scope | Ship peak steady ms | Ship setup frame 1 ms | O3 peak ms | Ship spilled | O3 spilled | FLASH code Δ | ITCM Δ | Ship board | Captured |
|---|---|--:|--:|--:|--:|--:|--:|--:|---|---|
| [RiverStones](shipping/profile_riverstones_teensy_2026-11-18.md) § ● | `fx_shader_draw`<br>worst: RiverStones 2, 58.167 ms / 16 fps | 🟢 58.167 (2) | 70.027 | — | 🟢 0/1667 (0.00%) | — | — | — | COM3 | ship 2026-11-18 23:32 |
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-11-17.md) / [O3](O3/profile_islamicstars_teensy_2026-11-17.md) § ● | `is_timeline_step`<br>worst: dodecahedron_hk35_ambo_hk62_ambo_relax_hk42, 58.059 ms / 16 fps | 🟢 58.059 (16) | 16.537 | 🟢 57.095 (16) | 🟢 0/4895 (0.00%) | 🟢 0/4895 (0.00%) | +28,528 B | +9,456 B | COM3 | ship 2026-11-17 09:09<br>O3 2026-11-17 10:34 |
| [Filigree](shipping/profile_filigree_teensy_2026-11-20.md) § ● | `fx_shader_draw`<br>worst: Jade Weave → Ember Lattice, 57.375 ms / 16 fps | 🟢 57.375 (2) | 76.754 | — | 🟢 0/1667 (0.00%) | — | — | — | COM3 | ship 2026-11-20 08:42 |
| [GridWarp](shipping/profile_gridwarp_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 56.343 | 80.625 | — | 🟢 0/479 (0.00%) | — | — | — | COM3 | ship 2026-11-17 09:41 |
| [AshCloud](shipping/profile_ashcloud_teensy_2026-11-17.md) / [O3](O3/profile_ashcloud_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 56.216 | 85.379 | 🔴 79.81 | 🟢 0/959 (0.00%) | 🔴 544/544 (100%) | +16,608 B | +12,112 B | COM4 | ship 2026-11-17 09:39<br>O3 2026-08-26 02:45 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-11-17.md) / [O3](O3/profile_shapeshifter_teensy_2026-09-25.md) § ● | `ss_draw_all`<br>worst: 8, 56.021 ms / 16 fps | 🟢 56.021 (9) | 76.786 | 🟢 60.05 (9) | 🟢 0/4303 (0.00%) | 🟢 0/2457 (0.00%) | +29,288 B | +24,496 B | COM3 | ship 2026-11-17 09:00<br>O3 2026-09-25 07:29 |
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-11-17.md) / [O3](O3/profile_meshfeedback_teensy_2026-08-26.md) § ● | `mf_feedback_flush`<br>worst: WigglingWormhole, 55.672 ms / 16 fps | 🟢 55.672 (12) | 6.470 | 🟢 58.32 (13) | 🟢 0/2895 (0.00%) | 🟢 0/6688 (0%) | +34,144 B | +21,776 B | COM3 | ship 2026-11-17 09:23<br>O3 2026-08-26 02:10 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-11-17.md) / [O3](O3/profile_hopffibration_teensy_2026-08-26.md) ● | `hf_render_trails` | 🟢 55.635 | 0.672 | 🟢 46.52 | 🟢 0/959 (0.00%) | 🟢 0/1088 (0%) | +19,616 B | +18,272 B | COM3 | ship 2026-11-17 08:45<br>O3 2026-08-26 01:27 |
| [HyperLattice](shipping/profile_hyperlattice_teensy_2026-11-17.md) / [O3](O3/profile_hyperlattice_architecture_teensy_2026-10-01.md) § ● | `hl_shader_draw`<br>worst: octet-4d-flight, 53.311 ms / 16 fps | 🟢 53.311 (9) | 45.039 | 🟢 36.975 (3) | 🟢 0/4799 (0.00%) | 🟢 0/2696 (0.00%) | +9,616 B | +4,912 B | COM4 | ship 2026-11-17 09:08<br>O3 2026-10-01 09:26 |
| [Millefiori](shipping/profile_millefiori_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 52.832 | 79.856 | — | 🟢 0/959 (0.00%) | — | — | — | COM3 | ship 2026-11-17 09:35 |
| [Raymarch](shipping/profile_raymarch_teensy_2026-11-17.md) / [O3](O3/profile_raymarch_teensy_2026-09-20.md) ● | `rm_shader_draw` | 🟢 52.759 | 94.040 | 🟢 56.04 | 🟢 0/479 (0.00%) | 🟢 0/1736 (0.00%) | +12,568 B | +7,904 B | COM3 | ship 2026-11-17 08:52<br>O3 2026-09-20 23:07 |
| [JewelMelt](shipping/profile_jewelmelt_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 52.547 | 79.906 | — | 🟢 0/959 (0.00%) | — | — | — | COM4 | ship 2026-11-17 09:45 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-11-17.md) / [O3](O3/profile_dreamballs_teensy_2026-08-26.md) § ● | `db_timeline_step`<br>worst: 9, 51.363 ms / 16 fps | 🟢 51.363 (10) | 49.416 | 🟢 34.32 (11) | 🟢 0/3199 (0.00%) | 🟢 0/3648 (0%) | +26,896 B | +12,352 B | COM3 | ship 2026-11-17 09:19<br>O3 2026-08-26 02:19 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-11-17.md) / [O3](O3/profile_ringspin_teensy_2026-09-24.md) ● | `rs_draw_rings` | 🟢 51.117 | 3.279 | 🟢 50.750 | 🟢 0/959 (0.00%) | 🟢 0/1087 (0.0%) | +15,400 B | +12,880 B | COM3 | ship 2026-11-17 08:54<br>O3 2026-09-24 20:18 |
| [MindSplatter](shipping/profile_mindsplatter_teensy_2026-11-17.md) / [O3](O3/profile_mindsplatter_architecture_teensy_2026-10-01.md) § ● | `msp_draw_particles`<br>worst: 7, 50.810 ms / 16 fps | 🟢 50.810 (8) | 0.372 | 🟢 54.539 (8) | 🟢 0/1279 (0.00%) | 🟢 0/1736 (0.00%) | +21,976 B | +17,200 B | COM3 | ship 2026-11-17 09:03<br>O3 2026-10-01 09:22 |
| [KaleidoscopeMelt](shipping/profile_kaleidoscopemelt_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 50.074 | 78.699 | — | 🟢 0/959 (0.00%) | — | — | — | COM3 | ship 2026-11-17 09:45 |
| [KaleidoscopeHexOil](shipping/profile_kaleidoscopehexoil_teensy_2026-11-17.md) / [O3](O3/profile_kaleidoscopehexoil_teensy_2026-08-26.md) § ● | `fx_shader_draw`<br>worst: Kaleidoscope Hex Oil 2, 48.354 ms / 16 fps | 🟢 48.354 (2) | 68.013 | 🟢 38.94 (3) | 🟢 0/1679 (0.00%) | 🟢 0/2208 (0%) | +13,456 B | +10,608 B | COM3 | ship 2026-11-17 09:26<br>O3 2026-08-26 02:49 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-11-17.md) / [O3](O3/profile_gsreactiondiffusion_teensy_2026-09-30.md) ● | `grd_render` | 🟢 47.331 | 24.651 | 🔴 279.686 | 🟢 0/479 (0.00%) | 🔴 538/538 (100.00%) | +16,024 B | +10,832 B | COM3 | ship 2026-11-17 08:43<br>O3 2026-09-30 19:15 |
| [MermaidSkin](shipping/profile_mermaidskin_teensy_2026-11-17.md) / [O3](O3/profile_mermaidskin_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 45.898 | 68.306 | 🟢 54.55 | 🟢 0/479 (0.00%) | 🟢 0/1088 (0%) | +16,544 B | +11,952 B | COM3 | ship 2026-11-17 09:37<br>O3 2026-08-26 02:43 |
| [KaleidoscopeMandala](shipping/profile_kaleidoscopemandala_teensy_2026-11-17.md) / [O3](O3/profile_kaleidoscopemandala_teensy_2026-08-26.md) § ● | `fx_shader_draw`<br>worst: Cup Hue, 45.881 ms / 16 fps | 🟢 45.881 (2) | 65.217 | 🟢 36.86 (3) | 🟢 0/1679 (0.00%) | 🟢 0/2368 (0%) | +13,656 B | +11,440 B | COM4 | ship 2026-11-17 09:25<br>O3 2026-08-26 03:17 |
| [ChromaticLichen](shipping/profile_chromaticlichen_teensy_2026-11-17.md) / [O3](O3/profile_chromaticlichen_architecture_teensy_2026-10-01.md) ● | `fx_shader_draw` | 🟢 45.730 | 67.214 | 🟢 35.218 | 🟢 0/479 (0.00%) | 🟢 0/1096 (0.00%) | +15,264 B | +8,672 B | COM4 | ship 2026-11-17 09:37<br>O3 2026-10-01 09:51 |
| [Iridescence](shipping/profile_iridescence_teensy_2026-11-18.md) ● | `fx_shader_draw` | 🟢 44.915 | 83.018 | — | 🟢 0/479 (0.00%) | — | — | — | COM3 | ship 2026-11-18 18:14 |
| [Involution](shipping/profile_involution_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 44.738 | 82.497 | — | 🟢 0/479 (0.00%) | — | — | — | COM3 | ship 2026-11-17 09:39 |
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-11-17.md) / [O3](O3/profile_displacementfield_teensy_2026-09-19.md) § ● | `df_timeline_step`<br>worst: noise, 44.616 ms / 16 fps | 🟢 44.616 (2) | 18.450 | 🟢 57.45 | 🟢 0/1935 (0.00%) | 🟢 0/1088 (0%) | +25,856 B | +22,144 B | COM3 | ship 2026-11-17 08:38<br>O3 2026-09-19 22:19 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-11-17.md) / [O3](O3/profile_bzreactiondiffusion_teensy_2026-08-26.md) ● | `bz_render` | 🟢 44.329 | 76.068 | 🟢 48.65 | 🟢 0/959 (0.00%) | 🟢 0/2048 (0%) | +12,760 B | +10,224 B | COM3 | ship 2026-11-17 08:32<br>O3 2026-08-26 01:16 |
| [Galaxies](shipping/profile_galaxies_teensy_2026-11-17.md) ● | `gx_particle_step` | 🟢 43.803 | 1.477 | — | 🟢 0/959 (0.00%) | — | — | — | COM3 | ship 2026-11-17 08:28 |
| [LiquidWindows](shipping/profile_liquidwindows_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 42.732 | 63.633 | — | 🟢 0/479 (0.00%) | — | — | — | COM4 | ship 2026-11-17 09:41 |
| [KaleidoscopeFlowers](shipping/profile_kaleidoscopeflowers_teensy_2026-11-17.md) / [O3](O3/profile_kaleidoscopeflowers_teensy_2026-08-26.md) § ● | `fx_shader_draw`<br>worst: Fine Grid, 41.644 ms / 16 fps | 🟢 41.644 (3) | 61.934 | 🟢 33.38 (4) | 🟢 0/2767 (0.00%) | 🟢 0/4128 (0%) | +14,576 B | +11,712 B | COM4 | ship 2026-11-17 09:18<br>O3 2026-08-26 03:03 |
| [KaleidoscopeHexBright](shipping/profile_kaleidoscopehexbright_teensy_2026-11-17.md) / [O3](O3/profile_kaleidoscopehexbright_teensy_2026-08-26.md) § ● | `fx_shader_draw`<br>worst: Hex Twin Wave, 41.633 ms / 16 fps | 🟢 41.633 (2) | 57.993 | 🟢 32.41 (3) | 🟢 0/1679 (0.00%) | 🟢 0/2368 (0%) | +14,520 B | +11,712 B | COM4 | ship 2026-11-17 09:29<br>O3 2026-08-26 03:14 |
| [KaleidoscopeSmooth](shipping/profile_kaleidoscopesmooth_teensy_2026-11-17.md) / [O3](O3/profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md) § ● | `fx_shader_draw`<br>worst: Direct Grid, 40.288 ms / 16 fps | 🟢 40.288 (4) | 63.892 | 🟢 30.314 (4) | 🟢 0/3839 (0.00%) | 🟢 0/4137 (0.00%) | +15,048 B | +8,624 B | COM3 | ship 2026-11-17 09:14<br>O3 2026-10-01 11:19 |
| [GridWorms](shipping/profile_gridworms_teensy_2026-11-17.md) ● | `fx_shader_draw` | 🟢 39.075 | 54.282 | — | 🟢 0/479 (0.00%) | — | — | — | COM4 | ship 2026-11-17 09:35 |
| [AlienBrain](shipping/profile_alienbrain_teensy_2026-11-17.md) / [O3](O3/profile_alienbrain_teensy_2026-08-26.md) § ● | `fx_shader_draw`<br>worst: Alien Brain 3, 37.447 ms / 16 fps | 🟢 37.447 (4) | 57.420 | 🟢 27.90 (5) | 🟢 0/3839 (0.00%) | 🟢 0/4768 (0%) | +14,560 B | +11,712 B | COM4 | ship 2026-11-17 09:13<br>O3 2026-08-26 02:24 |
| [MobiusGrid](shipping/profile_mobiusgrid_teensy_2026-11-17.md) / [O3](O3/profile_mobiusgrid_teensy_2026-08-26.md) § ● | `fx_shader_draw`<br>worst: Mobius Grid 2, 33.522 ms / 16 fps | 🟢 33.522 (2) | 55.942 | 🟢 21.44 (3) | 🟢 0/1679 (0.00%) | 🟢 0/2688 (0%) | +13,448 B | +10,560 B | COM3 | ship 2026-11-17 08:48<br>O3 2026-08-26 01:30 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-11-17.md) / [O3](O3/profile_hankinsolids_teensy_2026-09-28.md) § ● | `hk_timeline_step`<br>worst: truncatedCube, 32.980 ms / 16 fps | 🟢 32.980 (18) | 18.945 | 🟢 32.865 (18) | 🟢 0/2719 (0.00%) | 🟢 0/4137 (0.00%) | +17,088 B | +2,144 B | COM4 | ship 2026-11-17 09:22<br>O3 2026-09-28 16:45 |
| [GridSpace](shipping/profile_gridspace_teensy_2026-11-17.md) / [O3](O3/profile_gridspace_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 32.794 | 46.539 | 🟢 25.39 | 🟢 0/479 (0.00%) | 🟢 0/1088 (0%) | +14,424 B | +11,712 B | COM3 | ship 2026-11-17 09:33<br>O3 2026-08-26 02:33 |
| [CosmicEyeball](shipping/profile_cosmiceyeball_teensy_2026-11-17.md) / [O3](O3/profile_cosmiceyeball_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 30.734 | 48.727 | 🟢 22.97 | 🟢 0/959 (0.00%) | 🟢 0/1088 (0%) | +12,552 B | +10,320 B | COM3 | ship 2026-11-17 09:43<br>O3 2026-08-26 03:05 |
| [Comets](shipping/profile_comets_teensy_2026-11-17.md) / [O3](O3/profile_comets_teensy_2026-08-26.md) § ● | `cm_draw_trail`<br>worst: 5, 29.869 ms / 16 fps | 🟢 29.869 (12) | 0.869 | 🟢 28.71 (13) | 🟢 0/1919 (0.00%) | 🟢 0/4128 (0%) | +15,432 B | +12,944 B | COM3 | ship 2026-11-17 09:29<br>O3 2026-08-26 02:02 |
| [Fishbowl](shipping/profile_fishbowl_teensy_2026-11-17.md) / [O3](O3/profile_fishbowl_teensy_2026-08-26.md) ● | `fish_build_vertices` | 🟢 28.256 | 1.681 | 🟢 21.29 | 🟢 0/959 (0.00%) | 🟢 0/1088 (0%) | +27,856 B | +26,400 B | COM3 | ship 2026-11-17 08:34<br>O3 2026-08-26 01:18 |
| [AlienCore](shipping/profile_aliencore_teensy_2026-11-17.md) / [O3](O3/profile_aliencore_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 27.530 | 42.245 | 🟢 20.10 | 🟢 0/479 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | COM4 | ship 2026-11-17 09:33<br>O3 2026-08-26 02:30 |
| [AlienOcean](shipping/profile_alienocean_teensy_2026-11-17.md) / [O3](O3/profile_alienocean_teensy_2026-08-26.md) ● | `fx_shader_draw` | 🟢 27.306 | 49.224 | 🟢 24.98 | 🟢 0/479 (0.00%) | 🟢 0/1088 (0%) | +13,648 B | +11,440 B | COM3 | ship 2026-11-17 09:31<br>O3 2026-08-26 02:28 |
| [GnomonicStars](shipping/profile_gnomonicstars_teensy_2026-11-17.md) / [O3](O3/profile_gnomonicstars_teensy_2026-08-26.md) ● | `gn_draw_stars` | 🟢 26.078 | 15.179 | 🟢 26.29 | 🟢 0/479 (0.00%) | 🟢 0/1088 (0%) | +12,184 B | +11,104 B | COM3 | ship 2026-11-17 08:40<br>O3 2026-08-26 01:22 |
| [PetalFlow](shipping/profile_petalflow_teensy_2026-11-17.md) / [O3](O3/profile_petalflow_teensy_2026-08-26.md) ● | `pf_draw_rings` | 🟢 25.420 | 19.712 | 🟢 10.77 | 🟢 0/479 (0.00%) | 🟢 0/1088 (0%) | +22,272 B | +21,008 B | COM3 | ship 2026-11-17 08:50<br>O3 2026-08-26 01:32 |
| [SphericalHarmonics](shipping/profile_sphericalharmonics_teensy_2026-11-17.md) / [O3](O3/profile_sphericalharmonics_teensy_2026-08-26.md) § ● | `sh_rasterize`<br>worst: l=4, m=4, 12.841 ms / 16 fps | 🟢 12.841 (24) | 21.305 | 🟢 11.70 (24) | 🟢 0/1535 (0.00%) | 🟢 0/3488 (0%) | +8,016 B | +6,400 B | COM4 | ship 2026-11-17 09:31<br>O3 2026-08-26 01:57 |

● captures updated 2026-11-17, 2026-11-18, 2026-11-19.

Shipping results use only the authored show budget, including in-show transitions.
Setup frame 1 remains separate. Post-show capture padding is excluded from every
peak, spill denominator, phase readout and per-preset table.
A full cycle means the scheduled show slot; an additional wrap is not required.

§ Shipping cyclers carry one aligned line per preset colour bucket, red first;
(N) counts presets. Buckets include transitions and exclude setup frame 1.
Per-preset cost spans and hold-window selection are in the linked reports.

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
[report](shipping/profile_hyperlattice_teensy_2026-11-17.md#per-preset-table) gives the observed entry buckets
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

[Historical IslamicStars TS4 comparison](O3/profile_islamicstars_historical_ts4_teensy_2026-09-28.md).

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

## Retired effect profiles

- [Jade Filigree](retired/profile_jadefiligree_teensy_2026-11-19.md)
- [Ember Tracery](retired/profile_embertracery_teensy_2026-11-19.md)

## Effects excluded from the Phantasm image

KaleidoscopeStainedGlass and RingShower remain in the full roster; their rows retain their archived measurements.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [KaleidoscopeStainedGlass selective o3](excluded/profile_kaleidoscopestainedglass_selective_o3_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 43.00 | 🟢 0/1087 (0.00%) | 2026-10-06 18:50 |
| [KaleidoscopeStainedGlass global o3](excluded/profile_kaleidoscopestainedglass_global_o3_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 46.99 | 🟢 0/1088 (0%) | 2026-08-26 02:51 |
| [RingShower selective o3](excluded/profile_ringshower_selective_o3_teensy_2026-10-06.md) | `rsh_draw_rings` | 🟢 4.32 | 🟢 0/1087 (0.00%) | 2026-10-06 18:20 |
| [RingShower global o3](excluded/profile_ringshower_global_o3_teensy_2026-08-26.md) | `rsh_draw_rings` | 🟢 3.86 | 🟢 0/1088 (0%) | 2026-08-26 01:36 |

Low-resolution-only (96×20) effects:

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [MobiusRings Hole A/B](excluded/profile_mobiusrings_teensy_2026-10-13.md) | `mr_draw_grid` | 🟢 12.30 | 🟢 0/1088 (0.00%) | 2026-10-13 |

Dynamo and Thrusters are also low-resolution only and have no report here.

## Temporal sampling comparison

[GridSpace matched A/B](shipping/profile_gridspace_teensy_2026-10-17.md) preserves the historical same-source temporal-on/off measurements. Current shipping timing is in the canonical report above.

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
