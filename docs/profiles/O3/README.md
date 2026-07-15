# On-device effect profiles — **-O3** (2026-07-14)

The **-O3** twins of the shipping `-Os` profiles in [`../Os/`](../Os/README.md).
Built with the `profile_o3` PlatformIO env — identical to the shipping profile
(newlib-nano, DMA LEDs, `phantasm.ld`, N=4, `HS_PROFILE`) **except `-O3`
(`-ffast-math`) replaces `-Os`**, so the only variable vs the `-Os` reports is
the optimization level. Each effect renders one quadrant ≈ 10,368 px; a display
window is 62.5 ms (16 fps = 1 window, 8 fps = 2, …).

**-O3 is not the shipping config.** A single-effect image fits FlexRAM at -O3,
but the full 26-effect roster overflows ITCM (~363 KB / 12 banks) — which is why
the shipped Phantasm image is `-Os`. These reports measure what `-Os` costs in
speed. The head-to-head speedup / cadence / size table is in
[`../README.md`](../README.md).

## Ranked by -O3 render cost (heaviest first)

Render/frame = frame − `*_buffer_wait` idle, at -O3, in **ms**. **●** = re-profiled after a
later landed batch (2026-07-14 plot column-cull; DisplacementField 2026-07-15
fused-stack scan). Ranked heaviest first.

| Effect | Render/frame | Cadence | Dominant cost (scope) | -O3 speedup vs -Os |
|---|--:|---|---|--:|
| ● [Dynamo](profile_dynamo_teensy_2026-07-14.md) | 112.9 ms (peak) | 8 fps | feedback flush (`dy_filter_flush`) | 1.47× (flush) |
| ● [FlowField](profile_flowfield_teensy_2026-07-14.md) | 111.1 ms (peak) | 8 fps (pulses) | particle raster (`ff_particle_draw`) | ~1.3× (pop-jitter) |
| [Raymarch](profile_raymarch_teensy_2026-07-14.md) | 98 ms | 8 fps | volume ray-march (`rm_shader_draw`) | 1.23× |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-07-14.md) | 85 ms | 8 fps | SSAA raster + simulate | 1.29× |
| ● [DreamBalls](profile_dreamballs_teensy_2026-07-14.md) | 72.5 ms (peak) | 8 fps | wireframe raster (`db_mesh_plot`) | 1.50× |
| ● [DisplacementField](profile_displacementfield_teensy_2026-07-15.md) | 59.5 ms | **16 fps** | fused stack raster (`df_fused_scan`) | 1.39× |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-07-14.md) | 71 ms | 8 fps | SSAA raster (`bz_raster`) | 1.45× |
| ● [MeshFeedback](profile_meshfeedback_teensy_2026-07-14.md) | 69.4 ms (peak) | **8→16 fps** | feedback flush (`mf_feedback_flush`) | 1.79× (flush) |
| ● [MindSplatter](profile_mindsplatter_teensy_2026-07-14.md) | 68.0 ms | 8 fps (occ. 16) | particle raster (`msp_particle_scan`) | 1.47× |
| [IslamicStars](profile_islamicstars_teensy_2026-07-14.md) | 33–110 ms (per shape) | 8–16 fps | per-face SDF (`scan_mesh_raster`) | 1.16–1.51× |
| [Voronoi](profile_voronoi_teensy_2026-07-14.md) | 57 ms | **16 fps** | per-pixel KD (`vo_shade`) | 1.35× |
| ● [RingSpin](profile_ringspin_teensy_2026-07-15.md) | 52 ms | **16 fps** | fused ring-group raster (`rs_ring_scan`) | 1.30× |
| ● [ShapeShifter](profile_shapeshifter_teensy_2026-07-14.md) | 50.9 ms (peak) | 8↔16 fps | SDF scan (`ss_scan_dispatch`) + plot | 1.33× (blend 1.7×) |
| ● [HopfFibration](profile_hopffibration_teensy_2026-07-14.md) | 44.1 ms (peak) | **8→16 fps** | trail raster (`hf_trail_raster`) | 1.59× |
| [Flyby](profile_flyby_teensy_2026-07-14.md) | 42 ms | **16 fps** | shader (`fly_shader_draw`) | 1.82× |
| [HankinSolids](profile_hankinsolids_teensy_2026-07-14.md) | 29 ms | 16 fps | per-face SDF (`scan_mesh_raster`) | 1.29–1.38× |
| [Liquid2D](profile_liquid2d_teensy_2026-07-14.md) | 26 ms | 16 fps | shader (`lq_shader_draw`) | 1.14× |
| [Comets](profile_comets_teensy_2026-07-14.md) | 22 ms | 16 fps | point raster (`cm_point_scan`) | 1.55× |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-07-14.md) | 14.5 ms | 16 fps | field raster (`sh_rasterize`) | 1.39× |
| ● [ChaoticStrings](profile_chaoticstrings_teensy_2026-07-14.md) | 14.2 ms | 16 fps | multiline raster (`cs_multiline_draw`) | 1.40× |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-07-14.md) | 13.5 ms | 16 fps | star raster (`gn_star_scan`) | 1.19× (blend 2.4×) |
| ● [MobiusGrid](profile_mobiusgrid_teensy_2026-07-14.md) | 12.4 ms (peak) | 16 fps | curve raster (`mg_lines_draw`) | 1.42× |
| ● [PetalFlow](profile_petalflow_teensy_2026-07-14.md) | 11.5 ms | 16 fps | ring raster (`pf_ring_scan`) | 1.45× |
| [DistortedRing](profile_distortedring_teensy_2026-07-14.md) | 3–6 ms | 16 fps | SDF ring raster (`dr_ring_scan`) | 1.29× |
| ● [Thrusters](profile_thrusters_teensy_2026-07-14.md) | 2.35 ms | 16 fps | ring raster (`th_ring_draw`) | 1.43× |
| ● [RingShower](profile_ringshower_teensy_2026-07-14.md) | 2.20 ms | 16 fps | ring raster (`rsh_ring_plot`) | 1.42× |

**Bold cadence** = a frame-rate tier that -O3 gained vs -Os (see the comparison
table in [`../README.md`](../README.md)).

The 11 ● effects were re-captured 2026-07-14 after the plot column-cull batch;
their -O3 speedup is the compiler, unchanged by the batch. DisplacementField
and RingSpin were re-captured 2026-07-15 at their fused-scan tips (RingSpin's
16 fps lock is the fusion + -O3 together). The other rows carry their
original 2026-07-14 captures. IslamicStars now has a full
[per-shape report](../profile_islamicstars_pershape_teensy_2026-07-14.md).

## Reproduce

```
$env:PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=<Effect> -D HS_PROFILE_WINDOW=64'
pio run -e profile_o3 -t upload
python tools/profile_capture.py --seconds 75 --out build/profile_<effect>_o3.log
```
