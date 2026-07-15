# On-device effect profiles — Teensy 4.0, segmented mode (2026-07-14)

On-device timing profiles for **all 26 roster effects**, captured on a
bench-attached Teensy 4.0 running the shipping Phantasm configuration
(`POVSegmented<288, 4, 480>`, board = segment 0 master, `-Os`, newlib-nano, DMA
LEDs, flywheel + DMA ISRs live), via the `HS_PROFILE` cycle-counter harness
(`targets/Profile/Profile.ino`, `just profile <Effect>`).

**Point-in-time snapshots** (regenerate with `just profile <Effect>`). These
reports and the `HS_PROFILE`
instrumentation in the effect headers are working-tree only. Each effect
renders one **quadrant** = this segment's 72-row band × the display window's
144-col half ≈ **10,368 px**. A display window (half-revolution at 480 RPM) is
**62.5 ms**; `draw_frame` wall time quantizes up to whole windows, so cadence is
16 fps (1 window), 8 fps (2), 5.3 fps (3), etc. Every scope absorbs live ISR
time (that is the point of segmented-mode profiling).

The five originally-requested effects (DisplacementField, DreamBalls, Flyby,
IslamicStars, RingSpin) used two-pass captures (128- + 32-frame windows); the
other 21 used a single ~75 s pass at 64-frame windows.

## Ranked by active render cost (heaviest first)

"Render/frame" is the effect's own work (frame minus the `*_buffer_wait`
display-sync idle). Cadence is the observed steady/representative fps.

**●** = re-profiled 2026-07-14 after the plot column-cull batch (peak/steady
render, whichever bounds cadence). Ranked by `-Os` render/frame (ms), heaviest
first.

| Effect | Render/frame | Cadence | Dominant cost (scope) | Bound by |
|---|--:|---|---|---|
| ● [Dynamo](profile_dynamo_teensy_2026-07-14.md) | 166 ms (peak) | 4–8 fps | Trails flush (`dy_filter_flush` 74%) | feedback composite |
| ● [MeshFeedback](profile_meshfeedback_teensy_2026-07-14.md) | 124 ms (peak) | 6–8 fps | feedback flush (`mf_feedback_flush` 75%) | feedback composite |
| ● [FlowField](profile_flowfield_teensy_2026-07-14.md) | 121 ms (peak) | 8 fps (pulses) | particle trail raster (`ff_particle_draw` 94%) | particle rasterizer |
| [Raymarch](profile_raymarch_teensy_2026-07-14.md) | 121 ms | 8 fps | volume ray-march (`rm_shader_draw` 86%) | per-pixel shader (11.7 µs/px) |
| ● [DreamBalls](profile_dreamballs_teensy_2026-07-14.md) | 109 ms (peak) | 5–8 fps | wireframe raster (`db_mesh_plot` 77%) | mesh wireframe raster |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-07-14.md) | 108 ms | 8 fps | SSAA raster 57% + simulate 28% | RD simulate + SSAA shader |
| [DisplacementField](profile_displacementfield_teensy_2026-07-14.md) | 105 ms | 8 fps | SDF ring raster (`df_ring_scan` 69%) | SDF ring rasterizer |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-07-14.md) | 104 ms | 8 fps | SSAA raster (`bz_raster` 75%) | SSAA shader (light simulate) |
| ● [MindSplatter](profile_mindsplatter_teensy_2026-07-14.md) | 100 ms | 8 fps | particle raster (`msp_particle_scan` 77%) | particle rasterizer |
| [IslamicStars](profile_islamicstars_teensy_2026-07-14.md) | 43–133 ms (per shape) | 5–16 fps | per-face SDF (`scan_mesh_raster`) | mesh SDF rasterizer |
| [Voronoi](profile_voronoi_teensy_2026-07-14.md) | 77 ms | 8 fps | KD-nearest shade (`vo_shade` 61%) | per-pixel query (7.4 µs/px) |
| [RingSpin](profile_ringspin_teensy_2026-07-14.md) | 77 ms | 8 fps | SDF ring raster (`rs_ring_scan` 61%) | SDF ring raster (4.5× overdraw) |
| [Flyby](profile_flyby_teensy_2026-07-14.md) | 50–77 ms | 8–16 fps | stereographic shader (`fly_shader_draw` 62%) | per-pixel shader |
| ● [HopfFibration](profile_hopffibration_teensy_2026-07-14.md) | 70 ms (peak) | 8–16 fps | trail raster (`hf_trail_raster` 62%) | trail rasterizer |
| ● [ShapeShifter](profile_shapeshifter_teensy_2026-07-14.md) | 68 ms (peak) | 8–16 fps | SDF scan (`ss_scan_dispatch` 86%) + plot 12% | SDF ring raster |
| [HankinSolids](profile_hankinsolids_teensy_2026-07-14.md) | 40 ms | 8–16 fps | per-face SDF (`scan_mesh_raster` 54%) | mesh SDF rasterizer |
| [Comets](profile_comets_teensy_2026-07-14.md) | 32 ms | 16 fps | point raster (`cm_point_scan` 52%) | point rasterizer |
| [Liquid2D](profile_liquid2d_teensy_2026-07-14.md) | 30 ms | 16 fps | stereographic shader (`lq_shader_draw` 48%) | per-pixel shader |
| ● [ChaoticStrings](profile_chaoticstrings_teensy_2026-07-14.md) | 20 ms | 16 fps | multiline raster (`cs_multiline_draw` 22%) | multiline rasterizer (idle-bound) |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-07-14.md) | 20 ms | 16 fps | full-quadrant field raster (`sh_rasterize` 32%) | field rasterizer |
| ● [MobiusGrid](profile_mobiusgrid_teensy_2026-07-14.md) | 17.7 ms (peak) | 16 fps | curve raster (`mg_lines_draw`/`mg_rings_draw` 27%) | curve rasterizer (idle-bound) |
| ● [PetalFlow](profile_petalflow_teensy_2026-07-14.md) | 16.7 ms | 16 fps | ring raster (`pf_ring_scan` 25%) | ring rasterizer (idle-bound) |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-07-14.md) | 16 ms | 16 fps | star raster (`gn_star_scan` 26%) | star rasterizer |
| [DistortedRing](profile_distortedring_teensy_2026-07-14.md) | 4–8 ms | 16 fps | SDF ring raster (`dr_ring_scan` 6–13%) | single SDF ring (idle-bound) |
| ● [Thrusters](profile_thrusters_teensy_2026-07-14.md) | 3.35 ms | 16 fps | ring raster (`th_ring_draw` + `th_thruster_draw` 5%) | idle-bound |
| ● [RingShower](profile_ringshower_teensy_2026-07-14.md) | 3.12 ms | 16 fps | ring raster (`rsh_ring_plot` 4%) | idle-bound (lightest) |

The 11 ● effects were re-captured 2026-07-14 after the plot column-cull batch
(9ac8cebd, 62450701, 708d4b9b). Only DreamBalls and MindSplatter moved
materially (wireframe/particle edges span columns); the idle-bound ring effects
and the flush-bound Dynamo/MeshFeedback are essentially unchanged by the batch.
See [`../README.md`](../README.md) for the batch's measured effect. The other 15
rows carry their original 2026-07-14 captures.

## By category

- **Per-pixel shader / field** (cost ∝ pixels × per-pixel math): Raymarch
  (heaviest, 11.7 µs/px), Voronoi (KD-nearest, 7.4 µs/px), Flyby, Liquid2D,
  SphericalHarmonics.
- **Reaction-diffusion** (simulate + 4× SSAA shader): GSReactionDiffusion
  (simulate-heavy, ~28% of render), BZReactionDiffusion (simulate ~8%). Both
  SSAA-raster-dominated.
- **Feedback composite** (the `filters.flush` replay dominates, NOT the
  rasterizer): Dynamo (Trails), MeshFeedback. The two outliers of the roster.
- **Particle rasterizer** (ramps as the pool fills): FlowField, MindSplatter.
- **Mesh SDF rasterizer** (per-face SDF via `Scan::Mesh`, exposes render-side
  `scan_mesh_raster`/`scan_face_setup`): IslamicStars, HankinSolids.
- **Wireframe / trail / ring / point rasterizers** (`Plot::*` / `Scan::*`):
  DreamBalls, DisplacementField, RingSpin, ShapeShifter, DistortedRing,
  HopfFibration, Comets, ChaoticStrings, PetalFlow, GnomonicStars, MobiusGrid,
  Thrusters, RingShower.

## Cross-cutting observations

- **The rasterizer/shader is the whole cost almost everywhere.** Geometry,
  transform, LUT-bake, and simulate stages are <5% of the frame for every
  effect except the two reaction-diffusion sims (simulate 8–28%) and the two
  feedback effects (flush 63–78%). Optimization leverage is per-pixel /
  per-primitive raster cost and overdraw.
- **Feedback replay is the only non-raster bottleneck.** Dynamo and MeshFeedback
  are dominated by `filters.flush` (Trails / feedback-buffer composite), not by
  drawing — a distinct optimization target from every other effect.
- **Half the roster is idle-bound at 16 fps** with large headroom (buffer_wait
  60–99%): RingShower and Thrusters use ~5% of the frame; DistortedRing,
  MobiusGrid, GnomonicStars, PetalFlow, SphericalHarmonics, ChaoticStrings,
  Comets, Liquid2D all clear one window comfortably.
- **ISR steal correlates with memory pressure.** The light/idle effects sit at
  ~2.7–3.3% chip ISR steal; the heavy memory-bound per-pixel effects (RingSpin,
  Voronoi) push `isr_pack` averages up to ~20–27 µs and ~6–8% steal as the
  OCRAM display-buffer reads contend — the `packPixel` `min` floor (~4.25 µs)
  is identical everywhere, so the inflation is contention + serial-dump
  preemption, not extra ISR work.
- **Image sizes** span FLASH code 33.5 KB (SphericalHarmonics/Voronoi) to
  69 KB (IslamicStars); ITCM 21 KB to 44.6 KB (ShapeShifter). RAM2 free is a
  constant 4,736 B (the shared arena partition).

## Reproduce

```
just profile <Effect> [seconds]           # build + flash + capture
# or, with knobs:
$env:PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=<Effect> -D HS_PROFILE_WINDOW=64'
pio run -e profile -t upload
python tools/profile_capture.py --seconds 75 --out build/profile_<effect>.log
```

Raw capture logs are in `build/profile_<effect>.log`.
