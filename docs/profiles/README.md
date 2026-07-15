# On-device effect profiles — Teensy 4.0, segmented mode (2026-07-14)

On-device timing profiles for **all 26 roster effects**, captured on a
bench-attached Teensy 4.0 running the shipping Phantasm configuration
(`POVSegmented<288, 4, 480>`, board = segment 0 master, newlib-nano, DMA LEDs,
flywheel + DMA ISRs live), via the `HS_PROFILE` cycle-counter harness. Each
effect renders one **quadrant** ≈ **10,368 px**; a display window is **62.5 ms**,
so cadence quantizes to 16 fps (1 window), 8 fps (2), 5.3 fps (3), etc.

Captured at **two optimization levels**:

- [**`Os/`**](Os/README.md) — `-Os`, the **shipping** config. The full 26-effect
  Phantasm image only fits FlexRAM at `-Os` (per `platformio.ini`), so this is
  what actually ships.
- [**`O3/`**](O3/README.md) — `-O3` (`-ffast-math`), via the `profile_o3` env
  (identical to the shipping profile except the optimization level). `-O3` is
  **not shippable at roster scale** — a single-effect image fits, but the full
  roster overflows ITCM. These runs isolate the optimization level as the only
  variable to measure what `-Os` costs in speed.
- [**`selectiveO3/`**](selectiveO3/README.md) — the shipping `-Os` image with
  the landed `HS_O3` hot-loop regions (docs/selective_o3_spec.md), 2026-07-15:
  RingSpin and DisplacementField three-way comparisons against the same-tip
  `-Os` and global-`-O3` captures. Since the regions landed, `just profile
  <Effect>` measures this configuration — the `Os/` reports below predate it.

**Point-in-time snapshots** (regenerate with `just profile <Effect>`; numbers
age as the render code moves). Reports + the `HS_PROFILE` instrumentation +
the `profile_o3` env are working-tree only.

## -O3 vs -Os — render time (ms) and speedup

**Render** = the effect's own per-frame work (frame minus the `*_buffer_wait`
display-sync idle), at the steady/representative window, in **ms**. `O3×` = -Os
render ÷ -O3 render. Size Δ is the single-effect image growth. Ordered heaviest
`-Os` render first. **●** marks the 11 effects re-profiled 2026-07-14 after the
plot column-cull batch (9ac8cebd/62450701/708d4b9b); the other 15 carry their
original 2026-07-14 captures (the batch does not touch their `Scan::` path).
DisplacementField was re-profiled 2026-07-15 after the fused-stack batch
(ec61faa7..7d50b672), which crosses its whole cycle to 16 fps at -O3.
RingSpin was re-profiled 2026-07-15 after the fused trail-scan landing plus
the slim RingGroup rework (730b6d2f): -Os per-blend cost dropped ~29%
(~14 fps effective, up from ~9.4), and -O3 hard-locks 16 fps.

| Effect | Dominant scope | Render Os (ms) | Render O3 (ms) | O3× | FLASH Δ | ITCM Δ |
|---|---|--:|--:|--:|--:|--:|
| ● [Dynamo](O3/profile_dynamo_teensy_2026-07-14.md) | feedback flush | 165.8 | 112.9 | 1.47× | +53% | +75% |
| ● [MeshFeedback](O3/profile_meshfeedback_teensy_2026-07-14.md) | feedback flush | 124.4 | 69.4 | 1.79× | +49% | +53% |
| [Raymarch](O3/profile_raymarch_teensy_2026-07-14.md) | volume ray-march | 121 | 98 | 1.23× | +42% | +33% |
| ● [FlowField](O3/profile_flowfield_teensy_2026-07-14.md) | particle raster | 120.7† | 111.1† | ~1.3‡ | +48% | +68% |
| ● [DreamBalls](O3/profile_dreamballs_teensy_2026-07-14.md) | wireframe raster | 108.8 | 72.5 | 1.50× | +43% | +37% |
| [GSReactionDiffusion](O3/profile_gsreactiondiffusion_teensy_2026-07-14.md) | SSAA raster + sim | 108 | 85 | 1.29× | +49% | +78% |
| [BZReactionDiffusion](O3/profile_bzreactiondiffusion_teensy_2026-07-14.md) | SSAA raster | 104 | 71 | 1.45× | +44% | +69% |
| ● [MindSplatter](O3/profile_mindsplatter_teensy_2026-07-14.md) | particle raster | 100.0 | 68.0 | 1.47× | +48% | +67% |
| [IslamicStars](O3/profile_islamicstars_teensy_2026-07-14.md) | per-face SDF | 43–133§ | 33–110§ | 1.16–1.51§ | +47% | +41% |
| [Voronoi](O3/profile_voronoi_teensy_2026-07-14.md) | per-pixel KD | 77 | 57 | 1.35× | +62% | +84% |
| ● [DisplacementField](O3/profile_displacementfield_teensy_2026-07-15.md) | fused ring-stack raster | 82.6 | 59.5 | 1.39× | +56% | +77% |
| ● [RingSpin](O3/profile_ringspin_teensy_2026-07-15.md) | slim ring-group raster | 49.6 | 39.1 | 1.27× | +89% | +124% |
| [Flyby](O3/profile_flyby_teensy_2026-07-14.md) | stereographic shader | 77 | 42 | 1.82× | +64% | +94% |
| ● [HopfFibration](O3/profile_hopffibration_teensy_2026-07-14.md) | trail raster | 70.2 | 44.1 | 1.59× | +70% | +98% |
| ● [ShapeShifter](O3/profile_shapeshifter_teensy_2026-07-14.md) | SDF scan (+plot) | 67.8 | 50.9 | 1.33× | +50% | +64% |
| [HankinSolids](O3/profile_hankinsolids_teensy_2026-07-14.md) | per-face SDF | 40 | 29 | 1.29–1.38× | +44% | +43% |
| [Comets](O3/profile_comets_teensy_2026-07-14.md) | point raster | 32 | 22 | 1.55× | +57% | +77% |
| [Liquid2D](O3/profile_liquid2d_teensy_2026-07-14.md) | stereographic shader | 30 | 26 | 1.14× | +75% | +108% |
| ● [ChaoticStrings](O3/profile_chaoticstrings_teensy_2026-07-14.md) | multiline raster | 20.0 | 14.2 | 1.40× | +65% | +89% |
| [SphericalHarmonics](O3/profile_sphericalharmonics_teensy_2026-07-14.md) | field raster | 20 | 14.5 | 1.39× | +42% | +61% |
| ● [MobiusGrid](O3/profile_mobiusgrid_teensy_2026-07-14.md) | curve raster | 17.7 | 12.4 | 1.42× | +46% | +65% |
| ● [PetalFlow](O3/profile_petalflow_teensy_2026-07-14.md) | ring raster | 16.7 | 11.5 | 1.45× | +68% | +99% |
| [GnomonicStars](O3/profile_gnomonicstars_teensy_2026-07-14.md) | star raster | 16 | 13.5 | 1.19× | +43% | +63% |
| [DistortedRing](O3/profile_distortedring_teensy_2026-07-14.md) | SDF ring raster | 8 | 6 | 1.29× | +61% | +88% |
| ● [Thrusters](O3/profile_thrusters_teensy_2026-07-14.md) | ring raster | 3.35 | 2.35 | 1.43× | +53% | +73% |
| ● [RingShower](O3/profile_ringshower_teensy_2026-07-14.md) | ring raster | 3.12 | 2.20 | 1.42× | +50% | +71% |

† FlowField pulses (particle bursts, no fixed steady state); the peak-burst
render is reported and the burst population is **not pinned across flashes**, so
its ratio is soft. ‡ per-fragment O3 is ~1.3–1.45× but the O3 run happened to
catch a heavier burst — trust the per-fragment estimate over the peak ratio.
§ IslamicStars sweeps 24 shapes with a 3× render spread; the range is the
per-shape span — see the [per-shape report](profile_islamicstars_pershape_teensy_2026-07-14.md).
The Dynamo/MeshFeedback speedups are their **feedback-flush** (`Pixel::Feedback::flush`),
which the plot-cull batch does not touch — see below.

## Plot column-cull batch — measured effect (2026-07-14)

The batch (9ac8cebd column-arc cull, 62450701 planar-edge cull, 708d4b9b
deferred trail shading) adds a per-edge **column** cull to `Plot::rasterize`
under the segment clip, on top of the existing row cull. Re-profiling the 11
`Plot::`-rasterized effects (● above) at both levels shows the win is **highly
uneven** — it only pays where a large share of the frame is plot-edge raster
*and* the edges span columns:

- **Real wins**: **DreamBalls** (wireframe — long geodesic/planar edges cross the
  quadrant; −Os render ~145 → 109 ms) and **MindSplatter** (particle trails; −Os
  scan 107.75 → 96.48 ms, −10.5 %). These are the effects the batch was built
  for.
- **Little to none**: **HopfFibration** — trails are short and near-meridian, so
  they hit `edge_col_span`'s near-meridian no-cull fallback; its median render is
  unchanged from the pre-batch capture. Idle-bound ring effects (**PetalFlow,
  MobiusGrid, ChaoticStrings, Thrusters, RingShower**) already render well under
  one 62.5 ms window, so a smaller edge cull only widens `*_buffer_wait` idle —
  no cadence change.
- **Not beneficiaries despite using `Plot::`**: **MeshFeedback** (plot draw is
  5–7 % of the frame) and **Dynamo** (13 %) are bottlenecked on
  `Pixel::Feedback::flush`, which the batch does not touch; **ShapeShifter** is
  `Scan::`-dominated (plot is 12–18 %). Their headline speedups are the compiler
  (-O3), not the cull.

The **-O3-vs-Os** per-effect speedup (~1.4–1.6×) is orthogonal to the batch — it
is the compiler, and is unchanged by the cull.

## What -O3 buys (and why -Os still ships)

- **Render is 1.14×–1.82× faster at -O3**, clustering around **~1.4×**. The
  biggest wins are the memory/loop-bound compositors and rasterizers (Flyby
  shader 1.82×, MeshFeedback feedback flush 1.79×, HopfFibration trail 1.59×,
  DreamBalls/Comets ~1.55×); the smallest are the branch-heavy or already-lean
  shaders (Liquid2D 1.14×, GnomonicStars star-transform 1.19×, Raymarch's
  divergent ray-march 1.23×).
- **Per-pixel inner loops benefit most.** Wherever a `filter_blend` leaf is
  visible, its per-blend cost roughly **halves** at -O3 (~1.7×–2.4×: GnomonicStars
  231→96 cyc, RingSpin/DistortedRing ~2×, SphericalHarmonics 1.95×) — -O3's
  vectorization/unrolling pays off most in the tight SDF distance + blend loops.
- **7 effects cross a cadence tier**: MeshFeedback, Flyby, Voronoi (all 8→16 fps),
  DreamBalls & FlowField (5→8 fps), HopfFibration (jitter→steady 16), IslamicStars
  (mid-shapes 8→16); HankinSolids & RingSpin cross partially. The rest keep their
  cadence — the saved time falls into `*_buffer_wait` display-sync idle because
  the render either was already well under a window (the light 16 fps effects) or
  still overruns one after the speedup (the heaviest raster/RD/shader effects).
- **-O3 costs ~+40–100% FLASH and ~+40–110% ITCM per effect.** ShapeShifter tops
  out at **73 KB ITCM** single-effect. The full 26-effect roster at -O3 needs
  ~363 KB ITCM (12 of 32 FlexRAM banks) and overflows the phantasm FlexRAM split
  (`platformio.ini` §4.1) — which is exactly why the shipping image is `-Os`. The
  -O3 numbers here are only reachable as single-effect profile images.

## Reproduce

```
# -Os (shipping): just profile <Effect> [seconds]
# -O3 twin:
$env:PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=<Effect> -D HS_PROFILE_WINDOW=64'
pio run -e profile_o3 -t upload      # (or -e profile for -Os)
python tools/profile_capture.py --seconds 75 --out build/profile_<effect>_o3.log
```

Per-effect detail, counter trees, ISR figures, and the `## -O3 vs -Os` sections
are in the individual reports under [`Os/`](Os/README.md) and [`O3/`](O3/README.md).
Raw logs: `build/profile_<effect>.log` (-Os) and `build/profile_<effect>_o3.log`.
