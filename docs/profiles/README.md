# On-device effect profiles — Teensy 4.0, segmented mode (2026-07-14)

On-device timing profiles for **all 26 roster effects**, captured on a
bench-attached Teensy 4.0 running the shipping Phantasm configuration
(`POVSegmented<288, 4, 480>`, board = segment 0 master, newlib-nano, DMA LEDs,
flywheel + DMA ISRs live), via the `HS_PROFILE` cycle-counter harness. Each
effect renders one **quadrant** ≈ **10,368 px**; a display window is **62.5 ms**,
so cadence quantizes to 16 fps (1 window), 8 fps (2), 5.3 fps (3), etc.

The shipping image is **selective -O3**: an `-Os` base with the landed `HS_O3`
hot-loop regions (docs/selective_o3_spec.md) compiled at `-O3 -ffast-math`, the
rest at `-Os` so the full 26-effect roster still fits FlexRAM. Profiles are
captured at **two build configurations**:

- [**`shipping/`**](shipping/README.md) — the **shipping** image, via the
  `-Os` `profile` env (the `HS_O3` regions activate on the `-Os` device build, so
  this env measures the shipping config by construction). Since the regions
  landed 2026-07-15, `just profile <Effect>` measures this configuration.
  RingSpin and DisplacementField — the only two effects a region reaches — have
  dedicated three-way (`-Os` / selective / global-`-O3`) reports here.
- [**`O3/`**](O3/README.md) — **global `-O3`** (`-ffast-math`), via the
  `profile_o3` env (identical to the shipping profile except `-O3` replaces
  `-Os` **everywhere**). Global `-O3` is **not shippable at roster scale** — a
  single-effect image fits, but the full roster overflows ITCM (~363 KB /
  12 banks). These runs isolate the optimization level as the only variable:
  they are the **ceiling** the shipping image is measured against.

For the 24 effects **no `HS_O3` region reaches**, the shipping selective-O3
image is byte-identical to the `-Os` base, so their `-Os` numbers below **are**
their shipping numbers. RingSpin and DisplacementField (◆) bank most of the
global-`-O3` ceiling through their regions — see [`shipping/`](shipping/README.md).

**Point-in-time snapshots** (regenerate with `just profile <Effect>`; numbers
age as the render code moves). Reports + the `HS_PROFILE` instrumentation +
the `profile_o3` env are working-tree only.

## Global -O3 ceiling vs the -Os base — render time (ms) and speedup

**Render** = the effect's own per-frame work (frame minus the `*_buffer_wait`
display-sync idle), at the steady/representative window, in **ms**. **Base -Os**
is the shipping render for every effect except the two ◆ rows (see above).
`O3×` = base -Os render ÷ global-`-O3` render — the compiler ceiling. Size Δ is
the single-effect image growth at global `-O3`. Ordered heaviest base render
first. **●** marks the 11 effects re-profiled 2026-07-14 after the plot
column-cull batch (9ac8cebd/62450701/708d4b9b); the other 15 carry their
original 2026-07-14 captures (the batch does not touch their `Scan::` path).
DisplacementField was re-profiled 2026-07-15 after the fused-stack batch
(ec61faa7..7d50b672), which crosses its whole cycle to 16 fps at `-O3`.
RingSpin was re-profiled 2026-07-15 after the fused trail-scan landing plus
the slim RingGroup rework (730b6d2f): base per-blend cost dropped ~29%
(~14 fps effective, up from ~9.4), and global `-O3` hard-locks 16 fps.

| Effect | Dominant scope | Render base -Os (ms) | Render -O3 ceiling (ms) | O3× | FLASH Δ | ITCM Δ |
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
| ●◆ [DisplacementField](shipping/profile_displacementfield_teensy_2026-07-15.md) | fused ring-stack raster | 82.6 | 59.5 | 1.39× | +56% | +77% |
| ●◆ [RingSpin](shipping/profile_ringspin_teensy_2026-07-15.md) | slim ring-group raster | 49.6 | 39.1 | 1.27× | +89% | +124% |
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

◆ **RingSpin** and **DisplacementField** ship with an `HS_O3` region, so their
shipping render is **below** the base `-Os` column above; the selective-O3
captures (linked) close ~90% and ~88% of the base→ceiling gap respectively:
RingSpin hot scope 46.3 → **34.6** ms (ceiling 33.3), DisplacementField
NOISE-dwell 58.5 → **45.7** ms (ceiling 44.0). Every other row has no region,
so its base `-Os` render **is** its shipping render.

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
`Plot::`-rasterized effects (● above) at both configs shows the win is **highly
uneven** — it only pays where a large share of the frame is plot-edge raster
*and* the edges span columns:

- **Real wins**: **DreamBalls** (wireframe — long geodesic/planar edges cross the
  quadrant; base render ~145 → 109 ms) and **MindSplatter** (particle trails; base
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
  (global `-O3`), not the cull.

The **ceiling** per-effect speedup (~1.4–1.6×) is orthogonal to the batch — it
is the compiler, and is unchanged by the cull.

## What the -O3 ceiling buys (and how selective -O3 banks it)

- **Render is 1.14×–1.82× faster at global `-O3`**, clustering around **~1.4×**.
  The biggest wins are the memory/loop-bound compositors and rasterizers (Flyby
  shader 1.82×, MeshFeedback feedback flush 1.79×, HopfFibration trail 1.59×,
  DreamBalls/Comets ~1.55×); the smallest are the branch-heavy or already-lean
  shaders (Liquid2D 1.14×, GnomonicStars star-transform 1.19×, Raymarch's
  divergent ray-march 1.23×).
- **Per-pixel inner loops benefit most.** Wherever a `filter_blend` leaf is
  visible, its per-blend cost roughly **halves** at `-O3` (~1.7×–2.4×: GnomonicStars
  231→96 cyc, RingSpin/DistortedRing ~2×, SphericalHarmonics 1.95×) — `-O3`'s
  vectorization/unrolling pays off most in the tight SDF distance + blend loops.
  This is exactly the class of loop selective -O3 wraps: the RingSpin/
  DisplacementField `HS_O3` regions and the shared `Pipeline` blend sink recover
  ~90% of the per-blend win for a fraction of the ITCM (see
  [`shipping/`](shipping/README.md)).
- **7 effects cross a cadence tier at the ceiling**: MeshFeedback, Flyby, Voronoi
  (all 8→16 fps), DreamBalls & FlowField (5→8 fps), HopfFibration (jitter→steady
  16), IslamicStars (mid-shapes 8→16); HankinSolids & RingSpin cross partially.
  The rest keep their cadence — the saved time falls into `*_buffer_wait`
  display-sync idle because the render either was already well under a window
  (the light 16 fps effects) or still overruns one after the speedup (the
  heaviest raster/RD/shader effects).
- **Global `-O3` costs ~+40–100% FLASH and ~+40–110% ITCM per effect.**
  ShapeShifter tops out at **73 KB ITCM** single-effect. The full 26-effect
  roster at global `-O3` needs ~363 KB ITCM (12 of 32 FlexRAM banks) and
  overflows the phantasm FlexRAM split (`platformio.ini` §4.1) — which is why the
  shipping image is `-Os`-based, wrapping only the measured hot loops at `-O3`
  (selective -O3). The global `-O3` numbers here are only reachable as
  single-effect profile images.

## Reproduce

```
# shipping (selective -O3): just profile <Effect> [seconds]
# global -O3 ceiling:
$env:PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=<Effect> -D HS_PROFILE_WINDOW=64'
pio run -e profile_o3 -t upload      # (or -e profile for the shipping selective-O3 config)
python tools/profile_capture.py --seconds 75 --out build/profile_<effect>_o3.log
```

Per-effect detail, counter trees, and ISR figures for the two region effects are
in [`shipping/`](shipping/README.md); the global-`-O3` ceiling reports and
their `## -O3 vs base` sections are under [`O3/`](O3/README.md).
Raw logs: `build/profile_<effect>.log` (shipping) and `build/profile_<effect>_o3.log`.
