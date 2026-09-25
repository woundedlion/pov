# Core Subsystems

Section 7 of the [Holosphere README](https://github.com/woundedlion/pov/blob/master/README.md).

The [pullback stage-families specification](specs/pullback_stage_families_spec.md)
describes the shipped stage model.

## Contents

- [7.0 The Shader Interface](#70-the-shader-interface)
  - [The Fragment](#the-fragment)
  - [Shader Signatures](#shader-signatures)
  - [Register Conventions by Rasterizer](#register-conventions-by-rasterizer)
- [7.1 SDF Shapes and the Scan Rasterizer](#71-sdf-shapes-sdfh-and-the-scan-rasterizer-scanh)
  - [SDF Shape Primitives](#sdf-shape-primitives-sdfh)
  - [Volumetric Shapes](#volumetric-shapes-sdfvolumeh)
  - [CSG Operations](#csg-operations-sdfh)
  - [Scan Rasterization Primitives](#scan-rasterization-primitives-scanh)
  - [Near-Pole Azimuthal LOD](#near-pole-azimuthal-lod)
- [7.2 The Curve Rasterizer](#72-the-curve-rasterizer-ploth)
  - [Sampling Policy](#sampling-policy)
  - [Plot Primitives](#plot-primitives)
- [7.3 The Animation System](#73-the-animation-system-animationh)
  - [Animation Types](#animation-types)
  - [Orientation and Motion Blur](#orientation-and-motion-blur)
  - [OrientationTrail](#orientationtrail)
  - [VectorTrail and QuantizedVectorTrail](#vectortrail-and-quantizedvectortrail)
  - [`tween` and `deep_tween`](#tween-and-deep_tween)
  - [Animations and Mutable State](#animations-and-mutable-state)
- [7.4 Geometry Transformers](#74-geometry-transformers-transformerh)
  - [Displacement Fields](#displacement-fields)
  - [Pool Lifecycle](#pool-lifecycle)
  - [Standalone Utilities](#standalone-utilities)
- [7.5 Memory Architecture](#75-memory-architecture-memoryh-memorycpp)
  - [Compaction with `Persist<T>`](#compaction-with-persist)
  - [Additional Data Structures](#additional-data-structures)
- [7.6 The Color System](#76-the-color-system-corecolorcolorh)
  - [Palette Types](#palette-types)
  - [OKLCH Perceptual Color](#oklch-perceptual-color)
  - [The Gamut Boundary Grid](#the-gamut-boundary-grid)
  - [Palette Modifiers](#palette-modifiers)
  - [Additional Palette Types](#additional-palette-types)
  - [Recipe-Compiled Palettes](#recipe-compiled-palettes)
  - [Palette Cycling](#palette-cycling)
- [7.7 The Mesh System](#77-the-mesh-system-coremesh)
  - [Core MeshOps](#core-meshops-coremeshmeshh)
  - [Conway Operators](#conway-operators-conwayh)
  - [Hankin Pattern System](#hankin-pattern-system-hankinh)
  - [Solids Library](#solids-library-solidsh-solid_generatorsh)
- [7.8 Generators](#78-generators-memoryh)
- [7.9 The Preset System](#79-the-preset-system-controlchoreographyh)
- [7.10 Hardware Drivers](#710-hardware-drivers-dma_ledh-pov_singleh-pov_segmentedh)
  - [DMA LED Controller](#dma-led-controller-dma_ledh-hd107s_frameh-dma_led_coreh-dma_led_controllerh)
  - [Single-Teensy POV Driver](#single-teensy-pov-driver-pov_singleh)
  - [Multi-Teensy Segmented POV Driver](#multi-teensy-segmented-pov-driver-pov_segmentedh)
  - [Frame Sync Protocol: 1-Wire Signal Datasheet](#frame-sync-protocol-1-wire-signal-datasheet)
- [7.11 Mathematical Kernels](#711-mathematical-kernels-coremath)
- [7.12 Spatial Queries](#712-spatial-queries-corespatial)

---

## 7.0 The Shader Interface

All rasterizers — SDF scanline, curve plotting, mesh, volumetric, and full-screen shader — share a common shading model based on the `Fragment` struct and two function signatures.

### The Fragment

A `Fragment` (`render/shading.h`) is the data packet exchanged between rasterizers and shaders. It carries the pixel position, four general-purpose float registers, and the output color:

```cpp
struct Fragment {
  math::Vector pos;              // Position (typically a unit vector on the sphere)
  float v0 = 0.0f;        // Register 0: normalized progress t (0–1)
  float v1 = 0.0f;        // Register 1: arc length / distance
  float v2 = 0.0f;        // Register 2: stroke coverage / face ID
  float v3 = 0.0f;        // Register 3: auxiliary
  float size = 1.0f;      // Size metric for normalization
  float age = 0.0f;       // Age (for trail decay / motion blur)
  Color4 color = Color4(0, 0, 0, 0); // Output: shader writes RGBA here; defaults to transparent black (alpha 0.0)
};
```

The registers are *inputs* — populated by the rasterizer before the shader runs. The shader reads them and writes `color`. The rasterizer then forwards the color through the filter pipeline to the canvas.

### Shader Signatures

Two shader types are defined as zero-allocation `FunctionRef` callables (`concepts.h`):

| Signature | Type | Role |
|---|---|---|
| `FragmentShaderFn` | `void(const math::Vector &, Fragment &)` | Per-pixel/per-sample shader. Receives the world position and a pre-populated Fragment; writes `color`. Called for every rasterized point. |
| `VertexShaderRef` | `void(Fragment &)` | Per-vertex or per-pixel-center shader. Runs once before sub-sampling to set up expensive shared state in the Fragment registers. Optional on the `Plot::` primitives (a null callable is skipped); required by `Scan::Shader::draw`'s split vertex/fragment overload, which traps on a null one. |

`FunctionRef` is a non-owning, non-allocating type-erased callable (similar to `std::function_ref` from C++26). It captures a pointer to any lambda, functor, or function pointer with zero heap allocation — critical for ISR-safe code on Teensy.

Effects pass lambdas that capture their state:

```cpp
auto shader = [&](const math::Vector &p, Fragment &f) {
    float t = f.v0;           // read: normalized progress from rasterizer
    f.color = palette.get(t); // write: color from palette lookup
};
Scan::Ring::draw<W, H>(pipeline, canvas, basis, radius, thickness, shader);
```

### Register Conventions by Rasterizer

Each rasterizer family populates the Fragment registers with a consistent convention. Shaders can rely on these semantics:

**SDF Scanline Path** (`Scan::Ring`, `Scan::Star`, `Scan::PlanarPolygon`, `Scan::Flower`, `Scan::Line`, `Scan::Mesh`):

| Register | Source | Meaning |
|---|---|---|
| `v0` | `DistanceResult.t` | Normalized azimuth (0–1) for `Scan::Ring` and `Scan::Star`; normalized radial position for `Scan::PlanarPolygon`, `Scan::SphericalPolygon` and `Scan::Flower` (polar angle over the shape radius, so it passes 1 outside the body); unused (0) for `Scan::Line` and `Scan::Mesh` faces |
| `v1` | `DistanceResult.raw_dist` | Producer-specific: unsigned centerline/arc distance for rings and lines, polar angle for polygons and stars, or antipode scan distance for flowers; `Scan::Mesh` faces carry the signed edge distance instead — negative inside the face, in gnomonic plane units on a small (`linear_dist`) face and in radians on a large one — which `fragment_edge_dist()` turns into normalized inward depth as `-v1 / size` |
| `v2` | Set by rasterizer | Stroke AA coverage (0–1, also applied by Scan at plot time), 0 for solid shapes, or face index for `Scan::Mesh` (but see the per-face setup note below) |
| `v3` | `DistanceResult.aux` | Auxiliary — shape-dependent secondary parameter (0 when unused, including faces) |
| `size` | `DistanceResult.size` | Stroke half-width for stroke shapes, or radius or apothem for filled shapes (mesh `Face` floors it at 0.25× the face circumradius, so on a sliver face — whose true inradius approaches zero — the reported size overstates it without bound) |

With a per-face setup callback, `Scan::Mesh::draw_specialized` runs its minimal-fragment loop: only `v1` is refreshed per pixel. Other fragment inputs are unavailable and are poisoned with NaN in debug builds; shaders must not read `v2` or `size` through `mesh_face_index()` or `fragment_edge_dist()`. The face index and face size reach the shader through the setup callback instead. `Scan::Mesh::draw` does not accept this callback.

The `DistanceResult` struct is returned by each SDF shape's `distance<ComputeUVs>()` method. Distances and `size` use radians except for small `SDF::Face` shapes (inradius < 0.2), which use gnomonic tangent-plane units. The per-producer register table in `core/render/sdf/common.h` defines `t` and `raw_dist`:

```cpp
struct DistanceResult {
  float dist;        // Signed distance (negative = inside)
  float t;           // Shape-dependent parameter or angle
  float raw_dist;    // Unsigned / supplementary distance
  float aux;         // Auxiliary (0 for every current producer)
  float size = 1.0f; // Size metric
};
```

**Curve Plot Path** (`Plot::Line`, `Plot::Multiline`, `Plot::Ring`, `Plot::Polygon`):

| Register | Meaning |
|---|---|
| `v0` | Path progress (0.0 → 1.0 along the full curve) |
| `v1` | Cumulative arc length in radians |
| `v2` | Index of the segment's start vertex, interpolated toward the next index across the segment — only the control points land on whole numbers, interior samples are fractional; `Plot::Ring` writes the strided grid index (stride up to W/8), `Plot::Line` has no interior vertices and writes 0 on every sample, and `Plot::Mesh` writes a constant edge index instead |
| `v3` | Inherited from control-point Fragment (user-defined) |

Plot primitives interpolate registers between control-point Fragments via `Fragment::lerp_registers()`. The vertex shader, if provided, runs once per control point before rasterization. For `Plot::Polygon<Plot::PlanarProjection>`, `Plot::Star<Plot::PlanarProjection>`, and `Plot::Flower`, the rasterizer re-derives `v0`/`v1` from the rendered azimuthal-equidistant arc — which bows longer than the great-circle chord between vertices — so both stay consistent with the drawn position.

**Full-Screen Shader Path** (`Scan::Shader`):

Registers are not pre-populated — the shader receives only `pos` (reconstructed from pixel coordinates). The single-callback overload provides a `Color4(const Vector &)` interface. The two-callback overload separates per-pixel vertex setup from per-subsample fragment evaluation.

**Volumetric Path** (`Scan::Volume`):

The fragment shader receives `pos` set to the closest local-space hit point (in the SDF's coordinate frame) and `size` set to the closest signed distance. No register convention — the shader computes lighting from the local-space position directly.

## 7.1 SDF Shapes (`sdf.h`) and the Scan Rasterizer (`scan.h`)

The rendering pipeline splits shape definitions from rasterization. `sdf.h` defines the SDF shape primitives, each implementing three methods:

1. **`get_vertical_bounds()`** — analytic tight bounding box in pixel-Y space (phi angle range). Only rows within this range are scanned.
2. **`get_horizontal_intervals(y, out)`** — analytic scanline intervals per row. Called per row to skip empty columns without evaluating the distance function.
3. **`distance<ComputeUVs>(p, result)`** — signed distance from a sphere-surface point `p` to the shape boundary, plus texture coordinate and auxiliary data in `DistanceResult`.

`scan.h` contains `Scan::rasterize()`, which drives the scanline loop and anti-aliasing, plus convenience wrappers that pair SDF shapes with the rasterizer.

`sdf.h` is an umbrella over the headers in `core/render/sdf/`: the substrate every shape shares (azimuth intervals, row bounds, `DistanceResult`, the cap/annular span-emission helpers) in `common.h`, the polygon, star, flower and line leaves in `core/render/sdf/shapes.h`, the ring leaves in `rings.h`, the CSG operators in `csg.h`, `SDF::Face` with its congruence-class LUT in `face.h`, and the volumetric family in `core/render/sdf/volume.h`. Including `sdf.h` pulls in all six, so nothing outside needs to name them.

The `process_pixel` function applies anti-aliasing based on shape type:
- **Solid shapes**: quintic smoothstep over a 2-pixel AA band centered on the edge (`-pixel_width <= d <= pixel_width`). Full interior pixels (`d < -pixel_width`) skip AA math entirely. `pixel_width` is the compile-time constant `2π/W` — the angular width of one *equatorial* pixel — so the band is a fixed angular thickness at every latitude, and near the poles (where columns converge) it spans more than two columns. At 288×144 the row and equatorial column arcs nearly match. At 96×20 a row spans about 2.18 column arcs, so a horizontal edge has a narrower AA band in row units than a vertical edge has in column units.
- **Strokes**: opacity falloff across the full stroke thickness.

### SDF Shape Primitives (`sdf.h`)

| Shape | Description |
|---|---|
| `SDF::Ring` | Geodesic circle at a given radius and thickness |
| `SDF::DistortedRing` | Ring with per-azimuth radius perturbation via a callback |
| `SDF::PlanarPolygon` | Regular N-gon in the tangent plane of a basis vector |
| `SDF::SphericalPolygon` | Regular N-gon with geodesic (great-circle) edges |
| `SDF::Star` | N-pointed star using the standard inradius/circumradius construction |
| `SDF::Flower` | Inverted star (N-petal flower shape from the antipodal perspective) |
| `SDF::Line` | Geodesic line segment between two sphere-surface points |
| `SDF::Face` | Planar polygon face (used for mesh rendering) |

The table covers the effect-facing shapes. `sdf.h` also holds internal specializations that only the matching `scan.h` wrapper constructs — `SDF::FlatDistortedRing`, an undisplaced `DistortedRing` with an exact polar centerline distance, is instantiated by `Scan::DistortedRing::draw_flat` and never named by an effect.

### Volumetric Shapes (`sdf/volume.h`)

The 3D family marched by `Scan::Volume` lives in its own header. It shares the `SDF` namespace but not the scanline contract above: these shapes return a plain `float` distance in Cartesian ray-space, have no vertical bounds or horizontal intervals, and are reached only from the march loop. Everything it declares is a struct or a template, so an effect that draws only 2D shapes emits none of it.

| Shape | Description |
|---|---|
| `SDF::Torus` | 3D volumetric torus SDF with configurable major/minor radii (Cartesian ray-space, not a 2D sphere-surface shape) |
| `SDF::Warp::Twist` | Domain warp composed with a volumetric SDF via `SDF::WarpedVolume<Shape, Warp>` — e.g. `WarpedVolume<Torus, Warp::Twist>` twists a torus by oscillating Y around the ring azimuth, with an analytic Lipschitz bound for safe sphere-tracing (used by Raymarch) |

### CSG Operations (`sdf.h`)

Shapes can be combined using Constructive Solid Geometry:

```cpp
SDF::Union<Ring, Line>           // min(d_A, d_B)
SDF::SmoothUnion<Line, Line>     // smooth minimum with blending radius
SDF::Subtract<Ring, PlanarPolygon> // max(d_A, -d_B)
SDF::Intersection<Ring, Line>    // max(d_A, d_B) with interval intersection
SDF::AngularRepeat<Shape>        // N-fold angular repetition around an axis
```

`Union`, `SmoothUnion` and `Intersection` require both children to share
`is_solid`; a solid+stroke mix routes one winner through the wrong AA branch and
is rejected at compile time. `Subtract` tracks the minuend's solidity, so a
solid carved by a stroke (or the reverse) is legal. `SmoothUnion` adds a second
compile-time rule, `blends_smoothly`: its weld term needs a real signed distance
outside the surface, so `Ring`, `DistortedRing`, `FlatDistortedRing` and `Face`
— which clamp to a far sentinel past their reject band — cannot be its children,
directly or nested inside a combinator.

### Scan Rasterization Primitives (`scan.h`)

Convenience structs that construct an SDF shape and rasterize in a single `draw()` call:

| Primitive | Description |
|---|---|
| `Scan::Ring` | Rasterizes a ring (from `SDF::Ring`) |
| `Scan::RingGroup` | Fused single-pass rasterizer for a small group of rings — one scan over the union band paints every member in slot order, so the per-row interval math runs once instead of per ring. Fragments carry position, stroke coverage and size only (no UVs, no raw distance) |
| `Scan::Circle` | Disc (ring with radius-wide thickness) — stroke coverage ramps quintically from center to rim, and the shader styles it from register 2 |
| `Scan::Point` | Dot at a sphere-surface position, with the same center-to-rim coverage ramp |
| `Scan::Line` | Geodesic line segment between two points |
| `Scan::Star` | N-pointed star shape |
| `Scan::Flower` | N-petal flower shape |
| `Scan::DistortedRing` | Ring with per-azimuth radius perturbation |
| `Scan::DistortedRingStack` | Fused single-pass rasterizer for an evenly spaced same-axis stack of distorted rings — the per-pixel frame every ring shares is computed once and the candidate rings fall out of its polar angle by arithmetic; the shader takes a ring slot alongside the fragment |
| `Scan::PlanarPolygon` | Regular N-gon in the tangent plane |
| `Scan::SphericalPolygon` | Regular N-gon with geodesic (great-circle) edges |
| `Scan::Mesh` | Rasterizes all faces of a `MeshState` |
| `Scan::Shader` | Full-screen per-pixel shaders with configurable SSAA (super-sample anti-aliasing), across four entry points. `draw(canvas, shader)` takes a single fragment shader; `draw_cached(canvas, shader)` provides the same typed draw with its traversal placed in cached flash and is the path composed effects use. `draw(canvas, fragment_shader, vertex_shader)` separates a per-pixel vertex shader (called once at pixel center) from a per-subsample fragment shader (called SAMPLES×), so expensive per-pixel work is computed once — both callables are required, and a null one traps. The SSAA reaction-diffusion effects use `draw_grid`. `draw_grid(canvas, vertex_shader, pixel_shader)` hands the seeded fragment and the row's sub-pixel grid to a templated pixel shader that owns the sampling and returns the finished pixel. |
| `Scan::TransformedVolume` | Wraps an SDF shape with a world-space position and orientation quaternion for volumetric rendering |
| `Scan::Volume` | Volumetric ray-marcher that steps along the view direction through a `TransformedVolume`, applying a fragment shader at the hit point with configurable step count and AA width |

### Near-Pole Azimuthal LOD

A row at colatitude φ has horizontal pixel pitch `sin(φ)` times the vertical, so `1/sin(φ)` columns share one physical LED footprint and need only one shade between them. The scan walk offers those columns as a block of `pole_lod_aggressiveness / sin(φ)` (`core/render/render_policy.h`, clamped to `POLE_LOD_MAX_RUN = 32`), and the sink settles the whole block from one probe wherever the probe can vouch for it. Only full canvas-aligned blocks are offered, so an offer never straddles two blocks and a settled column always takes its shade from its own block's anchor. A block truncated by a clip or span edge goes per column instead, so the columns beside a segment seam shade at full resolution rather than from the anchor the neighbouring segment would have used.

`pole_lod_aggressiveness` is a hardware-calibrated knob, not a derived constant: the true masking width depends on the LED's angular size and the per-column exposure. 1.0 tracks the footprint exactly; smaller values stay inside it; 0 makes every offer one column and the walk bit-identical to an undecimated one. It defaults to 0 (`HS_POLE_LOD_DEFAULT`). Firmware compiles it in as a `constexpr` with no setter — at the default, the decimation branches fold away entirely — while host and WASM builds keep it mutable so it can be tuned live (§10.2 `setPoleLod`).

The knob reaches the walk, not every primitive. `Scan::RingGroup` and `Scan::DistortedRingStack` replace the per-ring walk with one fused scan over the group's union band and shade every column of it, so raising the knob leaves them undecimated. Their equivalence to rasterizing the members one by one is stated at aggressiveness 0 for exactly that reason.

## 7.2 The Curve Rasterizer (`plot.h`)

For drawing lines, curves, and paths, the `Plot` namespace provides a geodesic/planar rasterizer with adaptive step size. Each sub-step is sized from the curve's full 2-D screen-space speed (`sqrt(vx² + vy²)`, combining longitudinal and latitudinal motion), so samples land roughly one pixel apart everywhere on the curve regardless of latitude. The step is clamped to keep the equator near one sample per column and floored near the poles — where screen speed diverges — so pole oversampling stays bounded.

```cpp
Plot::Line::draw<W, H>(pipeline, canvas, start, end, fragment_shader);
Plot::Multiline::draw<W, H>(pipeline, canvas, vertices, fragment_shader);
```

`Plot::Multiline` accepts a `Fragments` array (an arena-backed `ArenaVector<Fragment>`), while `Plot::Line` accepts its two `Fragment` endpoints. The parametric primitives (`Ring`, `Polygon`, `DistortedRing`, `Star`, and `Flower`) take their geometric parameters and sample into a `Fragments` array internally. Each fragment carries position, texture registers (v0–v3), age, and color.

- **Edge interpolation** — how consecutive fragments are joined. *Geodesic* (the default) walks the great-circle arc between endpoints; *planar* interpolates along an azimuthal-equidistant straight line in a basis's tangent plane (for effects that live in a 2D local space). This is selected by whether a **planar basis** is supplied to the draw call (`null` ⇒ geodesic).

### Sampling Policy

`rasterize` takes its compile-time behavior as one `RasterConfig` NTTP (`rasterize<W, H, RasterConfig{.single_pass = true}>`), whose `sampling_policy` field sets the adaptive sample density. `DEFAULT` targets `SCREEN_STEP_PX` (0.9 px) and compiles the alternative away; `BALANCED` always trades samples for speed; `SELECTABLE` defers the choice to `RasterOptions::balanced_sampling`, so one instantiation serves both and the policy is picked per draw call. Only the single-pass rasterizer reads it — the cached-replay path always samples at the default density.

`RasterOptions` remains a by-value aggregate. `RasterLoop::closed(seam)` attaches optional seam registers to a closed path. `RasterProjection::planar(basis)` selects a planar chart; `RasterProjection::geodesic(flags)` accepts a bounded visibility span. `PointProjections` pairs equal-sized row and column arrays, with `paired(rows, cols)` for dynamic spans. The rasterizer checks span lengths against the actual polyline.

Balanced sampling stretches each adaptive step by `BALANCED_SCREEN_STEP_PX / SCREEN_STEP_PX` (1.25×), clamped to one base step (2π/W) and left exact below the pole floor (`MIN_POLE_SCALE * BALANCED_POLE_GUARD_SCALE` base steps), where spacing is already at its minimum. Two consequences:

- **Emitted alpha changes.** Sparser samples lay down less coverage per unit arc, so each fragment's alpha is scaled by `balanced_sample_alpha()` — gain `1 + (ratio - 1) * (0.88 - 0.20 * alpha)`, saturating at 1, with `ratio` the balanced step over the default step. The gain shrinks as alpha rises because opaque samples compound less; it is a linear fit to source-over accumulation, close to exact below alpha 0.4 and over-boosting above it, so a stroke past alpha ~0.85 saturates at 1 and loses its soft edge. A balanced draw is not pixel-identical to a default one.
- **Step evaluation is reused, on planar edges only.** Where the walk is locally straight and clear of the poles (tangent dot > 0.995, step change under 10%, step under 0.9 base steps), the next sample recomputes position only and carries the previous step forward, skipping the tangent and the screen-velocity step. The reuse needs the monotonic position-only entry that only `PlanarEdgeSampler` exposes; a geodesic edge takes the sparser steps and the alpha gain without it.

`ShapeShifter` is the sole caller: `SELECTABLE`, on for policy-selected stars at 32 or more contours and off for every other primitive. Of those, only the dense planar star's pole-crossing edges collect the step reuse; the spherical star's edges are geodesic.

### Plot Primitives

| Primitive | Description |
|---|---|
| `Plot::Line` | Geodesic line segment between two points |
| `Plot::Multiline` | Connected line strip from a sequence of fragments |
| `Plot::Ring` | Circle rasterized as a plotted polyline |
| `Plot::Polygon<Plot::PlanarProjection>` | Regular N-gon in the tangent plane |
| `Plot::Polygon<Plot::GeodesicProjection>` | Regular N-gon with geodesic (great-circle) edges |
| `Plot::DistortedRing` | Ring with per-azimuth radius perturbation via callback |
| `Plot::Star<Projection>` | N-pointed star with planar or geodesic edges |
| `Plot::Flower` | N-petal flower shape |
| `Plot::Mesh` | Wireframe mesh rendering with edge deduplication. The dedup bitset holds `DEDUP_CAPACITY = 128` vertices; a larger mesh traps as its faces are walked — at setup for `extract_edges()`, but every frame at render time for `draw()`. Conway operators pass 128 vertices within two or three ops, so a wireframe fed from an `OpLeg` chain must keep its vertex count under the cap |
| `Plot::ParticleSystem` | Particle trail rendering from `QuantizedVectorTrail` history |

## 7.3 The Animation System (`animation.h`)

The `Timeline` class manages a list of running `IAnimation` objects. Each frame, `timeline.step(canvas)` advances all active animations. Finished animations are removed; repeating animations are rewound. All animation types inherit from `AnimationBase` and support method chaining via `.then()` for sequencing.

Animation pause is opt-in per timeline event, not a global stop. Effects schedule parameter drivers and preset choreography with `timeline.add_pausable(..., &anims_paused)`; while the flag is set, both the animation step and any pending start delay are frozen. Passing a pause pointer to an animation constructor freezes only its `step()` call, so the timeline event's start delay still elapses. Events added with `add()` and motion advanced directly by `draw_frame()` continue to run, which lets the GUI pause animated controls without stopping ambient motion.

`animation.h` defines the contract every animation implements — `IAnimation`, the CRTP `AnimationBase`, and `Animation::Space` — and then includes nine fragment headers grouped by what they animate:

| Header | Subject | Contents |
|---|---|---|
| `timers.h` | Callbacks on a clock | `RandomTimer`, `PeriodicTimer` |
| [params.h](../core/animation/params.h) | A caller-owned parameter, written each frame | `Transition`, `Mutation`, `Progress`, `Driver`, `Lerp`, `ColorWipe`, the `Mobius*` family, `Ripple`, `Noise`, `BallDrop`, `NoiseProduct` |
| `motion.h` | An `Orientation` driven through space | `Path`/`ProceduralPath`, `Motion`, `Rotation`, `RandomWalk` |
| `trails.h` | Recorded history | `Trail` and its `OrientationTrail`/`VectorTrail` aliases — index 0 is the oldest snapshot and `length()-1` the newest, the ordering the JS simulator mirrors — plus `QuantizedVectorTrail`, the `TrailBody` per-body aggregate, and the `tween`/`deep_tween` traversals |
| `sprites.h` | Visible things | `Sprite`, `Particle`/`ParticleSystem` |
| `timeline.h` | Scheduling | `TimelineEvent`, `Timeline` |
| `opleg.h` | One Conway-chain morph leg, swept per frame | `OpLeg` |
| `segue.h` | How one mesh hands the sphere to the next | the `Segue` policies |
| `carousel.h` | Two persistent mesh slots + arena compaction | `MeshCarousel` |

The fragments compile only inside `animation.h` (a direct include fails with an `#error`); consumers include `animation.h` alone. Their types span three scopes rather than one: the animations themselves are in `namespace Animation`; the transition policies are in `namespace Segue`; and `TimelineEvent`/`Timeline`, `MeshCarousel`, `Path`/`ProceduralPath`, and the `tween`/`deep_tween` traversals sit at global scope.

### Animation Types

| Type | Description |
|---|---|
| `Rotation<W>` | Quaternion rotation of an `Orientation` around an axis, with optional repeat. Supports World and Local coordinate spaces. |
| `RandomWalk<W>` | Continuously perturbs an `Orientation` with smoothly changing random angular velocity driven by OpenSimplex2 noise. Configurable via `Options` presets (Languid, Energetic). |
| `Motion<W, CAP>` | Moves an `Orientation` along a `Path` or `ProceduralPath` (the path is a constructor argument; `CAP` is the orientation sub-frame capacity, default 4) |
| `Sprite` | Calls a draw function over a duration with fade-in and fade-out envelopes |
| `PeriodicTimer` | Fires a callback at regular intervals (once or repeatedly) |
| `RandomTimer` | Fires a callback after a random delay within a min/max range |
| `Transition` | Smoothly interpolates a float variable from its current value to a target over a duration with easing |
| `Mutation` | Applies a custom scalar function to a float variable over time with easing |
| `Progress` | Invokes a caller-supplied `void(float)` callback once per frame with eased progress, leaving the caller to write whatever state it drives |
| `Driver` | Continuously increments a float variable each frame (optionally wraps at 0..1) |
| `Lerp` | Type-erased interpolation between any `T` that implements `lerp(start, target, t)`. The caller owns start, subject, and target data; Lerp holds pointers and a type-erased lerp function. |
| `ColorWipe` | Smoothly interpolates a `GenerativePalette` between caller-owned start and target snapshots that must remain unchanged and outlive the animation |
| `ParticleSystem<W, CAPACITY>` | Physics simulation with emitters, attractors, friction, gravity. Particles have `QuantizedVectorTrail` history for trail rendering. |
| `Ripple` | Animates a `RippleParams` to expand a Ricker wavelet across the sphere |
| `MobiusWarp` | Animates `MobiusParams` to apply and release a Möbius transformation |
| `MobiusWarpCircular` | Animates `MobiusParams` for a circular warp that stays warped throughout, suitable for repeating effects |
| `MobiusWarpEvolving` | Continuously modulates `MobiusParams` over multiple frequencies for a non-repeating, evolving warp |
| `MobiusFlow` | Animates `MobiusParams` for a continuous loxodromic flow |
| `Noise` | Animates `NoiseParams` over time for flowing distortion fields |
| `BallDrop` | Animates a `BumpParams` to drop one spherical-cap bump from the north pole to the south along a fixed meridian, ramping the footprint envelope so the bump emerges from and vanishes into the poles |
| `NoiseProduct` | Integrates the time axis of a `NoiseProductParams` (`time += speed` per frame) to flow a two-octave product-noise field; perpetual |
| `OpLeg` | Animates one leg of a mesh operator chain: a Conway-operator parameter sweep along a `ConwayGraph` edge or recipe step, a hankin contact-angle sweep on a fixed seed, a relax or medial slerp, or a gated partition swap. Each frame it rebuilds the swept mesh in scratch, compiles it, attaches the leg's hoisted congruence classification, pre-blends the (from, to) palette ramps at the leg's crossfade weight, and hands the mesh to a draw callback — exactly one mesh drawn per frame. Each leg's `.then()` completion handler schedules the next, so a run of legs walks a whole morph path. |
| `MeshCarousel<SegueT>` | Double-buffered mesh transition system, parameterized on a compile-time segue policy (`namespace Segue`) that owns the transition's animation scheduling via `schedule_segue()` and shapes its rendering through phase-driven hooks (`opacity`/`fill`/`grade`, plus optional `warp`, per-face sweep ordering, and `retarget` for per-transition anchors). Manages a pair of `MeshState` buffers and exposes the front index (`front_index`/`set_front`); the flip itself belongs to the effect, which builds the incoming mesh into the back slot, calls `set_front` on it, and only then calls `schedule_segue` — so the sprite's captured slot index already names the new shape. `schedule_segue` enforces that ordering with an always-on `HS_CHECK` that the passed slot is already the front index, and forwards an optional pause gate to the policy's `schedule()` hook — a signature the `Segue::Schedulable` concept pins so no policy can shadow it with a shorter one. The concept sees the arity only; a policy taking the gate and never passing it on is caught by test (`test_segue_policies_forward_pause_gate`), which steps each policy's sprite under a set flag and requires the envelope to hold. `Segue::Crossfade` schedules one fading `Animation::Sprite` per transition and returns a next-transition delay that makes consecutive sprites **overlap**: each transition fades only its own incoming shape in (and back out), while the previous transition's sprite — still alive in its fade-out tail — keeps drawing the outgoing shape; no single call ever draws both meshes, but two are rasterized per overlap frame. The overlap length is configurable via the policy's `overlap` member (frames, clamped to the fade window; negative — the default — selects the full window, and `0` makes the schedule sequential so a single mesh renders per frame). `Segue::Dissolve` overlaps the same way but hands the two draws complementary `DissolveMask`s (same threshold and salt, opposite `invert`), so they partition the wireframe's edges and the overlap frame still costs one mesh's scan; it is the only policy that partitions rasterizer work rather than fragments in the shader, so effects pass its masks to `Plot::Mesh::draw`'s edge-list overload themselves — the only entry point that takes a mask, which is why a solid-mesh pair cannot dissolve. Every other segue is **sequential** (one mesh per frame): `IrisBloom` (faces contract to glowing center points, then the new tessellation blooms out), `Lace` (fill drains to a glowing edge band and floods back), `TerminatorSweep` (a day/night line pinned to the mesh sweeps across it at constant speed; once the line reaches a face, that face fades over a per-face random length in the live `fade_frames_min`/`fade_frames_max` range, fraying the front), `Shockwave` (an expanding wave erases outward from a point; its echo redraws), `Breakdown` (the pattern breaks down one topology class at a time — every face of a color family fades together, classes in a random order reshuffled per swap, each fully gone before the next starts), `SpinFlip` (rigid spin-up, swap hidden in POV motion blur), and `GoldConvergence` (palettes converge to molten gold around the swap). Used by IslamicStars (fixed to `Segue::TerminatorSweep`); HankinSolids keeps a single mesh and drives `OpLeg` directly instead, and DreamBalls uses `Segue::Crossfade` standalone with no carousel, at `overlap = 0` — the sequential case above, so its sprite hand-off never has two meshes on screen at once. |

### Orientation and Motion Blur

`Orientation<CAP>` stores a history of up to `CAP` quaternions (default 4) accumulated during one frame step. The template parameter is the history *capacity*, not the display width — the `World::Orient` and `World::OrientSlice` filters use `Orientation<>`, never `Orientation<288>`. The `World::Orient` filter iterates over this history to distribute motion blur: each point is plotted once per orientation step, with the `age` field increasing backward in time. This means fast-rotating effects naturally show streak-like motion blur with no extra code.

```cpp
timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600, ease_linear, true));
// orientation.length() grows by 1 per sub-step
// World::Orient distributes all steps → motion blur
```

`Orientation::upsample(count)` resamples the orientation history to a higher resolution via SLERP. This is used to rewrite the history when combining multiple animations for accurate parallel sub-frame path tracing — ensuring that concurrent rotations, motions, and walks all contribute to a single coherent set of intermediate orientations.

### OrientationTrail

`OrientationTrail<OrientationType, CAPACITY>` maintains a circular buffer of past `Orientation` snapshots, allowing effects to recall where an object was over previous frames. Each snapshot is a full `Orientation` (with its own sub-frame history).

### VectorTrail and QuantizedVectorTrail

`VectorTrail<CAPACITY>` maintains a circular buffer of past world-space `Vector` positions — 12 B per sample, stored exactly. Used by HopfFibration for its per-fiber trails.

`QuantizedVectorTrail<CAPACITY>` is the unit-sphere variant: each sample is three snorm16 components (6 B, half of `Vector`), clamped to [-1, 1] on record and decoded by value on `get()`, so the round-trip error is at most 1/65534 per component. Used by `ParticleSystem` to record per-particle trajectories for trail rendering.

### `tween` and `deep_tween`

Two traversal helpers linearize multi-level orientation history into a single callback loop:

| Function | Input | Description |
|---|---|---|
| `tween(orientation, callback)` | `Orientation<CAP>` | Iterates over the sub-frame quaternion history of a single orientation, calling `callback(quaternion, t)` for each step with `t ∈ (0, 1]`. Sub-frame 0 is the pose carried over from the previous frame's end and is skipped unless it is the only snapshot (which reads `t = 1`, age-neutral). Used by `World::Orient` to distribute motion blur. |
| `deep_tween(trail, callback)` | `OrientationTrail` (any `Tweenable`) | Flattens a trail of orientations into a single continuous traversal, calling `callback(quaternion, t)` with a global `t` spanning all frames and sub-frames. Used by the orientation-trail effects (Comets, Fishbowl) for rendering trails with full sub-frame accuracy. A bare `Orientation` has no per-frame structure to flatten and is rejected by the `Tweenable` concept — use `tween` for that. |
| `deep_tween_frames(trail, callback)` | `OrientationTrail` (any `Tweenable`) | Public frame-aware traversal that supplies the frame value, sub-frame index, global age, and normalized time; RingSpin uses it to preserve frame boundaries while rendering its trail. |

### Animations and Mutable State

Animations do not render directly — they mutate external state that the rendering pipeline reads. Each animation type targets a specific kind of mutable variable:

| Animation | Target State | What It Mutates |
|---|---|---|
| `Rotation`, `RandomWalk`, `Motion` | `Orientation<CAP>` | Quaternion orientation — pushes sub-frame steps into the orientation history, which `World::Orient` reads for motion blur |
| `Transition` | `float*` | Smoothly interpolates any float parameter (e.g. `speed`, `alpha`, `twist`) from current value to target with easing |
| `Mutation` | `float*` | Applies an arbitrary scalar function `f(t)` to a float over time (more general than `Transition`) |
| `Progress` | `void(float)` callback | Hands the caller eased progress each frame and writes nothing itself; every composed preset transition uses it to blend the authored parameter states |
| `Driver` | `float*` | Continuously increments a float each frame, optionally wrapping at 0..1 — used for phase accumulators |
| `Lerp` | `T*` (type-erased) | Interpolates any type with a `lerp()` function — `MeshState`, params structs, etc. The caller owns start, subject, and target; Lerp holds pointers |
| `ColorWipe` | `GenerativePalette*` | Interpolates palette keys between caller-owned start and target snapshots in OKLCH; both snapshots must remain unchanged and outlive the animation |
| `Ripple`, `MobiusWarp`, `Noise` | `RippleParams`, `MobiusParams`, `NoiseParams` | Animate transformer parameters (expansion radius, warp strength, noise scale) which the transformer pool reads during `MeshOps::transform()` |
| `BallDrop` | `BumpParams` | Walks the bump center down a meridian and re-derives the push axis from the stack's orientation, ramping the footprint envelope; the field pool sums the caps during `field()` |
| `NoiseProduct` | `NoiseProductParams` | Advances the field time axis so the two-octave product noise keeps flowing under live speed edits; the field pool reads it during `field()` |
| `ParticleSystem` | `Vector[]` positions | Physics simulation updates particle positions; `QuantizedVectorTrail` records history for trail rendering |

This separation means effects declare *what state exists* (orientations, floats, palettes) and *what animations drive that state* (rotations, transitions, drivers), but never manually interpolate or update values per-frame. The `Timeline` handles all timing, easing, sequencing, and cleanup:

```cpp
// Effect declares mutable state:
Orientation<> orientation;   // CAP is the sub-frame capacity, not the display width
float twist = 0.0f;
GenerativePalette palette;
GenerativePalette target_palette;
const auto palette_start = palette.snapshot();
const auto palette_target = target_palette.snapshot();

// Timeline drives state via animations:
timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600, ease_linear, true));
timeline.add(0, Animation::Transition(twist, 2.5f, 1000, ease_in_out_cubic));
timeline.add(0, Animation::ColorWipe(palette, palette_start, palette_target, 2000, ease_linear));

// Rendering reads state — no manual updates needed:
void draw_frame() {
    Canvas canvas(*this);
    timeline.step(canvas);  // all state updated automatically
    // orientation, twist, palette are now current-frame values
    filters.plot(canvas, v, palette.get(t), ...);
}
```

## 7.4 Geometry Transformers (`transformer.h`)

Transformers deform the sphere geometry before rendering. The `Transformer<ParamsT, AnimT, TransformFunc, CAPACITY>` class manages a pool of active transform instances, each with its own animated parameters:

```cpp
template <int CAPACITY>
using RippleTransformer = Transformer<Animation::RippleParams, Animation::Ripple,
                                      ripple_transform, CAPACITY>;
```

Available transformers:

| Transformer | Effect |
|---|---|
| `RippleTransformer` | Expands Ricker wavelets from a point, bending the sphere surface radially. Uses fast-reject dot-product heuristic — ~90-95% of vertices skip the `fast_acos` path. |
| `MobiusWarpTransformer` | Applies and releases a Möbius transformation |
| `MobiusWarpCircularTransformer` | Loops a Möbius warp continuously |
| `MobiusWarpGnomonicTransformer` | Möbius via gnomonic projection (preserves straight lines in hemisphere) |
| `NoiseTransformer` | Distorts surface positions with 3D simplex noise |

Transformers integrate with the `MeshOps::transform()` pipeline and can be chained: `MeshOps::transform(input, output, arena, ripple_transformer, orient_transformer)`. `transform()` takes any callable with `operator()(Vector)`, so a pool specialization and a plain adapter compose in the same call.

### Displacement Fields

`FieldTransformer<ParamsT, AnimT, FieldFunc, CAPACITY>` is the scalar counterpart: entities superpose by summation instead of composing as warps, so an effect can feed the summed field into a displacement path (e.g. a `DistortedRing` shift LUT). `field(p)` sums the active entities; `field_bound()` returns a per-frame upper bound on `|field()|` for sizing conservative culls. Where overlapping bodies must dominate instead of stacking, feed the per-entity values through `DominantFieldAccumulator`, which blends them magnitude-weighted (`sum(s³)/sum(s²)`).

| Field | Effect |
|---|---|
| `BallDropTransformer` | Spherical-cap bumps that fall pole-to-pole through a frame, bowing the surface away from each cap |
| `NoiseProductTransformer` | Two-octave product noise where octave 1 envelopes octave 2, so perturbations bunch where the envelope runs strong |

### Pool Lifecycle

Both classes derive from `TransformerPool`, which fixes the call order:

1. `init_storage(Arena&)` — from the effect's `init()`, after any `configure_arenas()` and before the first spawn. It also claims one of the shared `Timeline`'s `MAX_CLEAR_HOOKS` (4) clear-hook slots, so one `Timeline` carries at most four live pools; a fifth traps at registration.
2. `spawn(in_frames, args...)` — the returned pointer is transient; use it at the call site, not across frames.
3. `spawn_pausable(paused, in_frames, args...)` — same as `spawn()`, but the whole timeline event, start delay included, freezes while `*paused` is set, the way `Timeline::add_pausable` does. It is the only pool entry point that honours a GUI pause — `spawn()` animates straight through one — and the flag must outlive the event.
4. `spawn_pinned(in_frames, args...)` — same as `spawn()`, but the pointer may be retained (e.g. registered as a live GUI param). Valid only for an animation that never completes on its own — infinite, or repeating (it rewinds rather than reaching `done()`) — and is added before any finite timeline event, so compaction cannot shift it.
5. `prepare_frame()` — each frame before `transform()` / `field()`, whenever active params changed through animation or live config. The composition reads that prepared state but cannot verify it is current. Its per-entity hooks (`refresh_from(const ParamsT&)` for live config, `sync()` for derived state) are found by detection, so every `ParamsT` must declare one bool per hook — `static constexpr bool NEEDS_REFRESH_FROM` and `static constexpr bool NEEDS_SYNC` — each true only when that hook is carried. Either mismatch is its own `static_assert`, which is what turns a renamed or signature-drifted hook into a compile error instead of a silently unrefreshed entity.
6. `reclaim_storage(Arena&)` — from the after-reset callback of an arena that is compacted mid-effect (e.g. a mesh carousel). Spawned animations hold `Params` references into the slots, so the caller must replay the same allocation order after the reset as after `init_storage()`; the re-claimed blocks must land at their original addresses (asserted). A reset only rewinds the offset, so the untouched bytes carry live entities through.

### Standalone Utilities

`OrientTransformer<CAP>` (`transformer.h`) is a plain adapter struct, not a `Transformer<>` specialization: it holds a reference to an `Orientation<CAP>` and applies `orientation.orient()` to each vertex. It has no pool, no params and no lifecycle — effects construct one on the stack at the call site (a deduction guide takes `CAP` from the orientation) and hand it straight to `MeshOps::transform()`.

## 7.5 Memory Architecture (`memory.h`, `memory.cpp`)

A single contiguous memory block (`GLOBAL_ARENA_SIZE`) is partitioned into three arena allocators. It is 298 KiB on Teensy; the WASM module widens it to 512 KiB (`HS_GLOBAL_ARENA_BYTES` in the root [CMakeLists.txt](https://github.com/woundedlion/pov/blob/master/CMakeLists.txt)), because the chain interpreter's two arenas plus `ShaderChain`'s shared color resources outgrow the device-mirrored default. Individual effects can call `configure_arenas()` to repartition the block at runtime.

| Arena | Default Size | Purpose |
|---|---|---|
| `persistent_arena` | 266 KiB | Long-lived compiled mesh data, persists across frames |
| `scratch_arena_a` | 16 KB | Short-lived intermediate geometry (RAII scoped) |
| `scratch_arena_b` | 16 KB | Secondary scratch for ping-pong subdivision passes |

The native effect harnesses use `HS_GLOBAL_ARENA_BYTES=8388608` — and with it the persistent default, which is whatever the two 16 KB scratch arenas leave — so they can render every effect without OOMing. The device budget stays reachable as `DEVICE_GLOBAL_ARENA_SIZE` / `DEVICE_PERSISTENT_BUDGET`, which the per-effect footprint `static_assert`s check against instead, so an effect that outgrows the device still fails in the host suite.

Effects that need more scratch memory can repartition at init time:

```cpp
// The three sizes must not exceed GLOBAL_ARENA_SIZE (298 KiB on device); an
// over-subscribed partition traps at init() via HS_CHECK rather than silently
// scaling down. Under-subscription is allowed (the surplus is just unused),
// but partitioning the full budget is the norm. Here scratch is doubled at the
// expense of persistent space:
configure_arenas(234 * 1024, 32 * 1024, 32 * 1024);  // 234 + 32 + 32 = 298 KiB
```

A global that caches a pointer into arena storage registers an `ArenaResetHook` beside itself, and drops the pointer from the callback. `configure_arenas()` and `reset_persistent_arena()` (used by the mesh carousel's compaction) run the whole list before handing the storage out again, so the owner never has to be named by the allocator:

```cpp
inline void release_gamut_lut() { g_gamut_lut = GamutLut{}; }
inline const ArenaResetHook GAMUT_LUT_RESET_HOOK(release_gamut_lut);
```

`ScratchScope` provides stack-like RAII lifetime:

```cpp
{
    ScratchScope scratch_a_guard(scratch_arena_a);  // save offset
    // ... allocate from scratch_arena_a ...
}                                        // restore offset — all allocations freed
```

All functions that require scratch memory take explicit `Arena&` parameters — there are no hidden arena references or implicit state, outside the exceptions listed under [Why Arena Allocation?](https://github.com/woundedlion/pov/blob/master/README.md#2-engineering-philosophies):

```cpp
scratch_arena_a.reset();
scratch_arena_b.reset();
ScratchScope _a(scratch_arena_a);
ScratchScope _b(scratch_arena_b);
PolyMesh result = MeshOps::kis(mesh, scratch_arena_a, scratch_arena_b);
```

Conway operators take `(Arena& target, Arena& temp)`, generator functions take `(Arena& a, Arena& b)`, and `classify_faces_by_topology` takes `(Arena& scratch_a, Arena& scratch_b, Arena& persistent)`. This purely functional approach makes the memory layout during heavy geometric operations explicit at every call site.

<a id="compaction-with-persist"></a>

### Compaction with `Persist<T>`

`Persist<T>` is an RAII class that safely evacuates live data from the persistent arena, allowing it to be reset and defragmented, then automatically restores the data on destruction:

```cpp
{
    Persist<MeshState> p(live_mesh, scratch_arena_a, persistent_arena);
    reset_persistent_arena();
    // ... allocate fresh data into persistent_arena ...
}   // ~Persist: clones backup back into persistent_arena
```

### Additional Data Structures

| Type | Description |
|---|---|
| `ArenaVector<T>` | Arena-backed vector whose capacity is fixed between `bind()` calls — appending never grows it (`push_back` traps at capacity). A `bind()` that asks for more re-reserves and abandons the old block until the arena is reset. Copy-disabled, move-enabled. Debug builds detect use-after-free via arena generation tracking. |
| `ArenaSpan<T>` | Non-owning read-only view into an `ArenaVector` (explicit borrow) |

## 7.6 The Color System (`core/color/color.h`)

All internal color data is **16-bit linear light** (`uint16_t r, g, b` in range 0–65535). This avoids the precision loss and incorrect blending that occurs with gamma-encoded 8-bit values.

The conversion pipeline:
```
Input (sRGB 8-bit) → sRGB→linear LUT → Pixel (linear 16-bit) → blend ops
                                                                      ↓
FastLED output ← CRGB(gamma encode) ← linear→sRGB ← Pixel
```

`Color4` wraps `Pixel` with a float alpha channel. The canvas sink composites with a single straight-alpha "over" operation — `blend_alpha(α)`, i.e. `dst = src * α + dst * (1-α)`, applied in 16-bit linear light (see `filter.h`). There is no selectable blend-mode tag.

### Palette Types

| Type | Description |
|---|---|
| `ProceduralPalette` | Cosine palette: `a + b*cos(2π*(c*t + d))` per channel. Defined by 4 vec3 coefficients. |
| `Gradient` | OKLCH interpolation between a sorted list of (position, color) stops. |
| `GenerativePalette` | Procedurally generated palette from harmony rules (triadic, analogous, etc.) combined with brightness/saturation profiles. Supports snapshot/lerp for animated transitions. |

Twenty-seven named `ProceduralPalette` instances are pre-defined in the `Palettes` namespace: `DARK_RAINBOW`, `BLOOD_STREAM`, `VINTAGE_SUNSET`, `RICH_SUNSET`, `UNDERSEA`, `LATE_SUNSET`, `MANGO_PEEL`, `ICE_MELT`, `LEMON_LIME`, `ALGAE`, `EMBERS`, `FIRE_GLOW`, `DARK_PRIMARY`, `MAUVE_FADE`, `LAVENDER_LAKE`, `DESERT_ROSE`, `BRUISED_MOSS`, `BRUISED_BANANA`, `BRIGHT_SUNRISE`, `FIRE_AND_ICE`, `PEACH_POP`, `POPPED_PEACH`, `BLUE_LAGOON`, `ORANGE_CRUSH`, `PLUM_SUNRISE`, `CORAL_BLUE`, and `BRUISED_MANGO`. Six of them — `EMBERS`, `RICH_SUNSET`, `BRIGHT_SUNRISE`, `BRUISED_MOSS`, `LAVENDER_LAKE`, `POPPED_PEACH` — are the slots of `MeshPaletteBank`, the shared baked bank the mesh effects draw from.

### OKLCH Perceptual Color

Palette interpolation is performed in the OKLCH perceptual color space: both `Gradient` (color-stop interpolation) and `GenerativePalette` (harmony-key interpolation and animated transitions) build their tables in OKLCH by default. A `PaletteRecipe` setting `color_path` to `ColorPath::OKLAB_CARTESIAN` moves `GenerativePalette`'s key interpolation onto the rectangular OKLab path instead. The cosine `ProceduralPalette` is the exception — it evaluates its per-channel waveform directly in sRGB. The pipeline:

```
Pixel (linear 16-bit) → linear RGB float → OKLab (L, a, b) → OKLCH (L, C, h)
                                                                  ↓
                                              shortest-arc hue interpolation
                                                                  ↓
                                    OKLCH → OKLab → linear RGB → Pixel
                                                       ↓ (only if out of gamut)
                                              reduce chroma, hold hue + L
```

| Function | Description |
|---|---|
| `linear_rgb_to_oklab()` | Convert linear RGB to the OKLab perceptual space |
| `oklab_to_oklch()` | Convert OKLab (rectangular) to OKLCH (polar: Lightness, Chroma, Hue) |
| `lerp_oklch()` | Interpolate two OKLCH values with shortest-arc hue (avoids the red→green→blue detour) |
| `gamut_clip_preserve_chroma()` | Maps an out-of-gamut OKLab color back into the sRGB cube by reducing chroma while holding hue and lightness (walk-then-bisect on the chroma scale). The hue-preserving alternative to a per-channel RGB clip. Gated behind an in-gamut test (`oklab_to_linear_rgb_gamut`), so in-gamut colors — the vast majority — pay only the test and skip the search. |
| `hue_rotate()` | Perceptual hue rotation — rotates the (a,b) chroma plane in OKLab, preserving lightness and chroma. Forward nonlinearity uses `fast_cbrt` (hot per-pixel path); inverse is exact. Out-of-gamut results are chroma-reduced rather than per-channel clipped, which holds hue and stabilizes the feedback loop against saturated-color drift. The feedback `hue_fade` transform uses `hue_rotate_lms_matrix`; sphere-space hue noise uses `hue_rotate_lut_gamut`. |

### The Gamut Boundary Grid

The chroma clip brackets the sRGB boundary from a generated table (`core/color/gamut_lut.h`, emitted by `tools/gen_gamut_lut.py`) indexed by the diamond angle of (b, a) and by L. Each cell stores the minimum and maximum boundary chroma over the region it covers, so the true boundary of every ray in the cell lies inside the stored bracket at any resolution; the per-pixel path walks that bracket in `GAMUT_SCAN_STEPS` and bisects the straddling step `GAMUT_BRACKET_STEPS` times. Grid resolution sets how wide the bracket starts, and the bisection sets how far it is narrowed; what the floor bounds is the chroma deficit (about 0.003 at 256 × 128, 0.007 at 32 × 16). At every grid, including the shipped master, a handful of rays stride over a disconnected in-gamut interval and land past the first exit — an oversaturation of up to about 0.05 chroma that no amount of bisection recovers, since the wrong step is already selected; a finer grid reduces how many rays do this, not how far they overshoot.

The clip reads the 256 × 128 flash master by default. An effect that clips per pixel can arm an arena copy at the master's resolution or coarser, which buys read latency alone (RAM rather than QSPI flash):

| Function | Description |
|---|---|
| `init_gamut_lut(arena, angle_steps, l_steps)` | Downsamples the flash master into `arena` and points the clip at the copy. Both step counts must divide the master's 256 × 128 and stay at or above `GAMUT_LUT_MIN_ANGLE_STEPS` × `GAMUT_LUT_MIN_L_STEPS` (128 × 64), the coarsest grid the walk resolves — both trapped. Costs `gamut_lut_bytes(angle_steps, l_steps)`. Call from the effect's `init()`, after any `configure_arenas()`. |
| `release_gamut_lut()` | Drops the copy and points the clip back at the flash master. Registered as an `ArenaResetHook`, so `configure_arenas()` and the mesh carousel's compaction both run it before handing the storage out again. |

`MeshFeedback::init()` is the only production call site that arms an arena copy. It takes the full 256 × 128 grid, `gamut_lut_bytes(256, 128)` = 131,074 B of persistent arena. The Shader workbench, `ShaderChain`, and composed effects clip against the flash master without allocating a gamut copy.

### Palette Modifiers

Modifiers compose around any palette source at compile time via
`StaticPalette<Source, Coords<...>, Colors<...>, Wrap, Shade>`. There are two axes: a
**coordinate** chain that remaps the lookup parameter `t` *before* the source is
sampled, and a **color** chain that reshapes the resulting sample *after*, with
the original coordinate in hand. Both chains are inlined by fold expression with
zero runtime overhead. `Wrap` (default `true`) wraps the final coordinate into
`[0,1)` before the lookup — leave it on for cycling modifiers that overflow the
range; set it `false` for bounded remaps that must reach the source endpoints.
Both directions are `static_assert`ed: an unbounded modifier rejects
`Wrap=false`, and a bounded final modifier rejects `Wrap=true` (wrapping would
fold its 1.0 output to 0.0 and destroy the top endpoint). Only a modifier that
re-bounds *arbitrary* input (`WrapModifier`'s fold, `FoldModifier`'s triangle
wave, `InsetModifier`'s clamp) clears an unbounded predecessor; `ReverseModifier`,
`MirrorModifier`, `QuantizeModifier`, and `PinchModifier` are bounded on `[0,1]` but pass an out-of-range coordinate
straight through — chaining one after a cycling modifier needs a `WrapModifier`
between them and `Wrap=false`.

`Shade` (a `ShadeCoord`, default `MATCH_WRAP`) selects which coordinate the
color chain receives:

| `ShadeCoord` | Coordinate handed to `shade()` |
|---|---|
| `MATCH_WRAP` | The lookup coordinate when `Wrap` is on, the raw input when it is off |
| `LOOKUP` | The coordinate the source was sampled at, whatever `Wrap` is |
| `RAW_INPUT` | The raw pre-modifier input, whatever `Wrap` is |

It is a separate knob rather than a second meaning of `Wrap` because the
coordinate chain can force `Wrap` through the `static_assert`s above: a shade
that must read the un-wrapped input still composes behind a cycling modifier.

Coordinate modifiers (`modify(float) -> float`):

| Modifier | Effect |
|---|---|
| `CycleModifier` | Shifts the lookup parameter by a continuously incrementing offset (palette scrolling) |
| `BreatheModifier` | Oscillates the lookup parameter with a sinusoidal "breathing" envelope |
| `RippleModifier` | Applies a wavelet distortion to the lookup parameter |
| `FoldModifier` | Folds the parameter space (mirror at edges) to create ping-pong patterns |
| `PinchModifier` | Non-linearly warps the lookup parameter toward a focal point; the per-sample `powf` suits bake-time sampling |
| `QuantizeModifier` | Posterizes the palette into discrete bands |
| `ScaleModifier` | Scales and offsets the lookup parameter |
| `ReverseModifier` | Mirrors the lookup parameter (1.0 - t) |
| `MirrorModifier` | Maps [0,1] to [0,1,0] for a seamless symmetric loop |
| `InsetModifier` | Compresses the source domain into an inset window, clamping outside |
| `WrapModifier` | Folds the lookup parameter into `[0,1)` mid-chain, so a bounded modifier can follow a cycling one |
| `NoiseWarpModifier` | Displaces the lookup parameter with smooth value noise — the aperiodic counterpart to `RippleModifier` |
| `DriftModifier` | Meanders the whole palette along a per-frame noise walk (wanders, hesitates, reverses) |

Color modifiers (`shade(Color4, float) -> Color4`):

| Modifier | Effect |
|---|---|
| `AlphaFalloffShade` | Scales alpha by a caller-supplied falloff curve over the coordinate |
| `EdgeFadeShade` | Fades the sample color to black near the edges (opaque vignette) |
| `EdgeAlphaShade` | Fades the sample alpha near the edges (transparent vignette) |
| `HueSpinShade` | Rotates every sample's hue in OKLab by a driver amount (continuous hue cycling); the rotation folds into a per-frame memoized 3×3 |
| `HueWobbleShade` | Rotates hue by an amount that varies along the domain (iridescent drift); per-sample cost suits bake-time sampling |
| `SparkleShade` | Ignites sparse traveling glints where an evolving noise field exceeds a threshold |
| `ChromaPulseShade` | Breathes OKLab chroma between pastel and vivid on a per-frame memoized pulse |
| `LightnessGrainShade` | Grains brightness with evolving noise; uniform linear-RGB gain, so hue is exact below the saturation point (a gain above 1 clips bright channels) |
| `IridescentShade` | Adds a thin-film cosine sheen with per-channel phase offsets, saturating at white |

The noise-driven modifiers sample the deterministic `value_noise_1d`/`value_noise_2d`
hash lattice (`3dmath.h`) with a per-instance seed, so two modifiers on the same
driver decorrelate by seed. Frame-constant work memoizes against the driver
value (`HueSpinShade`'s rotation matrix, `ChromaPulseShade`'s pulse factor,
`DriftModifier`'s walk offset). The OKLab shades still pay a per-sample
conversion, so they pair well with `BakedPaletteStorage::rebake`, which re-samples a
256-entry LUT once per frame; the noise and cosine shades are cheap enough for
live per-pixel paths.

```cpp
// Compose a baked palette with a breathing coordinate modifier
StaticPalette<BakedPalette, Coords<BreatheModifier>> palette;

// A transparent vignette: inset the source, fade alpha at the edges
StaticPalette<ProceduralPalette, Coords<InsetModifier>,
              Colors<EdgeAlphaShade>, /*Wrap=*/false> vignette;

// Psychedelic composite: noise-warped coordinate, continuously spinning hue,
// glints riding on top
StaticPalette<ProceduralPalette, Coords<NoiseWarpModifier>,
              Colors<HueSpinShade, SparkleShade>> lava;
```

### Additional Palette Types

| Type | Description |
|---|---|
| `MutatingPalette` | Extends `ProceduralPalette` with continuous coefficient mutation between two procedural palettes |
| `SolidColorPalette` | Returns a single fixed color for every coordinate |
| `PaletteFacade<SP>` | Exposes a compile-time `StaticPalette` composition through the polymorphic `Palette` API, for preset tables and baking |
| `BakedPalette` | Read-only view of an arena-backed 256-entry color/alpha LUT. |
| `BakedPaletteStorage` | Owns mutation rights to a palette LUT and rebakes a `Palette` or `StaticPalette` source into it. |
| `NoiseHuePalette<Source>` | Applies a sphere-domain noise field as a spatial OKLab hue rotation over any palette source. Its shared hue-rotation and cube-map noise LUT preparation is used by ordinary effects, composed shader effects, and the Shader workbench. Call `hue_shift(direction, amount)` once when a whole primitive shares a noise coordinate, `noise_uv(cos_u, sin_u, cos_v, sin_v)` for a seamless two-axis surface field, or `get(t, direction, amount)` directly per sample. |

### Recipe-Compiled Palettes

`GenerativePalette` is not configured field by field — it is *compiled* from a
`PaletteRecipe`, a flat POD of authoring controls that canonicalizes into the
control keys the palette evaluates. The recipe is the persisted form, so it
carries a `schema_version` pinned to `PaletteRecipe::SCHEMA_VERSION`; a stored
recipe from another schema is rejected rather than silently misread.

| Field | Description |
|---|---|
| `input` | `PaletteInputWindow{offset, span}` — the sub-window of the source domain the recipe maps across |
| `domain` | `PaletteDomain`: `STRAIGHT`, `MIRROR` (ping-pong), `VIGNETTE` (fades in and out at both ends), `FALLOFF` (holds, then fades out), `LOOP` (seamless, whole-turn hue winding) |
| `easing` | `SegmentEase` between control keys: `LINEAR`, `COSINE`, `SMOOTHSTEP` |
| `color_path` | `ColorPath::OKLCH_ARC` (polar, shortest-arc hue) or `OKLAB_CARTESIAN` (rectangular, straight through the neutral axis) |
| `hue` | `HueControls`: `mode` (`HARMONY`, `SWEEP`, `CUSTOM`), `harmony` (`PaletteHarmony`), `direction` (`HueDirection`), `base_turns`, `spread_turns`, `sweep_turns`, `custom_turns[PALETTE_MAX_KEYS]` |
| `lightness` | `AxisControls`: `curve` (`AxisCurve`: `CONSTANT`, `ASCENDING`, `DESCENDING`, `BELL`, `CUP`, `CUSTOM`), `center`, `range`, `custom[]` |
| `chroma` | `ChromaControls`: the same curve/center/range/custom, plus `basis` (`ChromaBasis::LOCAL_GAMUT` or `ABSOLUTE`; `PATH_MINIMUM` holds its ordinal but is unimplemented and fails compilation) and `headroom` |
| `hue_torsion` | Shifts each key's hue by `hue_torsion * (L - 0.5)`, so the light and dark ends drift apart |
| `falloff_start` | Where the `FALLOFF` domain's fade reaches zero; must lie in `(2/3, 1)` under that domain, and is canonicalized back to its default under any other |

Compilation both validates and normalizes.
`GenerativePalette::try_compile(input, output, canonical, status)` returns false
on rejection and leaves `output` untouched; the constructor takes the same path
and fail-fast traps instead. `PaletteCompileStatus` carries the verdict:

| Member | Description |
|---|---|
| `code` | `PaletteCompileCode`: `OK`, `INVALID_SCHEMA`, `NON_FINITE`, `INVALID_ENUM`, `HUE_LIMIT`, `NON_INTEGER_LOOP_SWEEP`, `INVALID_FALLOFF_START` (`INCOMPATIBLE_OPTIONS` holds its ordinal but is never produced) |
| `field` | The `PaletteRecipeField` naming the offending control, so an authoring tool can point at it |
| `adjustments` | `PaletteAdjustments`: three `PaletteRecipeField` bitmasks — `wrapped_fields`, `clamped_fields`, `canonicalized_fields` — recording every silent normalization the compile applied |

A successful compile also hands back the `canonical` recipe: the normalized form
the palette was actually built from, which is what an authoring tool should
persist rather than the raw input.

The `PaletteRecipes` namespace collects the stock builders — `hue_turns()`,
`harmony()`, `balanced_analogous()`, `profile()`, `random_profile()`,
`random_base_turns()`,
`from_oklch_keys()`, `from_colors()`, `isolight_spectral_loop()` and
`tonal_monochrome()` — and `core/color/effect_palette_recipes.h` holds the
per-effect recipes the roster renders.

### Palette Cycling

`PaletteCycler` (`core/color/palette_cycler.h`) drives a display LUT through a
sequence of palettes over time: it dwells on an entry for `dwell_frames`, then
fades into the next over `fade_frames` (a zero dwell chains fades back to back).
Effects call `step()` once per frame and shade from `palette()`, a
`BakedPalette` — outside a fade the display is a bit-exact bake of the current
entry.

The fade mechanism is chosen per adjacent pair at `init()`. Two morph-compatible
`GenerativePalette`s fade by **key-space morph**, interpolating control keys for
perceptually coherent hue travel; every other pair — composed, prebaked, or
morph-incompatible — falls back to a **baked-LUT crossfade**. An `Entry`
accordingly accepts a `GenerativePalette`, any `Palette`, or a `BakedPalette`,
so a mixed sequence of up to `MAX_ENTRIES` cycles correctly. Entries are
caller-owned and must outlive the cycler; `Entry`'s rvalue overloads are deleted
so a temporary cannot bind.

`init_generated()` replaces the fixed entry array with a `NextPaletteFn`
provider that fills palette *n* on demand, for an endless non-repeating cycle;
successive palettes must stay morph-compatible and every retarget is fail-fast
checked. Arena cost is declared up front — `display_arena_bytes()`,
`crossfade_arena_bytes()`, `morph_arena_bytes()`, and the
`required_arena_bytes()` / `generated_arena_bytes()` worst cases — so an effect
sizes its arena without probing. `advance_without_display()` keeps the timeline
and the provider moving without rebuilding the LUT; a later `step()` rebuilds it
at the current phase.

## 7.7 The Mesh System (`core/mesh/`)

The mesh system uses these headers in `core/mesh/` and `core/render/sdf/`:

- **`base_mesh.h`** — Base mesh identities, bounds, and authoring labels
- **`relax_bake.h`** — Relax payload and source identity checks
- **`core/mesh/mesh.h`** — Core data structures (`PolyMesh`, `HalfEdgeMesh`) and fundamental `MeshOps` (compile, clone, classify)
- **`conway.h`** — Conway mesh operators and vertex transformations
- **`conway_graph.h`** — Constexpr 23-edge morph graph over the 18 simple-registry solids: per-edge operator/seed/reseed specs, bridge-aware walk weighting, and the closed `ORDERED_TOUR`
- **`recipe_types.h`** — The authored op-chain model: the `Op` operator set, one `OpStep`, and the `Recipe` chain a registry generator mirrors, split out so the model is not read out of the registry tables written in it
- **`recipe.h`** — Lowers an authored recipe to primitive steps (`expand_to_primitives`), sizes that lowering at compile time (`lowered_step_count`, `max_lowered_step_count`), replays either form through `SolidBuilder` (`build_recipe`, `build_steps`), and decides which lowered steps a morph leg can sweep (`is_morphable_step`)
- **`hankin.h`** — Hankin pattern compilation and dynamic update
- **`core/render/sdf/face_class_bake.h`** — Congruence-class clustering plus one canonical distance-LUT bake per class, allocated by descending face count under an 18 KB per-mesh budget
- **`core/render/sdf/face_classes.h`** — The class id space and the three record types the rasterizer binds per frame, split out so the clustering and bake machinery stays out of every rasterizer translation unit
- **`mesh_state.h`** — `MeshState`, the flat-array renderer format, split out so mesh, Conway, Hankin and solids code can share the renderer-facing representation without the construction machinery
- **`solid_generators.h`** — Hardcoded Platonic vertex/face tables, the `SolidBuilder` operator chain, and the named Archimedean / Catalan / Islamic Star Pattern generators
- **`solids.h`** — The three solid registries, the authored `Recipe` mirrors of the generators, and the name/index lookups over them
- **`relax_bake_specs.h`** — Authored bake names and iteration budgets. The extraction harness reads these inputs without loading generated payloads.
- **`relax_bakes_generated.h`** — Baked relaxed-mesh vertices behind `MeshOps::relax_baked`, generated by `tools/relax_bakes.py`; never hand-edited, regenerate with `<build>/relax_bake_gen | python tools/relax_bakes.py emit --stdin`. Grid constants come from the harness dump.

`PolyMesh` stores vertices and face connectivity via `ArenaVector` arrays. `MeshState` (in `mesh_state.h`) is the flat compiled format consumed by the renderer. `HalfEdgeMesh` provides a half-edge traversal structure built from either a `PolyMesh` or `MeshState`.

`MeshOps::require_closed_manifold()` temporarily allocates a `uint16_t` fan count for every index from zero through the largest referenced vertex. Budget `2 * (largest_index + 1)` bytes plus alignment in its scratch arena, up to 65,536 bytes at `MeshLimits::MAX_VERTICES`. This scratch is additional to the resident `HalfEdgeMesh` and is rewound before return.

### Core MeshOps (`core/mesh/mesh.h`)

| Operation | Description |
|---|---|
| `MeshOps::compile` | Convert a `PolyMesh` to the flat-array `MeshState` format used by the renderer |
| `MeshOps::clone` | Arena-safe deep copy |
| `MeshOps::classify_faces_by_topology` | Group faces by vertex count, sorted whole-degree-rounded interior angles, and neighbor topology for palette assignment. The angle vector is the discriminator that separates faces of equal side count, so a solid whose interior angles straddle a rounding boundary classifies differently either side of it |

### Conway Operators (`conway.h`)

All Conway *geometry* operators (`dual` through `bevel` below) take `(const PolyMesh& mesh, Arena& target, Arena& temp)` and return a `PolyMesh`; `truncate`, `expand`, `chamfer`, `snub` and `bevel` take a defaulted shape parameter `float t` after the arenas, and `snub` a `float twist` after that. `medial` takes the same two arenas but writes its two outputs through reference parameters — `(const PolyMesh& mesh, PolyMesh& out_a, ArenaVector<Vector>& out_b, Arena& target, Arena& temp)`. `ambo`, `truncate`, `expand`, `chamfer` and `snub` each carry a second overload taking a caller-built `const HalfEdgeMesh&` right after `mesh`: the output topology is the same at every `t`, so one connectivity build serves a whole parameter sweep instead of one per call. The one exception is `truncate` at exactly `t == 0.5`, which short-circuits to `ambo` and its different face census — a sweep that must hold one topology stops short of it (`ConwayGraph::T_EPS_AMBO`). `transform`, `transform_in_place`, `relax`, `relax_baked`, and `normalize` are listed in the same table but are mesh utilities with their own signatures. Every operator, primitive or composed, produces its `PolyMesh` into `target` and uses `temp` for intermediate computation; composed operators (`gyro`, `meta`, `needle`, `zip`, `bevel`) reuse the same internal ping-pong as their constituent ops, with an even-length composition starting its first step in `temp` so the last one lands in `target` (see the COMPOSITION POLARITY note in `conway.h`):

| Operation | Description |
|---|---|
| `MeshOps::transform` | Apply a chain of vertex transformers to produce a new `MeshState` |
| `MeshOps::transform_in_place` | Apply a chain of vertex transformers to a `MeshState`'s own vertices, leaving topology untouched |
| `MeshOps::dual` | Dual mesh (faces ↔ vertices) |
| `MeshOps::kis` | Raise a pyramid on each face |
| `MeshOps::ambo` | Truncate vertices to edge midpoints |
| `MeshOps::truncate` | Cut corners off the polyhedron (configurable depth) |
| `MeshOps::expand` | Separate faces (ambo of ambo) |
| `MeshOps::chamfer` | Bevel edges (hexagonal expansion) |
| `MeshOps::snub` | Chiral semi-regular polyhedron with twist (Newell-method face normals) |
| `MeshOps::gyro` | Gyro operator (= dual ∘ snub) |
| `MeshOps::meta` | Meta operator = kis ∘ dual ∘ ambo |
| `MeshOps::needle` | Needle operator = kis ∘ dual |
| `MeshOps::zip` | Zip operator = dual ∘ kis |
| `MeshOps::bevel` | Bevel operator = truncate ∘ ambo |
| `MeshOps::medial` | Both endpoint vertex sets of the dual morph on one shared medial (rectified) connectivity: `out_a` is `ambo(mesh)`, `out_b` the matching `ambo(dual(mesh))` positions |
| `MeshOps::relax` | Edge-length relaxation by spring forces on the unit sphere. |
| `MeshOps::relax_baked` | Substitute a flash-baked relax result for the pass. Its runtime checks catch a source/bake mismatch (dimensions, topology hash) and payload corruption (output hash, re-derived from the bake's own vertex bits) — they say nothing about freshness, since a positional retune that leaves connectivity intact passes every one of them. Freshness is the `unit_relax_bake_verify` ctest's job: it re-runs the live `relax` and asserts bit-exact equality with the committed payload |
| `MeshOps::normalize` | Project all vertices onto the unit sphere |

### Hankin Pattern System (`hankin.h`)

| Operation | Description |
|---|---|
| `MeshOps::compile_hankin` | Pre-compute topological data for fast Hankin pattern updates |
| `MeshOps::update_hankin` | Update dynamic vertices based on angle parameter (no new allocation when reusing a sufficiently-sized output mesh) |
| `MeshOps::hankin` | One-shot Hankin pattern generation (compile + update) |

`compile_hankin` produces a `CompiledHankin` struct containing base vertices, static midpoints, and dynamic instructions. `update_hankin` evaluates the dynamic vertices by sweeping the Hankin angle, producing the star polygon line intersections for each face. It re-binds the output mesh's vectors on every call, so it avoids new allocation only in the steady state — reusing the same output mesh against the same arena, already sized large enough.

### Solids Library (`solids.h`, `solid_generators.h`)

`solid_generators.h` provides constexpr vertex/face data for all Platonic solids plus procedural generators for Archimedean, Catalan, and Islamic Star Pattern families; `solids.h` organizes them into three registries. Firmware builds by name via `Solids::get_by_name(arena, a, b, name)`, which fails fast on an unknown name; the WASM bridge validates with `Solids::find_entry(name)` first and generates from the entry; the WASM geometry tools enumerate the registries by index with `Solids::get_entry(index)` to populate the picker, then build the selected solid by name:

| Registry | Count | Description |
|---|---|---|
| `simple_registry` | `PLATONIC_COUNT + ARCHIMEDEAN_COUNT` | 5 Platonic (tetrahedron through icosahedron) + 13 Archimedean solids |
| `catalan_registry` | `CATALAN_COUNT` | Duals of the Archimedean solids (triakisTetrahedron, rhombicDodecahedron, pentakisDodecahedron, etc.) |
| `islamic_registry` | `ISLAMIC_COUNT` | Complex multi-operator recipes producing Islamic star patterns from base solids |

Total: `Solids::NUM_ENTRIES`. Each count is a named constant in `solids.h`, `static_assert`ed against its registry's size (and their sum against `NUM_ENTRIES`), so the values live there rather than being restated here.

`Collections` namespace provides typed spans for iterating subsets: `get_platonic_solids()`, `get_archimedean_solids()`, `get_simple_solids()`, `get_catalan_solids()`, `get_islamic_solids()`.

`SolidBuilder` provides a fluent interface for chaining Conway operators with automatic arena swapping:

```cpp
return SolidBuilder(to_polymesh<Icosahedron>(a), a, b)
    .truncate()
    .dual()
    .build();
```

Islamic Star Pattern recipes chain multiple operators with Hankin pattern generation:

```cpp
SolidBuilder(dodecahedron(a, b), a, b)
    .hankin(54.0f * D2R).ambo().hankin(72.0f * D2R).build();
```

## 7.8 Generators (`memory.h`)

`memory.h` provides a single universal generation wrapper that manages arena lifecycle for all procedural geometry creation:

```cpp
namespace hs {
template <typename GenerateFn, typename... Args>
auto generate(Arena &target, GenerateFn &&fn, Args &&...args);
}
```

It resets both scratch arenas only at the outermost call (depth zero), scopes them on every call, then invokes `fn(target, scratch_a, scratch_b, args...)`. Nested calls preserve the caller's live scratch allocations. Direct registry lookups and effect geometry creation go through this wrapper for a deterministic arena lifecycle:

```cpp
auto mesh = hs::generate(persistent_arena, Solids::get_by_name, std::string_view("icosahedron"));
```

One deliberate exception: `SolidBuilder`'s fluent Conway chain (`solid_generators.h`) owns its own two-arena ping-pong, swapping the scratch arenas between operators, so it manages arena lifecycle directly rather than through `generate()`.

## 7.9 The Preset System (`control/choreography.h`)

`ChoreographedEffect<Derived, Params>` provides the runtime preset lifecycle: an `Effect` base owning the effect's live parameter set, its preset table, and the choreography that moves between presets. `Derived` declares its presets — a `PRESETS` table (`std::array<PresetEntry<Params>, N>`) and/or `PRESET_IDS` naming them — plus a `Segue` preset policy, a dwell, a parameter schema version, and a validity predicate:

```cpp
static constexpr Segue::Preset::Snap PRESET_SEGUE{};
static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
static constexpr uint16_t PRESET_DWELL_FRAMES = 241;
static constexpr std::array<PresetEntry<Params>, 12> PRESETS = {{ ... }};
static bool valid_params(const Params &p);
```

The effect calls `begin_choreography()` once from `init()` — it configures the engine preset controller from the table (or `PRESET_IDS`) and, under a `Segue::Preset::Fade` policy, arms the through-black envelope loop — and `step_choreography()` every frame, which retires the dwell and starts the next automatic transition. A single-preset effect compiles the dwell countdown out. An effect that clocks its own advancement (DreamBalls' sprite hand-off chain) skips the per-frame tick and calls `advance_preset()` from its own schedule instead.

Automatic transitions follow the policy — `Segue::Preset::Lerp` crossfades the live parameters into the target, `Segue::Preset::Snap` adopts immediately, `Segue::Preset::Fade` snaps inside the envelope's dark frame — while manual and synchronized selections always snap. Hooks specialize the mechanics: `preset_params(index)` (static, or a member when the effect patches entries at runtime, e.g. re-binding a noise pointer) overrides the `PRESETS[index]` lookup, and in its static form supplies the startup default too; `initial_params()` overrides whichever of `preset_params(0)` and `PRESETS[0]` would otherwise start the effect; shadowing `adopt_params(target)` re-derives dependent state after a snap; `transition_armed(target)` fires once as a crossfade arms, capturing the endpoint state the blend interpolates alongside the parameters; `blend_params(progress)` writes an in-flight Lerp; `set_preset_opacity(value)` receives a Fade envelope. The base also carries schema-versioned parameter snapshots: `serialize_parameters()` tags the live set with `PARAMETER_SCHEMA_VERSION`, and `restore_parameters()` rejects a snapshot taken under a different schema or failing `valid_params()`.

`control/preset_host.h` holds the controller the choreography drives: the committed index, the vetoable `apply_preset()` hook, and the manual `selectPreset`/`nextPreset`/`previousPreset` surface the WASM bridge calls. `control/presets.h` holds the table vocabulary: `PresetEntry<Params>` (the row type) and the free `constexpr` helper `all_presets_in_ranges(entries, in_ranges)`, which folds a slider-range predicate over an entry table so an effect can `static_assert` its whole preset table against its registered parameter ranges — a loop rather than an unrolled conjunction, so appended entries are covered automatically.

Promoted composed effects also have authored shader documents under `patterns/`. Their descriptor and preset-bank digests are checked against their effect headers by `scripts/promoted_digests.test.mjs`; the runtime uses `ChoreographedEffect` to play the resulting presets.

## 7.10 Hardware Drivers (`dma_led.h`, `pov_single.h`, `pov_segmented.h`)

Three hardware drivers form a layered stack. The DMA LED layer handles the SPI wire protocol across four headers: `hd107s_frame.h` (wire format and inline color correction), `dma_led_core.h` (framing, transfer length, stale-transfer predicate), `dma_led_controller.h` (double-buffer orchestration, templated on its transport) and `dma_led.h` (the Teensy SPI/DMA peripheral driver, the only Arduino-only piece). `pov_single.h` and `pov_segmented.h` sit above it and manage the POV column sweep, differing only in how many Teensys share the work; the segmented driver's ISR decisions are split into two further host-tested headers, `pov_handoff.h` and `pov_submit_gate.h`.

### DMA LED Controller (`dma_led.h`, `hd107s_frame.h`, `dma_led_core.h`, `dma_led_controller.h`)

Non-blocking DMA-based LED output for HD107S (APA102-compatible) LEDs on Teensy 4.x.  Enabled by `#define USE_DMA_LEDS` in the target's boilerplate header (`targets/Phantasm/phantasm_target.h`) before it includes the driver; `led.h` stays neutral and the default FastLED/WS2801 path remains as fallback. The FastLED fallback applies only to the single-board `POVDisplay`; the segmented `POVSegmented` driver `#error`s without `USE_DMA_LEDS` (FastLED's bit-bang `show()` masks IRQs for windows that break the sync symbol margins, which are derived from a mask window M ≈ 0), so DMA LEDs are mandatory on Phantasm.

| Class | Header | Role |
|---|---|---|
| `HD107SFrame<N>` | `hd107s_frame.h` | Pre-formatted DMA buffer for the HD107S protocol. `pack_pixel()` writes `Pixel` values directly into the frame buffer with inline color correction (color correction → temperature → brightness, then `linear_to_srgb8` to 8-bit sRGB), bypassing the CRGB intermediate. The buffer is 32-byte-aligned (`__attribute__((aligned(32)))`) and cleaned with `arm_dcache_flush()` (clean, no invalidate — the buffer is TX-only) for cache coherency. |
| `TeensySPIDMA` | `dma_led.h` | Low-level DMA+SPI driver wired to LPSPI4. Configures a `DMAChannel` with completion interrupt for fully async byte-stream transmission. |
| `DMALEDController<N>` | `dma_led_controller.h` | Double-buffered high-level controller. The ISR packs pixels into `back_frame()`, then `submit_frame()` flushes it and triggers async DMA, returning immediately. If the previous transfer is still in flight, `submit_frame()` **drops** the new frame (bumping `get_overrun_count()`) and returns false rather than spinning; a transfer that never completes is surfaced as a wedged-channel fault. The drop returns before the buffers swap, so `back_frame()` still holds the dropped pixels and a caller can re-submit them without repacking. |
| `next_buffer()`, `transfer_len()`, `transfer_us()`, `transfer_stale()` | `dma_led_core.h` | Free `constexpr` framing and watchdog math the controller's decisions derive from, host-tested without the peripherals. |

The canvas and color-correction pipeline use 16-bit linear values. `pack_pixel()`
converts them with `linear_to_srgb8` to the 8-bit sRGB channels sent on the SPI wire:

```cpp
// ISR path (per column): fetch the display buffer once, index it directly
const Pixel* buf = effect->display_buffer();               // 16-bit linear pixels
// Physical LED index comes from the single-source-of-truth map (pov_single_map.h),
// which applies the top-arm reversal / bottom-arm offset — never the raw row index.
frame.pack_pixel(pov::strip_top_led(y, S), buf[y * width + x]); // Pixel → HD107S frame
(void)ledController.submit_frame();                          // non-blocking DMA; returns false on overrun
```

### Single-Teensy POV Driver (`pov_single.h`)

`POVDisplay<S, RPM>` drives the Holosphere — one Teensy owns the entire LED strip.  An `IntervalTimer` ISR fires at `1,000,000 / (RPM/60) / width` µs intervals to advance one column:

```
Main Loop                              ISR (IntervalTimer)
──────────                             ───────────────────
effect->draw_frame()                   show_col() fires every N µs
  Canvas canvas(*this)                   for y in 0..S/2:
    render to bufs[cur]                    pack_pixel(strip_top_led(y,S),    get_pixel(x, y))                      // top arm
  ~Canvas → queue_frame()                  pack_pixel(strip_bottom_led(y,S), get_pixel(strip_opposite_col(x,W), y)) // bottom arm
                                         submit_frame() → async DMA
                                         x = (x+1) % width
                                         if x==0 || x==width/2: advance_display()
```

The top arm's physical LED ordering is reversed (LED 0 at the junction end, descending in Y), and the bottom arm shows the opposite half of the image (x offset by W/2).

`show_col()` discards `submit_frame()`'s overrun verdict: this driver carries no retry latch and no dark fallback, so a dropped column leaves the previous column lit for one extra period — the drop returns before the buffer flip. `run()` fail-fast-checks that one composite transfer fits inside a column period, which rules out the systematic overrun that would hold the strip on one frame and is what makes discarding the verdict sound.

The shipping `holosphere` sketch runs `RingSpin<96, 20>` with column strobing disabled. Its 1302 µs column period clears the 1160 µs FastLED transfer bound; the roster’s default RingSpin configuration keeps strobing enabled for the DMA targets.

| Parameter | Value (Holosphere) |
|---|---|
| S (total pixels) | 40 |
| RPM | 480 |
| Column interval | ~1302 µs (= 125 ms / 96 columns) |
| ISR duration | ~20 µs on the DMA path (`holosphere_dma`), which packs the column and starts an asynchronous transfer. The FastLED path clocks 40×24 bits at 6 MHz (~160 µs per show). A nonstrobed column normally takes ~160 µs after sufficient idle time; `CMinWait<1000>` can add up to 1000 µs. A strobed column adds an immediate blank show and its 1000 µs wait: ~1320 µs even when the first show does not wait, exceeding the ~1302 µs period. The current transfer guard conservatively budgets 1160 µs per show (2320 µs for strobing) and refuses a column period at or below that bound. |


### Multi-Teensy Segmented POV Driver (`pov_segmented.h`)

`POVSegmented<S, N, RPM>` drives Phantasm — N Teensys (4 by default, up to 8) each control a contiguous Y-segment on one arm. IDs `[0, N/2)` map to arm A and `[N/2, N)` to arm B. Within an arm, northern bands advance in +Y and southern bands run from the S pole toward the junction in -Y.

**Physical strip layout (N=4, S=288):**

```
Arm A                               Arm B (x offset by W/2)
┌──────────────────────────┐        ┌──────────────────────────┐
│ Seg 0 (top)              │        │ Seg 2 (top)              │
│ LED 0 at N pole (y=0)    │        │ LED 0 at N pole (y=0)    │
│ → LED 71 at junction     │        │ → LED 71 at junction     │
├────────── junction ──────┤        ├────────── junction ──────┤
│ Seg 1 (bottom, reversed) │        │ Seg 3 (bottom, reversed) │
│ LED 0 at S pole (y=143)  │        │ LED 0 at S pole (y=143)  │
│ → LED 71 at junction     │        │ → LED 71 at junction     │
└──────────────────────────┘        └──────────────────────────┘
```

At N=8, each arm has four 36-pixel bands: north outer `0–35`, north inner
`36–71`, south outer `143–108`, and south inner `107–72`. Arm A uses IDs 0–3;
arm B uses IDs 4–7. The firmware profile is compile-tested, but an eight-board
480-RPM rotor requires a separately qualified mounting, balance, cable, and
swept-envelope design before operation.

**Hardware ID detection**: Each Teensy reads `log2(N)` active-low GPIO straps: pin 21 (ID0), pin 22 (ID1), and pin 23 (ID2, N=8 only). The ID is `(~raw) & (N-1)`, so grounding a strap sets its bit and all-floating selects ID 0 (sync master). The header supports power-of-two `N ≤ 8`.

**Precomputed ISR indexing**: All per-segment mapping decisions are resolved at boot time into three precomputed values:

| Value | Description |
|---|---|
| `y_base` | Starting Y index for this segment's row band |
| `y_step` | +1 for northern bands, -1 for reversed southern bands |
| `arm_b` | Whether this segment is on arm B (x offset by W/2) |

The ISR accumulates the precomputed pixel offset and handles a unity envelope without per-pixel scaling:

```cpp
const Pixel* buf = effect->display_buffer();    // fast path: no per-pixel virtual dispatch
const int stride = segment_row_stride(segment, width);
int off = segment_pixel_base(segment, x_col, width);
if (effect->output_envelope_u16() == 65535u)
    for (int i = 0; i < PPS; ++i, off += stride)
        frame.pack_pixel(i, buf[off]);
else
    for (int i = 0; i < PPS; ++i, off += stride)
        frame.pack_pixel(i, effect->apply_output_envelope(buf[off]));
```

**ISR state machines**: the wake's non-trivial decisions are split out of the Arduino-only driver into two host-tested headers, which `run_wake_sequence()` drives in ISR order:

| State machine | Header | Role |
|---|---|---|
| `EffectHandoff<T>` | `pov_handoff.h` | Foreground↔ISR effect ownership: the teardown counter handshake, the acquire/release publish and adopt of a pending effect, the consumed-generation gate that keeps the ISR off a deleted instance, and the display-window (clip) alternation. The foreground constructs and deletes instances; the ISR only ever dereferences what `live()` handed it. |
| `SubmitGate`, `SyncPulseGate` | `pov_submit_gate.h` | The LED transport's accept/drop verdict and the sync pin's pulse width. Both submit paths — the fail-dark black frame and the image column — clear their pending state only on an accepted submit, and a dropped column latches a retry that the next wake re-submits without repacking, so a drop costs one flywheel wake (~54 µs) rather than a dark column. A wake that renders nothing has too short a body to carry a scheduled sync pulse, so the pin is held HIGH across the ISR boundary and dropped at the head of the next wake. |

**Effect transparency**: Effects are written against the full 288×144 canvas with no per-segment code. Each board clips rendering to its half-width segment band for the current display window (`clip_to_segment`), except stateful effects (`needs_full_frame()` / `persists_pixels()`), which render the full canvas; the ISR then packs this board's LEDs. Every board reseeds the shared `Pcg32` at every effect build from `HS_PHANTASM_EFFECT_SEEDS[]`, which `targets/Phantasm/phantasm_playlist.h` builds as `hs::stable_effect_seed(hs::stable_effect_id<name<CANVAS_W, CANVAS_H>>(#name))` (`core/platform/rng.h`) so an entry's stream follows its persisted effect ID (or class name when no ID is declared) rather than its roster position; `hs::epoch_seed(effect index)` (epoch 0 is the identity seed `1337`) is the fallback for a board that supplies no seed table. Either way a board's canvas depends only on the beacon-synchronized index — a mid-show joiner renders bit-identically to boards that have been up for hours.

| Parameter | Value (qualified N=4 default unless noted) |
|---|---|
| S (total pixels) | 288 |
| N (segments) | 4 qualified default; 8 compile-tested firmware profile |
| PPS (pixels per segment) | 72 at N=4; 36 at N=8 |
| RPM | 480 (8 rev/s, ~125 ms/rev) |
| Frame rate | 16 FPS (2 frames/rev — one per side; each side draws W/2 = 144 cols per 62.5 ms frame) |
| Column frequency | 2304 Hz |
| Column interval | ~434 µs (= 125 ms / 288 = 62.5 ms / 144) |
| Flywheel wake period | ~54.25 µs (= column interval / `OVERSAMPLE` = 8) |
| ISR duration (N=4 72px pack + DMA trigger) | Read per build from the `HS_ISR_PROFILE` accumulators (`g_flywheel_wake_cycles`, `g_column_pack_cycles`, `g_dma_submit_cycles`) that the `Profile` target dumps. A column-boundary ISR outrunning the wake period coalesces the wakes it overruns; the column index is derived from the cycle counter, so the next wake resumes at the time-correct column |

### Frame Sync Protocol: 1-Wire Signal Datasheet

Phantasm's boards stay coherent over **one wire**. Segment 0 (the **master**) is the conductor; segments 1 through N−1 (**downstream**) listen. Each board generates its own columns from a local **flywheel timebase** and snaps that timebase to count-coded pulse bursts the master broadcasts on the wire. Full design: `docs/specs/phantasm_frame_sync_spec.md`; host-tested protocol core: `hardware/pov_sync.h`.

The flywheel derives the column index from the free-running CPU cycle counter, never from counting timer interrupts:

```
x = ( x_boundary + (now − epoch) · (W/2) / cycles_per_half_rev )  mod W
                                                    └─ 64-bit intermediate
```

`epoch` is folded forward by exactly one half-revolution at every boundary crossing, so the 32-bit cycle counter's ~7.16 s wrap is structurally unobservable.  An interrupt-masked window (e.g. `FastLED.show()`) cannot drop columns — the ISR that runs after the mask reads the clock and resumes at the *time-correct* column.

**Pin / signal description.**

| Signal | Master (seg 0) | Downstream (segments 1–N−1) | Electrical | Direction |
|---|---|---|---|---|
| `SYNC` | pin 3 (GPIO out) | pin 3 (ext. interrupt in) | 3.3 V CMOS, active-high pulses, idle LOW | master → all |
| `MASTER_EN` | pin 5 (drive LOW) | pin 5 (drive HIGH) | 3.3 V CMOS, gates the external sync-out buffer | per-board strap out |
| `ID` straps | pins 21–23 (ID0–ID2 as required) | pins 21–23 (ID0–ID2 as required) | active-low, internal pull-ups; `ID_STRAPS = log2(N)` read | strap (board identity) |

The ID straps select the board: the build reads `ID_STRAPS = log2(N)` active-low bits and decodes `(~raw) & (N-1)`. N=4 reads ID0/pin 21 and ID1/pin 22; N=8 also reads ID2/pin 23. All-floating selects segment 0/master. `SYNC` is one shared pin 3 — the master drives it and downstream boards receive on its rising edge; `MASTER_EN` (pin 5) gates an external level shifter so only the master drives the shared bus. `SYNC` is the only inter-board connection; pin 4 is unused. It is assumed physically reliable (a hard, soldered line); a severed wire is out of scope (boards free-run and precess apart at crystal rate, a slow smear, never an instant break).

**Signal levels & symbol waveforms.** The wire idles LOW.  A **symbol** is a burst of short active-high pulses at a fixed pitch; **the meaning is the count of rising edges — pulse width carries no information.**  The pin is driven HIGH by the flywheel tick that schedules the pulse, not at ISR entry, and the rising edge is the only timed event.  A wake that also renders drops the pin before it returns, for a ~8–13 µs pulse — the path every scheduled pulse takes, since pulses fall on column boundaries.  A wake that renders nothing holds the pin across the ISR boundary and drops it at the head of the next wake, ~54 µs later.  Pulses are drawn narrow, to scale against the ~868 µs pitch:

```
 HALF — 1 pulse — marks boundary x = W/2 (144)
            ┌┐
 ───────────┘└──────────────────────────────────────────────────  idle LOW
            ▲
            └ boundary instant (x = W/2)

 ZERO — 3 pulses — marks boundary x = 0
            ┌┐      ┌┐      ┌┐
 ───────────┘└──────┘└──────┘└───────────────────────────────────
            ▲   └ 2-col pitch ┘
            └ boundary instant (x = 0)

 ZERO+EPOCH — 5 pulses — marks x = 0 AND advances the playlist
            ┌┐      ┌┐      ┌┐      ┌┐      ┌┐
 ───────────┘└──────┘└──────┘└──────┘└──────┘└───────────────────
            ▲
            └ boundary instant (x = 0)
```

A burst terminates when the wire stays quiet past the **gap timeout** (4 columns).  The consumer counts rising edges and classifies:

| Symbol | Edges | Marks | Carries | Rate |
|---|---|---|---|---|
| `HALF` | **1** | boundary `x = W/2` | half-rev phase + flip | 1 / rev |
| `ZERO` | **3** | boundary `x = 0` | half-rev phase + flip | 1 / rev |
| `ZERO+EPOCH` | **5** | boundary `x = 0` | phase + flip + **playlist advance** | 1 / effect (×R repeats) |
| `BEACON` | 5 base-8 digits @ `x ≈ W/4` | — (data channel) | absolute effect index + rev count, checksummed | rev ≡ 1 (mod 16) + first revs of an effect |
| *invalid* | any **even** count, or > 5 | — | discarded whole: no snap, no flip, no advance | — |

**Why count, not width:** on the i.MX RT each pin has a single latched interrupt flag, so an IRQ-mask window *delays* an edge's ISR but cannot lose the edge unless two edges fall inside one mask window.  With pulse pitch chosen **greater than the worst-case mask window M**, the edge *count* is exact even when `FastLED.show()` masks IRQs mid-symbol (on Phantasm's DMA LED path, M ≈ 0).  The alphabet is **odd-only, distance 2** — a single lost or spurious edge lands on an even (invalid) count and is discarded.  A glitch degrades to a *missed* symbol (covered by the local boundary crossing), **never** a *misclassified* one: *fail to "missed," never to "wrong."*

**AC timing characteristics.** At 480 RPM / 600 MHz / W = 288 (1 column = 434.03 µs = 260,417 cycles):

| Parameter | Symbol | Columns | Time | Cycles | Rule |
|---|---|---|---|---|---|
| Column period | T0 | 1 | 434.0 µs | 260,417 | `cycles_per_half_rev / (W/2)` |
| Boundary pulse pitch | t_PB | 2 | 868.1 µs | 520,833 | **pitch > M** ⇒ no edge lost to the latch |
| Beacon digit pitch | t_PD | 1 | 434.0 µs | 260,417 | checksum tolerates tighter pitch |
| Burst gap timeout | t_GAP | 4 | 1.736 ms | 1,041,667 | **> pitch + M** ⇒ a mask can't split one burst |
| Glitch filter (min edge spacing) | t_GF | — | 100 µs | 60,000 | edges closer than this are EMI — rejected |
| Master late-censor budget | t_LATE | ½ | 217 µs | 130,208 | first pulse later than this ⇒ skip whole symbol |
| ACQUIRE quiet-before guard | t_QB | 16 | 6.94 ms | 4,166,667 | a hard snap requires this much prior silence |
| Beacon interdigit timeout | t_BID | 24 | 10.4 ms | 6,250,000 | stale partial beacon frame dropped after this |
| Half-revolution | — | 144 | 62.5 ms | 37,500,000 | one image / one flip interval |
| Revolution | — | 288 | 125 ms | 75,000,000 | two flips, two boundary symbols |

All bursts are ≪ the 62.5 ms half-rev, so consecutive symbols never overlap. The table is nominal: cycle counts are the exact rational rounded, while `Config::cycles_per_column()` truncates its division (~2.6 ppm low). Only pitches and thresholds derive from it — flywheel position divides by the exact `cycles_per_half_rev`.

**One-revolution signal map.** Where each symbol lands across a single 125 ms revolution (beacon only on scheduled revolutions):

```
 column x →    0           72(W/4)        144(W/2)        216           288 ≡ 0
               │             │              │              │              │
 SYNC wire   ██ZERO     ░░░BEACON░░░      ██HALF                       ██ZERO
             (3 edges)  (5 digits, data)  (1 edge)                     (3 edges)
               │←─────────── half-rev = 62.5 ms ──────────→│
 display     flip A                       flip B                       flip A
 layer 1     snap φ                       snap φ                       snap φ
               │←──────────────────── revolution = 125 ms ────────────────────→│
```

Boundary symbols (`ZERO`/`HALF`) serve **two** layers at once: they snap the flywheel's column phase (Layer 1) *and* act as the exactly-once flip backstop (Layer 2).  The beacon rides the otherwise-quiet stretch at `x ≈ W/4`, separating the timing channel and the data channel **in time** on the same wire.

**The three disciplined layers.** Every layer reads the same flywheel timebase, so one snap corrects all three coherently; each also has an absolute reference on the wire that pulls it back if it drifts:

* **Layer 1 — Column phase.** Boundary symbols snap each flywheel twice per revolution; inter-snap crystal drift is **~0.006 column** at 40 ppm. The larger downstream phase term is the flywheel wake grid: up to **~54.25 microseconds (0.125 column)** per receiver before it processes a boundary. The master has no receive-wake delay. Both terms fit below one column; crystal drift alone does not bound the seam error.  In **LOCKED** a symbol is accepted only if its implied correction is **≤ G = 4 columns** and its boundary identity matches the flywheel's prediction (the plausibility gate).
* **Layer 2 — Buffer flip.** The local boundary crossing flips the display buffer; the symbol is a deduplicated backstop.  `try_flip`, keyed on boundary identity (boundaries strictly alternate `ZERO, HALF, …`), makes the flip **exactly-once** even when both the crossing and the symbol fire.  Losing both paths in one half-rev is the only glitch, and it self-heals the next half-rev.
* **Layer 3 — Content.** The playlist is **epoch-counted**, not `millis()`-gated.  Duration is **per roster entry**, not uniform: `HS_PHANTASM_EFFECT_LIST` carries a seconds column beside each name and `targets/Phantasm/Phantasm.ino` converts it to `EFFECT_REVOLUTIONS[]` at `seconds · RPM / 60`, spanning 38 s (304 revolutions) to 240 s (1,920 revolutions) across the playlist.  The master emits the `EPOCH` mark (plus R = 3 redundancy repeats) when the current entry's revolutions elapse; every board counts down to the same **absolute** commit boundary regardless of which copy it heard, constructs the next roster entry during the final K = 2-revolution **construction window** (display black on all boards simultaneously), and all swap to its frame 0 at the same boundary.  The beacon broadcasts the absolute effect index so a board that missed every epoch repeat corrects within ~2 s, and a rebooted board rejoins at the correct effect — **fail-dark, never fail-wrong** (a board with no established identity shows black rather than a guessed effect).  Every one of those revolution budgets is absolute, so on a 304-revolution entry the 25-revolution rejoin bound still costs 8% of the effect's airtime.

**Index beacon frame format.** The beacon is a **data** symbol (integrity by *rejection*, not by exactness).  Five base-8 digits at 1-column pitch, each digit a burst of `digit + 1` pulses, digits separated by 5 quiet columns (one past the gap timeout, so the decoder reliably terminates each digit):

```
 Frame = [ idx_hi  idx_lo  rev_hi  rev_lo  checksum ]   (5 digits, base-8)
           └── effect index 0–63 ──┘ └ rev mod 64 ┘  └ Σ(i+1)·dᵢ mod 8

 digit Dk transmitted as (Dk+1) pulses @ 1-col pitch, then a 5-col quiet gap:

         D0          D1               D2          …        D4
        ┌┐          ┌┐┌┐┌┐           ┌┐┌┐                 ┌┐┌┐┌┐
 ───────┘└──/ /─────┘└┘└┘└──/ /──────┘└┘└──/ /───────────┘└┘└┘└──────────
        │←Dk+1 pulses→│   │←5-col quiet (terminates digit)→│
        │←──────────── frame ≈ 26 ms worst case (≪ half-rev) ───────────→│
```

Any checksum mismatch, wrong digit count, out-of-range digit, or stale partial frame **drops the whole frame** — the next beacon is ≤ 2 s away.  Schedule: revolution 1 of every 16 (`rev ≡ 1 mod 16` — never rev 0, so a just-powered board meets clean isolated boundary symbols first), plus the first revs of a fresh effect; silent during a pending commit.

**Receiver state machine.** Each downstream board is in one of two states.  The master is born `LOCKED` with identity (effect 0, rev 0) — it *is* the reference and never snaps:

```
                  first accepted snap
   ┌──────────────┐ ──────────────────▶ ┌──────────────┐
   │   ACQUIRE    │                      │    LOCKED    │
   │  (display    │                      │ (disciplined,│
   │   black)     │ ◀────────────────── │  rendering)  │
   └──────────────┘  R = 4 consecutive   └──────────────┘
                     gate rejections
                     (~2 revolutions)

 ACQUIRE : accept any *valid* symbol unconditionally (hard snap), but only
           on a burst preceded by ≥ t_QB (16 col) of wire silence — so a
           beacon digit train can't capture a just-rebooted board mid-frame.
           The train's FIRST digit is preceded by silence exactly as a
           boundary symbol is, so it can still be mistaken once; the
           R-rejection fallback bounds the recovery (spec §9.1 mis-snap row).
           Renders black until it has BOTH phase (a snap) AND identity
           (epoch/beacon).
 LOCKED  : accept a valid symbol only if implied correction ≤ G (4 col) AND
           boundary identity matches the prediction. Else reject (telemetry,
           no snap, no flip). After R rejections the board concludes its OWN
           timebase is at fault and falls back to ACQUIRE (the escape hatch
           that stops a genuinely-lost board from rejecting good symbols
           forever).
```

**Epoch commit sequence.** `EPOCH` at ZERO boundary **B** schedules an absolute commit at **B + R + K** (R = 3 repeats, K = 2 construction revolutions).  A board hearing any repeat infers its position in the train from its own revolution count and lands on the *same* boundary:

```
 ZERO boundary:   B        B+1      B+2      B+3      B+4      B+5
 master emits:   ●EPOCH    ○rpt     ○rpt     ○rpt      —        —
                 (5 edges) (5)      (5)      (5)
 commit_in_revs:   5        4        3        2        1        0
                 │←─── announce: OLD effect still renders ────│
                                            │←── construct ──→│ swap → NEW
                 ░░░░░ display BLACK: envelope zero from B ░░░░░  frame 0
 all boards:     ░░ outgoing renders, but output is dark ░░░░░░── new effect
```

The construction window is identical (K revolutions) on every board because construction can't begin before B+R — only then is the window's start common knowledge regardless of which copy each board heard.  An effect that can't construct inside K revolutions trips `HS_CHECK` (fail-fast).  All boards reseed `hs::random()` per effect build from `HS_PHANTASM_EFFECT_SEEDS[]` in `targets/Phantasm/phantasm_playlist.h` — `hs::stable_effect_seed(hs::stable_effect_id<name<CANVAS_W, CANVAS_H>>(#name))` in `core/platform/rng.h`, with `hs::epoch_seed(effect index)` the fallback when no seed table is supplied.  That seed is identical on every visit, so an entry replays the same stream each time it comes round, and the new instance is bit-identical across boards no matter what each board rendered — or whether it even existed — before the epoch.

**Output envelope.** Nothing on the strip cuts at the commit: the ISR scales every packed column by `effect_output_envelope` (`pov_sync_content.h`), a fade-through-clear driven from the synchronized `rev_in_effect` and the column index alone.  An entry fades up over its first two revolutions and back down over its last two, and once `rev_in_effect` reaches the entry's configured duration the envelope is **zero**.  That is exactly the revolution at which the master starts the EPOCH train, so the LEDs are already dark at **B** and stay dark through the announce phase as well as the construction window — R + K = 5 revolutions, ~0.63 s at 480 RPM — before the incoming entry fades up from its own frame 0.  The window itself, not the counter, is what forces the zero: the beacon carries only six bits of revolution, so a board that joined an entry longer than 64 revolutions counts congruent to the master rather than equal to it and would otherwise read full brightness at **B**.  It steps to black with everyone else and misses only the ramp into it.  Otherwise a pure function of already-synchronized state, it needs nothing extra on the wire and every board computes the same value for the same column; the pack loop tests for full brightness first, so a mid-effect revolution pays no multiply.

**Live-takeover join grid.** The epoch commit only aligns a *swap* between running boards; a board with nothing live yet — at boot, or after a reboot — reaches its first constructed effect whenever its own identity arrives, which is later downstream than on the master. It therefore does not take that effect live on arrival: it waits for a ZERO crossing whose revolution-in-effect is a multiple of `join_grid_revs` (4), marked by `TickActions::join_boundary`. The master sits on the same grid, so at boot every board goes live at the *same* crossing with aligned frame counters instead of the master leading by however long downstream identity took; a mid-show rejoin waits ≤ 4 revolutions, well inside the enforced 25-revolution rejoin budget. That budget is a `Config::valid()` relation, not prose: `rejoin_bound_revs()` is the widest beacon-to-beacon gap plus the join-grid wait — beacon period 16 + EPOCH repeats 3 + construction window 2 (beacons are suppressed for the whole commit window) + join grid 4 = 25 revolutions, ~3.1 s at 480 RPM — and `valid()` rejects any configuration whose bound exceeds `rejoin_budget_revs`. The 16-revolution figure is the beacon cadence alone. The grid must divide 64 so a beacon's mod-64 revolution count lands on the master's grid. Where the epoch commit traps on a missing effect, the join is conditional: `EffectHandoff::joinable` also requires the pending build generation to be unconsumed and to match the one the wire advertises, so a visibility lag just joins one grid step later. No join boundary is marked while a commit is pending — the epoch path owns that swap.

**Concurrency & failure modes.** Two ISRs per board, **single-writer** by construction.  The sync-wire RISING ISR is a pure *publisher* — glitch filter, edge count, first-edge timestamp into a small mailbox, nothing else.  The flywheel ISR (waking ~8× per column) is the sole *consumer/owner* of all sync state: it claims terminated bursts, classifies, gates, snaps, flips, and runs epoch scheduling.  The hot path is ~7-of-8 wakes doing one cycle-counter read and a 64-bit position compute (≈1 % CPU at 600 MHz); only a column change packs pixels and submits DMA.

| Event | Layer 1 (column) | Layer 2 (flip) | Layer 3 (content) |
|---|---|---|---|
| Masked-IRQ window (`FastLED.show()`) | resumes at time-correct column | unaffected | unaffected |
| 1 dropped boundary symbol | coasts ≤ 1 rev (~0.01 col); re-snaps next | local crossing still flips | unaffected |
| 1 spurious / EMI edge | even count discarded, or gate rejects | identity dedup no-ops it | epoch refractory + gate guard it |
| Late-emitted symbol | master self-censors; residual gate-rejected | crossing flips on time regardless | unaffected |
| 1 board renders slow (drops a frame) | — | shows prior frame for 1 period | stateless: heals next frame/beacon; stateful: heals next epoch |
| 1 dropped epoch symbol | — | — | R repeats; missed-all-R corrected by next beacon (~2 s) |
| Board reboots mid-show | re-acquires phase from next valid symbol | resumes flipping once LOCKED | rejoins correct effect via beacon, ≤ 25 revs (~3.1 s); dark until then |
| Sync wire severed (out of scope) | free-runs on own crystal; precesses ≥ 1 col in ~10–20 s | keeps flipping locally | holds last effect; slow drift, never an instant break |

The flywheel ISR maintains telemetry counters (symbols accepted / gate-rejected / discarded, beacons ok / rejected, index corrections, epochs refractory-ignored, lock transitions, flips, emissions censored / aborted, longest coast) that the foreground reports behind `hs::debug` at ≤ 1 Hz — so any degradation the protocol absorbs silently is still visible at a glance.


## 7.11 Mathematical Kernels (`core/math/`)

The math headers provide coordinate transforms and scalar fields used by the
plot, scan, and pullback paths. They do not own a canvas or effect lifecycle.

| Headers | Surface |
|---|---|
| `core/math/3dmath.h`, `core/math/4dmath.h` | Vectors, matrices, complex arithmetic, fast scalar approximations, and four-dimensional rotations. |
| `core/math/geometry.h` | Umbrella for periodic, pixel-mapping, and spherical helpers. |
| `core/math/periodic.h` | Periodic coordinates and index wrapping. |
| `core/math/pixel_mapping.h` | Sphere/pixel coordinates and display latitude conventions, using `H_OFFSET` defined in `core/platform/platform.h`. |
| `core/math/spherical.h` | Canonical axes, tangent frames, spherical distances, and parallel transport. |
| `core/math/projections.h` | Bonne, Peirce quincuncial, Airocean, folded sinusoidal, and equirectangular sphere-to-plane kernels. |
| `core/math/stereographic.h`, `core/math/mobius.h`, `core/math/lenses.h` | Stereographic projections, fractional-linear transforms, and sphere-domain lens kernels. |
| `core/math/rotate.h`, `core/math/projection_patterns.h` | Angle wrapping, canvas-to-sphere projection, and shared projected-pattern coordinates. |
| `core/math/noise_field.h`, `core/math/spherical_field.h`, `core/math/spherical_harmonics.h` | Noise sampling, spherical fields, and harmonic evaluation. |
| `core/math/easing.h`, `core/math/waves.h`, `core/math/interpolate.h` | Scalar easing curves, periodic waves, and interpolation helpers. |

Projection results carry their coordinate and validity contracts in the headers.
The pullback policies in `core/render/pullback/` bind these kernels to effect
parameters and frame state.

## 7.12 Spatial Queries (`core/spatial/`)

`kd_tree.h` provides arena-backed three-dimensional nearest-neighbor queries.
`reaction_graph.h` defines the Fibonacci reaction lattice and node lookup;
`reaction_graph.cpp` stores its generated neighbor table. The generator of record
is `scripts/generate_reaction_graph.py`. Reaction-diffusion effects use these
neighbors to advance their fields without rebuilding adjacency per frame.
