# Core Subsystems

The engine subsystem reference for the [Holosphere README](https://github.com/woundedlion/pov/blob/master/README.md).

## Contents

- [The Shader Interface](#the-shader-interface)
  - [The Fragment](#the-fragment)
  - [Shader Signatures](#shader-signatures)
  - [Register Conventions by Rasterizer](#register-conventions-by-rasterizer)
  - [Composed pullback metadata](#composed-pullback-metadata)
  - [Column-time sampling](#space-time-rendering)
  - [Spherical ray rendering](#spherical-ray-rendering)
- [SDF Shapes and the Scan Rasterizer](#sdf-shapes-sdfh-and-the-scan-rasterizer-scanh)
  - [SDF Shape Primitives](#sdf-shape-primitives-sdfh)
  - [Volumetric Shapes](#volumetric-shapes-sdfvolumeh)
  - [CSG Operations](#csg-operations-sdfh)
  - [Scan Rasterization Primitives](#scan-rasterization-primitives-scanh)
  - [Near-Pole Azimuthal LOD](#near-pole-azimuthal-lod)
  - [Continuous Radial Surfaces](#continuous-radial-surfaces)
- [The Curve Rasterizer](#the-curve-rasterizer-ploth)
  - [Sampling Policy](#sampling-policy)
  - [Plot Primitives](#plot-primitives)
- [The Animation System](#the-animation-system-animationh)
  - [Animation Types](#animation-types)
  - [Orientation history](#orientation-history)
  - [OrientationTrail](#orientationtrail)
  - [VectorTrail and QuantizedVectorTrail](#vectortrail-and-quantizedvectortrail)
  - [`tween` and `deep_tween`](#tween-and-deep_tween)
  - [Animations and Mutable State](#animations-and-mutable-state)
- [Geometry Transformers](#geometry-transformers-transformerh)
  - [Displacement Fields](#displacement-fields)
  - [Pool Lifecycle](#pool-lifecycle)
  - [Standalone Utilities](#standalone-utilities)
- [Memory Architecture](#memory-architecture-memoryh-memorycpp)
  - [Compaction with `Persist<T>`](#compaction-with-persist)
  - [Additional Data Structures](#additional-data-structures)
- [The Color System](#the-color-system-corecolorcolorh)
  - [Palette Types](#palette-types)
  - [OKLCH Perceptual Color](#oklch-perceptual-color)
  - [The Gamut Boundary Grid](#the-gamut-boundary-grid)
  - [Palette Modifiers](#palette-modifiers)
  - [Additional Palette Types](#additional-palette-types)
  - [Recipe-Compiled Palettes](#recipe-compiled-palettes)
  - [Palette Cycling](#palette-cycling)
- [The Mesh System](#the-mesh-system-coremesh)
  - [Core MeshOps](#core-meshops-coremeshmeshh)
  - [Conway Operators](#conway-operators-conwayh)
  - [Hankin Pattern System](#hankin-pattern-system-coremeshhankinh)
  - [Solids Library](#solids-library-solidsh-solid_generatorsh)
- [Generators](#generators-memoryh)
- [The Preset System](#the-preset-system-controlchoreographyh)
- [Hardware Drivers](#hardware-drivers-dma_ledh-pov_singleh-pov_segmentedh)
  - [DMA LED Controller](#dma-led-controller-dma_ledh-hd107s_frameh-dma_led_coreh-dma_led_controllerh)
  - [Single-Teensy POV Driver](#single-teensy-pov-driver-pov_singleh)
  - [Multi-Teensy Segmented POV Driver](#multi-teensy-segmented-pov-driver-pov_segmentedh)
  - [Frame Sync Protocol: 1-Wire Signal Datasheet](#frame-sync-protocol-1-wire-signal-datasheet)
- [Mathematical Kernels](#mathematical-kernels-coremath)
- [Spatial Queries](#spatial-queries-corespatial)

---

## The Shader Interface

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

The registers are *inputs* — populated by the rasterizer before the shader runs. The shader reads them and writes `color`. Plot rasterizers forward the color through the filter pipeline to the canvas. `Scan::Shader` instead overwrites the canvas directly with premultiplied color, bypassing filters and destination blending.

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

**SDF Scanline Path** (`Scan::Ring`, `Scan::Star`, `Scan::PlanarPolygon`, `Scan::SphericalPolygon`, `Scan::Flower`, `Scan::Line`, `Scan::Mesh`):

| Register | Source | Meaning |
|---|---|---|
| `v0` | `DistanceResult.t` | Normalized azimuth (0–1) for `Scan::Ring` and `Scan::Star`; normalized radial position for `Scan::PlanarPolygon`, `Scan::SphericalPolygon` and `Scan::Flower` (polar angle over the shape radius, so it passes 1 outside the body); unused (0) for `Scan::Line` and `Scan::Mesh` faces |
| `v1` | `DistanceResult.raw_dist` | Producer-specific: unsigned centerline/arc distance for rings and lines, polar angle for polygons and stars, or antipode scan distance for flowers; `Scan::Mesh` faces carry the signed edge distance instead — negative inside the face, in gnomonic plane units on a small (`linear_dist`) face and in radians on a large one — which `fragment_edge_dist()` turns into normalized inward depth as `-v1 / size` |
| `v2` | Set by rasterizer | Stroke AA coverage (0–1, also applied by Scan at plot time), 0 for solid shapes, or face index for `Scan::Mesh` (but see the per-face setup note below) |
| `v3` | `DistanceResult.aux` | Auxiliary pass-through register; every current producer writes 0 |
| `size` | `DistanceResult.size` | Stroke half-width for stroke shapes, or radius or apothem for filled shapes (mesh `Face` floors it at 0.25× the face circumradius, so on a sliver face — whose true inradius approaches zero — the reported size overstates it without bound) |

With a per-face setup callback, `Scan::Mesh::draw_specialized` runs its minimal-fragment loop: only `v1` is refreshed per pixel. `age` is an in/out register carried forward unless refreshed. Other fragment inputs are unavailable and are poisoned with NaN in debug builds; shaders must not read `v2` or `size` through `mesh_face_index()` or `fragment_edge_dist()`. The face index and face size reach the shader through the setup callback instead. `Scan::Mesh::draw` does not accept this callback.

Each SDF shape's `distance<ComputeUVs>()` method writes a `DistanceResult`. Rings, lines and spherical polygons report angular distances; planar polygons, stars and flowers use chart-plane distances, and small `SDF::Face` shapes (reference size `max(inradius, 0.25 × circumradius)` < 0.2) use gnomonic tangent-plane units. `size` follows the producer's distance units. The per-producer register table in `core/render/sdf/common.h` defines `t` and `raw_dist`:

```cpp
struct DistanceResult {
  float dist;        // Signed distance (negative = inside)
  float t;           // Per-shape register; see the producer table
  float raw_dist;    // Per-shape register (signed for mesh faces)
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

Plot primitives interpolate registers between control-point Fragments via `Fragment::lerp_registers()`. The vertex shader, if provided, runs once per control point before rasterization. For `Plot::Polygon<Plot::PlanarProjection>`, `Plot::Star<Plot::PlanarProjection>`, and `Plot::Flower`, the rasterizer re-derives `v0`/`v1` from the rendered arc: azimuthal-equidistant edges, with a geodesic fallback near the basis antipode. Both registers stay consistent with the drawn position.

**Full-Screen Shader Path** (`Scan::Shader`):

The finished premultiplied color replaces each canvas pixel directly. No filter pipeline or destination blending is applied; alpha below one darkens the output.

Registers are not pre-populated — the shader receives only `pos` (reconstructed from pixel coordinates). The single-callback overload accepts `shader(const Vector &)` or `shader(const Vector &, int x, int y)` and returns `Color4` or an already premultiplied `Pixel`. The two-callback overload separates per-pixel vertex setup from per-subsample fragment evaluation.

**Volumetric Path** (`Scan::Volume`):

The fragment shader receives `pos` set to the closest local-space hit point (in the SDF's coordinate frame) and `size` set to the closest signed distance. No register convention — the shader computes lighting from the local-space position directly.

### Composed pullback metadata

`Pullback::Spec` supplies a typed `Pipeline<B>` plus the `HARMONY`, `ANIMATED_PROJECTION`, `TEMPORAL_SAMPLING`, `FACTORED_COLUMNS`, `MODULATION` and `MODULATION_EVALUATION` controls. `Pullback::MetadataFor<Spec>` inspects the normalized semantic leaves, including nested placement groups and optional `void` entries. It derives projection, sampling and field coverage, transfer kind, lens policy, surface placement, path tracking and color modes without binding an executable frame. `ComposedEffect::Metadata` exposes that topology to slider registration and optional storage. Authored scalar topology fields are rejected.

### Space-time rendering

A POV display lights each column at its own moment. With two opposing arms,
columns x and x + W/2 are lit together at surface time

```
t(x) = (x mod (W/2)) / (W/2)
```

For the 288-column canvas there are 144 distinct times in a frame. A frame is
the interval between two effect ticks. The arm cuts lie at 143|144 and 287|0:
physically neighbouring columns there belong to consecutive frame buffers.
Drawing each buffer from one instant makes moving geometry jump by a frame's
motion at a cut. Sampling geometry and shading across the frame lets those
columns meet at their actual emission times instead.

`Spacetime::SamplingSurface<W>` describes the relation between columns, arms
and time indices; `Spacetime::surface_time` provides the runtime-width map.
The still selection is `Spacetime::flat_surface()`: drawing uses the current
endpoint, t = 1. The simulator retains the `setTemporalSampling` binding,
whose Boolean selects the complementary timed mode. The sampling selection
does not advance animation: effects still step their mutable state once per
frame and retain the endpoints needed to evaluate intermediate times.

[Continuous radial surfaces](#continuous-radial-surfaces) evaluate their
star family at each column's sampling time through the same column walk.

### Motion, worldlines and poses

A moving point traces a worldline; an edge joins two endpoint worldlines into
a worldsheet. `Spacetime::Worldline` evaluates a chord, a quadratic through a
middle capture, or a circular arc about an axis. `Spacetime::Worldlines`
retains per-vertex start, middle and end captures in arena storage, carries a
frame's end into the next frame's start, and selects an interpolant using
`Spacetime::PathBound`. Effects and draw helpers supply the captures; a new
sprite or topology change can hold the start at the end. The interpolant's
`deviation()` supplies a pad for swept bounds.

`Spacetime::PoseTable` samples an orientation's sub-frame history at each
surface time and at the frame end. Both arms share the same entries. World
stages compose these poses with geometry motion; `Filter::World::Orient`
prepares its table from the orientation it references. `sweep_arc()` bounds
the world stage's frame motion for culling.

#### Points and plotted edges

`Plot::PointRoute` prepares a clip gate widened by the point motion bound.
Preparation retains the render clip's exact content key; plotting checks it
against the key cached by `Effect` when clip bounds or margins change.
Its `plot` takes a world-space end, a lazy start callable, and a shader. The
start and shading are evaluated only for points that can reach the clip.
A projected chord gives a screen crossing through
`Spacetime::screen_crossings`; equating the point's linear longitude with
an arm's lit longitude gives one division per arm. A projected path that
bends too far is divided into pieces sampled from the same world chord.
`Spacetime::crossings` and `piecewise_crossings` provide the corresponding
longitude-track operations.

A crossing is splatted through the sink's anti-alias stage. When adjacent tap
columns require different positions, including at an arm cut, the route
reads the point's position at each tap column's own time. A shader without a
time argument is shaded once for that point; a shader taking t is sampled
at the emitted crossing's time.

The column walk draws an edge at each column's known surface time. It
combines the endpoint worldlines with the world's pose there, intersects
the resulting geodesic with that column's meridian, and emits its samples
through the filter chain. `Plot::EdgeWorldsheet` caches geodesic bases at
lattice times; `Plot::SweptEdgeCull` rejects edges whose swept bounds miss the
clip before their geodesic bases and sampling pitch are built. Longitude
half-planes are tested before the height bound. `Plot::EdgeBatch` shares
captured vertex data across edges. History-free walks whose taps need no
position extrapolation also bound samples by the physical northern AA cap.
Two-point hulls specialize at their flash callers. Signed row bounds skip
normalization only when the padded pairwise dot products stay positive.
A directional spacing bound certifies edges whose midpoint lattice stays
within the sample footprint; other edges retain the detailed spacing check
and use column-local sampling where needed.
`Plot::Mesh::draw_timed` routes table-compatible world stages through
`draw_columns` and `rasterize_column_edges`, with an optional
`Plot::GeometryCallback` for vertex motion. Pipelines with source history use
`rasterize_recording_column_edges`; other pipelines compile out recording.
The still selection uses the ordinary mesh raster. Other timed plotted
configurations retain their
existing raster dispatch.

`Plot::Ring::draw_timed` and `Plot::DistortedRing::draw_timed` share
`draw_ring_columns`. Their circular vertex paths take two endpoint captures
in `Spacetime::Worldlines`, selected through `Spacetime::PathBound`, and an
optional `Spacetime::PoseTable` supplies world rotation. Tiny steady edges
use `Plot::PointRoute`; the remaining polyline edges enter the column walk
as one `Plot::EdgeBatch`. Their still selection uses the ordinary shape draw.
Deforming rings opt into `Plot::RasterConfig::column_edge_lengths`: the same
walker derives their edge lengths from each column's vertices and uses signed
first-arm meridians for both arms.

#### Turning planar contours

`Plot::PlanarChords` strokes closed polylines that are straight in an
azimuthal-equidistant chart. Its timed `draw_closed` draws a contour turning
rigidly over the frame. It takes the frame-end vertices, the rotations from the
frame end to the frame start and to mid-frame, and a bound on the contour's
angular speed. `prepare` first records the frame's motion bound. Screen anchors
along each edge carry quadratic motion through the frame start, middle and end.
Coherent chord samples solve their column crossing directly; other chords emit
each column they reach at that column's time. An anchor interval takes a
display pole's chart instead when it lies within `POLE_PIECE_ROWS` of a pole,
or when its motion would bend more than `POLAR_TIME_ERROR_PX` row pitches over
the frame at its latitude. There `Plot::PolarWalk` (`polar_walk.h`) fits each
anchor's chart position as a `Plot::PolarAnchor`, a quadratic through the three
captures. It then intersects each chart-straight piece with the diameter that
is each crossed column's meridian, at that column's time.

#### Moving mesh faces

`Scan::deform_mesh` prepares a face scan from the effect's motion plan,
start vertices and indexed vertex evaluator. `Scan::detail::scan_swept_faces`
classifies each active face and reuses captured vertices between passes.
Faces clear of an arm cut share a scan with vertices captured at their
seeing-column times. Faces near a cut draw in arm-clipped passes; faces that
reach beyond the shared pass use a `Scan::FaceArmCopy` with its own vertices
on each arm. Across-cut vertex times continue the incident edges' times.

A face reaching a pole is placed in each half of the arm under a
`TimeWindow`. A polar or wide face whose motion would make those placements
stray too far instead uses `Scan::PoleFaceSweep`: it retains vertex
worldlines and is scanned one column pair at a time. The swept face scanner
prepares edge-plane polynomials and evaluates them at that pair's time.
This classification accounts for the face's azimuth span, motion and reach
toward a pole.

#### Scan bases and fields

Forward scans keep shapes built at the frame end and carry camera motion in
the rays. `Scan::ColumnBasisTable<W,H>` prepares two rotated vectors per
surface time from the orientation's sub-frame history. Its
`Scan::ColumnBasis` view describes the fixed-meridian direction rotated
from t(x) into the current endpoint's world. `Scan::ColumnBasisScope`
installs the view on the canvas; scans read `canvas.column_basis()`.
A pixel direction is built from azimuth[k] and up[k], shared by both arms.
The view's `sweep_arc` widens conservative bounds. A stationary frame or the
flat selection uses the fixed meridian.

`Scan::rasterize`, `Scan::Mesh`, `Scan::RingGroup` and
`Scan::DistortedRingStack` honour the basis. Compact rasterized shapes use
`Scan::lit_cap` to bound the columns that can see them. A birth fraction
skips columns lit before the shape exists. A drift adds a shape's independent
rotation still to come as `drift * (1 - t)`; `Scan::settle_drift` can place a
compact shape once at the middle of the times that see it. Ring-group slots
also interpolate their individual axes. `Scan::Volume` uses
`Scan::VolumeColumns`, with a centre and rotation per surface time.

Shaders working in a moving object's coordinates use
`Scan::ObjectBasisTable` and `Scan::Shader::walk_object_grid`; each sample
receives its object-space centre and column frame. A field shader with a
trailing float t is detected by `Spacetime::reads_time`.
`Spacetime::shade_at_time` and `shade_at_column` dispatch that argument.
`TimedShaderFn` is the corresponding interface for rasterize, point and
volume scans. Their flat selection shades at t = 1. `FacePaletteShader`
samples face alpha, palette crossfade and counterpart weight per column
while those values move; `Animation::OpLeg::Shading` supplies its retained
crossfade endpoint. `state_moves` requests the timed scan even when the
orientation is stationary.

`Scan::Shader::draw_column_times` prepares a shader state once for each
surface time and reuses it for both arms and across rows. HyperLattice
instead retains its render state and prepares trace contexts through
`draw_columns`, interpolating rotation phases, camera flight and continuous
parameters. Voronoi retains site endpoints, prepares moved sites once per
surface time, and shades each column from the candidates classified at its
block corners. Both hold geometry at discrete changes and retain their
instantaneous drawing selection.

#### Retained fields and pullbacks

Composed effects opt into `Pullback::Spec::TEMPORAL_SAMPLING` when their state
supports continuous interpolation; the whole roster enables it. Each column
samples between the previous and current animation
states at `(x % (W / 2)) / float(W / 2)`. Both opposing arms use the same
half-turn phase, independent of the rendering clip. This introduces up to one
frame of animation delay while making the endpoint of one sweep meet the start
of the next. The first draw uses one state for both endpoints.

`Scan::Shader::draw_column_times` prepares the pipeline once per sweep time,
for both columns lit at it, and reuses it across rows. Orientations use spherical interpolation; parameters use their
declared interpolation curves; clocks follow their shortest periodic arc.
Retained palette and hue tables supply two colorizer samples at the same
interpolated geometry, with their results blended in time. Noise tables are
baked once per frame. Effects without this opt-in retain the frame-wide path
and allocate no temporal history.


Further opt-ins need continuous parameter curves, clock advances below half
their period, and spatially periodic phase wraps. Discrete preset changes do
not become continuous through this option. GridSpace's whole-cell affine
winding satisfies the wrap condition.

A Mobius lens is sampled separately from its preset interpolation, which holds
the previous coefficients. `Animation::MobiusWarpCircular` publishes the angle
and magnitude it writes into b, and each column regenerates b from the shortest
arc of that phase over the current a, c and d. A lens whose b was replaced by a
snap or write, or whose a, c or d changed, takes the current lens in every
column. A column lens that fails `Pullback::valid` takes the current lens.

The chain interpreter retains the previous endpoint as a mirror of the
committed program's param and state blocks. Each column samples every
operator's params along their transition curves and its state through the
state's `temporal_sample` overload, which follows clocks and quaternion arcs
with `Spacetime::sample_clock` and `Spacetime::TemporalRotation`; states without clocks are read at
the current endpoint. `ShaderChain` keeps the previous `FrameContext` in a
`FrameHistory` whose palettes and hue tables point at owned copies, and
blends the colorize operator's colours from both endpoints' resources at the
column's interpolated controls. The first frame, a recompile or snapshot
restore, and changes to topology or snapping fields draw the current
endpoint in every column, except fields admitted by the operator's
`temporal_hold_field` policy. Modulation recipe changes keep blending endpoint
colours with the current modulation recipe; matching recipe topology allows
mechanic amount interpolation.

Manual shaders reuse the same facility, from
[field.h](../core/render/spacetime/field.h), over their
own frame struct.
`Spacetime::FrameHistory<Frame>` retains the previous endpoint of a frame
that owns or shares immutable resources, and `Spacetime::sample_frame` calls the
frame's static `sample(from, to, fraction)` for every sweep time before the
current endpoint. SphericalHarmonics samples its spin through `Spacetime::TemporalRotation` and its
amplitude linearly. A morph rollover continues the outgoing harmonic pair at
full blend into the incoming pair at zero blend; harmonic indices are never
interpolated. A manual preset changes the morph epoch, so that frame takes the
current morph in every column while the spin keeps its sweep time.


`Pullback::FrameInterval` prepares the composed field's endpoint intervals
and evaluates them with `at(t)`. `Spacetime::ClockInterval` and
`Spacetime::sample_clock` interpolate periodic clocks along their shortest
arc; `Spacetime::sample_frame` delegates a manual frame's interpolation to
its static `sample` function.

#### Stepped simulations

A `TimeWindow` ([canvas.h](../core/render/canvas.h)) restricts window-aware
column loops to the half-turn range [begin, end) on both arms, optionally
restricted to one arm, inside the driver's fixed clip. Windows do not nest;
closing one admits every column again. `Scan::Shader::walk_grid`,
`draw_grid` and `walk_object_grid` honour the window through
`Scan::walk_window_columns`. `Scan::Mesh` face and deformation passes and
`Scan::DistortedRingStack` also honour it. Generic `Scan::rasterize` writes
its whole clip.
`Spacetime::time_slab` divides the half-turn into ordered time slabs.
`Spacetime::simulate_in_slabs` advances the substeps that end each slab's
share of the frame, then shades under its window. The flat selection runs
all substeps before one frame-wide shade.

#### Pixel feedback

`Filter::Pixel::Feedback` reads the previous frame's buffer, in which every
column was lit at its own time. With column-time sampling,
`Filter::Pixel::FeedbackColumnTime`
([feedback_column_time.h](../core/render/filter/feedback_column_time.h)) gives
the frames between a target pixel and the source it reads, `1 + t(x) - t(s)`,
blending the lit times of a read between columns as its taps are blended.

- **Within an arm:** a Style warp displacing D columns per frame reads the
  source lit `HALF / (HALF + D)` frames earlier (`elapsed_scale`, at most 2).
- **Past the arm's forward cut:** it reads the other arm, lit a frame before
  that.
- **Behind the backward cut:** no stored frame holds the content, so the read
  slides from the cut column lit last to the arm's first column.

Lattice offsets carry the elapsed scale of their displacement, and the stored
column offset is clamped to half a turn before its 16-bit store. Cells whose
sources stay within their arm fade by the cell's elapsed time through
`Filter::Pixel::FeedbackFadeTable`, and the rest resolve per pixel. The flat
selection reads every source one frame back.

#### Verification

Native tests compare a timed draw, column by column, with the same geometry
frozen at that column's surface time
([temporal_oracle.h](../tests/temporal_oracle.h)), and a draw passes when both
of these hold:

- **Per column:** the summed luminance error is no larger than the error the
  frozen frame shows when its geometry is displaced 0.2 row along its meridian.
- **Whole frame:** the summed error is at most one twentieth of the endpoint
  render's error against the same frozen references.

Anti-alias taps reuse their sample's position except across arm cuts, or where
the error bound needs each tap column's own position. Frame-end parity hashes
pin the still selection.

### Spherical ray rendering

The reusable ray entry is [ray.h](../core/render/ray.h). `Raycast::PreparedCamera`
embeds unit display directions in 3D or a 4D slice with world-distance near/far
limits and a radial start. Pattern transforms are separate from the orthonormal
camera embedding. `VolumeQuery`, `PlacedQuery`, and `DomainQuery3/4` adapt existing
volume shapes and ambient fields without introducing another framebuffer traversal.

`surface_search` finds the first verified boundary, with independent exterior
and interior clearance guarantees and explicit unresolved/exhaustion outcomes.
Its bounded driver, `Raycast::march_steps`, is shared with
`Raycast::closest_approach`, the closest-approach search `Scan::Volume` uses;
that caller keeps its approximate silhouette and occluder handling. Surface
proximity alone never contributes opacity.

`trace_events` merges a compile-time bounded set of monotone analytic streams.
It preserves distinct material/merge identities, coalesces duplicate coverage
within a relative tolerance, and reports candidate/layer exhaustion. Adapters may declare a merge-group capacity; it defaults to their stream count.
The lattice uses one slot for a junction. Exceeding that bounded storage reports
exhaustion after consuming the valid pending contributions.
Geometry-specific candidates live in the lattice and triangular-framework SDF
headers. `Appearance`, `shade_events`, `shade_surface`, and `verified_filter`
reuse straight-alpha `LayerComposite` composition and shade verified subrays
independently. Unresolved subrays retain diagnostic failure status.
An `Appearance` needs finite, nonnegative `inv_far` and must keep
`t * inv_far <= 1` across its nonnegative trace interval;
`LatticeTrace::prepare` checks this when preparing a frame.

[RayStage](../core/render/pullback/ray.h) binds a frame-prepared renderer to the
existing `SphereSample -> Color4` contract. HyperLattice owns presets, palettes,
wrapped legacy flight animation, and its admitted tuple table; reusable geometry,
event traversal, shading, and composition reside in core. Curved periodic fields
and marched lattice/volume adapters remain explicit experiments pending resource
admission; merely including their templates emits no registered firmware effect.

## SDF Shapes (`sdf.h`) and the Scan Rasterizer (`scan.h`)

The rendering pipeline splits shape definitions from rasterization. `sdf.h` defines the SDF shape primitives, each implementing three methods:

1. **`get_vertical_bounds()`** — conservative row bounds, including the anti-aliasing fringe. Only rows within this range are scanned.
2. **`get_horizontal_intervals(y, out)`** — scanline intervals per row, or `false` to request a full-row scan. Handled rows skip empty columns without evaluating the distance function.
3. **`distance<ComputeUVs>(p, result)`** — signed shape-distance value at a sphere-surface point `p`, plus texture coordinate and auxiliary data in `DistanceResult`. Chart-plane values are not exact spherical boundary distances.

`scan.h` is an umbrella over `core/render/scan/`: `core/render/scan/raster.h` defines `Scan::rasterize()`, which drives the scanline loop and anti-aliasing; `core/render/scan/shapes.h`, `core/render/scan/mesh.h`, `core/render/scan/shader.h` and `core/render/scan/volume.h` hold the convenience wrappers that pair SDF shapes with the rasterizer.

`sdf.h` is an umbrella over six of the headers in `core/render/sdf/`: the substrate every shape shares (azimuth intervals, row bounds, `DistanceResult`, the cap/annular span-emission helpers) in `core/render/sdf/common.h`, the polygon, star, flower and line leaves in `core/render/sdf/shapes.h`, the ring leaves in `rings.h`, the CSG operators in `csg.h`, `SDF::Face` with its congruence-class LUT in `face.h`, and the volumetric family in `core/render/sdf/volume.h`. Including `sdf.h` pulls in all six; the lattice, trace, framework, periodic-surface and face-class headers are included directly by their consumers.

The `process_pixel` function applies anti-aliasing based on shape type:
- **Solid shapes**: quintic smoothstep over a 2-pixel AA band centered on the edge (`-pixel_width <= d <= pixel_width`). Full interior pixels (`d < -pixel_width`) skip AA math entirely. `pixel_width` is the compile-time constant `2π/W` — the angular width of one *equatorial* pixel — so the band is a fixed angular thickness at every latitude, and near the poles (where columns converge) it spans more than two columns. At 288×144 the row and equatorial column arcs nearly match. At 96×20 a row spans about 2.4 column arcs (2.5 on the ideal profile the host build uses), so a horizontal edge has a narrower AA band in row units than a vertical edge has in column units.
- **Strokes**: opacity falloff across the full stroke thickness.

### SDF Shape Primitives (`sdf.h`)

| Shape | Description |
|---|---|
| `SDF::Ring` | Geodesic circle at a given radius and thickness |
| `SDF::DistortedRing` | Ring whose centerline is shifted per azimuth, either by a callback (polar-separation distance) or by a knot polyline (exact polyline distance) |
| `SDF::PlanarPolygon` | Regular N-gon in the tangent plane of a basis vector |
| `SDF::SphericalPolygon` | Regular N-gon with geodesic (great-circle) edges |
| `SDF::Star` | N-pointed star using the standard inradius/circumradius construction |
| `SDF::Flower` | Inverted star (N-petal flower shape from the antipodal perspective) |
| `SDF::Line` | Geodesic line segment between two sphere-surface points |
| `SDF::Face` | Planar polygon face (used for mesh rendering) |

The table covers the effect-facing shapes. `sdf.h` also holds internal specializations that only the matching `scan.h` wrapper constructs — `SDF::FlatDistortedRing`, an undisplaced `DistortedRing` with an exact polar centerline distance, is instantiated by `Scan::DistortedRing::draw_flat` and never named by an effect.

### Volumetric Shapes (`sdf/volume.h`)

The 3D family marched by `Scan::Volume` lives in its own header. It shares the `SDF` namespace but not the scanline contract above: these shapes return a plain `float` distance in Cartesian ray-space, have no vertical bounds or horizontal intervals, and are reached from the `Scan::Volume` march loop and `Raycast::VolumeQuery` adapters. Everything it declares is a struct or a template, so an effect that draws only 2D shapes emits none of it.

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
| `Scan::RingGroup` | Fused single-pass rasterizer for a small group of rings — one scan over the union band paints every member in slot order, so the per-row interval math runs once instead of per ring. Fragments carry position, stroke coverage, size and age only (no UVs, no raw distance) |
| `Scan::Circle` | Disc (ring with radius-wide thickness) — stroke coverage ramps quintically from center to rim, and the shader styles it from register 2 |
| `Scan::Point` | Dot at a sphere-surface position, with the same center-to-rim coverage ramp |
| `Scan::Line` | Geodesic line segment between two points |
| `Scan::Star` | N-pointed star shape |
| `Scan::Flower` | N-petal flower shape |
| `Scan::DistortedRing` | Ring with per-azimuth radius perturbation |
| `Scan::DistortedRingStack` | Fused single-pass rasterizer for a same-axis stack of knot-mode distorted rings (shared basis, zero phase). A caller-owned `CandidateTable`, rebuilt each draw, maps the pixel's polar and azimuth cell to candidate rings; the per-pixel frame is computed once and the shader takes a ring slot alongside the fragment |
| `Scan::PlanarPolygon` | Regular N-gon in the tangent plane |
| `Scan::SphericalPolygon` | Regular N-gon with geodesic (great-circle) edges |
| `Scan::Mesh` | Rasterizes all faces of a `MeshState` |
| `Scan::Shader` | Full-screen per-pixel shaders with configurable SSAA (super-sample anti-aliasing), across typed entry points. `draw(canvas, shader)` takes a single fragment shader; `draw_cached(canvas, shader)` provides the same typed draw with its traversal placed in cached flash and is the path composed effects use. `draw(canvas, fragment_shader, vertex_shader)` separates a per-pixel vertex shader (called once at pixel center) from a per-subsample fragment shader (called SAMPLES×), so expensive per-pixel work is computed once — both callables are required, and a null one traps. BZReactionDiffusion uses `draw_grid`; GSReactionDiffusion uses `walk_grid(canvas, pixel_shader)`, which hands the centre vector and the row's `SsaaGrid` to the pixel shader. Under temporal sampling both use `walk_object_grid(canvas, basis, pixel_shader)`, which hands over the object-space centre and column frame instead. `draw_block_coherent(canvas, block, positions, scratch, classify, shade)` supplies Voronoi's block classification and shading path, and `draw_block_columns(canvas, block, scratch, classify, prepare_time, prepare_column)` the same classification column by column, preparing each sweep time once and each block of a column from its corner cells and their site union. `draw_columns(canvas, prepare)` prepares one shader per column and accepts Color4 or premultiplied Pixel results; `draw_column_times(canvas, prepare_time, prepare)` prepares state once per sweep time and a shader from it for each of the two columns lit at that time. `draw_grid(canvas, vertex_shader, pixel_shader)` hands the seeded fragment and the row's sub-pixel grid to a templated pixel shader that owns the sampling and returns the finished pixel. `draw_cached(canvas, shader, begin_row)` also takes a per-row setup callback while preserving wrapped clip margins. |
| `Scan::TransformedVolume` | Wraps an SDF shape with a world-space position and orientation quaternion for volumetric rendering |
| `Scan::Volume` | Volumetric ray-marcher that steps along the view direction through a `TransformedVolume`, applying a fragment shader at the hit point with configurable step count and AA width |

### Near-Pole Azimuthal LOD

A row at colatitude φ has horizontal pixel pitch that scales with `sin(φ)`; its ratio to row spacing also depends on the display aspect and latitude span. The scan walk offers those columns as a block of `pole_lod_aggressiveness / sin(φ)` (`core/render/render_policy.h`, clamped to `POLE_LOD_MAX_RUN = 32`), and the sink settles the whole block from one probe wherever the probe can vouch for it. Only full canvas-aligned blocks are offered, so an offer never straddles two blocks and a settled column always takes its shade from its own block's anchor. A block truncated by a clip or span edge goes per column instead, so the columns beside a segment seam shade at full resolution rather than from the anchor the neighbouring segment would have used.

`pole_lod_aggressiveness` is a hardware-calibrated knob, not a derived constant: the true masking width depends on the display aspect, LED angular size and per-column exposure. Smaller values reduce the run; 0 makes every offer one column and the walk bit-identical to an undecimated one. It defaults to 0 (`HS_POLE_LOD_DEFAULT`). Firmware compiles it in as a `constexpr` with no setter — at the default, the decimation branches fold away entirely — while host and WASM builds keep it mutable so it can be tuned live ([README §7.2](../README.md#72-the-wasm-bridge) `setPoleLod`).

The knob reaches the walk, not every primitive. `Scan::RingGroup` and `Scan::DistortedRingStack` replace the per-ring walk with one fused scan over the group's union band and shade every column of it, so raising the knob leaves them undecimated. Their equivalence to rasterizing the members one by one is stated at aggressiveness 0 for exactly that reason.

### Continuous Radial Surfaces

`core/render/scan/radial_surface*.h` renders a continuous family of nested
planar stars as a single surface. Drawing a dense stack of contours one by one
approaches this surface as the spacing shrinks; the radial surface evaluates
the limit directly, so its cost scales with pixels and sweep intervals, not
with contour count.

**Model.** One unit star is scaled and rotated by a complex multiplier
`c(s) = r(s)·e^{iφ(s)}` along a sweep parameter `s`. Each chart point shows
its *first arrival*: the smallest `s` whose star contains it. The query
returns that arrival as an interval and fraction, and the caller maps it to a
colour. When colour varies continuously with `s`, the interior is a smooth
gradient that a centre sample reproduces exactly. Colour jumps only where
first arrival jumps:

- **Terminal boundaries**: the outline of the last star, and of the first.
- **Tip traces**: the path of an outer tip that sweeps sideways into space no
  earlier star reached. Points on either side arrive at very different `s`,
  which shows as a fold where a twisting family overlaps itself.
- **Knot-edge minima**: segments of a knot's star where the piecewise-linear
  sweep changes direction, so first arrival has a local minimum along the
  edge.

**Supported family.** `RadialSurface<MAX_INTERVALS, MAX_SIDES>::bind` takes:

- 3 to `MAX_SIDES` outer tips;
- an inner-to-outer vertex ratio in `(0, cos(π/sides)]`, where the upper limit
  is the regular `sides`-gon;
- an interval count up to `MAX_INTERVALS`.

`set_knots` then supplies `intervals + 1` knots, each a multiplier and its
unwrapped argument. Only the first knot may be zero, and adjacent arguments
must differ by less than one edge sector, `π/sides`. Between knots the
multiplier is linear in the complex plane. Any radius and rotation profile
within these limits is accepted.

The family is restricted to this shape because the closed-form query depends
on three properties:

- Every star is a similarity of the unit star, so each edge's inside test is
  a fixed quadratic in `c` and therefore in the interval fraction.
- One shared edge constant serves every edge, because the star is regular.
- Each edge owns an equal angular sector, so the edge facing a point is one
  multiply.

Irregular, morphing or curved profiles need a different kernel.

**Point query.** `RadialSurface::sample(p, theta, first_interval)` tests the
facing edge `e` with `cross(p, c·e) − h·|c|² ≤ 0`, where `h` is the unit
star's edge constant. With `c` linear in the interval fraction `t`, that test
is quadratic in `t`. The facing edge changes at most once per interval, when a
corner crosses the point's argument; the query solves for that event linearly
and splits the interval there. Each piece has a closed-form descending root, so
there is no iteration. Intervals are visited outward from the caller's lower
bound and the walk stops at the first hit. Blocks of `BLOCK_INTERVALS` whose
swept disk and angular range cannot reach the point are skipped. The walk is
bounded by `MAX_INTERVALS`.

**Edge area coverage.** `RadialSurfaceEdgeAA::sample` takes the pixel as an
affine square in the chart, centred at `p` with column and row derivatives
`dx` and `dy`. For each pixel it:

1. Runs the centre point query.
2. Walks intervals from the footprint's lower bound and records the candidate
   boundaries clipped to the footprint, keeping the nearest. Its trace point
   and birth stay within that clipped segment. A strict cross-product
   certificate skips clipping when the projected point is already contained;
   uncertain points retain the full clipping path. The walk ends
   once a star at or beyond the centre hit contains the whole footprint,
   because nothing later can show.
3. Returns the centre colour when every intersecting candidate starts after a
   star that already covers the footprint. Such boundaries are hidden.
4. Otherwise treats the nearest boundary as a straight line, computes the
   exact area of the pixel square on each side, and blends the centre colour
   with the colour across the line. A second point query just past the
   boundary supplies that colour. Samples stay within the affine footprint,
   including a line footprint at a display pole. A tip on a strictly growing
   family uses the clipped tip's own arrival: earlier stars have smaller outer
   radii and cannot cover that point. Where two compatible boundaries meet at
   a tip, their coverages multiply and their front sample lies inside both
   half-planes. A last-star outline behind another boundary queries its inside
   and outside on that boundary's back side. Its colour and background coverage
   are used only when those contained queries find a sheet and background,
   respectively.

Each pixel resolves at most one jump, against a straight-line boundary. Where
more sheets overlap within one pixel the result is a bounded local
approximation, not a full partition of the pixel. Every pixel costs at least
one point query; an AA pixel adds the boundary walk and one more query.

**Candidate mask.** Most pixels are nowhere near a jump, but finding that out
through the full AA path means building the footprint and walking the
boundaries. `RadialSurfaceCandidateMask<64, 32>` answers it from a table:

- **Layout.** Each hemisphere gets 64 rows in squared chart radius by 32
  angular cells. The cells span one tip period, `2π/sides`, and repeat around
  the circle. That is 8 KiB per hemisphere.
- **Contents.** Each cell holds 32 bits of source-interval groups, where a
  group is the interval index shifted right by `interval_shift()`.
- **Preparation.** `prepare(surface, footprint_radius, complex_motion,
  adjacent_phase_motion)` enumerates every candidate boundary in one period
  through `RadialSurfaceBoundaryCandidates`, pads it by the footprint radius
  and the motion bounds, and marks the cells it can reach.
- **Motion bounds.** `complex_motion` bounds how far any multiplier moves
  during the mask's lifetime. A family with fixed radii may also pass
  `adjacent_phase_motion`, the drift of adjacent knot phase differences, for a
  tighter bound.
- **Hidden boundaries stay marked.** The mask records geometry only;
  visibility is decided by the AA step.

`interval_bits(radius, theta)` returns zero when no boundary can reach the
pixel for as long as the mask is valid. Such a pixel takes the point query
alone, without constructing a footprint. Nonzero bits restrict the AA walk to
the marked groups.

**Drawing on the sphere.** `draw_radial_surface<W, H>(canvas, near, far, area,
basis, columns, prepare, ensure, lower_bound, shade)` splits the sphere at
the equator into two azimuthal-equidistant charts. They are centred on
`basis.v` (`near`) and `-basis.v` (`far`), so chart radius is angle from the
pole in radians and a family reaching `π/2` fills its hemisphere. Pixels are
visited by `Shader::draw_column_times`. `columns` is the per-column pullback
from `prepare_column_basis`, or empty for a frozen frame. The caller supplies:

- **`prepare(column)`**: called once per sweep time, to evaluate
  time-dependent state.
- **`ensure(far)`**: called before the first query in a hemisphere at that
  time, to rebuild that family's knots and to refresh its mask when the mask
  has expired.
- **`lower_bound(radius)`**: the first interval that could contain a point at
  that chart radius. It must never skip an interval holding an earlier hit.
- **`area(surface, p, theta, radius, footprint, lower_bound, shade, stats)`**:
  the standard `Scan::RadialSurfaceArea<Surface, COUNT, MASK>` integrator
  checks the hemisphere candidate mask and calls `surface.sample` for zero
  cells without building a footprint. Marked cells build `footprint()` and
  call `RadialSurfaceEdgeAA::sample` with the candidate bits. `MASK=false`
  sends every pixel through edge AA; `COUNT` enables the coverage census.
- **`shade(hit)`**: maps a first arrival to a linear `SurfaceColor` on the
  16-bit channel scale. The renderer clamps that colour and writes it as a
  pixel.

Surfaces and masks allocate from caller arenas. `RadialSurface::storage_bytes`
gives the bound for a family, and the kernel is emitted in flash.

The Phantasm linker groups the radial query, edge AA, candidate-mask preparation
and ShapeShifter surface column walk into a contiguous flash block before the
generic hot sections. Its start has a fixed offset within a 16 KiB cache way.
Name patterns select the input sections; a linker assertion
bounds the block to the 32 KiB instruction cache. Other effects' hot code follows
this block.

**ShapeShifter P1.** ShapeShifter routes its authored dense planar-star preset
through this kernel and keeps every other preset on its contour paths. The
route applies only at P1's exact parameters: 288 contours of a 7-tip sine
sweep with screen-balanced spacing. Each hemisphere family has 24 intervals
per radian of amplitude, capped at 240. Its radii follow the screen-balanced
spacing out to the equator and its rotation is the travelling sine twist. Its
colour runs from the palette's end at the pole to the palette's midpoint at
the equator. Knots are rebuilt for each hemisphere
at every column time. A hemisphere's mask is rebuilt once the column time has
moved more than eight column steps, and never across an arm cut. Both
families, both masks and the per-interval tables fit in two 32 KiB scratch
arenas.

`tests/test_radial_surface.h` checks the query and coverage against
independent swept-edge oracles, `tests/test_radial_surface_mask.h` checks that
static and moving masks cover every boundary, and
`tests/test_shapeshifter_surface.h` pins P1's column-time seams and image
parity.

## The Curve Rasterizer (`plot.h`)

`plot.h` is an umbrella over `core/render/plot/` (`core/render/plot/cull.h`, `core/render/plot/raster.h` with its sample and planar-arc helpers in `core/render/plot/raster_walk.h`, `core/render/plot/shapes.h`, `core/render/plot/mesh.h`, `core/render/plot/particles.h`, `core/render/plot/chords.h`, and `core/render/plot/point_route.h`). For drawing lines, curves, and paths, the `Plot` namespace provides a geodesic/planar rasterizer with adaptive step size. Each sub-step is sized from the curve's full 2-D screen-space speed (`sqrt(vx² + vy²)`, combining longitudinal and latitudinal motion), so samples land roughly one pixel apart everywhere on the curve regardless of latitude. The step is clamped to keep the equator near one sample per column and floored near the poles — where screen speed diverges — so pole oversampling stays bounded.

```cpp
Plot::Line::draw<W, H>(pipeline, canvas, start, end, fragment_shader);
Plot::Multiline::draw<W, H>(pipeline, canvas, vertices, fragment_shader);
```

`Plot::Multiline` accepts a `Fragments` array (an arena-backed `ArenaVector<Fragment>`), while `Plot::Line` accepts its two `Fragment` endpoints. The parametric primitives (`Ring`, `Polygon`, `DistortedRing`, `Star`, and `Flower`) take their geometric parameters and sample into a `Fragments` array internally. Each fragment carries position, texture registers (v0–v3), age, and color.

- **Edge interpolation** — how consecutive fragments are joined. *Geodesic* (the default) walks the great-circle arc between endpoints; *planar* interpolates along an azimuthal-equidistant straight line in a basis's tangent plane (for effects that live in a 2D local space). `Plot::Polygon` and `Plot::Star` select it with their `Projection` template parameter, the other parametric primitives fix it, and `rasterize` takes it as `RasterProjection::planar(basis)` or `RasterProjection::geodesic()`.

### Sampling Policy

`rasterize` takes its compile-time behavior as one `RasterConfig` NTTP (`rasterize<W, H, RasterConfig{.single_pass = true}>`), whose `sampling_policy` field sets the adaptive sample density. `DEFAULT` targets `SCREEN_STEP_PX` (0.9 px) and compiles the alternative away; `BALANCED` always trades samples for speed; `SELECTABLE` defers the choice to `RasterOptions::balanced_sampling`, so one instantiation serves both and the policy is picked per draw call. Only the single-pass rasterizer reads it — the cached-replay path always samples at the default density.

`RasterOptions` remains a by-value aggregate. `RasterLoop::closed(seam)` attaches optional seam registers to a closed path. `RasterProjection::planar(basis)` selects a planar chart; `RasterProjection::geodesic(flags)` accepts a bounded visibility span. `PointProjections` pairs equal-sized row and column arrays, with `paired(rows, cols)` for dynamic spans. The rasterizer checks span lengths against the actual polyline.

Balanced sampling stretches each adaptive step by `BALANCED_SCREEN_STEP_PX / SCREEN_STEP_PX` (1.25×), clamped to one base step (2π/W) and left exact below the pole floor (`MIN_POLE_SCALE * BALANCED_POLE_GUARD_SCALE` base steps), where spacing is already at its minimum. Two consequences:

- **Emitted alpha changes.** Sparser samples lay down less coverage per unit arc, so each fragment's alpha is scaled by `balanced_sample_alpha()` — gain `1 + (ratio - 1) * (0.88 - 0.20 * alpha)`, saturating at 1, with `ratio` the balanced step over the default step. The gain shrinks as alpha rises because opaque samples compound less; it is a linear fit to source-over accumulation, close to exact below alpha 0.4 and over-boosting above it, so a stroke past alpha ~0.85 saturates at 1 and loses its soft edge. A balanced draw is not pixel-identical to a default one.
- **Step evaluation is reused, on planar edges only.** Where the walk is locally straight and clear of the poles (tangent dot > 0.995, step change under 10%, step under 0.9 base steps), the next sample recomputes position only and carries the previous step forward, skipping the tangent and the screen-velocity step. The reuse needs the monotonic position-only entry that only `PlanarEdgeSampler` exposes; a geodesic edge takes the sparser steps and the alpha gain without it.

`ShapeShifter` reaches balanced adaptive sampling on its stationary route: its own sampled draw enables it for spherical stars at 32 or more contours, and `Plot::PlanarChords` hands stationary dense planar-star pole intervals to `rasterize` with balanced sampling. The planar pole intervals can reuse their step; spherical-star edges are geodesic.

Authored sine-driven presets P2 through P9 instead carry contour twist and orientation to column time. Moving planar contours capture start, middle and end anchors, cull intervals from their motion bounds, and use coherent screen chords where their physical curvature permits. Other intervals use `Plot::PolarWalk`: quadratic trajectories in the display pole's equidistant chart intersect each column at its own time, with radial AA taps and the adaptive rasterizer's minimum physical pitch. Spherical contours use the shared column edge walk with their local lattice pitch. P1 retains stationary sampling. Custom non-sine waveforms retain frame-end twist under column-time orientation.

`Plot::GeometryCallback::uniform_color` declares a color-only shader whose color is constant along each edge throughout the frame. The column walk shades that color once per edge, using the shader's world or source coordinates at the frame end. Age and material registers retain their per-sample interpolation, including taps across an arm cut.

### Plot Primitives

| Primitive | Description |
|---|---|
| `Plot::PlanarChords` | Stationary and column-time chart-straight contour stroker |
| `Plot::PolarWalk` | Column-time AA strokes in a display pole chart |
| `Plot::PlanarBandSplit` | ClipBand edge splitter with per-edge visibility flags |
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

## The Animation System (`animation.h`)

The `Timeline` class manages a list of running `IAnimation` objects. Each frame, `timeline.step(canvas)` advances all active animations. Finished animations are removed; repeating animations are rewound. All animation types inherit from `AnimationBase` and support method chaining via `.then()` for sequencing.

Preset pause holds automatic preset and phase progression while rendering and ambient motion continue. Parameter-driving timeline events opt in to the same pause gate. Effects schedule parameter drivers with `timeline.add_pausable(..., &anims_paused)`, and preset choreography passes the same flag as its transition event's pause gate (including every in-flight `Lerp` departure); while the flag is set, both the animation step and any pending start delay are frozen. Passing a pause pointer to an animation constructor freezes only its `step()` call, so the timeline event's start delay still elapses. Events added with `add()` and motion advanced directly by `draw_frame()` continue to run, which lets the GUI hold parameter tuning without stopping ambient motion.

`animation.h` defines the contract every animation implements — `IAnimation`, the CRTP `AnimationBase`, and `Animation::Space` — and then includes nine fragment headers grouped by what they animate:

| Header | Subject | Contents |
|---|---|---|
| `timers.h` | Callbacks on a clock | `RandomTimer`, `PeriodicTimer` |
| [params.h](../core/animation/params.h) | A caller-owned parameter, written each frame | `Transition`, `Mutation`, `Progress`, `Driver`, `Lerp`, `ColorWipe`, the `Mobius*` family, `Ripple`, `Noise`, `BallDrop`, `NoiseProduct` |
| [motion.h](../core/animation/motion.h) | An `Orientation` driven through space | `Path`/`ProceduralPath`, `Motion`, `Rotation`, `RandomWalk` |
| `core/animation/trails.h` | Recorded history | `Trail` and its `OrientationTrail`/`VectorTrail` aliases — index 0 is the oldest snapshot and `length()-1` the newest, the ordering the JS simulator mirrors — plus `QuantizedVectorTrail`, the `TrailBody` per-body aggregate, and the `tween`/`deep_tween` traversals |
| `core/animation/sprites.h` | Visible things | `Sprite`, `Particle`/`ParticleSystem` |
| [timeline.h](../core/animation/timeline.h) | Scheduling | `TimelineEvent`, `Timeline` |
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
| `Mutation` | Applies a custom scalar function to a float variable over time with easing; `Mutation::Options` names repeat and step-only pause settings |
| `Progress` | Invokes a caller-supplied `void(float)` callback once per frame with eased progress, leaving the caller to write whatever state it drives |
| `Driver` | Continuously increments a float variable each frame (optionally wraps at 0..1) |
| `Lerp` | Type-erased interpolation between any `T` that implements `lerp(start, target, t)`. The caller owns start, subject, and target data; Lerp holds pointers and a type-erased lerp function. |
| `ColorWipe` | Smoothly interpolates a `GenerativePalette` between caller-owned start and target snapshots that must remain unchanged and outlive the animation |
| `ParticleSystem<W, CAPACITY>` | Physics simulation with emitters, attractors, friction, gravity. The default `Particle<TRAIL_LEN>` storage has `QuantizedVectorTrail` history for trail rendering. The trailing storage policy can select `PointParticle` with `TRAIL_LEN=0`: 24 bytes per star, two-float octahedral direction, float velocity, color seed and lifetime, without history. |
| `Ripple` | Animates a `RippleParams` to expand a Ricker wavelet across the sphere |
| `MobiusWarp` | Animates `MobiusParams` to apply and release a Möbius transformation |
| `MobiusWarpCircular` | Animates `MobiusParams` for a circular warp that stays warped throughout, suitable for repeating effects |
| `MobiusWarpEvolving` | Modulates `MobiusParams` through independently wrapped double-precision phases at multiple frequencies |
| `MobiusFlow` | Animates `MobiusParams` for a continuous loxodromic flow |
| `Noise` | Animates `NoiseParams` along a bounded, smoothly reversing time loop for flowing distortion fields |
| `BallDrop` | Animates a `BumpParams` to drop one spherical-cap bump from the north pole to the south along a fixed meridian, ramping the footprint envelope so the bump emerges from and vanishes into the poles |
| `NoiseProduct` | Animates `NoiseProductParams` along the same bounded time loop to flow a two-octave product-noise field; perpetual |
| `OpLeg` | Animates one leg of a mesh operator chain: a Conway-operator parameter sweep along a `ConwayGraph` edge or recipe step, a hankin contact-angle sweep on a fixed seed, a relax or medial slerp, or a gated partition swap. Each frame it rebuilds the swept mesh in scratch, compiles it, checks the compiled face count against the hoisted classification, maps faces to pre-blended ramps via the leg's hoisted topology classification, and hands the mesh and per-face ramp table (`Shading`) to a draw callback — exactly one mesh drawn per frame. Each leg's `.then()` completion handler schedules the next, so a run of legs walks a whole morph path. |
| `MeshCarousel<SegueT>` | Owns two mesh states and schedules the selected segue. The effect builds the back slot, calls `set_front`, then `schedule_segue`; the schedule forwards the preset pause gate. Shipping rendering supports `Segue::TerminatorSweep` in IslamicStars and standalone `Segue::Crossfade` in DreamBalls (with zero overlap). The other seven policies are library prototypes without production render consumers: IrisBloom, Lace, Shockwave, Breakdown, SpinFlip, GoldConvergence and Dissolve. Their `warp`, `reorder`, `mask_pair`, and fragment hooks require a matching effect draw path; changing the policy alone does not wire those hooks. |

Noise and NoiseProduct use a double-precision phase with a 4096-unit period.
Their published time coordinate stays linear over the initial phase range
[-960, 960], then reverses through C2 turns at coordinate extrema ±1000.
The field eventually retraces its trajectory. Constructors interpret the
current time as a phase seed modulo the period; a non-finite seed initializes
the internal phase to zero. Zero and non-finite speeds preserve the published
coordinate. Live configuration refresh and later writes to the published time
do not reseed the clock. Copies keep independent clock values.

MobiusWarpEvolving allocates 64 bytes of channel phases in the persistent arena.
Copy construction allocates an independent block; copy assignment reuses its
destination block or allocates one when the destination was moved from. Moves
transfer the block without allocating and leave the source empty. These
blocks remain reserved until the arena is reset and must outlive their
animations. GnomonicStars spawns one pinned warp per effect initialization.

### Orientation history

`Orientation<CAP>` stores a history of up to `CAP` quaternions (default 4) accumulated during one frame step. The template parameter is the history *capacity*, not the display width — the `Filter::World::Orient` and `Filter::World::OrientSlice` filters use `Orientation<>`, never `Orientation<288>`. Index 0 is the previous frame's end pose and the last index the current pose; `get(float t)` samples that sub-frame path at frame fraction `t` by slerping between the bracketing entries. The explicit-time rotation in `Filter::World::Orient` reads a `Spacetime::PoseTable`, the path's rotation matrices at every column time and the frame end, built once per frame in storage `init_storage()` takes from the persistent arena for the effect's canvas width. Ordinary plots use the current endpoint. Both routes preserve the `age` register. Edge clip culling uses the end pose with `sweep_arc()` padding; still-raster sample spacing visits the history poses independently.

```cpp
timeline.add(0, Animation::Rotation<W>(orientation, math::Y_AXIS, 2 * PI_F, 600, ease_linear, {.repeat = true}));
// orientation.length() upsamples to min(substeps + 1, CAP) when larger
// Filter::World::Orient samples the steps at each point's capture time
```

`Orientation::upsample(count)` resamples the orientation history to a higher resolution via SLERP. This is used to rewrite the history when combining multiple animations for accurate parallel sub-frame path tracing — ensuring that concurrent rotations, motions, and walks all contribute to a single coherent set of intermediate orientations.

### OrientationTrail

`OrientationTrail<OrientationType, CAPACITY>` maintains a circular buffer of past `Orientation` snapshots, allowing effects to recall where an object was over previous frames. Each snapshot is a full `Orientation` (with its own sub-frame history).

### VectorTrail and QuantizedVectorTrail

`VectorTrail<CAPACITY>` maintains a circular buffer of past world-space `Vector` positions — 12 B per sample, stored exactly. Used by HopfFibration for its per-fiber trails.

`QuantizedVectorTrail<CAPACITY>` is the unit-sphere variant: each sample is three snorm16 components (6 B, half of `Vector`), clamped to [-1, 1] on record and decoded by value on `get()`, with a quantization step of 1/32767 and additional float rounding during encoding and decoding. Used by `ParticleSystem` to record per-particle trajectories for trail rendering.

### `tween` and `deep_tween`

Two traversal helpers linearize multi-level orientation history into a single callback loop:

| Function | Input | Description |
|---|---|---|
| `tween(orientation, callback)` | `Orientation<CAP>` | Iterates over the sub-frame quaternion history of a single orientation, calling `callback(quaternion, t)` for each step with `t ∈ (0, 1]`. Sub-frame 0 is the pose carried over from the previous frame's end and is skipped unless it is the only snapshot (which reads `t = 1`, age-neutral). |
| `deep_tween(trail, callback)` | `OrientationTrail` (any `Tweenable`) | Flattens a trail of orientations into a single continuous traversal, calling `callback(quaternion, t)` with a global `t` spanning all frames and sub-frames. Used by the orientation-trail effects (Comets, Fishbowl) for rendering trails with full sub-frame accuracy. A bare `Orientation` has no per-frame structure to flatten and is rejected by the `Tweenable` concept — use `tween` for that. |
| `deep_tween_frames(trail, callback)` | `OrientationTrail` (any `Tweenable`) | One callback per contributing trail frame, `callback(qs, ts, count)`; concatenated `(qs, ts)` follows `deep_tween` emission order. `qs` points into frame storage and `ts` into a buffer reused across callbacks; neither pointer outlives the call. RingSpin preserves frame boundaries with this traversal. |

### Animations and Mutable State

Parameter and motion animations mutate external state that the rendering pipeline reads. Sprite animations call their draw function directly. The state-driving animation types target these variables:

| Animation | Target State | What It Mutates |
|---|---|---|
| `Rotation`, `RandomWalk`, `Motion` | `Orientation<CAP>` | Quaternion orientation — upsamples the history to the frame's sub-step count (capped at CAP) and writes each sub-step's rotation, which `Filter::World::Orient` samples at capture time |
| `Transition` | `float*` | Smoothly interpolates any float parameter (e.g. `speed`, `alpha`, `twist`) from current value to target with easing |
| `Mutation` | `float*` | Applies an arbitrary scalar function `f(t)` to a float over time (more general than `Transition`) |
| `Progress` | `void(float)` callback | Hands the caller eased progress each frame and writes nothing itself; every composed preset transition uses it to blend the authored parameter states |
| `Driver` | `float*` | Continuously increments a float each frame, optionally wrapping at 0..1 — used for phase accumulators |
| `Lerp` | `T*` (type-erased) | Interpolates any type with a `lerp()` member — params structs, styles, etc. The caller owns start, subject, and target; Lerp holds pointers |
| `ColorWipe` | `GenerativePalette*` | Interpolates palette keys between caller-owned start and target snapshots in OKLCH; both snapshots must remain unchanged and outlive the animation |
| `Ripple`, `MobiusWarp`, `Noise` | `RippleParams`, `MobiusParams`, `NoiseParams` | Animate transformer parameters (expansion radius, warp strength, noise time axis) which the transformer pool reads during `MeshOps::transform()` |
| `BallDrop` | `BumpParams` | Walks the bump center down a meridian and re-derives the push axis from the stack's orientation, ramping the footprint envelope; the field pool sums the caps during `field()` |
| `NoiseProduct` | `NoiseProductParams` | Advances the bounded field time loop under live speed edits; the field pool reads it during `field()` |
| `ParticleSystem` | `Vector[]` positions | Physics simulation updates particle positions; the default storage records `QuantizedVectorTrail` history for trail rendering |

Effects can declare *what state exists* (orientations, floats, palettes) and schedule animations to drive it. The `Timeline` handles timing, easing, sequencing, and cleanup for those scheduled animations; effects can also advance state directly in `draw_frame()`:

```cpp
// Effect declares mutable state:
Orientation<> orientation;   // CAP is the sub-frame capacity, not the display width
float twist = 0.0f;
GenerativePalette palette;
GenerativePalette target_palette;
const GenerativePalette::Snapshot palette_start = palette.snapshot();
const GenerativePalette::Snapshot palette_target = target_palette.snapshot();

// Timeline drives state via animations:
timeline.add(0, Animation::Rotation<W>(orientation, math::Y_AXIS, 2 * PI_F, 600, ease_linear, {.repeat = true}));
timeline.add(0, Animation::Transition(twist, 2.5f, 1000, ease_in_out_cubic));
timeline.add(0, Animation::ColorWipe(palette, palette_start, palette_target, 2000, ease_linear));

// Rendering reads state — no manual updates needed:
void draw_frame() {
    Canvas canvas(*this);
    timeline.step(canvas);  // all state updated automatically
    // orientation, twist, palette are now current-frame values
    pipeline.plot(canvas, v, palette.get(t).color, 0.0f, 1.0f);
}
```

## Geometry Transformers (`transformer.h`)

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

Transformers integrate with the `MeshOps::transform()` pipeline and can be chained: `MeshOps::transform(input, output, arena, ripple_transformer, orient_transformer)`. `transform()` takes any callable with `operator()(Vector)`, so a derived pool and a plain adapter compose in the same call.

### Displacement Fields

`FieldTransformer<ParamsT, AnimT, FieldFunc, CAPACITY>` is the scalar counterpart: entities superpose by summation instead of composing as warps, so an effect can feed the summed field into a displacement path (e.g. a `DistortedRing` shift LUT). `field(p)` sums the active entities; `field_bound()` returns a per-frame upper bound on `|field()|` for sizing conservative culls. Where overlapping bodies must dominate instead of stacking, feed the per-entity values through `DominantFieldAccumulator`, which blends them magnitude-weighted (`sum(s³)/sum(s²)`).

| Field | Effect |
|---|---|
| `BallDropTransformer` | Spherical-cap bumps that fall pole-to-pole through a frame, bowing the surface away from each cap |
| `NoiseProductTransformer` | Two-octave product noise where octave 1 envelopes octave 2, so perturbations bunch where the envelope runs strong |

### Pool Lifecycle

Both classes derive from `TransformerPool`, which fixes the call order:

1. `init_storage(Arena&)` — from the effect's `init()`, after any `configure_arenas()` and before the first spawn. It also claims one of the shared `Timeline`'s `MAX_CLEAR_HOOKS` (4) clear-hook slots, and each `ChoreographedEffect` claims one hook, so a choreographed effect fits at most three pools; exhausting the hook table traps at registration.
2. `spawn(in_frames, args...)` — the returned pointer is transient; use it at the call site, not across frames. The animation must be finite and non-repeating; anything else traps (use `spawn_pinned()`).
3. `spawn_pausable(paused, in_frames, args...)` — same as `spawn()`, but the whole timeline event, start delay included, freezes while `*paused` is set, the way `Timeline::add_pausable` does. It is the only pool entry point that honours a GUI pause — `spawn()` animates straight through one — and the flag must outlive the event.
4. `spawn_pinned(in_frames, args...)` — same as `spawn()`, but the pointer may be retained (e.g. registered as a live GUI param). Valid only for an animation that never completes on its own — infinite, or repeating (it rewinds rather than reaching `done()`) — and is added before any finite, non-repeating timeline event, so compaction cannot shift it.
5. `prepare_frame()` — each frame before `transform()` / `field()`, whenever active params changed through animation or live config. The composition reads that prepared state but cannot verify it is current. Its per-entity hooks (`refresh_from(const ParamsT&)` for live config, `sync()` for derived state) are found by detection, so every `ParamsT` must declare one bool per hook (or supply them through a `transformer_detail::ExternalParamsHooks<ParamsT>` specialization) — `static constexpr bool NEEDS_REFRESH_FROM` and `static constexpr bool NEEDS_SYNC` — each true only when that hook is carried. Either mismatch is its own `static_assert`, which is what turns a renamed or signature-drifted hook into a compile error instead of a silently unrefreshed entity.
6. `reclaim_storage(Arena&)` — from the after-reset callback of an arena that is compacted mid-effect (e.g. a mesh carousel). Spawned animations hold `Params` references into the slots, so the caller must replay the same allocation order after the reset as after `init_storage()`; the re-claimed blocks must land at their original addresses (asserted). A reset only rewinds the offset, so the untouched bytes carry live entities through.

### Standalone Utilities

`OrientTransformer<CAP>` (`transformer.h`) is a plain adapter struct, not a `TransformerPool`-derived pool: it holds a reference to an `Orientation<CAP>` and applies `orientation.orient()` to each vertex. It has no pool, no params and no lifecycle — effects construct one on the stack at the call site (a deduction guide takes `CAP` from the orientation) and hand it straight to `MeshOps::transform()`.

## Memory Architecture (`memory.h`, `memory.cpp`)

A single contiguous memory block (`GLOBAL_ARENA_SIZE`) is partitioned into three arena allocators. It is 298 KiB on Teensy; the WASM module widens it to 512 KiB (`HS_GLOBAL_ARENA_BYTES` in the root [CMakeLists.txt](https://github.com/woundedlion/pov/blob/master/CMakeLists.txt)), because the chain interpreter's two arenas, `ShaderChain`'s two column-time blocks of the same size and its shared color resources outgrow the device-mirrored default. Individual effects can call `configure_arenas()` to repartition the block at runtime.

| Arena | Default Size (Teensy) | Purpose |
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

Geometry operators take explicit `Arena&` parameters. Global scratch consumers include `MeshCarousel::compact_*`, `OpLeg`, `Filter::Pixel::Feedback::flush()`, `Plot::rasterize`, `gate_trail_edges`, `Plot::Mesh`, `Plot::ParticleSystem`, the cull projection helpers, and `RecipeBuild`; see [Why Arena Allocation?](../README.md#33-why-arena-allocation). Scratch lifetimes are LIFO: callbacks allocate only under their own `ScratchScope` and never call `reset()`, which would invalidate their caller's live buffers.

```cpp
ScratchScope scope_a(scratch_arena_a);
ScratchScope scope_b(scratch_arena_b);
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

## The Color System (`core/color/color.h`)

`Pixel` stores RGB channels in **16-bit linear light** (`uint16_t r, g, b` in range 0–65535). This avoids the precision loss and incorrect blending that occurs with gamma-encoded 8-bit values.

The conversion pipelines:
```
Input (sRGB 8-bit) → sRGB→linear LUT → Pixel (linear 16-bit) → blend ops
                                                                      ↓
FastLED output ← CRGB(linear PWM) ← linear16_to_linear8 ← Pixel
DMA HD107S output <- adaptive 5-bit brightness + 8-bit BGR pack
                  <- brightness, temperature and color gains <- Pixel
```

`Color4` wraps `Pixel` with a float alpha channel. The canvas sink composites with a single straight-alpha "over" operation — `blend_alpha(α)`, i.e. `dst = src * α + dst * (1-α)`, applied in 16-bit linear light (see `core/color/pixel.h`). There is no selectable blend-mode tag.

### Palette Types

| Type | Description |
|---|---|
| `ProceduralPalette` | Cosine palette: `a + b*cos(2π*(c*t + d))` per channel. Defined by 4 vec3 coefficients. |
| `Gradient` | OKLCH interpolation between a sorted list of (position, color) stops. |
| `GenerativePalette` | Procedurally generated palette from harmony rules (triadic, analogous, etc.) combined with brightness/saturation profiles. Supports snapshot/lerp for animated transitions. |

Twenty-seven named `ProceduralPalette` instances are pre-defined in the `Palettes` namespace: `DARK_RAINBOW`, `BLOOD_STREAM`, `VINTAGE_SUNSET`, `RICH_SUNSET`, `UNDERSEA`, `LATE_SUNSET`, `MANGO_PEEL`, `ICE_MELT`, `LEMON_LIME`, `ALGAE`, `EMBERS`, `FIRE_GLOW`, `DARK_PRIMARY`, `MAUVE_FADE`, `LAVENDER_LAKE`, `DESERT_ROSE`, `BRUISED_MOSS`, `BRUISED_BANANA`, `BRIGHT_SUNRISE`, `FIRE_AND_ICE`, `PEACH_POP`, `TIDAL_JADE`, `BLUE_LAGOON`, `ORANGE_CRUSH`, `PLUM_SUNRISE`, `CORAL_BLUE`, and `BRUISED_MANGO`. Six of them — `EMBERS`, `RICH_SUNSET`, `BRIGHT_SUNRISE`, `BRUISED_MOSS`, `LAVENDER_LAKE`, `TIDAL_JADE` — are the slots of `MeshPaletteBank`, the shared baked bank the mesh effects draw from.

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
| `hue_rotate()` | Perceptual hue rotation — rotates the (a,b) chroma plane in OKLab, preserving lightness and chroma before gamut mapping. Forward nonlinearity uses `fast_cbrt` (hot per-pixel path); inverse is exact. Out-of-gamut results are chroma-reduced rather than per-channel clipped, which holds hue and stabilizes the feedback loop against saturated-color drift. The feedback `hue_fade` transform uses `hue_rotate_lms_matrix`; sphere-space hue noise uses `hue_rotate_lut_gamut`. |

### The Gamut Boundary Grid

The chroma clip brackets the sRGB boundary from a generated table (`core/color/gamut_lut.h`, emitted by `tools/gen_gamut_lut.py`) indexed by the diamond angle of (b, a) and by L. Each cell stores guarded extrema of sampled first-exit chroma estimates, which do not certify enclosure of every ray's exact first exit; the per-pixel path walks that bracket in `GAMUT_SCAN_STEPS` and bisects the straddling step `GAMUT_BRACKET_STEPS` times. Grid resolution sets how wide the bracket starts, and the bisection sets how far it is narrowed; coarser grids increase the chroma deficit (about 0.003 at 256 × 128, 0.0043 at the minimum accepted 128 × 64). At every grid, including the shipped master, a handful of rays stride over a narrow out-of-gamut gap and land past the first exit — an oversaturation of up to about 0.05 chroma that no amount of bisection recovers, since the wrong step is already selected; a finer grid reduces how many rays do this, not how far they overshoot.

The clip reads the 256 × 128 flash master by default. An effect that clips per pixel can arm an arena copy at the master's resolution or coarser, which improves read latency (RAM rather than QSPI flash). Coarser grids also change the brackets and clipping results:

| Function | Description |
|---|---|
| `init_gamut_lut(arena, angle_steps, l_steps)` | Downsamples the flash master into `arena` and points the clip at the copy. Both step counts must divide the master's 256 × 128 and stay at or above `GAMUT_LUT_MIN_ANGLE_STEPS` × `GAMUT_LUT_MIN_L_STEPS` (128 × 64), the coarsest accepted grid (measured deficit 0.0043 and 9 overshoot rays) — both trapped. Costs `gamut_lut_bytes(angle_steps, l_steps)`. Call from the effect's `init()`, after any `configure_arenas()`. |
| `release_gamut_lut()` | Drops the copy and points the clip back at the flash master. Registered as an `ArenaResetHook`, so `configure_arenas()` and the mesh carousel's compaction both run it before handing the storage out again. |

`MeshFeedback::init()` is the only production call site that arms an arena copy. It takes the full 256 × 128 grid, `gamut_lut_bytes(256, 128)` = 131,074 B of persistent arena. The Shader workbench, `ShaderChain`, and composed effects clip against the flash master without allocating a gamut copy.

The feedback flush's colour path (`hue_fade_apply`, `hue_fade_apply2`) does not walk the bracket: `lms_cbrt_transform_rgb_lut` (two pixels at a time in `lms_cbrt_transform_rgb2_lut`) rescales an out-of-gamut pixel's chroma straight onto the armed grid's cell minimum through `lms_cbrt_scale_to_gamut_lut`, so a clipped trail pixel costs one grid read. The sampled cell minimum has no certified deficit bound against the exact boundary. The feedback loop re-clips every bright trail pixel each frame, which is where the per-pixel walk was spent.

### Palette Modifiers

Modifiers compose around any palette source at compile time via
`StaticPalette<Source, Coords<...>, Colors<...>, Wrap, Shade>`. There are two axes: a
**coordinate** chain that remaps the lookup parameter `t` *before* the source is
sampled, and a **color** chain that reshapes the resulting sample *after*, with
the coordinate selected by `Shade`. Both chains are inlined by fold expression with
zero runtime overhead. `Wrap` (default `true`) wraps the final coordinate into
`[0,1)` before the lookup — leave it on for cycling modifiers that overflow the
range; set it `false` for bounded remaps that must reach the source endpoints.
Both directions are `static_assert`ed: an unbounded modifier rejects
`Wrap=false`, and a bounded final modifier rejects `Wrap=true` (wrapping would
fold its 1.0 output to 0.0 and destroy the top endpoint). Only a modifier that
re-bounds *arbitrary* input (`WrapModifier`'s fold, `FoldModifier`'s triangle
wave, `InsetModifier`'s clamp) clears an unbounded predecessor; `ReverseModifier`
and `PinchModifier` are bounded on `[0,1]` but pass an out-of-range coordinate
straight through, and `MirrorModifier` maps one to another out-of-range value. `QuantizeModifier` caps values above 1 and quantizes negative
coordinates without rebounding them. Chaining any of these after a cycling
modifier needs a `WrapModifier` between them and `Wrap=false`.

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
| `RippleModifier` | Applies a travelling sine distortion to the lookup parameter |
| `FoldModifier` | Folds the parameter space (mirror at edges) to create ping-pong patterns |
| `PinchModifier` | Non-linearly warps the lookup parameter toward a focal point; the per-sample `powf` suits bake-time sampling |
| `QuantizeModifier` | Posterizes the palette into discrete bands |
| `ScaleModifier` | Multiplies the lookup parameter by a scale factor |
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
| `LightnessGrainShade` | Grains brightness with evolving noise; uniform linear-RGB gain before channel rounding and clipping (a gain above 1 clips bright channels) |
| `IridescentShade` | Adds a thin-film cosine sheen with per-channel phase offsets, saturating at white |

The noise-driven modifiers sample the deterministic `value_noise_1d`/`value_noise_2d`
hash lattice (`3dmath.h`) with a per-instance seed, so two modifiers on the same
driver decorrelate by seed. Frame-constant work memoizes against the driver
value (`HueSpinShade`'s rotation matrix, `ChromaPulseShade`'s pulse factor,
`DriftModifier`'s walk offset). The OKLab shades still pay a per-sample
conversion, so they pair well with `BakedPaletteStorage::rebake`, which re-samples a
256-entry LUT once per frame (the source must be a `Wrap=false` composition;
`rebake` rejects wrapping sources); the noise and cosine shades are cheap enough for
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
| `ColorModulation::PositiveLightnessLift` | Raises OKLab lightness over positive scalar lobes, preserving alpha and exact zero identity. |
| `LayerComposite` | Composes straight-alpha color layers. |
| `ColorModulation::PatternModulatedPalette<PaletteSource, Field, Mechanic>` | Composes a palette with an immutable typed scalar field and an independent color operation. Sphere fields use admitted directions and explicit direct or baked-cube evaluation. |
| `ColorModulation::HueRotation` | Applies signed periodic hue rotation, preserving alpha and exact zero identity. |
| `ColorModulation::Opacity` | Multiplies straight alpha by a scalar mask while preserving RGB. |

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
| `chroma` | `ChromaControls`: an `axis` (`AxisControls`: curve/center/range/custom), plus `basis` (`ChromaBasis::LOCAL_GAMUT` or `ABSOLUTE`; `PATH_MINIMUM` holds its ordinal but is unimplemented and fails compilation) and `headroom` |
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
`tonal_monochrome()` — and `effects/common/palette_recipes.h` holds the
per-effect recipes the roster renders.
`effects/generated/mind_splatter_palette_bank.h` is the MindSplatter recipe
baked at 256 base hues by `tools/mindsplatter_palette_gen.cpp`.

### Palette Cycling

`PaletteCycler` (`core/color/palette_cycler.h`) drives a display LUT through a
sequence of palettes over time: it dwells on an entry for `dwell_frames`, then
fades into the next over `fade_frames` (a dwell of 0 behaves as 1: the next fade begins on the following step).
Effects call `step()` once per frame and shade from `palette()`, a
`BakedPalette` — outside a fade, after `step()`, the display is a bit-exact
bake of the current entry; after `advance_without_display()` it holds the last
displayed bake until the next `step()`.

The fade mechanism is chosen per adjacent pair at `init()`. Two morph-compatible
`GenerativePalette`s fade by **key-space morph**, interpolating control keys for
perceptually coherent hue travel; every other pair — composed, prebaked, or
morph-incompatible — falls back to a **baked-LUT crossfade**. An `Entry`
accordingly accepts a `GenerativePalette`, any `Palette`, or a `BakedPalette`,
so a mixed sequence of up to `MAX_ENTRIES` cycles correctly. Entries are
caller-owned and must outlive the cycler; `Entry`'s rvalue overloads are deleted
so a temporary cannot bind.

`init_generated()` replaces the fixed entry array with a `NextPaletteFn`
provider that fills palette *n* on demand, without a fixed entry-count limit;
successive palettes must stay morph-compatible and every retarget is fail-fast
checked. Its `display_interval_frames` sets the fade display-bake cadence
(1–255 frames, at most the fade length), with fade weight sampled at each
interval start; the clock and provider still advance every frame, and fade
endpoints bake exactly. Composed effects use an interval of 2 frames.
Arena cost is declared up front — `display_arena_bytes()`,
`crossfade_arena_bytes()`, `morph_arena_bytes()`, and the
`required_arena_bytes()` / `generated_arena_bytes()` worst cases — so an effect
sizes its arena without probing. `advance_without_display()` keeps the timeline
and the provider moving without rebuilding the LUT; a later `step()` rebuilds it
at the current sampled phase.

## The Mesh System (`core/mesh/`)

The mesh system uses these headers in `core/mesh/` and `core/render/sdf/`:

- **`base_mesh.h`** — Base mesh identities, bounds, and authoring labels
- **`relax_bake.h`** — Relax payload and source identity checks
- **`core/mesh/mesh.h`** — Core data structures (`PolyMesh`, `HalfEdgeMesh`) and fundamental `MeshOps` (compile, clone, classify)
- **`conway.h`** — Conway mesh operators and vertex transformations
- **`conway_graph.h`** — Constexpr 23-edge morph graph over the 18 simple-registry solids: per-edge operator/seed/reseed specs, shortest-path routing to uniformly random destinations, and the closed `ORDERED_TOUR`
- **`recipe_types.h`** — The authored op-chain model: the `Op` operator set, one `OpStep`, and the `Recipe` chain a registry generator mirrors, split out so the model is not read out of the registry tables written in it
- **`recipe.h`** — Lowers an authored recipe to primitive steps (`expand_to_primitives`), sizes that lowering at compile time (`lowered_step_count`, `max_lowered_step_count`), replays either form through `SolidBuilder` (`build_recipe`, `build_steps`), and decides which lowered steps a morph leg can sweep (`is_morphable_step`)
- **`core/mesh/hankin.h`** — Hankin pattern compilation and dynamic update
- **`core/render/sdf/face_class_bake.h`** — Congruence-class clustering plus one canonical distance-LUT bake per class, allocated by descending face count under an 18 KB per-mesh budget
- **`core/render/sdf/face_classes.h`** — The class id space and the three record types the rasterizer binds per frame, split out so the clustering and bake machinery stays out of every rasterizer translation unit
- **`mesh_state.h`** — `MeshState`, the flat-array renderer format, split out so mesh, Conway, Hankin and solids code can share the renderer-facing representation without the construction machinery
- **`solid_generators.h`** — Umbrella over `solid_tables.h` (hardcoded Platonic vertex/face tables), `solid_builder.h` (the `SolidBuilder` operator chain) and `procedural_solids.h` (the named Archimedean / Catalan / Islamic Star Pattern generators)
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

All Conway *geometry* operators (`dual` through `bevel` below) take `(const PolyMesh& mesh, Arena& target, Arena& temp)` and return a `PolyMesh`; `truncate`, `expand`, `chamfer`, `snub` and `bevel` take a defaulted shape parameter `float t` after the arenas, and `snub` a `float twist` after that. `medial` takes the same two arenas but writes its two outputs through reference parameters — `(const PolyMesh& mesh, PolyMesh& out_a, ArenaVector<Vector>& out_b, Arena& target, Arena& temp)`. `ambo`, `truncate`, `expand`, `chamfer` and `snub` each carry a second overload taking a caller-built `const HalfEdgeMesh&` right after `mesh`: the output topology is the same at every `t`, so one connectivity build serves a whole parameter sweep instead of one per call. The one exception is `truncate` at exactly `t == 0.5`, which short-circuits to `ambo` and its different face census — a sweep that must hold one topology stops short of it (`ConwayGraph::T_EPS_AMBO`). `transform`, `transform_in_place`, `relax`, `relax_baked`, `normalize`, and `reconcile_vertices` are listed in the same table but are mesh utilities with their own signatures. Every operator, primitive or composed, produces its `PolyMesh` into `target` and may use `temp` for intermediate computation; composed operators (`gyro`, `meta`, `needle`, `zip`, `bevel`) reuse the same internal ping-pong as their constituent ops, with an even-length composition starting its first step in `temp` so the last one lands in `target` (see the COMPOSITION POLARITY note in `conway.h`):

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
| `MeshOps::relax_baked` | Substitute a flash-baked relax result for the pass. Its runtime checks catch a source/bake mismatch (dimensions, topology hash, quantized source-vertex hash) and payload corruption (output hash, re-derived from the bake's own vertex bits) — they say nothing about freshness, since a change to `relax` itself that leaves its input intact passes these checks. Freshness is the `unit_relax_bake_verify` ctest's job: it re-runs the live `relax` and asserts bit-exact equality with the committed payload |
| `MeshOps::normalize` | Project all vertices onto the unit sphere |
| `MeshOps::reconcile_vertices` | Certify and reorder a nearest-vertex bijection into a non-scratch target arena using 4%/8%/10% z bands |

### Hankin Pattern System (`core/mesh/hankin.h`)

| Operation | Description |
|---|---|
| `MeshOps::compile_hankin` | Pre-compute topological data for fast Hankin pattern updates |
| `MeshOps::update_hankin` | Update dynamic vertices based on angle parameter (no new allocation when reusing a sufficiently-sized output mesh) |
| `MeshOps::hankin` | One-shot Hankin pattern generation (compile + update) |

`compile_hankin` produces a `CompiledHankin` struct containing base vertices, static midpoints, and dynamic instructions. `update_hankin` evaluates the dynamic vertices by sweeping the Hankin angle, producing the star polygon line intersections for each face. It re-binds the output mesh's vectors on every call, so it avoids new allocation only in the steady state — reusing the same output mesh against the same arena, already sized large enough.

### Solids Library (`solids.h`, `solid_generators.h`)

`solid_generators.h` provides constexpr vertex/face data for all Platonic solids plus procedural generators for Archimedean, Catalan, and Islamic Star Pattern families; `solids.h` organizes them into three registries. Firmware effects build by registry index via `Solids::get_entry(index)`; `Solids::get_by_name(arena, a, b, name)` builds by name and fails fast on an unknown name; the WASM bridge validates with `Solids::find_entry(name)` first and generates from the entry; the WASM geometry tools enumerate the registries by index with `Solids::get_entry(index)` to populate the picker, then build the selected solid by name:

| Registry | Count | Description |
|---|---|---|
| `simple_registry` | `PLATONIC_COUNT + ARCHIMEDEAN_COUNT` | 5 Platonic (tetrahedron through icosahedron) + 13 Archimedean solids |
| `catalan_registry` | `CATALAN_COUNT` | Duals of the Archimedean solids (triakis tetrahedron, rhombic dodecahedron, pentakis dodecahedron, etc.) |
| `islamic_registry` | `ISLAMIC_COUNT` | Complex multi-operator recipes producing Islamic star patterns from base solids |

Total: `Solids::NUM_ENTRIES`. Counts are derived in `solids.h` from registry sizes and the `BaseMesh` Platonic boundary; assertions pin the registry slice boundaries and `BaseMesh` census.

`Collections` namespace provides typed spans for iterating subsets: `get_platonic_solids()`, `get_archimedean_solids()`, `get_simple_solids()`, `get_catalan_solids()`, `get_islamic_solids()`.

IslamicStars selects its playback presets by `IslamicStarsDetail::SHAPE_NAMES`,
in that list's order. The geometry registry also retains shapes outside the
playback roster for the solids tools.

`SolidBuilder` provides a fluent interface for chaining Conway operators with automatic arena swapping:

```cpp
return SolidBuilder(to_polymesh<Icosahedron>(a), a, b)
    .truncate()
    .dual()
    .build();
```

Islamic Star Pattern recipes chain multiple operators with Hankin pattern generation:

```cpp
SolidBuilder(Platonic::dodecahedron(a, b), a, b)
    .hankin(54.0f * D2R).ambo().hankin(72.0f * D2R).build();
```

## Generators (`memory.h`)

`memory.h` provides (through `core/memory/generate.h`) a generation wrapper that manages arena lifecycle around a supplied callable:

```cpp
namespace hs {
template <typename GenerateFn, typename... Args>
auto generate(Arena &target, GenerateFn &&fn, Args &&...args);
}
```

It resets both scratch arenas only at the outermost call (depth zero), scopes them on every call, then invokes `fn(target, scratch_a, scratch_b, args...)`. Nested calls preserve the caller's live scratch allocations. Callers can use the wrapper for registry lookups and effect geometry creation:

```cpp
auto mesh = hs::generate(persistent_arena, Solids::get_by_name, std::string_view("icosahedron"));
```

`SolidBuilder`'s fluent Conway chain (`solid_builder.h`) owns its own two-arena ping-pong, swapping the scratch arenas between operators, so it manages arena lifecycle directly rather than through `generate()`.

## The Preset System (`control/choreography.h`)

`ChoreographedEffect<Derived, Params>` provides the runtime preset lifecycle: an `Effect` base owning the effect's live parameter set, its preset table, and the choreography that moves between presets. `Derived` declares its presets — a `PRESETS` table (`std::array<PresetEntry<Params>, N>`) or a bound static `preset(index)` returning one entry with `PRESET_IDS` naming its rows; a `PRESETS` table may also supply `PRESET_IDS` — plus a dwell, a parameter schema version, and either `parameter_fields()` (from which the base derives `valid_params()` and `blend_params()`; a `Params::FIELDS` table supplies it) or its own `valid_params()`. Each entry pairs a parameter set with its departure, the `Segue::Preset::Departure` that carries the effect from that preset to the next:

```cpp
static constexpr Segue::Preset::Fade DEPARTURE{16};
static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
static constexpr uint16_t PRESET_DWELL_FRAMES = 224;
static constexpr std::array<PresetEntry<Params>, 12> PRESETS = {{{..., DEPARTURE}, ...}};
static bool valid_params(const Params &p);
```

The effect calls `begin_choreography()` once from `init()`, which configures the engine preset controller from the presets (or `PRESET_IDS`), and `step_choreography()` every frame, which retires the dwell and starts the departure to the next preset. `cycle_frames()` conservatively sums every `PRESET_DWELL_FRAMES` dwell and every departure but the last preset's. A Lerp or Fade completing frame also consumes the next dwell's first frame, so the actual departure-to-departure period is dwell + departure - 1. The show budget carries one slack frame per included Lerp or Fade before rounding up to seconds; Snap adds no slack. A departure must be shorter than the dwell. A single-preset effect compiles the dwell countdown out. An effect that clocks its own advancement (DreamBalls' sprite hand-off chain, IslamicStars' per-shape segue) skips the per-frame tick and calls `advance_preset()` from its own schedule instead.

Automatic transitions follow the departing preset's policy — `Segue::Preset::Lerp` crossfades the live parameters into the target, `Segue::Preset::Snap` adopts immediately, `Segue::Preset::Fade` dims and brightens, adopting the target at the first step at or past half progress (an exact black frame requires a positive even frame count) — while manual and synchronized selections always snap. Departures are chosen per preset, so an effect can morph between related presets and fade between unrelated ones. Pausing freezes a fade at full opacity. `choreography_hooks()` declares a type-only `Control::HookSet` of `Control::HookBinding` entries. Each entry binds a customization point to its static or member function pointer; exact signatures, duplicate entries, fallback bindings, and declaration ownership are checked at construction. Composed effects extend their base list with `extend<Derived, ...>()`, or explicitly select `Control::HookOwnership::INHERITED` when using only composed defaults. Bound hooks specialize the mechanics: a member `preset_params(index)` overrides the entry's parameters when the effect patches them at runtime, e.g. re-binding a noise pointer; `initial_params()` overrides the first preset as the startup default; binding `adopt_params(target)` re-derives dependent state after a snap; `transition_armed(target)` fires once as a crossfade arms, capturing the endpoint state the blend interpolates alongside the parameters; `blend_params(progress)` writes an in-flight Lerp; `finish_blend(target)` can specialize automatic target adoption at a Lerp endpoint or Fade midpoint, retaining fields owned by independent animations; `set_preset_opacity(value)` receives a Fade's opacity. A Lerp departure requires `blend_params` and a Fade requires `set_preset_opacity`; both are checked at compile time. The base also carries schema-versioned parameter snapshots: `serialize_parameters()` tags the live set with `PARAMETER_SCHEMA_VERSION`, and `restore_parameters()` rejects a snapshot taken under a different schema or failing `valid_params()`.

`control/preset_host.h` holds the controller the choreography drives: the committed index, the vetoable `apply_preset()` hook, and the manual `selectPreset`/`nextPreset`/`previousPreset` surface the WASM bridge calls. `control/presets.h` holds the table vocabulary: `PresetEntry<Params>` (the row type) and the free `constexpr` helper `all_presets_in_ranges(entries, in_ranges)`, which folds a slider-range predicate over an entry table so an effect can `static_assert` its whole preset table against its registered parameter ranges — a loop rather than an unrolled conjunction, so appended entries are covered automatically.

Promoted composed effects also have authored shader documents under `patterns/`. `scripts/promoted_digests.test.mjs` checks the header/document/product-group roster; `scripts/generate_composed_presets.test.mjs` pins descriptor and preset-bank digests against effect headers; the runtime uses `ChoreographedEffect` to play the resulting presets.

## Hardware Drivers (`dma_led.h`, `pov_single.h`, `pov_segmented.h`)

Three hardware drivers form a layered stack. The DMA LED layer handles the SPI wire protocol across four headers: `hd107s_frame.h` (wire format and inline color correction), `dma_led_core.h` (framing, transfer length, stale-transfer predicate), `dma_led_controller.h` (double-buffer orchestration, templated on its transport) and `dma_led.h` (the Teensy SPI/DMA peripheral driver, the only Arduino-only piece). `pov_single.h` and `pov_segmented.h` sit above it and manage the POV column sweep, differing only in how many Teensys share the work; the segmented driver's ISR decisions are split into two further host-tested headers, `pov_handoff.h` and `pov_submit_gate.h`. The host-tested protocol is split across `pov_sync.h`, `pov_sync_protocol.h`, `pov_sync_flywheel.h`, `pov_sync_content.h`, and `pov_sync_emitter.h`; `pov_segment_frame.h` preserves the inactive frame half; `POVSegmented::pack_column` packs columns, and `pov_segment_map.h` owns segment-ID and physical-strip index arithmetic.

### DMA LED Controller (`dma_led.h`, `hd107s_frame.h`, `dma_led_core.h`, `dma_led_controller.h`)

Non-blocking DMA-based LED output for HD107S (APA102-compatible) LEDs on Teensy 4.x.  Enabled by `#define USE_DMA_LEDS` in the target's boilerplate header (`targets/Phantasm/phantasm_target.h`) before it includes the driver; `led.h` stays neutral and the default FastLED/WS2801 path remains as fallback. The FastLED fallback applies only to the single-board `POVDisplay`; the segmented `POVSegmented` driver `#error`s without `USE_DMA_LEDS` (FastLED's bit-bang `show()` masks IRQs for windows that break the sync symbol margins, which are derived from a mask window M ≈ 0), so DMA LEDs are mandatory on Phantasm.

| Class | Header | Role |
|---|---|---|
| `HD107SFrame<N>` | `hd107s_frame.h` | Pre-formatted DMA buffer for the HD107S protocol. `pack_pixel()` writes `Pixel` values directly into the frame buffer with inline color correction (color correction → temperature → brightness, then adaptive 5-bit pixel brightness and 8-bit linear PWM), bypassing the CRGB intermediate. The buffer is 32-byte-aligned (`__attribute__((aligned(32)))`) and cleaned with `arm_dcache_flush()` (clean, no invalidate — the buffer is TX-only) for cache coherency. |
| `TeensySPIDMA` | `dma_led.h` | Low-level DMA+SPI driver wired to LPSPI4. Configures a `DMAChannel` with completion interrupt for fully async byte-stream transmission. |
| `DMALEDController<N, Transport>` | `dma_led_controller.h` | Double-buffered high-level controller. The ISR packs pixels into `back_frame()`, then `submit_frame()` flushes it and triggers async DMA, returning immediately. If the previous transfer is still in flight, `submit_frame()` **drops** the new frame (bumping `get_overrun_count()`) and returns false rather than spinning; a transfer that never completes is surfaced as a wedged-channel fault. The drop returns before the buffers swap, so `back_frame()` still holds the dropped pixels and a caller can re-submit them without repacking. |
| `next_buffer()`, `transfer_len()`, `transfer_us()`, `transfer_stale()` | `dma_led_core.h` | Free `constexpr` framing and watchdog math the controller's decisions derive from, host-tested without the peripherals. |

The canvas and color-correction pipeline use 16-bit linear values. `pack_pixel()`
packs them into shared 5-bit pixel brightness and 8-bit linear PWM channels on the SPI wire. Brightness changes at the dominant channel's rounding boundary to keep gray ramps monotonic; the channels are rescaled and rounded at that brightness. Black uses brightness 31 and zero PWM. This preserves dim detail without another gamma conversion:

```cpp
// ISR path (per column): fetch the display buffer once, index it directly
const Pixel* buf = effect->display_buffer();               // 16-bit linear pixels
// Physical LED index comes from the single-source-of-truth map (pov_single_map.h),
// which applies the top-arm reversal / bottom-arm offset — never the raw row index.
frame.pack_pixel(pov::strip_top_led(y, S), buf[y * width + x]); // Pixel → HD107S frame
(void)led_controller.submit_frame();                          // non-blocking DMA; returns false on overrun
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

The top arm's physical LED ordering is reversed (LED 0 at the S pole, row S/2-1; midpoint LEDs S/2-1 and S/2 at the N pole), and the bottom arm shows the opposite half of the image (x offset by W/2).

`show_col()` discards `submit_frame()`'s overrun verdict: this driver carries no retry latch and no dark fallback, so a dropped column leaves the previous column lit for one extra period — the drop returns before the buffer flip. `run()` fail-fast-checks that one composite transfer fits inside a column period, which rules out the systematic overrun that would hold the strip on one frame and is what makes discarding the verdict sound.

The shipping `holosphere` sketch runs `RingSpin<96, 20>` with column strobing disabled. Its 1302 µs column period clears the 1181 µs FastLED transfer bound; the roster’s default RingSpin configuration keeps strobing enabled for the DMA targets.

| Parameter | Value (Holosphere) |
|---|---|
| S (total pixels) | 40 |
| RPM | 480 |
| Column interval | ~1302 µs (= 125 ms / 96 columns) |
| ISR duration | ~20 µs on the DMA path (`holosphere_dma`), which packs the column and starts an asynchronous transfer. The FastLED path clocks 40×24 bits at 6 MHz (at least 181 µs per show). A nonstrobed column normally takes at least 181 µs after sufficient idle time; `CMinWait<1000>` can add up to 1000 µs. A strobed column adds an immediate blank show and its 1000 µs wait: at least 1362 µs even when the first show does not wait, exceeding the ~1302 µs period. The current transfer guard budgets 1181 µs per show (2362 µs for strobing) and refuses a column period at or below that bound. |


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

**ISR state machines**: the wake's non-trivial decisions are split out of the Arduino-only driver into two host-tested headers, which `run_wake_sequence()` in `pov_submit_gate.h` drives in ISR order:

| State machine | Header | Role |
|---|---|---|
| `EffectHandoff<T>` | `pov_handoff.h` | Foreground↔ISR effect ownership: the teardown counter handshake, the acquire/release publish and adopt of a pending effect, the consumed-generation gate that keeps the ISR off a deleted instance, and the display-window (clip) alternation. The foreground constructs and deletes instances; the ISR only ever dereferences what `live()` handed it. |
| `SubmitGate`, `SyncPulseGate` | `pov_submit_gate.h` | The LED transport's accept/drop verdict and the sync pin's pulse width. Both submit paths — the fail-dark black frame and the image column — clear their pending state only on an accepted submit, and a dropped column latches retries on subsequent wakes, without repacking, until the transport accepts it or a new column or dark state replaces it. A wake that renders nothing has too short a body to carry a scheduled sync pulse, so the pin is held HIGH across the ISR boundary and dropped at the head of the next wake. |

**Effect transparency**: Effects are written against the full 288×144 canvas with no per-segment code. Each board clips rendering to its half-width segment band for the display window the frame will be shown in, one ahead of the open one (`clip_to_segment`), except stateful effects (`needs_full_frame()` / `persists_pixels()`), which render the full canvas; the ISR then packs this board's LEDs. Every board reseeds the shared `Pcg32` at every effect build from `HS_PHANTASM_EFFECT_SEEDS[]`, which `targets/Phantasm/phantasm_playlist.h` builds as `hs::stable_effect_seed(hs::effect_seed_key<name<CANVAS_W, CANVAS_H>>(#name))` (`core/platform/effect_identity.h`) so an entry's stream follows its immutable RNG key, independent of canonical ID and roster position; `hs::epoch_seed(effect index)` (epoch 0 is the identity seed `1337`) is the fallback for a board that supplies no seed table. Aligned epoch reconstruction starts the same content stream on every board. A mid-show join can have a frame offset until a subsequent aligned epoch; the beacon index alone does not align rendered frames.

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

`epoch` is folded forward by exactly one half-revolution at every boundary crossing. Servicing the flywheel within its signed-safe coast window preserves phase across the 32-bit cycle counter's ~7.16 s wrap; an unattended full wrap loses revolution history. Within that window an interrupt-masked interval does not accumulate column-phase drift, though intervening column outputs can be skipped — the ISR that runs after the mask reads the clock and resumes at the *time-correct* column.

**Pin / signal description.**

| Signal | Master (seg 0) | Downstream (segments 1–N−1) | Electrical | Direction |
|---|---|---|---|---|
| `SYNC` (rev 1.1) | pin 3 (GPIO out) | pin 3 (ext. interrupt in) | 3.3 V CMOS, active-high pulses, idle LOW | master → all |
| `SYNC` (rev 1.2) | pin 4 (TX through U1), pin 3 (RX) | pin 3 (ext. interrupt in), pin 4 held LOW | 3.3 V GPIO / 5 V bus, active-high pulses, idle LOW | master → all |
| `MASTER_EN` | pin 5 (drive LOW) | pin 5 (drive HIGH) | 3.3 V CMOS, gates the external sync-out buffer | per-board strap out |
| `ID` straps | pins 21–23 (ID0–ID2 as required) | pins 21–23 (ID0–ID2 as required) | active-low, internal pull-ups; `ID_STRAPS = log2(N)` read | strap (board identity) |

The ID straps select the board: the build reads `ID_STRAPS = log2(N)` active-low bits and decodes `(~raw) & (N-1)`. N=4 reads ID0/pin 21 and ID1/pin 22; N=8 also reads ID2/pin 23. All-floating selects segment 0/master. On rev 1.1, `SYNC` uses pin 3 — the master drives it and downstream boards receive on its rising edge; `MASTER_EN` (pin 5) gates an external level shifter so only the master drives the shared bus. `SYNC` is the only inter-board connection. Rev 1.1 leaves pin 4 unused. Rev 1.2 transmits from pin 4 through U1 and receives on pin 3 on every board; downstream boards hold pin 4 LOW. It is assumed physically reliable; a severed wire is out of scope (boards free-run and precess apart at crystal rate, a slow smear, never an instant break).

**Signal levels & symbol waveforms.** The wire idles LOW.  A **symbol** is a burst of short active-high pulses at a fixed pitch; **the meaning is the count of rising edges — pulse width carries no information.**  The pin is driven HIGH by the flywheel tick that schedules the pulse, not at ISR entry, and the rising edge is the only timed event.  A wake that also renders drops the pin before it returns, for a ~8–13 µs pulse — the usual lit-column path, since pulses fall on column boundaries.  A wake that renders nothing holds the pin across the ISR boundary and drops it at the head of the next wake, ~54 µs later.  Pulses are drawn narrow, to scale against the ~868 µs pitch:

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

**Why count, not width:** on the i.MX RT each pin has a single latched interrupt flag, so an IRQ-mask window *delays* an edge's ISR but cannot lose the edge unless two edges fall inside one mask window.  With pulse pitch chosen **greater than the worst-case mask window M**, the edge *count* is exact across interrupt-masked windows (Phantasm requires DMA LEDs).  The alphabet is **odd-only, distance 2** — a single lost or spurious edge lands on an even (invalid) count and is discarded.  A single lost or spurious edge degrades to a *missed* symbol (covered by the local boundary crossing). Multiple glitches can instead produce another accepted odd count.

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
* **Layer 3 — Content.** The playlist is **epoch-counted**, not `millis()`-gated.  Duration is **per roster entry**, not uniform: `HS_PHANTASM_EFFECT_LIST` carries a seconds column beside each name and `targets/Phantasm/Phantasm.ino` converts it to `EFFECT_REVOLUTIONS[]` at `seconds · RPM / 60`, spanning 30 s (240 revolutions) to 306 s (2,448 revolutions) across the playlist.  The master emits the `EPOCH` mark (plus R = 3 redundancy repeats) when the current entry's revolutions elapse; every board counts down to the same **absolute** commit boundary regardless of which copy it heard, constructs the next roster entry during the final K = 2-revolution **construction window** (display black on all boards simultaneously), and all swap to its frame 0 at the same boundary.  The beacon broadcasts the absolute effect index so a board that missed every epoch repeat corrects after two confirming beacons (up to 21 revolutions between beacons; recovery within 4 s, spec section 9.1), and a rebooted board rejoins at the correct effect — **fail-dark, never fail-wrong** (a board with no established identity shows black rather than a guessed effect).  Every one of those revolution budgets is absolute, so on a 240-revolution entry the 25-revolution rejoin bound still costs 10% of the effect's airtime.

**Index beacon frame format.** The beacon is a **data** symbol (integrity by *rejection*, not by exactness).  Five base-8 digits at 1-column pitch, each digit a burst of `digit + 1` pulses, digits separated by 5 quiet columns (one past the gap timeout, so the decoder reliably terminates each digit):

```
 Frame = [ idx_hi  idx_lo  rev_hi  rev_lo  checksum ]   (5 digits, base-8)
           └── effect index 0–63 ──┘ └ rev mod 64 ┘  └ d0⊕d1⊕d2⊕d3

 digit Dk transmitted as (Dk+1) pulses @ 1-col pitch, then a 5-col quiet gap:

         D0          D1               D2          …        D4
        ┌┐          ┌┐┌┐┌┐           ┌┐┌┐                 ┌┐┌┐┌┐
 ───────┘└──/ /─────┘└┘└┘└──/ /──────┘└┘└──/ /───────────┘└┘└┘└──────────
        │←Dk+1 pulses→│   │←5-col quiet (terminates digit)→│
        │←─────────── frame ≤ 55 columns / 23.9 ms (≪ half-rev) ──────────→│
```

Any checksum mismatch, wrong digit count, out-of-range digit, or stale partial frame **drops the whole frame** — the next beacon is ≤ 16 revolutions (2 s) away, or ≤ 21 (~2.6 s) across a commit window.  Schedule: revolution 1 of every 16 (`rev ≡ 1 mod 16` — never rev 0, so a just-powered board meets clean isolated boundary symbols first), plus the first revs of a fresh effect; silent during a pending commit.

**Receiver state machine.** Each downstream board is in one of two states.  The master is born `LOCKED` with identity (effect 0, rev 0) — it *is* the reference and never snaps:

```
                  first accepted snap
   ┌──────────────┐ ──────────────────▶ ┌──────────────┐
   │   ACQUIRE    │                      │    LOCKED    │
   │  (display    │                      │ (disciplined,│
   │   black)     │ ◀────────────────── │  rendering)  │
   └──────────────┘  reject_fallback = 4   └──────────────┘
                     gate rejections
                     (~2 revolutions)

 ACQUIRE : accept any *valid* symbol unconditionally (hard snap), but only
           on a burst preceded by ≥ t_QB (16 col) of wire silence — so a
           beacon digit train can't capture a just-rebooted board mid-frame.
           The train's FIRST digit is preceded by silence exactly as a
           boundary symbol is, so it can still be mistaken once; the
           reject_fallback rejection threshold bounds the recovery (spec §9.1 mis-snap row).
           Renders black until it has BOTH phase (a snap) AND identity
           (beacon).
 LOCKED  : accept a valid symbol only if implied correction ≤ G (4 col) AND
           boundary identity matches the prediction. Else reject (telemetry,
           no snap, no flip). After reject_fallback rejections the board concludes its OWN
           timebase is at fault and falls back to ACQUIRE (the escape hatch
           that stops a genuinely-lost board from rejecting good symbols
           forever).
```

**Epoch commit sequence.** `EPOCH` at ZERO boundary **B** schedules an absolute commit at **B + R + K** (R = 3 repeats, K = 2 construction revolutions). A board with an absolute synchronized revolution count infers a repeat's position in the train and lands on the *same* boundary. A beacon-joined modulo-64 count can instead commit late on its first epoch (spec §6.3.1):

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

With absolute synchronized revolution counts, the construction window is identical (K revolutions) on every board because construction can't begin before B+R — only then is the window's start common knowledge regardless of which copy each board heard. An effect that can't construct inside K revolutions trips `HS_CHECK` (fail-fast). All boards reseed `hs::random()` per effect build from `HS_PHANTASM_EFFECT_SEEDS[]` in `targets/Phantasm/phantasm_playlist.h` — `hs::stable_effect_seed(hs::effect_seed_key<name<CANVAS_W, CANVAS_H>>(#name))` using the identity selector in `core/platform/effect_identity.h`, with `hs::epoch_seed(effect index)` the fallback when no seed table is supplied. That seed is identical on every visit, so an entry replays the same stream each time it comes round, and the new instance is bit-identical across boards no matter what each board rendered — or whether it even existed — before the epoch.

**Output envelope.** Nothing on the strip cuts at the commit: the ISR scales every packed column by `effect_output_envelope` (`pov_sync_content.h`), a fade-through-clear driven from the synchronized `rev_in_effect` and the column index alone.  An entry fades up over its first two revolutions and back down over its last two, and once `rev_in_effect` reaches the entry's configured duration the envelope is **zero**.  That is exactly the revolution at which the master starts the EPOCH train, so the LEDs are already dark at **B** and stay dark through the announce phase as well as the construction window — R + K = 5 revolutions, ~0.63 s at 480 RPM — before the incoming entry fades up from its own frame 0.  The window itself, not the counter, is what forces the zero: the beacon carries only six bits of revolution, so a board that joined an entry longer than 64 revolutions counts congruent to the master rather than equal to it and would otherwise read full brightness at **B**.  It steps to black with everyone else and misses only the ramp into it.  Otherwise a pure function of already-synchronized state, it needs nothing extra on the wire and every board computes the same value for the same column; the pack loop tests for full brightness first, so a mid-effect revolution pays no multiply.

**Live-takeover join grid.** The epoch commit only aligns a *swap* between running boards; a board with nothing live yet — at boot, or after a reboot — reaches its first constructed effect whenever its own identity arrives, which is later downstream than on the master. It therefore does not take that effect live on arrival: it waits for a ZERO crossing whose revolution-in-effect is a multiple of `join_grid_revs` (4), marked by `TickActions::join_boundary`. The master sits on the same grid, so at boot every board goes live at the *same* crossing with aligned frame counters instead of the master leading by however long downstream identity took; a mid-show rejoin waits ≤ 4 revolutions, well inside the enforced 25-revolution rejoin budget. That budget is a `Config::valid()` relation, not prose: `rejoin_bound_revs()` is the widest beacon-to-beacon gap plus the join-grid wait. Beacons are suppressed for the whole commit window, so `commit_beacon_gap_revs()` computes each roster entry's gap across its commit from its duration; for the shipped roster that is at most 15 + EPOCH repeats 3 + construction window 2 + 1, and with join grid 4 the bound is 25 revolutions, ~3.1 s at 480 RPM. `valid()` rejects any configuration whose bound exceeds `rejoin_budget_revs`, including an entry one revolution past a multiple of the beacon period. The 16-revolution figure is the beacon cadence alone. The grid must divide 64 so a beacon's mod-64 revolution count lands on the master's grid. Where the epoch commit traps on a missing effect, the join is conditional: `EffectHandoff::joinable` also requires the pending build generation to be unconsumed and to match the one the wire advertises, so a visibility lag just joins one grid step later. No join boundary is marked while a commit is pending — the epoch path owns that swap.

**Concurrency & failure modes.** Two ISRs per board, **single-writer** by construction.  The sync-wire RISING ISR is a pure *publisher* — glitch filter, edge count, first-edge timestamp into a small mailbox, nothing else.  The flywheel ISR (waking ~8× per column) is the sole *consumer/owner* of all sync state: it claims terminated bursts, classifies, gates, snaps, flips, and runs epoch scheduling.  The hot path is ~7-of-8 wakes doing one cycle-counter read and a 64-bit position compute (≈1 % CPU at 600 MHz); column changes pack pixels, and a pending DMA submission can retry on later wakes.

| Event | Layer 1 (column) | Layer 2 (flip) | Layer 3 (content) |
|---|---|---|---|
| Masked-IRQ window | resumes at time-correct column | unaffected | unaffected |
| 1 dropped boundary symbol | coasts ≤ 1 rev (~0.01 col); re-snaps next | local crossing still flips | unaffected |
| 1 spurious / EMI edge | even count discarded, or gate rejects | identity dedup no-ops it | epoch refractory + gate guard it |
| Late-emitted symbol | master self-censors; residual gate-rejected | crossing flips on time regardless | unaffected |
| 1 board renders slow (drops a frame) | — | shows prior frame for 1 period | frame counters can remain offset until the next epoch; beacons resync only `rev_in_effect` |
| 1 dropped epoch symbol | — | — | R repeats; missed-all-R requires two confirming beacons (within 4 s; spec section 9.1) |
| Board reboots mid-show | re-acquires phase from next valid symbol | resumes flipping once LOCKED | rejoins correct effect via beacon, ≤ 25 revs (~3.1 s); dark until then |
| Sync wire severed (out of scope) | free-runs on own crystal; precesses ≥ 1 col in ~10–20 s | keeps flipping locally | holds last effect; slow drift, never an instant break |

The flywheel ISR maintains telemetry counters (symbols accepted / gate-rejected / discarded, beacons ok / rejected, index corrections, epochs refractory-ignored, lock transitions, flips, emissions censored / aborted, longest coast) that the foreground reports behind `hs::debug` at ≤ 1 Hz — so any degradation the protocol absorbs silently is still visible at a glance.


## Mathematical Kernels (`core/math/`)

The math headers provide coordinate transforms and scalar fields used by the
plot, scan, and pullback paths. They do not own a canvas or effect lifecycle.

| Headers | Surface |
|---|---|
| `core/math/3dmath.h`, `core/math/4dmath.h` | Vectors, matrices, complex arithmetic, fast scalar approximations, and four-dimensional rotations. |
| `core/math/geometry.h` | Umbrella for periodic, pixel-mapping, and spherical helpers. |
| `core/math/periodic.h` | Periodic coordinates and index wrapping. |
| `core/math/pixel_mapping.h` | Sphere/pixel coordinates using the calibrated LED-center latitude profile. |
| `core/math/display_geometry.h` | Display latitude profiles, row spacing, and pole-cap geometry. |
| `core/math/spherical.h` | Canonical axes, sphere generators, tangent frames and bases, and parallel transport. |
| `core/math/projections.h` | Bonne, Peirce quincuncial, Airocean, folded sinusoidal, and equirectangular sphere-to-plane kernels. |
| `core/math/stereographic.h`, `core/math/mobius.h`, `core/math/lenses.h` | Stereographic projections, fractional-linear transforms, and sphere-domain lens kernels. |
| `core/math/rotate.h`, `core/math/projection_patterns.h` | Angle wrapping, canvas-to-sphere projection, and shared projected-pattern coordinates. |
| `core/math/noise_field.h`, `core/math/spherical_field.h`, `core/math/spherical_harmonics.h` | Noise sampling, spherical fields, and harmonic evaluation. |
| `core/math/easing.h`, `core/math/waves.h`, `core/math/interpolate.h` | Scalar easing curves, periodic waves, and interpolation helpers. |

Projection results carry their coordinate and validity contracts in the headers.
The pullback policies in `core/render/pullback/` bind these kernels to effect
parameters and frame state.

## Spatial Queries (`core/spatial/`)

`kd_tree.h` provides arena-backed three-dimensional nearest-neighbor queries.
`KDTree::nearest_moved` answers the same query for points that have moved by a
bounded distance since the build, ranking them by a caller-supplied distance.
`reaction_graph.h` defines the Fibonacci reaction lattice and node lookup;
`reaction_graph.cpp` stores the generated neighbor table, its run-length
encoding (`neighbor_runs`/`neighbor_run_index`) and immutable float node positions.
The generator of record is `scripts/generate_reaction_graph.py`; its provenance
gate diffs the whole generated file. The analytic `ReactionGraph::node()`
remains the reference, with a native comparison across all 7,680 positions.
Reaction-diffusion effects use these neighbors without rebuilding adjacency.

| Interface | Description |
|---|---|
| `ReactionGraph::node_positions` | Generated immutable float positions for the Fibonacci lattice. |
| `ReactionGraph::neighbor_runs` / `ReactionGraph::neighbor_run_index` | Run-length encoding of the generated neighbor table. |
| `ReactionGraph::CubemapLUT` | Arena-backed 64 × 64 cube-face lookup of nearest-node seeds; `BYTES` is 48 KiB. |

`ReactionDiffusionBase` (`effects/common/reaction_diffusion.h`) can either
build canonical positions in its persistent arena or bind the flash array
through `init_lattice<true>()`. Its arena-budget check takes the matching `FLASH_NODES` flag; with temporal sampling disabled the oriented positions remain a
per-frame scratch allocation. GS uses flash positions, releasing 92,160 bytes of
persistent storage for its palette cache. BZ retains the arena-backed default.

With temporal sampling enabled, BZ prepares the frame's lattice basis, a
`Scan::ObjectBasisTable` of its view orientation, and shades from the
unrotated lattice through `Scan::Shader::walk_object_grid`, so each column sees
the lattice at its sweep time without an oriented copy. It advances and shades
through `Spacetime::simulate_in_slabs`, one sweep slab per substep, and the
basis serves every slab. With temporal sampling disabled it runs both substeps
and shades the lattice oriented at the current endpoint once.

With temporal sampling enabled, GS shades from its flash lattice through the
same basis in six slabs, one per substep. Its float field stays in scratch
across the slabs; each slab quantizes B into the persistent Q16 field through
the frame's dissolve mask, marking the nodes at the render floor as it goes,
dilates the marks into its two-ring cull bits, and blends the colour
noise between the previous and current frame's tables at the slab's end,
keeping the previous table in a second persistent noise buffer. The pigment
step and the lifecycle outcome are settled before the first substep: a frame
that ends a dissolve stays black, and a reseed frame shades the seeded field
whole. The Q16 state ends every frame as the frame-wide order, used with
temporal sampling disabled, leaves it.

GS reuses the lattice cubemap projection for noise sampling, reversing the
horizontal coordinate on the X faces and the vertical coordinate on the Y faces
to match the noise cubemap orientation.

GS caches 16 concentration coordinates by 15 shared-noise coordinates for each
seed palette in float linear RGB, then interpolates those colors. Cache validity follows palette
reseeding and changes to Hue Shift or Shimmer. Hue Shift magnitudes above 0.25,
Shimmer above 0.4, and the unmodified palette use the procedural color path.
Four antialiasing samples determine coverage from signed kernel mass relative to
the visible B threshold. Two diagonal samples supply reciprocal total weights;
their mean estimates the other two reciprocals for color concentration. A
conservative center-support certificate falls back to the complete four-weight
kernel when the refined center may lie outside kernel support, and a signed-mass
guard does so near numerical coverage boundaries.
The refined center supplies its packed two-palette pigment mixture; each
contributing palette is evaluated once at the shared concentration. Bilinear
RGB weights accumulate directly into the pixel, with one final quantization.
The fidelity oracle bounds this color approximation against four procedural
color samples through growth and reseeding, and a dense aggregate pigment
reference bounds nearest-node color selection separately. The scalar nearest-
pigment reference bounds arithmetic rounding; the stencil oracle compares both
geometries under the nearest-pigment color rule.
