# Spherical perspective and independent patterns

**Status: PROPOSED (2026-09-26).** This document specifies a generalization of
[HyperLattice](../../effects/HyperLattice.h); it does not describe shipped APIs.
Implementation is a separate change. Dimensional rift is removed from the
proposed design.

## 1. Decisions

Spherical perspective is a camera and rendering facility. A pattern defines
geometry in an ambient space. A sampling domain embeds the camera's rays into
that space. A tracing backend finds contributions along those rays. A renderer
shades and composites those contributions onto the display sphere.

These are independent choices, subject to declared compatibility:

| Choice | Responsibility | Initial values |
| --- | --- | --- |
| Pattern | Geometry, repetition, feature identity, material identity | Cubic wire lattice; additional patterns below |
| Sampling domain | Ambient dimension and camera embedding | Spatial 3D; Slice through 4D |
| Tracing backend | Intersection search, ordering, termination | Ordered analytic events; bounded distance-field marching |
| Appearance | Palette, fog, opacity, feature coloring | Depth or pattern-feature coloring |

3D and 4D may be called **sampling modes** in the UI, but the architectural
name is **sampling domain**. They do not select a different spherical
projection or tracing algorithm. The proposed labels are **3D** and **4D slice**.
The latter requires a pattern with an actual 4D definition.

There is no dimensional-rift mode, blend coefficient, or interpolation between
3D and 4D distance metrics. Pattern geometry, sampling domain, and backend are
discrete choices. Ordinary parameters may interpolate within a compatible
configuration; changes between configurations use the effect's transition
policy rather than blending unrelated geometry queries.

## 2. Camera and dimensional semantics

For a unit display direction n, define the ambient ray by:

```text
d = E n
p(t) = c + (r + t) d,     near <= t <= far
```

Here c is the camera center in ambient world units, r is the radial start
offset, and t is distance beyond that start. E has three orthonormal columns.
Consequently d is unit length, and t has the same world-distance meaning in
both sampling domains. Rotation and translation are frame-prepared values.

In 3D, E is a 3-by-3 rotation and c is a 3D position. In 4D, E is a 4-by-3
orthonormal embedding and c is a 4D position. All rays lie in the affine
three-dimensional slice c + image(E). Rotations involving the fourth axis
change that slice's orientation; translation perpendicular to it changes its
offset. Other rotations turn the view within the slice.

The radial offset skips the interior of a sphere in the viewing slice. Rays
remain collinear with rays from c. It is neither lens magnification nor a
four-dimensional perspective divide. At r = 0, sampling starts at the center.
Fog and near/far controls use t; a pattern query receives p(t). This convention
matches the current effect's distance beyond the display sphere.

The 4D domain samples the intersection of ambient geometry with the viewing
slice. It does not project all 4D objects into 3D and does not integrate along
the slice normal. Such operations would require separately specified modes.

### 2.1 Geometry dimension matters

A 4D pattern is not automatically the 3D pattern with a fourth coordinate
ignored. Every pattern declares supported ambient dimensions and its geometry
in each. A cubic lattice family can define edges of cubic cells in 3D and
edges of hypercubic cells in 4D, with a positive wire radius in ambient units.

An infinitesimally thin 4D edge generally meets a generic 3D slice only at
isolated points, if at all. Finite wire radius produces visible sections;
continuous 3D-looking wireframes are not promised. Other choices, such as
thickened 4D faces or a 4D hypersurface, produce different slice structures
and must be named as different geometry, not silently substituted for edges.

A pattern defined only in 3D remains unavailable in the 4D domain unless it
explicitly provides a 4D extension. Extruding it unchanged along the fourth
axis is one possible extension, but must be labeled as such.

### 2.2 Lattices as volume distance fields

Both wire lattices admit a volume distance query. For ambient dimension N,
let q[i] be the distance of coordinate i to the nearest integer, in cell units:

```text
q[i] = abs(p[i] - round(p[i]))
edge_distance = sqrt(sum(q[i]^2) - max(q[i]^2))
field = cell_size * edge_distance - wire_radius_world
```

The largest component is discarded because that coordinate is free along the
nearest axis-aligned edge. N = 3 gives cubic edges; N = 4 gives hypercubic
edges. The unsigned edge distance is exact. Subtracting the radius defines
the union of thickened edges and gives exact exterior distance; inside
overlapping tubes it need not be exact distance to the union boundary.
An inside-start or exit-search algorithm must respect that distinction.

The current HyperLattice evaluates related distances only at selected grid
crossings. The same geometry can instead be supplied to a volume marcher.
These backends need not produce identical coverage or visibility: the current
crossing approximation and its finite layer limit are explicit compatibility
choices, not part of the lattice's mathematical definition.

## 3. Integration with existing rendering code

Build on the existing volume and pullback facilities rather than introducing
a parallel SDF library or another framebuffer traversal:

| Existing component | Integration decision |
| --- | --- |
| [Volume shapes](../../core/render/sdf/volume.h) | Reuse `SDF::Torus`, `SDF::WarpedVolume`, their distance queries, warp bounds, and precondition checks |
| [Volume renderer](../../core/render/scan/volume.h) | Extract a reusable ray kernel from `Scan::Volume::trace_closest`; retain `Scan::Volume::draw` as the existing orthographic caller |
| `Scan::TransformedVolume` | Preserve placement semantics; reuse its world-to-local ray transform for placed 3D shapes |
| [HyperLattice](../../effects/HyperLattice.h) | Retain sphere-sample acquisition and prepared pullback-stage integration; delegate tracing and pattern queries |
| [Layer compositor](../../core/color/layer_composite.h) | Reuse `LayerComposite` for ordered contributions |
| [Raymarch effect](../../effects/Raymarch.h) | Keep its placement, shading, and orthographic behavior as a regression consumer |

The existing draw path uses parallel rays for each placed volume. Its
`trace_closest` helper already accepts a ray origin and direction, but its
termination is coupled to a bounding radius and it reports closest clearance
rather than a complete hit/status record. Extract the shared numerical kernel
with an explicit ray interval, footprint policy, and result/status contract.
Keep a compatibility wrapper for the existing caller. Do not copy its loop
into the generalized effect or route spherical rays through its orthographic
draw method. Preserve overrelaxation fallback and first-silhouette ownership.

Existing shape queries returning float distance remain valid. The current
`SDF::VolumeShape` concept also requires normal and fragment population for
warping; retain that public contract. A trace-only capability can require just
distance without forcing new shapes to invent unused fragment registers.
Appearance adapters reuse existing normals and fragment helpers when useful.
Preserve `WarpedVolume::precision` and trace preconditions: cheap clearance
bounds must not masquerade as accurate near-surface samples.

### 3.1 Use 3D ray coordinates for both sampling domains

The shared marcher can continue using `math::Vector`. In camera-relative slice
coordinates x, a prepared pattern adapter evaluates:

```text
query3(x) = pattern3.distance(c3 + E3 x)
query4(x) = pattern4.distance(c4 + E4 x)
ray_in_slice(t) = (r + t) n
```

Additional pattern-local transforms are applied inside the adapter with metric
correction. The 4D definition uses `math::Vec4` and existing 4D math internally.
The marcher receives a 3D callable distance query in either case; no second
4D march loop or mandatory conversion of existing volume shapes to Vec4 is
needed. An ambient distance bound restricted to the slice stays conservative,
although it is generally not the exact SDF of the resulting 3D intersection.

For slice-surface lighting, transform the ambient gradient by the transpose
of E and normalize in 3D. Handle a vanishing projected gradient explicitly;
do not pass a 4D normal to a 3D lighting function. Finite differences of the
restricted query are an optional alternative with separately measured cost.

Here an SDF volume means solid geometry described throughout a spatial domain.
It is already compatible with spherical perspective. It does not imply fog
or participating-medium integration, which uses a different accumulation rule.

## 4. Component boundaries

The integration follows the existing
[pullback pipeline](pullback_pipeline_spec.md): a prepared stage accepts a
sphere sample and returns composited color. It delegates to the following
components rather than owning pattern-specific branches.

```text
Sphere sample -> Camera/domain -> Ray in slice coordinates
                                      |
Pattern definition -> Prepared query adapter -> Tracing backend
                                                   |
                                      Ordered geometric contributions
                                                   |
                                      Appearance -> Compositor -> Color
```

### 4.1 Pattern definition

A pattern owns its topology, repetition rules, geometry parameters, supported
dimensions, stable feature/material identifiers, and local coordinate mapping.
Its preparation step may build bounded lookup data or acceleration data.
Its geometry is independent of camera pose, palette, fog, framebuffer layout,
and the ray termination policy.

Pattern-local units may be convenient for repetition. The query adapter must
convert results to the shared world-distance convention. Anisotropic scaling
or shear requires correct metric conversion or a conservative distance bound;
normalizing a transformed direction must not silently change the returned t.

### 4.2 Query adapter

An adapter exposes the pattern's geometry to a supported backend. Analytic
intersection formulas necessarily depend on geometry; they belong here, not
in the generic ray loop. Adapters return geometry and material data, never
final display colors or fogged opacity.

A pattern may support either or both of these initial query capabilities:

| Capability | Adapter supplies | Backend supplies |
| --- | --- | --- |
| Analytic events | Bounded cursor initialization and next candidate intersections or coverage events | Event ordering, grouping, range checks, budgets, consumption |
| Distance field | Field evaluation, conservative distance bound, material/feature lookup and optional gradient | Step selection, root refinement, hit continuation, budgets |

Capability selection occurs during frame preparation. The renderer must not
require every pattern to implement a signed distance function, and the analytic
backend must not assume integer planes, orthogonal axes, or cubic cells.

### 4.3 Tracing backend

The backend consumes a ray in slice coordinates, prepared query adapter,
footprint, and trace limits. Analytic adapters map that ray to ambient geometry
when needed. It emits contributions in nondecreasing t. It owns near/far
clipping, traversal state, progress checks, event grouping, and early exit.
Pattern adapters do not own a second hidden compositing or traversal loop.

Analytic adapters may expose several monotone event streams, for example one
per plane family. They must advance strictly after consuming an event and
signal exhaustion. The backend merges the streams; a fixed-capacity heap or
small linear scan is an implementation choice. Capacity is declared per
compiled adapter, not hard-coded to three or four coordinate axes.

Distance-field adapters must document their guarantees. An exact signed
distance or conservative bound permits safe marching. An arbitrary implicit
value is not a distance: its adapter must provide a valid bound, such as a
Lipschitz bound, or select a separately documented approximate search. Missing
thin surfaces must not be presented as an exact intersection guarantee.

For 4D slicing, ambient distance to ambient geometry is a conservative bound
on distance to an intersection reachable within the slice. It can support
safe progress but may converge slowly or approach geometry outside the slice.
Hit acceptance must verify the defined geometry at the sampled position;
exhausting the step budget is not a hit.

### 4.4 Contribution contract

Each contribution contains finite world-distance t, geometric coverage in
[0, 1], material identity, and optional feature identity and ambient normal.
It also declares whether it is an exact surface event or an approximate
coverage event. Pattern-local information needed for material evaluation may
be carried in a bounded payload. No payload contains owning allocations.

Events carry a merge identity. Coincident duplicate reports of the same
feature are coalesced; distinct surfaces or materials are not merged merely
because their distances are close. A lattice adapter may explicitly identify
a junction as one coverage layer. Equal-distance ordering is deterministic.

An angular pixel footprint travels with the ray. At t its world footprint is
derived from r + t, then transformed conservatively into pattern coordinates.
Geometric filtering uses this footprint; appearance applies distance fading
afterward. This is distinct from the legacy lattice's particular filter
approximation, which may be retained in its compatibility adapter initially.

### 4.5 Appearance and compositing

Appearance maps geometric contributions to color and opacity using depth,
material, or feature identity. It owns palettes, near fading, fog, and any
lighting model. A normal is optional; unlit patterns need not compute one.
Axis coloring becomes an adapter-provided feature palette mapping rather than
an assumption that every pattern has four axes.

The existing front-to-back layer compositing semantics remain shared. The
compositor stops at the configured opacity threshold. A surface backend must
define entry/exit handling so a solid object is not shaded twice by accident;
an explicitly translucent two-sided material may request both crossings.

Continuous participating media need interval integration and are outside the
initial contribution contract. They require a future backend and contribution
extension, not fake surface events at arbitrary march steps.

## 5. Configuration and execution

Configuration has separate pattern parameters, sampling-domain pose, camera
range, appearance parameters, and trace-quality limits. A registry declares
supported pattern/domain/backend combinations and their bounded resources.
Unsupported combinations fail validation before rendering; the UI presents
only supported combinations. Backend selection defaults from capabilities
and quality requirements and need not be an ordinary user control.

Frame preparation resolves the concrete combination once. Firmware uses
templates or equivalent static dispatch inside the pixel loop; no per-pixel
registry lookup, heap allocation, or virtual geometry dispatch is required.
Compile only admitted combinations to control flash and ITCM growth.
Prepared state is immutable during rendering; cursors and compositors are
per-ray state. Existing segmented rendering must see one consistent frame.

Trace limits separately bound analytic candidates, march steps, refinement
iterations, emitted layers, and geometric range. Every backend terminates on
exhaustion. Exhaustion returns accumulated valid contributions and a diagnostic
status; it does not invent a hit or turn an unfinished search into a confirmed
miss. Diagnostics are aggregated for tests/profiling, not shown as extra color.

Current lattice shells mean a plane-count limit per axis, not spherical shells
or a universal depth control. Keep that behavior as a lattice-specific event
limit during migration. New patterns use their own declared traversal limits
under the renderer's global work budget. Do not expose a generic shells knob.

## 6. Pattern scope and initial implementations

Arbitrary patterns means an extensible geometry contract, not a promise to
render every mathematical field within a fixed Teensy frame budget. A new
pattern adds a definition and adapters without modifying the camera,
compositor, or unrelated backend implementations.

The first implementation should demonstrate all of these cases:

1. **Cubic wire lattice**, with explicit 3D and 4D definitions and an analytic
   event adapter preserving the current bounded plane-crossing renderer.
   That renderer is a coverage approximation, not exact ray/cylinder tracing.
2. **Triangular-prism wire framework**, initially 3D only, using nonorthogonal
   plane families and analytic events. This verifies that the analytic loop
   contains no cubic-axis or three-stream assumptions.
3. **Existing volume SDF**, using `SDF::Torus` and then `SDF::WarpedVolume`
   through the extracted volume kernel with spherical rays. This verifies
   direct reuse of existing shapes without requiring plane crossings.
4. **Lattice volume query**, using the formula above in 3D and 4D-slice
   adapters through the same marcher. Compare geometry with the analytic
   adapter while documenting expected coverage differences.

Decorated cells, repeated analytic objects, and other implicit surfaces can
then reuse these boundaries. A 4D extension of a surface is a separate
authored definition with declared slice semantics. Adding further ambient
dimensions, arbitrary runtime shader code, 4D-to-3D projection, volume
integration, and a general scene graph are outside this initial scope.

## 7. Migration and compatibility

1. Extract camera/domain preparation, appearance, and layer consumption from
   HyperLattice behind the shared stage. Keep its current two endpoint domains
   and preserve their output before adding other patterns.
2. Move cubic/hypercubic geometry and plane candidate production into a lattice
   adapter. Keep its grouping, filtering, and horizon behavior explicit as
   compatibility behavior, with reference tests rather than assumptions that
   it is an exact surface renderer.
3. Remove dimensional rift from enum metadata, presets, transitions, parameter
   validation, tests, and documentation when implementation lands. Preserve
   stable serialized identities for surviving modes. If legacy numeric IDs
   are retained, reserve the removed value rather than renumbering 4D into it.
   Decode legacy rift selections as 3D with ordinary 3D geometry; never retain
   an invisible metric-blending path. The decoder must explicitly recognize
   the former value; unrelated invalid values still fail validation.
4. Extract and validate the existing volume kernel with Raymarch unchanged,
   then connect existing SDF shapes and the lattice queries to spherical rays.
   Introduce the noncubic analytic pattern and evaluate quality and device
   cost independently. Gyroids and other new fields follow after this reuse
   is demonstrated, with bounds specified for each field.

Retain the public HyperLattice effect identifier during the initial migration
so playlists and saved effect selection keep working. Its controls can expose
the generalized pattern/domain choices. A public rename to SphericalPerspective
is a separate compatibility decision, not necessary to establish the boundary.

## 8. Acceptance criteria

- Camera tests establish unit ray directions, correct radial starts, common
  t units, and equivalence of 3D sampling to a fourth-axis extrusion viewed
  in the corresponding 4D slice.
- Slice tests establish that rotations involving the fourth axis change the
  embedding and that no perspective divide or slice-normal integration occurs.
- Adapter/backend tests cover empty geometry, parallel and tangent rays,
  coincident features, negative cell coordinates, near/far boundaries,
  finite payloads, strict cursor progress, and deterministic event ordering.
- Marching tests use geometry with known intersections to check conservative
  progress, roots, continuation after a hit, and budget-exhaustion status.
- Existing Raymarch captures and volume tests remain valid after kernel
  extraction. A torus and warped torus render through both camera paths using
  the same shape types and numerical kernel. A 4D lattice slice uses that
  kernel through a 3D query adapter, without a second marching implementation.
- Reference captures preserve the surviving lattice modes during extraction;
  intentional filtering changes receive separate comparisons. Rift disappears
  from new controls and legacy selections follow the documented migration.
- Adding the triangular framework and implicit surface requires no pattern
  branches in camera, appearance infrastructure, or compositor. Generic
  backend code depends only on its query capability contract.
- Relevant native tests, firmware builds, and size/layout gates pass. Record
  RAM1 code, RAM1 variables, FLASH data, and changes from baseline. Profile
  admitted configurations on the Teensy with fixed provenance and report
  traversal exhaustion alongside frame cost. No performance equivalence is
  assumed between analytic and marching backends.

This spec-only change requires documentation validation; firmware validation
and captures belong to the implementation changes above.
