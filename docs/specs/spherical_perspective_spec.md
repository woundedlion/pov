# Spherical perspective and independent patterns

**Status: IMPLEMENTED architecture, revision 5 (2026-09-27).** The ray core,
legacy lattice migration, noncubic framework, volume/lattice query adapters,
and periodic-surface experiments are implemented. Firmware admission is separate
from architectural availability; experimental configurations are available in
the simulator and explicitly opted-in device builds, while standard firmware
keeps only the admitted analytic presets. See the
implementation validation and admission report (supporting artifact removed)
for measured limits. Dimensional rift is removed.

The experimental preset guide (supporting artifact removed)
records the preview presets, build switch, and bounded failure behavior. Their
availability does not change their experimental admission status.

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
| Appearance | Palette, fog, opacity | Depth coloring |

3D and 4D may be called **sampling modes** in the UI, but the architectural
name is **sampling domain**. They do not select a different spherical
projection or tracing algorithm. The proposed labels are **3D** and **4D slice**.
The latter requires a pattern with an actual 4D definition.

There is no dimensional-rift mode, blend coefficient, or interpolation between
3D and 4D distance metrics. Pattern geometry, sampling domain, and backend are
discrete choices. Ordinary parameters may interpolate within a compatible
configuration; changes between configurations switch one complete parameter
block at the transition boundary specified in section 5.

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
All inputs are finite, r >= 0, and 0 <= near < far. Scale belongs to the
pattern transform, not E; degenerate transforms fail preparation.

In 3D, E is a 3-by-3 rotation and c is a 3D position. In 4D, E is a 4-by-3
orthonormal embedding and c is a 4D position. All rays lie in the affine
three-dimensional slice c + image(E). Rotations involving the fourth axis
change that slice's orientation; translation perpendicular to it changes its
offset. Other rotations turn the view within the slice.

The radial offset skips the interior of a sphere in the viewing slice. Rays
remain collinear with rays from c. It is neither lens magnification nor a
four-dimensional perspective divide. At r = near = 0, sampling starts at the center.
Fog and near/far controls use t; a pattern query receives p(t). This convention
matches the current effect's distance beyond the display sphere.

The 4D domain samples the intersection of ambient geometry with the viewing
slice. It does not project all 4D objects into 3D and does not integrate along
the slice normal. Such operations would require separately specified modes.

Camera animation remains effect-owned. Generic camera position uses world
units and is never reduced modulo a cell. HyperLattice's current origin wraps
in cell coordinates every frame; retain that only in the legacy lattice
animation adapter and convert its pose into the world convention. Other
periodic patterns may fold coordinates by their declared periods inside their
queries. A placed torus has no such wrapping. Changing cell size must not
silently reinterpret a nonlattice camera position.

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
first map the ambient world point p into lattice cell coordinates u. For a
rigid lattice pose with origin a, rotation R, and positive uniform cell size s:

```text
u = transpose(R) (p - a) / s
q[i] = abs(u[i] - round(u[i]))
edge_distance = sqrt(sum(q[i]^2) - max(q[i]^2))
field = s * edge_distance - wire_radius_world
```

The largest component is discarded because that coordinate is free along the
nearest axis-aligned edge. N = 3 gives cubic edges; N = 4 gives hypercubic
edges. The unsigned edge distance is exact. Subtracting the radius defines
the union of thickened edges and gives exact exterior distance; inside
overlapping tubes it need not be exact distance to the union boundary.
An inside-start or exit-search algorithm must respect that distinction.
In floating-point code, sum the retained squared components directly instead
of subtracting the largest from their total, to avoid cancellation near edges.

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
The existing `probe_occluder` and `volume_edge_coverage` behavior also remains
under the orthographic compatibility policy. Its coarse probe is not a safe
general-purpose exit finder and is not exposed as generic multi-hit tracing.

Two policies share the extracted stepping machinery: **legacy closest
approach** preserves Raymarch's output, minimum-step floor, and approximate
coverage; **surface search** uses declared query guarantees and reports
unresolved progress instead of forcing a step across unknown geometry.
Extraction does not turn the legacy closest-approach test into a certified
intersection test. Keep existing profile scope attribution and optimization
boundaries measurable during extraction.

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
These equations specify semantics, not mandatory per-step matrix multiplies.
Compose rigid transforms per frame/ray and evaluate a prepared query along
`origin + t * direction` in native 3D or 4D coordinates when cheaper. The
scalar progress loop is shared; no second 4D marching algorithm or mandatory
conversion of existing shapes to Vec4 is needed. A 3D callable wrapper is
also valid. An ambient distance bound restricted to the slice stays conservative,
although it is generally not the exact SDF of the resulting 3D intersection.

For slice-surface lighting, transform the ambient gradient by the transpose
of E and normalize in 3D. Handle a vanishing projected gradient explicitly;
do not pass a 4D normal to a 3D lighting function. Finite differences of the
restricted query are an optional alternative with separately measured cost.

Here an SDF volume means solid geometry described throughout a spatial domain.
It is already compatible with spherical perspective. It does not imply fog
or participating-medium integration, which uses a different accumulation rule.

### 3.2 File layout and ownership

All reusable geometry and rendering machinery belongs in the core engine,
including concrete reusable patterns and their query adapters. A helper does
not remain effect-private merely because HyperLattice is its first consumer.
The effect owns the choice of geometry and settings, not their implementation.

The following paths are the implementation layout. The header-based templates
do not introduce a separate library build or runtime service.

| Path | State | Responsibility |
| --- | --- | --- |
| core/render/ray.h | Implemented | Public umbrella for reusable ray rendering; namespace `Raycast` |
| core/render/ray/contract.h | Implemented | Ray interval, footprint, trace limits/status, geometric contribution, and query capability contracts; no effect types or framebuffer traversal |
| core/render/ray/camera.h | Implemented | Spherical ray construction, validated 3D/4D slice embedding, world-distance conventions, and prepared camera transforms |
| core/render/ray/query.h | Implemented | Generic volume-query adaptation, placement/domain composition, prepared per-ray evaluation, and projected normals; geometry-specific formulas stay with their pattern |
| core/render/ray/march.h | Implemented | Shared scalar stepping/refinement kernel, first-boundary policy, and extracted legacy closest-approach policy |
| core/render/ray/events.h | Implemented | Bounded merge of analytic candidate streams, event grouping, ordered emission, and traversal budgets |
| core/render/ray/shade.h | Implemented | Reusable depth/feature appearance policies, fog/near fade, verified subray filtering, and contribution consumption through `LayerComposite`; receives palettes and settings from callers |
| core/render/pullback/ray.h | Implemented | Reusable prepared `SphereSample -> Color4` stage binding camera, query, backend, and appearance; no new canonical carriers |
| core/render/sdf/lattice.h | Implemented | Prepared cubic/hypercubic lattice crossings, analytic plane-event adapters, feature identities and shading |
| core/render/sdf/framework.h | Implemented | Triangular-prism and octet 3D/4D framework geometry with analytic query/event adapters |
| core/render/sdf/octet_trace.h | Implemented | Octet cell traversal and bounded analytic events |
| core/render/sdf/lattice_trace.h | Implemented | Cubic and hypercubic lattice cell traversal |
| core/render/sdf/cellular_wire.h | Implemented | Diamond, hexagonal and rhombic cellular wires |
| core/render/sdf/affine_lattice.h | Implemented | Affine lattice geometry and traversal |
| core/render/sdf/periodic_shells.h | Implemented | Periodic shell geometry and candidate traversal |
| core/render/sdf/lattice_field.h | Implemented | WireLattice world-unit distance and gradient queries |
| core/render/sdf/periodic_surface.h | Implemented | Cosine and gyroid level-set definitions, period/isovalue parameters, bounds, gradients, and surface-query adapters |
| core/render/sdf/volume.h | Existing | Existing torus, warps, and warped-volume geometry; retain their public shape contracts |
| core/render/scan/volume.h | Existing | Orthographic scan/cull/draw front end, `TransformedVolume` compatibility API, and legacy occluder-probe/plot behavior; delegates extracted stepping to the ray core |
| core/render/scan/shader.h | Existing | Sphere-sample traversal and cached scan behavior; does not learn about pattern choices |
| core/color/layer_composite.h | Existing | Existing straight-alpha front-to-back compositor; no duplicate accumulator in the new subsystem |
| core/math/3dmath.h and core/math/4dmath.h | Existing | Vector, matrix, and rotation algebra reused by camera/query preparation; any generally useful algebra extension stays here |
| effects/HyperLattice.h | Existing | Effect identity, registered parameters, admitted tuple table, presets, choreography, camera/object animation, palette selection, frame state, and call into the shared stage |
| effects/Raymarch.h | Existing | Existing placement animation, effect-specific material choices, and orthographic draw invocation |

Reusable appearance operations go in the ray core; authored palette recipes
continue using the existing color facilities. Generic repetition/transformation
code goes in the query layer or existing math helpers. A pattern-specific
distance bound, intersection formula, or feature-color coordinate belongs
with its core geometry adapter. Legacy visual compatibility alone is not a
reason to keep reusable lattice geometry in the effect header. The effect's
legacy wrapped flight animation remains effect-owned because it chooses a
motion path rather than defining geometry.

Keep the existing `Scan::TransformedVolume` API. Generic placement/domain
composition is implemented in the ray query layer and may be delegated to by
that wrapper; the new core must not include the scan renderer to obtain a
transform. Likewise, the shared marcher cannot include the effect header to
obtain a distance evaluator or a tracing constant.

### 3.3 Dependency direction and extraction boundaries

The lower-level contracts depend only on existing math/platform/value types.
Camera and generic query preparation depend on those contracts. Core pattern
headers implement query capabilities and may include the contracts, but not
the tracing backends. March/event backends depend on capabilities, not a
concrete pattern catalog. Appearance depends on contribution contracts and
existing color/shading primitives. The prepared pullback stage composes these
helpers using the existing pullback contract.

Effects select concrete core patterns and backends and supply frame-owned
state. Neither the ray core nor any pattern header includes files under
effects, workbench, or targets, knows the HyperLattice configuration enum, or
consults the effect registry. Low-level core geometry must compile and be
testable without instantiating an effect or a canvas. New includes must not
create a cycle through the existing SDF or scan umbrellas; leaf headers use
leaf dependencies, and consumers include the narrowest public entry needed.

Extract in dependency order: contracts and queries; camera and shared stepping;
core pattern adapters and event traversal; shared appearance and pullback
stage; effect wiring. Move reusable numerical logic out of
`HyperLatticeDetail`, rather than retaining effect-private implementations
behind core forwarding wrappers. The orthographic wrapper and spherical
stage must instantiate the same march implementation. Unused template pattern
inventory must not emit executable code, lookup data, or static registration.

### 3.4 Test and documentation locations

| Path | State | Responsibility |
| --- | --- | --- |
| tests/test_ray.h | Implemented | Camera/domain, query capabilities, shared marcher, analytic-stream merger, filtering, statuses, and synthetic-geometry tests without effect dependencies |
| tests/test_sdf_patterns.h | Implemented | Cubic/hypercubic fields, framework geometry, periodic surfaces, bounds, scale/translation invariance, and native reference comparisons |
| tests/test_hyper_lattice.h | Existing | Preset/configuration/schema behavior, choreography, frame preparation, and effect-level compatibility captures |
| tests/test_scan.h | Existing | Orthographic volume regressions, including first-graze and background-graze ownership after extraction |
| tests/test_pullback.h | Existing | New shared stage's prepared-state and `SphereSample -> Color4` contract integration |
| docs/specs/spherical_perspective_spec.md | Existing | This subsystem's architecture and numerical contracts |
| docs/subsystems.md | Existing | Shipped ray/query facilities and their public entry points, updated when implementation lands |
| docs/effects.md | Existing | User-facing pattern/configuration behavior, updated when implementation lands |
| docs/profiles/ | Existing | Device reports and source/build provenance through the current profiling workflow |

Register new suites through the existing native test harness and its coverage
gates. Core tests must directly instantiate core helpers; effect white-box
access is reserved for effect-owned behavior. Update tracked file maps and
documentation references when the proposed files are created, not by adding
empty placeholders in this spec change.

## 4. Component boundaries

The integration follows the existing
[ranked pullback pipeline](pullback_stage_families_spec.md): a prepared
`Stage::Contract` accepts `Pullback::SphereSample` and returns `Color4`, as
HyperLattice does today. Rays, queries, and contributions are internal helper
types, not new canonical pullback carriers or stage families. This proposal
does not expand the workbench interpreter or its runtime schema.
`Scan::Shader::draw_cached` continues to own HyperLattice's traversal. Shared
helpers belong beside the existing volume rendering/query facilities; effect
controls, animation, and the admitted-configuration table stay effect-owned.
Conceptual responsibilities need not become separately allocated objects or
require a new general-purpose registry.

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
| Distance field | Membership, declared clearance guarantees, surface validation, material/feature lookup and optional gradient | Step selection, bounded refinement, first-surface search, budgets |

Capability selection occurs during frame preparation. The renderer must not
require every pattern to implement a signed distance function, and the analytic
backend must not assume integer planes, orthogonal axes, or cubic cells.

Distance adapters declare these guarantees independently:

- Sign or membership evaluation, with a stated zero set.
- Conservative nonnegative exterior clearance in world units.
- Optional conservative interior clearance; it is never inferred by taking
  the absolute value of a signed result.
- Surface verification/refinement and optional slice-space filtering support.

The existing float-returning shape methods stay unchanged; an adapter supplies
these declarations and any extra operations. `WarpedVolume` currently corrects
positive distances but leaves negative raw values uncorrected, so its initial
adapter is exterior-only. The lattice field is 1-Lipschitz in world units:
its absolute value supplies a conservative boundary clearance on either side,
even though its interior value is not an exact boundary distance.

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

Distance-field adapters must document their guarantees. A conservative bound
for the current side permits safe marching. An arbitrary implicit
value is not a distance: its adapter must provide a valid bound, such as a
Lipschitz bound, or select a separately documented approximate search. Missing
thin surfaces must not be presented as an exact intersection guarantee.

For 4D slicing, ambient distance to ambient geometry is a conservative bound
on distance to an intersection reachable within the slice. It can support
safe progress but may converge slowly or approach geometry outside the slice.
Hit acceptance must verify the defined geometry at the sampled position;
exhausting the step budget is not a hit.

Surface verification uses an adapter-provided intersection certificate or a
root bracket with a verified membership change and bounded refinement.
A sign bracket proves a crossing exists, not that it is the first one:
emitting it also requires the preceding ray interval to be cleared or the
adapter to isolate its earliest root. Proximity alone is an approximate
candidate. If progress or root ordering cannot be resolved within the budget,
return `UNRESOLVED`; tangency without a certificate is not a confirmed miss.
An overrelaxed step is accepted only under the query's clearance guarantees.
The conservative policy never substitutes a minimum step larger than its
certified clearance to escape a stall.

Refinement evaluations may probe beyond that clearance inside a bounded
candidate interval without declaring the intervening interval empty. A
membership-changing bracket no wider than the position tolerance, with the
preceding interval safely cleared, locates the first boundary to that tolerance.
This permits existing torus queries without an analytic quartic solver. Failed
verification returns `UNRESOLVED`; a probe is never a forced traversal step.

Numerical position/root tolerances are positive and scale-aware, independent
of the angular footprint (which is zero at the camera center). Adapters state
their error allowance. Nonfinite query results terminate with `INVALID_QUERY`.
Mathematical certification here is always subject to the declared numerical
tolerance; results do not claim floating-point exactness.

### 4.4 Initial surface policy and trace outcome

The initial SDF backend emits the **first forward boundary**, with opaque
surface material. Starting outside seeks entry. Starting inside seeks exit
only when the adapter declares interior clearance and surface verification;
otherwise return `UNSUPPORTED_START` without fabricating a surface. Lattice
adapters support both sides; the initial warped-volume demonstration keeps
the entire radial start surface outside its bounds. An explicitly verified
boundary at the start is emitted once and terminates that ray.

There is no generic SDF transparency, exit skipping, or repeated-surface
continuation in the first implementation. Analytic event adapters retain
layered rendering. SDF multi-hit rendering is a later capability requiring
certified departure from the consumed boundary and tests for thin adjacent
components; restarting at `t + epsilon` is insufficient. This limitation is
visible in configuration admission, not hidden inside individual patterns.

Every trace returns bounded work counters and a terminal status: `SURFACE`,
`RANGE_COMPLETE`, `SATURATED`, `BUDGET_EXHAUSTED`, `UNRESOLVED`,
`UNSUPPORTED_START`, or `INVALID_QUERY`. `RANGE_COMPLETE` means the chosen
backend completed its declared search, not that approximate lattice events
prove absence of geometry between crossings. Partial valid contributions
survive exhaustion or failure; a nearest unresolved candidate is not emitted
as a certified surface. Legacy coverage results remain explicitly approximate.

### 4.5 Contribution contract

Each contribution contains finite world-distance t, geometric coverage in
[0, 1], material identity, and optional feature identity and slice-space normal.
It also declares whether it is a verified surface event or an approximate
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

In a 4D slice, ambient clearance must not feed a surface AA ramp: a ball lying
just outside the slice has small ambient clearance but no slice intersection.
The initial reference filter uses a fixed bounded set of angular subrays with
individually verified intersections; coverage is their hit fraction and each
subray shades independently. Unresolved samples remain diagnostic failures,
not certified misses. A cheaper slice-space filter requires its own validated
estimator. The legacy analytic lattice's approximate filtering remains labeled
as such and is compared separately from this reference.

### 4.6 Appearance and compositing

Appearance maps geometric contributions to color and opacity by depth. It
owns the palette, near fading, fog, and any lighting model. A normal is
optional; unlit patterns need not compute one.

The existing front-to-back `LayerComposite` semantics and saturation threshold
remain shared; do not add a second alpha convention or a new threshold control.
Feed it straight color and coverage times material opacity and appearance
fading; its result remains straight-alpha `Color4`. Surface normals exposed
to 3D lighting are in slice coordinates; any ambient normal stays adapter-local.
For the initial SDF policy, fog/near fading changes appearance but does not
request surfaces behind the first boundary. Analytic lattice layers retain
their existing fading and front-to-back composition.

Continuous participating media need interval integration and are outside the
initial contribution contract. They require a future backend and contribution
extension, not fake surface events at arbitrary march steps.

## 5. Configuration and execution

Configuration has separate pattern parameters, sampling-domain pose, camera
range, appearance parameters, and trace-quality limits. A bounded effect-local
constexpr table declares admitted `(pattern, domain, backend, policy)` tuples,
their defaults, parameter ranges, and resource limits. This is not a new
engine-wide registry.

Expose separate **Pattern** and **View** enums through the existing
`ParamHost` registration and admission hooks. Pattern chooses the geometry;
View chooses 3D perspective or a three-dimensional slice through 4D geometry.
The patterns are cubic, octet, diamond, hexagonal-prism honeycomb,
rhombic-dodecahedral cell edges, affine cubic wires, and periodic spherical
shells. Cubic, octet, affine cubic, and shells support both ambient dimensions;
the remaining cellular graphs admit only 3D. Configuration IDs index the
admitted table, not an arithmetic product of pattern and domain. Selecting a
3D-only pattern from a 4D view adopts its 3D defaults. Unsupported restores
are rejected. Non-cubic patterns and the octet 4D flight preset require
`HS_ENABLE_HYPERLATTICE_EXPERIMENTS`: they are available in the simulator and
opted-in device builds, not standard firmware.

Cellular wires use finite-strut closest-approach contributions and rectangular
cell traversal with bounded neighbors and world-space antialiasing. Sheared
cubic wires transform the lattice basis while evaluating coverage in the
ambient Euclidean metric. Shells use analytic sphere boundary roots; a 4D
slice intersects the actual hypersurface. Shear, stretch, and shell radius
are continuous geometry parameters and interpolate during transitions.
Camera translation wraps by each geometry's translation lattice. The shell
radius range stays below half a cell. Stretch applies only to affine cubic wires.
Traversal exhaustion preserves previously composited layers and is reflected
in Unfinished Rays.

Use one bounded, trivially copyable effect `Params` with stable storage for
registered fields. A configuration change adopts that row's geometry defaults
as one block; it never reinterprets the previous shape's parameters. Common
appearance controls may persist. Unsupported snapshots or external tuples
fail validation without mutating live state. Hidden/inapplicable parameters
must not alter the current configuration; schema refresh uses existing hooks.

Keep `ChoreographedEffect` and its existing transition cancellation semantics.
`Params::lerp` blends all continuous geometry and appearance fields. Discrete
configuration and enum members switch at progress 0.5. Manual presets and valid restores snap as they do today;
a manual parameter edit cancels an in-flight transition. New prepared state
becomes visible together at a frame boundary, never midway through a segment.

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

### 5.1 Resource admission

Architectural demonstrators may run in native tests or the simulator without
shipping on the Teensy. Every firmware-admitted tuple must record its tested
parameter range, work limits, persistent/per-ray memory, and fixed-source
device profile. Use the repository's current size/layout budgets and actual
driver display deadline, not constants copied into this spec.

Production admission requires passing those gates, zero observed deadline
spills in the declared motion/preset/transition sweep, and zero unresolved,
invalid, or exhausted surface searches in its reference validation set.
Analytic adapter horizon truncation is a declared approximation, not hidden
budget exhaustion. An experimental configuration may use a documented degraded
fallback and measured error threshold, but cannot be listed as fully admitted.
Record exclusion when a configuration fails; compiling every demonstrator
into shipping firmware is not an acceptance condition.

## 6. Pattern scope and initial implementations

Arbitrary patterns means an extensible geometry contract, not a promise to
render every mathematical field within a fixed Teensy frame budget. A new
pattern adds a definition and adapters without modifying the camera,
compositor, or unrelated backend implementations.

The implementation stages should demonstrate these cases, with device
admission assessed separately:

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

### 6.1 Periodic curved surfaces: performance experiment

Prioritize a 3D periodic curved surface after the existing torus adapter.
Start with the cosine level set commonly used to approximate the Schwarz P
surface; then evaluate the trigonometric gyroid approximation. These are nodal
approximations, not exact minimal-surface equations or exact SDFs; see the
[level-set definitions in the original research](https://pmc.ncbi.nlm.nih.gov/articles/PMC2709202/).

```text
X = k*x, Y = k*y, Z = k*z, k = 2*pi / period
P = cos(X) + cos(Y) + cos(Z) - iso
G = sin(X)*cos(Y) + sin(Y)*cos(Z) + sin(Z)*cos(X) - iso
```

Choose the boundary of `P <= 0` or `G <= 0` for the first experiment, rendered
as a two-sided first boundary so the camera can travel through either region.
A thick sheet `abs(G) <= thickness` is a separate geometry option; its field
thickness is not constant world-space wall thickness.

The component derivative bounds give conservative global Lipschitz constants
`sqrt(3)*k` for P and `2*sqrt(3)*k` for G. Thus `abs(field)/L` is a safe
clearance lower bound on either side of the zero set, not a hit-distance or
coverage estimate. Fast trig substitutions must account for approximation
error in their bounds or explicitly select an approximate policy. Dividing by
the local gradient magnitude alone is not a global safe-step guarantee.

Benchmark 8, 12, 16, and 24 total field-query budgets per ray, including
verification/refinement queries, at the real segmented resolution. Start
with depth coloring and a single sample per pixel; add analytic-gradient
lighting and reference supersampling as separately measured costs. Compare
against a high-budget native reference over camera translation, cell period,
isovalue, grazing rays, and starts in both regions. Report image differences,
unresolved fraction, mean/peak evaluations, render time, spills, and memory.

The historical HyperLattice capture (supporting artifact removed)
uses a 10,368-pixel live quadrant and a 62.5 ms display window at 600 MHz. That is a
gross budget of about 3,617 cycles per live sample before other work. Sixteen
queries per sample would leave at most about 226 cycles per query if nothing
else ran; the useful allowance is smaller. Its 50.095 ms observed peak is a
baseline, not an additive cost for a replacement backend. These observations
make a bounded experiment worthwhile but do not establish curved-surface
performance. Reprofile under the active driver and use section 5.1 to decide
admission. Optimize proven hot paths only after quality and query guarantees
are established.

## 7. Migration and compatibility

1. Extract camera/domain preparation, appearance, and layer consumption from
   HyperLattice behind the shared stage. Keep its current two endpoint domains
   and preserve their output before adding other patterns.
2. Move cubic/hypercubic geometry and plane candidate production into a lattice
   adapter. Keep its grouping, filtering, and horizon behavior explicit as
   compatibility behavior, with reference tests rather than assumptions that
   it is an exact surface renderer.
3. Remove dimensional rift from enum metadata, transitions, parameter
   validation, tests, and documentation. There are already just two authored
   presets, `cubic-flight` and `hypercube-flight`; retain those string IDs.
   Use dense current selection indices and bump `PARAMETER_SCHEMA_VERSION`
   from version 9 when the layout/semantics change. Existing
   `restore_parameters` rejects older schema snapshots without mutation;
   preserve that behavior. Do not reserve a hole in the option array or
   silently reinterpret old raw numeric selections. There is no existing
   historical snapshot decoder to extend. Any future persisted-format import
   requires its own versioned mapping; it is outside this implementation.
4. Extract and validate the existing volume kernel with Raymarch unchanged,
   then connect existing SDF shapes and the lattice queries to spherical rays.
   Introduce the noncubic analytic pattern and evaluate quality and device
   cost independently. Run the periodic-surface experiment in section 6.1
   after this reuse is demonstrated.

For endpoint compatibility, the old lattice origin is measured in cells,
`sphere_radius` and `wire_radius` scale by cell size to obtain world lengths,
and the old ray t/far distance already uses world units. Preserve old motion,
near-fade, softness, and AA behavior in the compatibility policy until a
separately compared change replaces them. Finite-volume demonstrations use
an unwrapped stationary camera with a bounded animated object placement and
an exterior radial start, rather than inheriting lattice flight coordinates.

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
- Marching tests use known intersections to check exterior and interior
  bounds, first-boundary ordering, starts on surfaces, tangencies, stalls,
  thin neighboring components, and each terminal status. No generic SDF
  continuation capability is implied by these first-surface tests.
- A 4D ball centered just beyond the slice, with separation less than the
  nominal AA width, produces no certified intersection or phantom coverage.
  Exact/reference comparisons include zero projected gradients and unresolved
  subrays. Center starts, tiny scales, and nonfinite queries are covered.
- Existing Raymarch captures and volume tests remain valid after kernel
  extraction. A torus and warped torus render through both camera paths using
  the same shape types and numerical kernel. A 4D lattice slice uses that
  kernel through a 3D query adapter, without a second marching implementation.
- Reference captures preserve the surviving lattice modes during extraction;
  intentional filtering changes receive separate comparisons. Preserve the
  volume first-graze and background-graze regressions. Rift disappears from
  controls; prior-schema snapshots are rejected without state mutation.
- Configuration tests cover every admitted tuple, invalid selection, automatic
  transitions across dimensions, same-tuple interpolation, manual edits,
  pause/resume, restores, and one consistent prepared state per segmented
  frame. A finite-volume camera never jumps at lattice-cell boundaries.
- Adding the triangular framework and implicit surface requires no pattern
  branches in camera, appearance infrastructure, or compositor. Generic
  backend code depends only on its query capability contract.
- The file ownership and dependency rules in sections 3.2-3.4 hold. Core
  geometry/query/trace tests compile without including effects. HyperLattice
  contains no duplicate generic distance, camera, tracing, or compositing
  implementation; reusable patterns and their analytic adapters live in core.
- Relevant native tests, firmware builds, and size/layout gates pass. Record
  RAM1 code, RAM1 variables, FLASH data, and changes from baseline. Profile
  admitted configurations on the Teensy with fixed provenance and report
  traversal exhaustion alongside frame cost. No performance equivalence is
  assumed between analytic and marching backends.
- Each demonstrator records a shipping or experimental/excluded decision
  under section 5.1. The periodic-surface report includes the quality/cost
  sweep in section 6.1; failure to meet the device budget is a reported result,
  not a reason to weaken surface guarantees silently.

Implementation validation covers native regression tests, firmware size/layout,
and fixed-source device captures. Experimental exclusions remain part of the
reported result rather than implicit surface approximations.
