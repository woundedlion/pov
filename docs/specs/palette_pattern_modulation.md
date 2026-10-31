# Pattern sources and palette modulation

**Status: IMPLEMENTED.** Shared pattern sources feed independently owned, ordered
color-modulation stages. The historical qualification section records the earlier
single-layer implementation; current stack validation accompanies each change.

## Problem and scope

Generalize the existing spatial palette noise modifier into independent pattern
sources and color operations, using its existing sphere coordinates.

The abstraction has two axes: **pattern source** and **modulation mechanic**.
Noise plus hue rotation is the existing noise modifier; noise plus positive
lightness lift is shimmer. Grid plus opacity uses the same protocol.
Keep ordinary palette lookup, `get(float)`, unchanged.

A palette colorizer maps Field to Color. Each subsequent modulation maps Color
to Color, so a chain can contain several ordered layers. Each layer independently
selects a noise, sine, interference, smooth-grid or displacement source and a
hue, lightness or opacity mechanic. Graphs, expression parsing, planar modulation
coordinates and remapping remain outside this contract.

## Existing behavior

[The scalar sphere storage](../../core/pattern/sphere_lut.h) contains six
24 by 24 signed-byte faces. The independent
[hue table](../../core/color/hue_rotation_lut.h) stores palette-by-hue colors. Dominant-axis ratios make lookup invariant to positive
direction scaling. Baking evaluates normalized directions through the existing
periodic noise loop.

[The modulation stage](../../core/render/pullback/color.h) reads the sphere and
path metadata retained with its incoming color. Sphere sources sample
`FieldSample::sphere`; displacement reads accumulated path length.

[The carrier contract](../../core/render/pullback/contract.h) defines the frame:
PLANE sampling retains the projection-frame direction; SPHERE sampling retains
the camera-frame direction. Plane endomorphisms do not rewrite the sphere point.
Inherit these semantics without another transform; this is not universal world space.

[Positive lightness lift](../../core/color/modulation.h) shares the table and lifts
OKLab lightness over positive lobes, preserving alpha: another source/mechanic pair.

## Shared sources with pullback

Build one canonical pattern library for both pullback rendering and modulation,
not a second catalogue. Extract reusable signed kernels, passive parameter records,
prepared math and source-specific phase math from
[pullback sources](../../core/render/pullback/source.h): planar rings, spiral,
interference, coupled grid, sphere rings, noise and other reusable raw kernels.
Reuse lower-level lattice/SDF kernels where appropriate without relocating those
subsystems. Keep genuinely different algorithms explicit: PlaneCoupledGrid and
SphereSmoothGrid are different domain-specific sources, not duplicate names for
one promised pattern.

Likewise shared noise must express basis/octave selection, contrast shaping and
loop construction as source parameters or meaningful algorithm variants. Current
pullback noise applies octave/contrast logic; palette noise samples one configured
generator value. Preserve those formulas explicitly in the common source library,
not consumer-specific copies or silent substitutions; mechanics do not select them.

Core sources know typed coordinates and pure values, not renderer carriers,
bindings or controls. Pullback retains provider/state adapters, coordinate
admission, Field metadata, projection transforms, weighting, coverage, conversion
to [0, 1], and carrier provenance/path length. Its Source policies and Sample
operators are real renderer bindings, not compatibility wrappers. Both those
policies and operator implementations call the shared kernels directly; remove
old copied free functions and forwarding aliases. UI descriptors wrap/visit each
source's single parameter record rather than duplicate its fields.

**A source does not inherently need a projection.** A planar source needs a
chosen mapping only when its consumer starts with sphere coordinates. Pullback
already supplies plane coordinates; sphere sources take directions directly.
Projection stays in the renderer's coordinate pipeline, not the shared kernel.
The library supports PlanePosition and SphereDirection now for existing consumers.
Directional modulation binds SphereDirection; displacement binds accumulated
path length. FieldSample does not retain warped plane coordinates. Supporting a planar modulation source
later needs an explicit coordinate binding, not a silent plane/sphere conversion
or universal context bag.

## Scalar and typed coordinate contracts

A spatial pattern takes its typed coordinate and immutable prepared state and
returns a finite scalar in [-1, 1]. The renderer-bound displacement source
instead returns finite, nonnegative accumulated path length without unit-range
clamping. Pattern output is independent of palette contents
and mechanic. It depends on the palette lookup coordinate only if the caller
explicitly supplies that coordinate as its ScalarCoordinate. Zero has no universal
identity meaning: a mechanic determines what it means. The source knows
coordinates, time and its own parameters; it does not know hue, lightness or
alpha. A mechanic receives color and a scalar, and knows neither coordinates
nor which source produced the scalar.

Unify the protocol across domains, not palettes into dimensional variants:

| Domain type | Meaning |
|---|---|
| ScalarCoordinate | A float coordinate for a one-dimensional pattern. |
| PlanePosition | A two-dimensional point, represented by a pair or Complex. |
| Position3D | A volumetric position; distance from the origin is meaningful. |
| SphereDirection | A unit direction represented by Vector3; radius is discarded. The sphere is intrinsically two-dimensional. |

A field exposes its Coordinate type, shared with its source. The wrapper accepts
that type; mechanics remain domain-agnostic. Avoid a universal context bag,
Vector3 for every domain, or per-pixel domain switching. The palette-by-hue LUT
is internal acceleration, not a user-facing 2D palette.

ScalarCoordinate and Position3D bindings remain future work. For SphereDirection,
the boundary admits vectors whose computed
float squared length is finite and at least `math::EPS_NORMALIZE_SQ` (currently
1e-12), matching [Vector normalization](../../core/math/3dmath.h). Explicitly
check finiteness before normalization; the existing lower-bound check alone
does not reject infinity. Apply this admission to directions and source axes.
Construct the unit SphereDirection tag only at this validated boundary.

Normalize once per sample and share that coordinate across retained endpoints.
Cube baking supplies unit directions; cube lookup may elide normalization using
scale-invariant ratios after admission. Migrate raw-vector callers to this
validated boundary. Test scale invariance only while both vectors remain
admitted. Neither path interprets direction length as radius or normalizes
Position3D; `FieldSample` gains no unit-length requirement.

Prepare validates parameters and axes. Debug assertions diagnose nonfinite
custom output; the scalar boundary clamps finite numerical overshoot, not
incorrect source normalization. Spatial sampling is immutable and bounded, with no allocation, random-state
advancement, I/O or clock reads. Time and seed are explicit.

## Generic protocol and initial patterns

Use a small policy protocol, consistent with the existing prepare/sample split
in [pullback sources](../../core/render/pullback/source.h). The following types
are schematic, not existing declarations:

```cpp
struct SphereDirection; // Constructed by the validated domain adapter.

struct DirectionalSine {
  using Coordinate = SphereDirection;
  using Params = SineParams;
  using State = SineState;
  struct Prepared { math::Vector axis; float frequency; float phase_turns; };
  static State init_state(const Params&);
  static State advance(const State&, const Params&);
  static Prepared prepare(const Params&, const State&);
  static float sample(const Coordinate&, const Prepared&);
};

template<class Pattern> struct DirectSphereField;
template<class Pattern> struct BakedSphereField;

struct Opacity {
  using Params = OpacityParams;
  struct Prepared { float amount; };
  static Prepared prepare(const Params&);
  static Color4 apply(Color4 color, float scalar, const Prepared&);
};

template<class PaletteSource, class Field, class Mechanic>
struct PatternModulatedPalette {
  const PaletteSource* palette;
  typename Field::Prepared field;
  typename Mechanic::Prepared mechanic;
  Color4 get(float value, const typename Field::Coordinate& coordinate) const;
};
```

Each source declares bounded passive State, possibly empty, and pure initialization,
advance and prepare operations; stateless operations may be no-ops. SineState has
one phase in turns. Existing multi-phase sources retain primary/secondary/angle
state with their declared radian periods, noise retains its loop phase, and other
inputs remain explicit. Preparation may take additional declared pure-value inputs,
such as a renderer-supplied sphere-ring axis; it never takes renderer resources.
Camera, spatial walks, projection and provider machinery stay renderer-owned.

The field adapter exposes prepare/sample, wrapping either direct evaluation or
cube storage. The mechanic exposes prepare/apply and returns a valid straight-alpha
`Color4`. The generic wrapper samples its base palette, samples the field, then
applies the mechanic. HueRotation and PositiveLightnessLift follow the same
mechanic protocol as Opacity. Nonzero hue/color LUT fusion is a qualified
approximation of this reference, not bitwise equivalence. Accelerated
routes must still honor the exact identity rules below.

The modulation binding admits these sphere sources alongside the displacement
source and the planar/sphere catalogue extracted for pullback:

| Source | Definition and controls |
|---|---|
| Noise | Adapter around the current configured noise generator, spatial scale and sphere-loop construction. Seed and generator configuration belong to source identity. |
| Directional sine | `sin(2*pi*(frequency*dot(axis,direction) + phase + offset))`; unit axis, nonnegative frequency and offset in turns. Frequency is cycles per unit axial coordinate, so the sphere spans twice that many cycles along the axis. |
| Sphere interference | Half the sum of two directional sines, with separate axes, frequencies and offsets. Both use the same loop phase; the fixed half gain bounds the result without dynamic normalization. |
| SphereSmoothGrid | Product of two directional sines over independently normalized axes, with independent frequencies and offsets and a shared loop phase. Admit finite, nondegenerate axes at any angle, including parallel axes. The product remains in [-1, 1]; this is not the existing planar coupled-grid algorithm. |
| Displacement | Direct accumulated path length from the incoming color carrier. No direction normalization, cube, frequency or moving source phase. |

New sphere sine/interference/grid phases use turns per engine update, including
negative speeds, and repeat after one turn. Existing source clocks keep their
own units, periods and rate coupling; do not force them into one float clock.
Modulation speed controls share noise limits of +/-0.001 turns per update and
analytic limits of +/-1 turn per update across workbench schemas and composed
sliders. Composed unit-mechanic amount sliders use [0, 1]; the workbench amount
schema remains [-4, 4], with each mechanic's validation enforcing its own domain.
Composed preset banks retain preset 0's source kind and mechanic kind across
every row; their registered modulation controls stay bound to those fields.
Prepare hoists frame-constant work. Admission qualifies upper frequency
and argument bounds compatible with the chosen trigonometry and the project's
fast-trig limits. Fold phase and keep arguments within that range; this proposal
does not invent a universal cutoff or promise detail free of spatial aliasing.

Custom C++ sources and mechanics follow these bounded, immutable protocols.
Compile-time extension needs a stable recipe ID and serializable parameter
schema for replay. Built-in selection resolves during preparation; callback
pointers are not recipes.

## Mechanics and valid combinations

| Pattern source | Hue rotation | Positive lightness lift | Opacity |
|---|---|---|---|
| Noise | Existing noise hue modifier | Existing shimmer | Organic alpha mask |
| Directional sine | Hue bands | Lightness bands | Smooth alpha bands |
| SphereSmoothGrid | Hue grid | Lightness grid | Grid alpha mask |
| Displacement | Hue turns per path-distance unit | Distance-driven lightness lift | Distance-driven alpha mask |

These combinations, and interference with any mechanic, are structurally valid.
Admit each recipe only after qualifying its resources, approximation error and
execution cost. Do not create an enum value per source/mechanic pair.

Hue rotation resolves `hue_turns = amount * scalar`, with signed amount in
turns, periodic hue lookup and existing gamut mapping. It preserves alpha.
PositiveLightnessLift resolves `lift = amount * clamp(scalar, 0, 1)` for amount in
[0, 1], then moves OKLab lightness toward white with the existing gamut-aware
chroma reduction. It preserves alpha; negative lobes do not darken. These two
mechanics leave color unchanged at scalar zero.

Opacity resolves `u = clamp((scalar + 1) / 2, 0, 1)`, then
`gain = lerp(1, u, amount)` for amount in [0, 1]. It returns unchanged RGB and
`alpha = input_alpha * gain`. Scalar zero is a half mask at full amount, not
identity. At full amount, scalar -1 removes alpha and scalar +1 preserves it.
For reference and accelerated routes, **amount zero is
exact Color4 identity for every mechanic**. A zero resolved hue rotation,
including zero after periodic wrapping, returns the exact base palette sample;
bypass the hue/color LUT. Positive lightness lift also returns its exact input
when lift is zero. These identities apply before downstream brightness/coverage.
No mechanic changes the source range or asks for a target-specific value.

Existing active LUT routes can resample palette color even at zero hue shift
because their palette axis is quantized. Remove that behavior in the canonical
implementation: zero resolved rotation samples the base palette exactly. Review
resulting render changes and update goldens intentionally.

Zero amount may skip field evaluation and preparing unused mechanic/color
resources, while source clocks continue. Keep any still-needed retained endpoint
resources valid. Configuration validation remains mandatory even when skipped.

The colorizer samples its palette. Color-to-Color layers accumulate their
selected mechanics before the final RGB brightness envelope. Hue layers add
turns across lightness and opacity layers. Positive lightness lifts compose as
`1 - product(1 - lift)` toward white. Temporal endpoints blend this deferred
state before the combined lift and hue rotation apply in OKLab, followed by
one final gamut mapping and RGB conversion per pixel. No intermediate layer
clips RGB. A palette hue lookup remains valid only when the blended base color,
palette coordinate and table are unchanged; otherwise finalization uses the
blended RGB base. Zero combined adjustments retain the exact blended base
color. Opacity changes only alpha and does not interrupt color composition.
The result remains straight alpha. Final premultiplication stays in
[the scan shader](../../core/render/scan/shader.h), after the chain; opacity
must not multiply RGB early. Qualified optimized kernels implement this same
contract; they do not become independent noise-hue or shimmer engines.

## Direct and baked sphere evaluation

SphereDirection recipes explicitly select evaluation strategy. Only sphere
directions normalize; future Position3D fields preserve radius.

* **Direct:** evaluate the prepared source at the admitted unit direction.
  Detail remains limited by scan resolution and floating-point math.
* **Baked cube:** evaluate the source into the existing six-face layout, clamp
  and round to signed bytes with the existing quantization, then bilinearly
  sample it. The cube table becomes reusable scalar storage independent of color.

Do not switch strategy automatically with frequency or load. Recipes select the
same path across native/WASM/firmware: migrated noise initially baked, analytic
sources initially direct. Baked analytic recipes require error qualification.

Cube storage cannot reproduce arbitrary frequencies or sharp boundaries.
Quantization contributes at most 1/254 scalar error at a bake point; spatial
interpolation adds source-dependent error. Each baked recipe declares admitted
parameter ranges and measured bounds; reject unqualified configurations. Preserve
cube face conventions and duplicated edge evaluation to avoid seams.

## Ownership, cache invalidation and replay

The effect/runtime arena owns prepared state and tables; views borrow through
the draw and retained interval. Direct endpoints own/copy complete source inputs
or borrow immutable retained generators. A live pointer to a mutated generator
is insufficient. Do not retain mutable controls or overwrite a still-read table.

The current noise cache keys scale and phase under a documented precondition
that its generator configuration and destination remain unchanged. Replace it
with one complete field cache identity covering source
kind/recipe version, all source parameters, seed and generator configuration,
source state/prepared inputs, and bake layout/version. A configuration generation may replace the
configuration portion only if every mutation increments it. Reset cache validity
on ownership changes and restore. Mechanic kind, amount and other mechanic
controls, and palette changes, do not invalidate the spatial table. Any
mechanic/color cache has a separate identity: palette changes invalidate the
hue-rotation color table; spatial pattern changes do not.

Snapshots store enabled state, source recipe/parameters, evaluation strategy,
seed where applicable, typed/versioned source state, and mechanic kind/parameters.
Derived tables and pointers are not serialized.
Restore invalidates derived state and rebuilds needed resources deterministically
before first use; invalid or unavailable recipes fail explicitly. Exact continuation of
a partially drawn interval also needs the previous endpoint recipe/resources
and interval position; restoring one current endpoint alone must not claim
that guarantee.

Update the independent owners and codecs together:
[composed runtime](../../core/render/pullback/composed_runtime.h),
[retained intervals](../../core/render/pullback/composed_interval.h),
[snapshot records](../../core/render/pullback/runtime_snapshot.h),
[operator snapshots](../../core/render/pullback/operators/snapshot.h),
[workbench host](../../workbench/shader/chain_host.h),
[WASM snapshot codec](../../targets/wasm/chain_snapshot_codec.h) and
[TypeScript declarations](../../targets/wasm/engine_types.d.ts).

## Temporal behavior and topology

Temporal colors blend deferred endpoint state before one perceptual finish.
Source phase and parameters are not additionally interpolated inside retained
endpoint evaluations.

For a new pattern interval with entirely unchanged modulation topology,
retain the two prepared source states or two baked tables. Evaluate each at
the column direction. Sample each endpoint palette with the shared column
mapping, apply the mechanic using shared column-interpolated continuous mechanic
controls. Accumulate hue turns and lightness lift without converting or clipping
either endpoint, and apply alpha coverage/value opacity to each endpoint.
Interpolate the base RGB with Q16 rounding and alpha with the original floating
fraction. Interpolate accumulated hue along the shortest arc on the unit circle;
interpolate the resolved lift and final brightness gain linearly. Apply the
combined hue and lift once to this blended base, gamut-map once, then apply
the final RGB brightness gain. Fraction zero and one preserve their endpoint
state exactly. Do not interpolate raw scalars before resolving mechanics: a
negative lobe must not cancel an endpoint's positive lightness lift. In particular,
equal base palette alpha does not imply equal modulated alpha: the common-alpha
fast path must be disabled for opacity unless equality is proven after modulation.
Alpha-preserving mechanics may use an equality fast path where that premise is
proven. Pin interpolation weights and rounding with tests; do not assume
algebraic rearrangements are bit-identical.
The direct perceptual finish clamps gamut-mapped floating RGB, applies brightness,
then quantizes once to RGB16. A unit gain preserves the unscaled quantization;
zero gain returns exact black with the original alpha. Hue-table samples retain
their existing RGB16 brightness scaling.

`Spec::TEMPORAL_COLOR_PHASE` opts eligible authored modulation providers into a
different interior color-time model. Blend the base palette RGB and alpha first,
sample eligible spatial sources at the shortest interpolated column phase, then
apply each hue/lift layer once. Smooth grids interpolate their axial phases.
Direct configured-single noise interpolates the loop clock and resolves its loop
offset once per column, then samples the retained configured generator once per
pixel. Noise admission requires unchanged noise parameters and generator seed;
baked and octave-contour noise retain endpoint-tail evaluation. Displacement
layers use the sampled geometric path length. All color-tail stages must declare and admit this protocol; an
unknown stage, unsupported source/mechanic, or changed source recipe keeps the
complete endpoint-tail model. Literal fractions zero and one remain canonical.
Grid preparation retains only the normalized axes, frequencies and phases;
it does not copy the noise-preparation variant into column storage. Noise
preparation borrows the current immutable generator and stores the column loop
offset; its endpoint owner survives the complete column draw. Millefiori and
AshCloud enable this model; other authored specs default to the endpoint-tail
model.
Both models compose hue and lift before one final gamut mapping. Alpha is
unchanged because the admitted shared mechanics preserve it.

`Spec::PHASE_NOISE_CACHE` optionally approximates AshCloud's direct configured-single
noise within a dodecahedral chamber. It leaves the authored source, evaluation
strategy, parameters and document topology unchanged. Each endpoint bakes a
floating scalar patch in the chamber's x/y, z/y chart. Baking rotates each chart
direction by that endpoint's projection and samples its immutable noise source.
A column recovers chart coordinates with the inverse sampled projection, samples
both endpoint patches bilinearly, interpolates their scalars at column time and
applies the admitted hue/lift mechanic once. The explicit chart bounds and sample
counts are declared by `ChamberNoisePatch`. Cache views retain actual column
time independently of the palette blend weight, which can collapse to one when
both source clocks and palettes are equal but the projection changes.

Cache admission is restricted to the qualified authored noise recipe: configured
single simplex, zero contrast, scale 2.4029531, speed 0.00015, seed 6047 and recipe
version 1. Other admitted controls retain direct evaluation without baking an
unused patch. Large clock or projection steps suspend caching; the bounds are
`ChamberNoisePatch::MAX_CLOCK_STEP` and `MAX_PROJECTION_STEP_SQ`. Step checks follow
the short clock arc and treat antipodal quaternions as one orientation. A rejected
endpoint is not rebaked by repeated preparation; subsequent continuous endpoints
prime the alternating storage again.

Two runtime-owned slots alternate only when the endpoint identity changes.
Their identities include the complete active noise recipe, seed, loop phase and
projection quaternion. Frame views retain the expected quaternion, and source
identity is checked against the retained immutable prepared endpoint. The previous
frame's slot remains immutable through the next draw. Direct source retention
never marks its ordinary cube active or requires history cube storage.

Literal endpoints and static frames retain canonical evaluation. Missing or stale
patches, changed projection controls and directions outside the chart retain the
existing direct column-phase sampler; out-of-chart coordinates are not clamped.
Changed noise recipes and unsupported variants retain the existing temporal
admission and endpoint-tail fallback. Cache storage and frame views exist only
for an opted-in direct-noise resource and remain inside the existing persistent
arena budget. The cache changes color sampling, leaving geometry, path length,
coverage and alpha intact. Qualification compares its RGB16 output against the
same effect with this cache alone disabled; its measured error is recorded
independently of the palette and radial approximations. `Color::Modulate`
exposes this additional internal approximation through
`PHASE_NOISE_APPROXIMATION`. The existing ordinary hue/noise LUT metrics do not
bound its additional chamber-cache error; qualification is finite capture
coverage, not a universal bound for all future frames and color controls.

`Spec::PALETTE_LAB_CACHE` optionally supplies the eligible shared tail with a
palette-aware base cache. Each node first blends the previous/current RGB16
palette samples with the canonical Q16 weight, then converts that RGB to OKLab.
The cache has 64 palette-coordinate nodes and three time rows (zero, half and
one); a column prepares its neighboring time rows, and a pixel interpolates
four Lab nodes. The cached value is copied into the color carrier before the
same composed hue, lift, gamut mapping and brightness finish.

Two effect-owned buffers preserve consecutive frame views. A frame never bakes
into the buffer borrowed by its previous endpoint; carriers copy Lab values,
not time-view pointers. Literal/static draws, changed color-family controls,
nonopaque palette alpha, and unsupported/custom color tails retain canonical
RGB sampling. Opaque palette alpha preserves the common coverage and opacity.
Derived color policies must explicitly own their
shared-base contract before the stage can bypass their endpoint operations.
Cached carriers enter after the sole base blend; interior blending of an
already cached carrier is outside that internal protocol.
Eligible phase tails adjust deferred state without reading the carrier's RGB.

The `PALETTE_RGB_BLEND_TO_LAB_LUT` oracle compares the cache to RGB-first palette
sampling, not interpolation of endpoint Lab. Qualification keeps existing hue
and framebuffer references unchanged and records the new cache's measured
RGB16 error independently.

A source kind/custom recipe version, evaluation strategy, mechanic kind or
enabled-state change is a topology change for that layer. Hold that layer
against its **entire current modulation**:
source parameters and current prepared state/resources (initialized on a source
kind/version switch), mechanic typed parameters/state
and enabled state, including fraction zero. Continuous mechanic controls interpolate
only when the entire topology is unchanged, even if the mechanic kind alone matches.
For example, sine-to-grid with opacity amount
0.2-to-0.8 uses the current grid and amount 0.8 throughout the switch interval.
Hue turns are never interpolated as opacity amount.
Apply that mechanic to both palette endpoints and interpolate their deferred
color states before the shared finalization, as above. The switch may therefore be
discontinuous; expose and test that policy rather than claiming cross-pattern
morphing. Ordinary unchanged-topology intervals retain exact endpoint behavior.
The fraction-zero result uses the previous palette with the held current
modulation, including its alpha effect; it need not equal the previous fully
modulated frame. Never reinterpret an old source's or mechanic's resources
under the new type.
Disabling therefore holds identity across that interval; enabling holds the new
prepared modulation. Other unchanged-topology layers retain their own endpoints.
At each temporal fraction, evaluate the complete ordered color suffix at each
endpoint without finalizing it, blend the deferred states, then finalize once.
Blending only the last color stage would resolve upstream mechanics against
the wrong endpoint.
The FIELD-to-COLOR stage may expose `run_pair` to evaluate retained and current
carriers together during an interior blend. Generated palettes share mapping,
brightness and opacity evaluation, then sample their separate immutable palette
and hue resources. Each result still traverses its complete color suffix.
Stages without this optional hook use individual endpoint evaluation; fractions
zero and one select the corresponding endpoint directly when the hook exists.

### Burst-ripple radial transfer

`Spec::RIPPLE_RADIAL_TABLE` defaults off; Millefiori and AshCloud enable it for temporal
column rendering. Each persistent resource owns radial curves of net rotation
and absolute travel for its burst centers. Frame preparation bakes current wave
ages and ages minus one, preserving wave order. Bitwise-identical ordered ages
and wave controls reuse a retained current curve for the previous endpoint. Column views borrow the current
resource's tables and blend those endpoints only across a one-frame interval;
table contents remain valid until the next draw. Initial identical frame
endpoints use the current curve. Static/direct evaluation and literal interval
fractions 0 and 1 use the canonical grouped kernel, preserving document parity
and endpoint identities.

The curve is a spatial and temporal approximation, with one final sphere
rotation per center and unchanged downstream pixel shading. Union support
bounds reject inactive rays before angular lookup. Unsupported controls,
parameter changes, skipped frame intervals, pole crossings and invalid lookup
cells fall back to grouped evaluation. The optional curves, retained keys and cell-validity bytes occupy 17,700 bytes
inside the existing persistent arena, counted by the composed footprint gate;
no arena or stack budget is increased.

The runtime opts immutable endpoint pairs into a derived cell-validity mask:
256 bytes per burst group, replacing the same four negative-path checks at
lookup. `prepare_radial_tables` defaults to ordinary validation for public or
manually edited curves. Its `immutable_view` opt-in clears mask readiness before
any update and publishes a new mask only after both endpoints are finalized.
Disabling interpolation also clears readiness. Curve edits require invalidating
that derived view before lookup; interpolation and fallback arithmetic stay
unchanged. Exhaustive cell fixtures and 600-frame update/reuse comparisons retain
bitwise direction and travel parity with ordinary validation.

Qualification covers geometry, unit norm, path continuity, burst births/deaths,
slot replacement and shared or random centers. The default preset's 1.536 million
sample qualification observed maximum direction difference 0.00463 radians and
squared-norm error below 3.6e-7. Sixteen parameter corners covering 24.576 million
samples observed maximum direction difference 0.00771 radians. Full-frame motion
and visual continuity qualify color changes; legacy RGB16 equality is not an
acceptance gate for this temporal approximation.

### Validated dodecahedral matrices

AshCloud's interior columns may select a reflection matrix from a coarse cube
cell. The constexpr tables hold the 120 dodecahedral reflection-group matrices
and six 64-by-64 cell maps in FLASH; they allocate no persistent arena or RAM
storage. A cell selects a candidate transform, rather than a quantized output
direction. The transformed input must lie strictly inside all three canonical
chamber walls, beyond the lookup's wall margin. A rejected candidate restarts
the ordered canonical fold from the original input. Literal column fractions
zero and one retain the canonical fold and document-built chain parity.

The accepted transform represents the same reflection-group action, with
floating-point differences from matrix multiplication. Canonical comparison,
unit-norm checks and transformed wall fixtures qualify those differences and
pin rejection near every chamber wall. Changing the cube cell classification
cannot bypass the chamber validation.

## Canonical configuration and migration

Each `color.modulate.v1` node owns one configuration containing source controls,
evaluation strategy, mechanic controls and enabled state. Disabled modulation is
identity. Validate recipes regardless of enabled state or amount. The
`colorize.generated-palette.v4` node owns palette-family selection, its family
controls, mapping, brightness and value opacity; it has no modulation controls
or hue-shift mode. Existing cycling-generated, generated, procedural, gradient and mutating
palette families remain independent of the modulation stack.

Each layer advances its own source clock once per update, including disabled and
zero-amount layers. Mechanic, evaluation and same-source parameter edits retain
that clock; source-kind/version changes initialize the selected source. Restored
clocks retain typed source identity. Disabled layers may skip sampling and bakes.

The displacement source reads accumulated path length from the field carrier.
It is direct, has no source clock motion, and preserves distances greater than
one. Opposing warps contribute both traveled distances. For hue rotation its
amount is turns per unit of traveled distance. Other mechanics clamp their
derived unit-range factors, keeping lightness and alpha within their supported
ranges.
Sphere patterns retain their bounded [-1, 1] scalar contract.

Shader documents keep document schema 2 and version their operator IDs. Import
converts a legacy `colorize.generated-palette.v3` node to the palette node,
followed by a displacement/hue layer and then its spatial modulation layer.
Both layers remain present when disabled, preserving dormant controls and
preset changes. Names are `displacement` and `modulation`, with `-2`, `-3` and
later suffixes when labels collide; original chain labels are reserved first.
Hue-shift mode becomes the displacement layer's enabled flag and hue-shift amount
becomes its amount. Existing spatial controls move to the spatial layer.
Serialization fields and staggered parameter groups follow the renamed controls.
Fresh documents begin with a palette node and add only requested layers.
Legacy standalone `modulate.field.v1` nodes become unified modulation nodes with
the same labels and recipes. Legacy `modulate.hue-shift.v1` nodes become
displacement/hue nodes in place; their default enabled mode is preserved.
Standalone-node conversion never inserts another generated pair of layers.

Snapshot schema 4 stores modulation recipes by modulation-node instance.
Palette clocks use `color-clock-v3`; layer clocks use `modulation-clock-v1`.
Legacy snapshot migration preserves the spatial recipe, seed and clock while
initializing the displacement source at phase zero. Derived resources are rebuilt.
Endpoint lookup uses exact instance labels, including after reorder or state
migration; a hash collision cannot merge two layers. Each active layer owns its
prepared source and caches, independently of the palette and other layers.

Migrate all noise-hue and shimmer callers, includes, controls, snapshots and
tests directly. Noise plus HueRotation and noise plus PositiveLightnessLift
are configurations of the generic implementation, not named compatibility
wrappers. Remove old palette headers after their final callers migrate and update
documentation references to the replacement modules.
Reuse existing noise-loop, generator and gamut math; review intentional render
changes from exact zero identity and unified rounding. Bump snapshot/preset
schema versions where needed. Supported legacy documents and snapshots convert
explicitly to the current representation; unsupported versions fail without
changing live state.

Source and mechanic selectors have independent stable identifiers in the new
schema. Custom policies are compiled in and replayed through registered recipes,
never persisted callback pointers. Existing one-dimensional NoiseWarpModifier
and RippleModifier remain prelookup remaps outside this Color4 mechanic contract.

## Code organization

Templates and small kernels match the repository's existing style; no per-platform backend copies or source/mechanic
pair files are needed. New APIs live in dedicated pattern/modulation namespaces,
not new global Noise, Grid or Opacity names.

Header layout:

```text
core/pattern/
  types.h                 Typed coordinates, source Params/State records/IDs
                          with explicit units and evaluation selector.
  field.h                 Source/field contracts and direct evaluation adapter.
  sources.h               Shared analytic kernels/preparation/phase math,
                          including existing planar and sphere sources.
  noise.h                 Noise policy; isolated FastNoiseLite dependency.
  sphere_lut.h            Scalar cube geometry, views, sampling, generic baking
                          and cache support; no noise-generator dependency.
core/color/
  modulation_types.h      Mechanic parameter records/IDs and canonical config;
                          imports passive pattern types, no samplers/UI/Color4.
  modulation.h            Three mechanics and generic palette composition.
  hue_rotation_lut.h      Palette-by-hue color table and qualified color fusion.
core/render/pullback/
  modulation.h            Built-in binding, parameter metadata/registration,
                          validation, scheduling and owned source state/resources.
  modulation_temporal.h   Retention/equality and completed-endpoint color kernels.
  color.h                 Palette crossing, ordered mechanics, deferred hue and finish.
  operators.h             Independent palette and modulation operator schemas.
tests/pattern/
  fields.h                Typed-domain admission and direct source cases.
  sphere_lut.h            Generic bake, seams, sampling and invalidation cases.
tests/color/modulation.h  Mechanics/composition, identity and color-LUT cases.
```

Source policies alias pattern records; mechanics alias color records. The
canonical config assembles them without implementation imports or duplicate
definitions. Codecs use this passive schema; UI metadata and serialization
implementations stay separate.

Dependencies flow integration -> color/pattern -> math/support. Pattern code
never imports color/render/controls/UI; analytic users never import FastNoiseLite.
Generic color composition and mechanics import no concrete sources, samplers,
cube storage or registry. Passive record imports are allowed, but mechanic apply
consumes only Color4, scalar and its own parameters. Composition templates the
supplied Field/Mechanic. Color fusion and scalar storage remain separate.

The renderer binding joins admitted built-ins and resolves dispatch before hot
loops. Shared helpers accept config/prepared views; they never import the
higher-level colorizer, composed runtime or workbench host. Entrypoints call
helpers, not vice versa. Ownership stays in arenas. Add no global registry state,
singleton or all-built-in umbrella facade.

The scalar cube geometry, sampler and generic bake/cache live in sphere storage;
configured generator evaluation lives in the noise policy. Palette-by-hue tables
and sampling live in the color LUT, and generic composition and positive-lightness
lift live in the color module. The retired combined modifier headers have no
forwarding headers or aliases.

The colorizer/composed runtime and workbench host use shared validation, clocks
and preparation. Existing snapshot/operator codecs carry canonical state; WASM
JSON and TypeScript definitions stay in their target files referenced above.
Schema and binding changes land together.

Migrate [composed source bindings](../../core/render/pullback/composed_policies.h),
[providers](../../core/render/pullback/composed_providers.h) and
[sample operators](../../core/render/pullback/operators/sample.h) to the shared
source library. The existing source header retains renderer bindings and field
metadata, not another kernel catalogue. Do not move camera walks, projections,
lenses or provider factories into the pattern module.

Place new source/domain/cube cases in the proposed pattern tests and mechanic
cases in the proposed color test. Register them through the existing native
module/CMake harness. Extend
[composed-effect temporal tests](../../tests/composed_effect/temporal_sampling.h)
and [shader-chain lifecycle tests](../../tests/shader_chain/lifecycle_cases.h)
for integration; retain migrated coverage from
[existing palette cases](../../tests/color/palette_layer.h).

Adding a source changes its pattern records/kernels, integration registration and
tests, not mechanics. Adding a mechanic changes its color records/implementation,
registration and tests, not source kernels. Extend
[include-direction checks](../../tools/include_paths.py) and
[their fixtures](../../tools/include_path_tests/test_include_paths.py) for these
boundaries. Add standalone
header compilation proving analytic, schema and generic composition independence.

## Budgets and acceptance

One current signed-byte cube costs 3,456 bytes; retaining both endpoints costs
6,912 bytes. Optional float copies cost another 13,824 bytes per endpoint.
Reuse existing tables for migrated noise users. Direct policies allocate no
cube, and effects with no pattern allocate no pattern resources. Do not enable
float mirrors or instantiate every built-in for every effect by default.
Account separately for hue/color LUTs, optional mirrors, retained prepared source
and mechanic state, immutable generator configuration, cache metadata and maximum
simultaneously live endpoint/transition resources; cube bytes are not the total.
Each additional layer consumes ordinary chain operator, parameter and arena
budgets. Direct displacement layers allocate no scalar cube; baked spatial
layers retain independent cubes. Limits come from the exported catalog.

Implementation must report incremental RAM1 code, RAM1 variables and FLASH
data against the baseline, passing the repository's current firmware gates.
Qualify the selected, bounded capability set across representative source by
mechanic by applicable strategy by temporal-mode combinations on Teensy under
the segmented driver. Measure preparation and scan/column timing, including both
endpoints, especially direct PositiveLightnessLift's OKLab/gamut conversions and
opacity. Admission must stay inside timing budgets; no universal cycle estimate
or exhaustive Cartesian instantiation is required. Convert domains once and share
coordinates, avoid duplicate source evaluation for identical endpoints, and skip
identity work. Resolve dispatch outside hot loops. Baking performs 3,456 source
evaluations per changed key and none for an unchanged key. Scan allocation is forbidden.

Acceptance requires:

* Migrated noise/shimmer caller coverage across seeds, positive/negative speeds,
  zero amount/shift and temporal paths. Assert exact Color4 identity at zero
  amount and resolved hue rotation, qualify nonzero approximation error, and
  review intentional render-golden changes instead of preserving old LUT numerics.
* Source known-value, bounds/periodicity, admitted scale-invariance and rejected
  invalid-axis, below-floor or nonfinite squared-length tests.
* Cross-axis tests using one source with all mechanics and one mechanic with
  different sources. Zero amount and mechanic-specific scalar-zero rules hold;
  verify positive-lightness gamut and alpha behavior.
* Compile-time rejection of mismatched domains and implicit Position3D-to-direction
  conversion; volumetric radius survives. No extra renderer plumbing is needed.
* Asymmetric orientation cases proving projection-frame versus camera-frame
  provenance, and cube edge/corner seam tests from both adjacent faces.
* Direct-versus-baked error qualification over admitted ranges, separating
  quantization and spatial error.
* Temporal tests whose expected endpoint RGB differs from scalar-first blending,
  including changed phase/scale, stationary resources and topology holds.
  Test source changes with same-kind amount changes, mechanic/strategy changes,
  enabled transitions and their fraction-zero current-state holds.
  Verify opacity's modulated alpha at both endpoints and interior fractions,
  including equal base alphas but unequal masks; scan premultiplies exactly once.
* Cache tests for every key component, palette-only and mechanic-only changes,
  with no spatial rebake when only the mechanic changes; snapshot replay
  after cache reset, and unknown-recipe/unsupported-version rejection.
* Source-state lifecycle tests: one scheduled owner, source-specific units/rate
  coupling, empty-state sources, disabled advancement, retention on same-source
  parameter/mechanic/strategy changes, initialization on kind/version switches,
  typed snapshot replay, rejection of invalid disabled recipes, and concurrent
  displacement and spatial layers, including hue cancellation and
  lightness/opacity composition.
* Raw-sample parity between pullback and modulation adapters given identical
  pattern/coordinate/parameters/prepared state, before LUT approximation.
  Preserve crossing value conversion, weighting, coverage, provenance and path
  semantics with independent renderer tests.
* Include-direction fixtures, standalone header checks and migrated caller builds;
  production callers/tests no longer import retired modifier headers/types, and
  documentation references move to the replacements when those files are deleted.
* Native and WASM validation plus firmware memory/timing evidence for the
  implementation.


## Historical single-layer qualification (2026-10-29)

The measurements below describe the earlier combined-colorizer implementation.
They qualify the shared pattern kernels and record its resource baseline; they
do not establish memory totals or timing limits for an arbitrary modulation stack.

The canonical records and pure sources live in `core/pattern`; independent
mechanics and the generic palette adapter live in `core/color`. Composed effects
and the workbench share the scheduled source clock, immutable preparation and
completed-endpoint interpolation. The former noise/shimmer modifier headers are
removed. The 19 composed effects and their shader documents use canonical recipes.

The historical snapshot schema 3 stored one modulation recipe per colorizer
and a typed, versioned source clock. The current format is schema 4, described
in the migration contract above. Derived cubes and color tables are rebuilt after restore.
This guarantees scheduled-state replay; a snapshot does not serialize a partially
consumed display interval. Composed parameter schema versions each advance once;
complete layout pins include inactive source records and scalar types.

Initial admission is deliberately bounded: sphere sine, interference and smooth
grid use direct evaluation through axial frequency 64. Noise admits scales
1/64 through 64 directly for configured-single simplex with zero contrast;
configured-single and octave-contour noise admit scales 1/64 through 8 with the
baked strategy. Native and WebAssembly evaluation admit every source/mechanic
combination, including direct octave-contour noise and octave-contour lightness.
The four historical spatial sources share the bounded [-1, 1] scalar interface;
changing a mechanic does not change the source recipe. `Pullback::Modulation::valid` checks semantic
validity. `Pullback::Modulation::device_qualified` separately checks the measured
Teensy timing envelope, and `Pullback::Modulation::admitted` applies that additional
policy on Arduino builds, including disabled recipes.

### Scalar approximation

`pattern_qualification` separates signed-byte quantization from interpolation.
The qualification matrix uses nine formula/contrast recipes, nine scales from
1/64 through 8, three seeds, three phases and 4,096 directions, with additional
20,000-direction configured-single sweeps. Strict and optimized runs pass 43,650
and 199,125 checks respectively. Maximum measured quantization error is
0.003937006 (approximately 1/254). Across the full matrix, maximum spatial error
is 1.88681579 and maximum combined error is 1.88824224. These are empirical
sampled maxima, not universal bounds over every seed and direction.

At scale 8, the dense configured-single sweep measures maximum combined error
1.42269; at scales through 0.125 it remains below 0.005. A 24 by 24 cube is a
coarse field at high spatial frequencies. Callers requiring raw-source fidelity
must use an admitted direct recipe. Cube layout version 2 uses symmetric face
coordinates: sampled edge/corner disagreement stays below 2.8e-6 in strict math
and 1.68e-6 with optimized math, including the scale-64 seam diagnostic.

Run the native `pattern_qualification` target to reproduce the CSV. Source,
cache, domain and mechanic tests cover strict and optimized configurations.
A mutation restoring the former asymmetric cube arithmetic fails 30 seam checks.

### Render compatibility

The strict ShaderChain render pin is `8218c82b269c111d`; the recording configuration
(`-Os -DNDEBUG -ffast-math -fno-finite-math-only`) pins `6a3bf6437900ea91`.
The strict fingerprint changes only two total red codes and one total blue code
for the reviewed frame corpus. RGB geometry and alpha contracts remain covered
independently. The migrated 488-row frozen snapshot corpus preserves non-color
parameters, palettes, seeds and source clocks.

The fused colorizer uses the spatial hue table whenever hue modulation is active.
Native and WebAssembly direct-oracle comparisons retain the existing limits:
PEIRCE maximum channel errors 132 or less against 256 allowed; HUE maximum
errors 2,820 or less against 3,066 allowed. No oracle tolerance is widened.


### Resource ownership

Measured with the pinned ARM GCC 15.2.1 Cortex-M7 hard-float ABI and shipping
`-Os` flags:

| Resource | Bytes |
|---|---:|
| Passive canonical recipe | 152 |
| Scheduled source clock | 12 |
| Immutable prepared endpoint | 260 |
| Owned FastNoiseLite generator within a noise endpoint | 52 |
| Source variant reservation within an endpoint | 80 |
| Cache identity / sphere bake cache storage | 152 / 160 |
| Current direct resources, including endpoint and cache | 428 |
| Current baked resources, including one byte cube | 3,880 |
| Hue rotation table / owner with generation counter | 6,144 / 6,148 |
| One byte cube / optional external float cube shape | 3,456 / 13,824 |
| Optional external float hue-table shape | 12,288 |

A temporal baked effect owns two cubes (6,912 bytes), two hue tables (12,288
bytes), two immutable source endpoints and one current cache. Direct effects own
no cubes; a non-temporal baked effect owns one. Float mirrors are not allocated
by these composed runtimes. Views and interval copies borrow owned storage.
`Pullback::Modulation::retain` copies immutable endpoints and active byte cubes,
then rebinds each retained cube view to its owner?s storage.
Effect transitions destroy the outgoing effect before constructing the incoming
one, so the maximum is one effect's current and retained resources.

For the GridSpace stage chain at 288 by 144, temporal runtime State sizes are
1,660 bytes with no modulation, 14,640 direct and 21,544 baked. Their complete
persistent arena budgets are 6,920 / 19,996 / 26,900 bytes, including palette
storage, registered controls, retained palette, names and alignment. The matching
non-temporal State sizes are 472 / 7,048 / 10,500 bytes. Capability storage stays
allocated while its runtime recipe is disabled. A no-modulation specialization
retains passive control layout but owns no pattern source, cube or hue resources.

Native 64-bit pointers produce different totals: Prepared is 296 bytes and direct
and baked resources are 480 / 3,928 bytes. Resource accounting uses actual target
sizes; these figures do not estimate the call stack. Native arena/stack gates and
WebAssembly high-water checks independently validate transient storage.


Temporal preparation retains one additional column context in the persistent
arena (396 bytes ARM, 448 bytes native). `FrameInterval::at_into` writes into
this context without returning another full frame through the scan stack.
Non-temporal specializations allocate no column context. The strict-FP
WebAssembly composed-effect stack peak is 1,744 bytes, below the unchanged
2,048-byte creep gate and 8,192-byte physical stack.

### Device timing and admission

The [recipe qualification report](../profiles/shipping/profile_gridspace_palette_patterns_teensy_2026-10-29.md)
records eleven captures with exact source manifests, ELF hashes, preparation,
column and scan scopes, live ISR measurements, and complete-window frame counts.
They use the shipping selective-O3 image, segmented 288 by 144 rendering at
480 RPM, amount 0.8, noise scale 8 and analytic frequencies 64.

The driver reserves its existing two-revolution 250 ms commit window for
construction and the first frame while the display is black. Subsequent frames
have a 62.5 ms render window. Startup-only spills remain visible in the report;
they are not counted as steady deadline failures or hidden by averaging.
Configured-single noise/direct hue and lightness peak at 54.302 and 59.847 ms
after startup. Octave noise/baked hue and opacity peak at 52.066 and 46.035 ms.
Direct grid/lightness peaks at 49.339 ms; non-temporal interference/hue peaks at
32.588 ms. Source and mechanic kernel files remain byte-identical through the
final frame-storage optimization.

Direct three-octave noise/lightness exceeds the steady deadline in every
captured frame (352/352 complete-window frames). Baked octave noise/lightness
still exceeds it in 270/543 post-startup frames. Admission therefore rejects
direct octave-contour noise and octave-contour noise/lightness under either
strategy on Teensy. This timing policy does not restrict native or WebAssembly
authoring. Disabled recipes undergo the same target admission. No memory, stack,
render-oracle or driver timing budget is raised.


A final capture after the retained-column change validates another 1,088 frames
of direct configured-single noise/lightness: 70.385 ms at startup and 59.831 ms
thereafter, with no steady spills. The timing report retains both campaign and
final-image provenance.


### Firmware memory

All ten PlatformIO environments pass their existing size/layout gates using
pinned Python 3.12.13, PlatformIO 6.1.19 and ARM GCC 15.2.1. Project translation
units additionally compile with `-Werror`; no first-party warnings are emitted.
This is incremental validation, not a cold compiler-warning census. Unchanged
vendored framework warnings remain outside the first-party policy.

| Phantasm region | Baseline `7e08beaa2` | Final storage implementation `0448cff3c` | Delta |
|---|---:|---:|---:|
| RAM1 code | 161,080 | 164,296 | +3,216 |
| RAM1 variables | 312,736 | 312,736 | 0 |
| FLASH data | 853,952 | 861,212 | +7,260 |

The code crosses an ITCM allocation-bank boundary. RAM1 space for locals is
14,944 bytes, above the unchanged 12,288-byte gate by 2,656 bytes; RAM2 free
space remains 4,224 bytes. The input tree, final core manifest, ELF identities,
complete region table and raw build output were retained with the qualification
artifacts. The final input-tree SHA-256 is
`716cf11239aaa3bae8122cfbc0db542d723612f226d69400895342884f8cd04a`.


### Final validation

* Native default and excluded targets build; all 125 CTests pass on the final
  implementation. The ShaderChain module passes 5,172,544 assertions.
* Release and strict-FP WebAssembly smoke tests each pass all 44 effects at
  both supported resolutions, including snapshot refusal/replay, stack, arena,
  parameter bindings and generated-document checks.
* JavaScript tooling passes 252 tests and 1,468 assertions across 16 nonempty
  files; one optional browser test is skipped.
* Include-direction, standalone headers, schema-layout mutations, source/cache
  cases, 498 pinned fail-fast cases and the guard census pass. The full native
  suite includes frozen-capture determinism and the unchanged render oracles.
* Assertions-off builds, optimized pattern/color tests, ten firmware environments,
  twelve segmented-driver captures and final ARM/native resource measurements
  retain their compiler flags, input hashes and raw evidence.

No compatibility modifier headers remain. GS retains the corrected symmetric
cube sampling; its reviewed FULL strict and fast render pins are updated, with
no widened color, gamut, alpha or determinism bounds.
