# Pattern sources and palette modulation

**Status: IMPLEMENTED.** Shared sources, independent mechanics, temporal endpoints, codecs and migrated callers are validated on native, WebAssembly and Teensy targets.

## Problem and scope

Generalize the existing spatial palette noise modifier into independent pattern
sources and color operations, using its existing sphere coordinates.

The abstraction has two axes: **pattern source** and **modulation mechanic**.
Noise plus hue rotation is the existing noise modifier; noise plus positive
lightness lift is shimmer. Grid plus opacity uses the same protocol.
Keep ordinary palette lookup, `get(float)`, unchanged.

First delivery provides one slot: noise/sine/interference/grid sources and
hue/lightness/opacity mechanics, for composed pullback and the workbench.
Graphs, expression parsing, additional modulation coordinate bindings, remapping and
multiple-mechanic ordering are deferred. Migrate existing callers directly,
with no compatibility wrappers, fallback configuration or parallel implementation.

## Existing behavior

[The scalar sphere storage](../../core/pattern/sphere_lut.h) contains six
24 by 24 signed-byte faces. The independent
[hue table](../../core/color/hue_rotation_lut.h) stores palette-by-hue colors. Dominant-axis ratios make lookup invariant to positive
direction scaling. Baking evaluates normalized directions through the existing
periodic noise loop.

[The pullback colorizer](../../core/render/pullback/color.h) samples this field
at `FieldSample::sphere`; the needed spatial input already exists.

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
The initial modulation binding still supplies only SphereDirection: FieldSample
does not retain warped plane coordinates. Supporting a planar modulation source
later needs an explicit coordinate binding, not a silent plane/sphere conversion
or universal context bag.

## Scalar and typed coordinate contracts

A pattern takes its typed coordinate and immutable prepared state and returns
a finite scalar in [-1, 1]. Pattern output is independent of palette contents
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
incorrect source normalization. Sampling is immutable and bounded, with no
allocation, random-state advancement, I/O or clock reads. Time and seed are explicit.

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

The modulation binding initially admits these four sphere sources alongside
the existing planar/sphere catalogue extracted for pullback:

| Source | Definition and controls |
|---|---|
| Noise | Adapter around the current configured noise generator, spatial scale and sphere-loop construction. Seed and generator configuration belong to source identity. |
| Directional sine | `sin(2*pi*(frequency*dot(axis,direction) + phase + offset))`; unit axis, nonnegative frequency and offset in turns. Frequency is cycles per unit axial coordinate, so the sphere spans twice that many cycles along the axis. |
| Sphere interference | Half the sum of two directional sines, with separate axes, frequencies and offsets. Both use the same loop phase; the fixed half gain bounds the result without dynamic normalization. |
| SphereSmoothGrid | Product of two directional sines over independently normalized axes, with independent frequencies and offsets and a shared loop phase. Admit finite, nondegenerate axes at any angle, including parallel axes. The product remains in [-1, 1]; this is not the existing planar coupled-grid algorithm. |

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

These combinations, and interference with any mechanic, are structurally valid.
Admit each recipe only after qualifying its resources, approximation error and
execution cost. Do not create an enum value per source/mechanic pair.

Hue rotation resolves `hue_turns = amount * scalar`, with signed amount in
turns, periodic hue lookup and existing gamut mapping. It preserves alpha.
PositiveLightnessLift resolves `lift = amount * max(0, scalar)` for amount in
[0, 1], then moves OKLab lightness toward white with the existing gamut-aware
chroma reduction. It preserves alpha; negative lobes do not darken. These two
mechanics leave color unchanged at scalar zero.

Opacity resolves `u = (scalar + 1) / 2`, then
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

The colorizer order is palette sampling, combined path/spatial hue rotation,
non-hue mechanic application, existing RGB brightness envelope, then alpha
multiplication by coverage and value opacity.
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

For modulation use endpoint-color interpolation, not scalar-first blending.
Pullback may retain its existing source-phase interpolation policy: consumer
temporal policies can differ while identical prepared source inputs produce
identical raw samples. Review color numeric changes instead of retaining a
parallel engine.

For a new pattern interval with entirely unchanged modulation topology,
retain the two prepared source states or two baked tables. Evaluate each at
the column direction. Sample each endpoint palette with the shared column
mapping, apply the mechanic using shared column-interpolated continuous mechanic
controls, then apply shared RGB brightness and alpha coverage/value opacity to
each endpoint. Interpolate these completed straight-alpha endpoint colors,
including **both RGB and alpha**. Do not additionally
interpolate source phase/parameters inside those endpoint evaluations. This
avoids double interpolation and aligns direct and baked resource semantics.
Do not interpolate the scalar first and apply the mechanic once. In particular,
equal base palette alpha does not imply equal modulated alpha: the common-alpha
fast path must be disabled for opacity unless equality is proven after modulation.
Alpha-preserving mechanics may use an equality fast path where that premise is
proven. Pin interpolation weights and rounding with tests; do not assume
algebraic rearrangements are bit-identical.

A source kind/custom recipe version, evaluation strategy, mechanic kind or
enabled-state change is a topology change. Hold the **entire current modulation**:
source parameters and current prepared state/resources (initialized on a source
kind/version switch), mechanic typed parameters/state
and enabled state, including fraction zero. Continuous mechanic controls interpolate
only when the entire topology is unchanged, even if the mechanic kind alone matches.
For example, sine-to-grid with opacity amount
0.2-to-0.8 uses the current grid and amount 0.8 throughout the switch interval.
Hue turns are never interpolated as opacity amount.
Apply that mechanic to both palette endpoints and interpolate their modulated
colors after shared downstream controls, as above. The switch may therefore be
discontinuous; expose and test that policy rather than claiming cross-pattern
morphing. Ordinary unchanged-topology intervals retain exact endpoint behavior.
The fraction-zero result uses the previous palette with the held current
modulation, including its alpha effect; it need not equal the previous fully
modulated frame. Never reinterpret an old source's or mechanic's resources
under the new type.
Disabling therefore holds identity across that interval; enabling holds the new
prepared modulation.

## Canonical configuration and migration

Use one modulation configuration containing source selector/parameters,
evaluation strategy, mechanic selector/parameters and enabled flag (default true).
It is the sole representation. Disabled modulation is identity. Validate recipe
and parameters regardless of enabled state or amount.

Each configured source instance has one state owner: the renderer stores and
schedules it using the source's pure initialization/advance operations. Advance
once per update even when modulation amount is zero or enabled is false. Changing
mechanic, strategy or same-source parameters retains State; new rates apply on
the next advance. Source kind/version changes initialize that source's State,
never reinterpret the old type. Snapshot restore restores typed state. Disabled
modulation may skip samples/bakes; reenable rebuilds needed resources from current
state. Palette-mapping clocks remain separate; no dormant source state is kept.

Keep the existing PATH_LENGTH hue feature explicitly configured in the renderer,
outside this spatial slot. Path-length and spatial hue rotations add in turns
before one wrapped hue lookup; a zero combined rotation returns the exact base
palette sample. Apply lightness or opacity modulation after hue rotation, before
the brightness envelope and coverage. Both features may be enabled together.
Remove old noise-mode ownership machinery; do not invent a path-length pattern
or duplicate hue API to retain it.

Migrate all noise-hue and shimmer callers, includes, controls, snapshots and
tests directly. Noise plus HueRotation and noise plus PositiveLightnessLift
are configurations of the generic implementation, not named compatibility
wrappers. Remove old palette headers after their final callers migrate and update
documentation references to the replacement modules.
Reuse existing noise-loop, generator and gamut math; review intentional render
changes from exact zero identity and unified rounding. Bump snapshot/preset
schema versions where needed and reject unsupported old versions explicitly;
do not emulate old wire shapes or maintain dual representations.

Source and mechanic selectors have independent stable identifiers in the new
schema. Custom policies are compiled in and replayed through registered recipes,
never persisted callback pointers. Existing one-dimensional NoiseWarpModifier
and RippleModifier remain prelookup remaps outside this Color4 mechanic contract.

## Code organization

Use ten cohesive header-only units. Templates and small kernels match the
repository's existing style; no per-platform backend copies or source/mechanic
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
new boundaries; the current broad layer gate does not enforce them. Add standalone
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
  PATH_LENGTH and generic modulation, including additive hue cancellation and
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


## Implementation qualification (2026-10-29)

The canonical records and pure sources live in `core/pattern`; independent
mechanics and the generic palette adapter live in `core/color`. Composed effects
and the workbench share the scheduled source clock, immutable preparation and
completed-endpoint interpolation. The former noise/shimmer modifier headers are
removed. The 19 composed effects and their shader documents use canonical recipes.

Snapshot schema 3 stores one modulation recipe per colorizer and a typed,
versioned source clock. Derived cubes and color tables are rebuilt after restore.
This guarantees scheduled-state replay; a snapshot does not serialize a partially
consumed display interval. Composed parameter schema versions each advance once;
complete layout pins include inactive source records and scalar types.

Initial admission is deliberately bounded: sphere sine, interference and smooth
grid use direct evaluation through axial frequency 64. Noise admits scales
1/64 through 64 directly for configured-single simplex with zero contrast;
configured-single and octave-contour noise admit scales 1/64 through 8 with the
baked strategy. Native and WebAssembly evaluation admit every source/mechanic
combination, including direct octave-contour noise and octave-contour lightness.
All source samples share the bounded [-1, 1] scalar interface; changing a mechanic
does not change the source recipe. `Pullback::Modulation::valid` checks semantic
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
