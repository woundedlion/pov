# Pattern sources and palette modulation

**Status: PROPOSED.** Specifies future behavior, not implemented controls.

## Problem and scope

The spatial palette noise modifier already has the coordinates needed for
analytic patterns. It samples a scalar field over sphere directions and uses
that scalar to rotate palette hue. Generalize both the field source and the
operation applied to color, independently.

The abstraction has two axes: **pattern source** and **modulation mechanic**.
Noise plus hue rotation is the existing noise modifier; noise plus positive
lightness lift is shimmer. Grid plus opacity uses the same protocol.
Keep ordinary palette lookup, `get(float)`, unchanged.

The first implementation delivers one generic modulation slot with noise,
directional sine, interference and smooth grid sources, and hue rotation,
positive lightness lift and opacity mechanics, for composed pullback effects
and the shader workbench. It does not add a modifier graph, user-entered
expressions, projected coordinate domains, palette-coordinate displacement or
multiple-mechanic ordering. Existing wrapper consumers remain compatible;
migrating every effect is unnecessary.

## Existing behavior

[The spatial palette utilities](../../core/color/noise_hue_palette.h) own two
different approximations: a palette-coordinate by hue-rotation color table,
and a six-face spatial noise table. The latter contains 24 by 24 signed-byte
samples per face. Its sampler uses dominant-axis ratios, so a direction need
not have unit length. Positive scaling leaves its mathematical lookup
coordinates unchanged. The bake evaluates noise at normalized face directions,
scaled and translated around the existing periodic noise loop.

[The pullback colorizer](../../core/render/pullback/color.h) samples this field
at `FieldSample::sphere`. This is already spatial input, distinct from the
one-dimensional coordinate modifiers in
[palette composition](../../core/color/composition.h).

[The carrier contract](../../core/render/pullback/contract.h) defines the frame:
PLANE sampling retains the projection-frame direction; SPHERE sampling retains
the camera-frame direction. Plane endomorphisms do not rewrite the sphere point.
Inherit these semantics without another transform; this is not universal world space.

[Shimmer](../../core/color/noise_shimmer_palette.h) already shares the spatial
table but consumes only positive noise lobes to lift OKLab lightness toward
white, preserving alpha. This is a special case of the generic source/mechanic
pair, not a separate permanent shimmer subsystem.

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
that type; mechanics remain domain-agnostic. Do not pass a universal context bag,
force every domain into Vector3, or switch coordinate domains per pixel. Ordinary
palettes still map a scalar lookup coordinate to color. The palette-coordinate
by hue-rotation LUT is internal acceleration, not a user-facing 2D palette.

Initial delivery supplies the sphere specialization; other types fit the protocol
but their renderer plumbing and presets are deferred. Its public compatibility
adapter accepts existing finite, nonzero vectors. Direct evaluation normalizes
before calling the pattern; cube baking supplies unit directions. Cube lookup
may elide normalization using its existing scale-invariant ratios, including
legacy fast paths. Neither interprets length as radius or normalizes Position3D.
Test scale invariance with float tolerance; `FieldSample` gains no unit requirement.

Prepare validates finite parameters, normalizes axes and rejects zero axes.
Custom sources must satisfy the output contract; debug assertions diagnose
nonfinite output, and the shared scalar boundary clamps finite numerical
overshoot. Clamping is not a way to conceal an incorrectly normalized custom
pattern. Sampling performs no allocation, random-state advancement, mutation,
I/O or clock reads. All time and seed inputs are explicit prepared values.
Custom sample cost must be bounded; unbounded iteration is outside the protocol.

## Generic protocol and initial patterns

Use a small policy protocol, consistent with the existing prepare/sample split
in [pullback sources](../../core/render/pullback/source.h). The following types
are schematic, not existing declarations:

```cpp
struct PatternClock {
  float phase_turns;
};

struct SphereDirection { math::Vector unit; };

struct DirectionalSine {
  using Coordinate = SphereDirection;
  struct Params { math::Vector axis; float frequency; float offset_turns; };
  struct Prepared { math::Vector axis; float frequency; float phase_turns; };
  static Prepared prepare(const Params&, PatternClock);
  static float sample(const Coordinate&, const Prepared&);
};

template<class Pattern> struct DirectSphereField;
template<class Pattern> struct BakedSphereField;

struct Opacity {
  struct Params { float amount; };
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

The field adapter exposes prepare/sample, wrapping either direct evaluation or
cube storage. The mechanic exposes prepare/apply and returns a valid straight-alpha
`Color4`. The generic wrapper samples its base palette, samples the field, then
applies the mechanic. HueRotation and PositiveLightnessLift follow the same
mechanic protocol as Opacity. Specialized palette-coordinate/hue LUT kernels may
fuse these operations for speed while preserving the observable contract.
PatternHuePalette and PatternShimmerPalette are not primary abstractions.

Provide four initial source policies:

| Source | Definition and controls |
|---|---|
| Noise | Adapter around the current configured noise generator, spatial scale and sphere-loop construction. Seed and generator configuration belong to source identity. |
| Directional sine | `sin(2*pi*(frequency*dot(axis,direction) + phase + offset))`; unit axis, nonnegative frequency and offset in turns. Frequency is cycles per unit axial coordinate, so the sphere spans twice that many cycles along the axis. |
| Two-wave interference | Half the sum of two directional sines, with separate axes, frequencies and offsets. Both use the same loop phase; the fixed half gain bounds the result without dynamic normalization. |
| Smooth grid | Product of two directional sines over orthogonal unit axes, with independent frequencies and offsets and a shared loop phase. Validate axis orthogonality. The product remains in [-1, 1]. This is a smooth grid over sphere directions, not a planar projection or a precomputed opacity mask. |

Phase is in turns and speed is turns per engine update, matching current
palette noise timing. Advance it once per update, wrap to [0, 1), and preserve
negative-speed behavior. Do not substitute wall-clock seconds. The analytic
sources return to the same field after one phase turn. Source preparation
hoists frame-constant work outside scanning. Admission qualifies upper frequency
and argument bounds compatible with the chosen trigonometry and the project's
fast-trig limits. Fold phase and keep arguments within that range; this proposal
does not invent a universal cutoff or promise detail free of spatial aliasing.

An effect may supply any custom C++ policy satisfying this protocol. Compile-time
instantiation provides extensibility; this is not an expression interpreter or
an unrestricted runtime callback registry. Built-ins may have a small tagged
selector, resolved during preparation into specialized sampling paths. Custom
sources need a stable recipe identifier and serializable parameter schema to
participate in workbench replay; a function pointer is not a snapshot recipe.
Custom mechanics follow the same bounded, immutable discipline and require
their own stable replay recipe and schema.

## Mechanics and valid combinations

| Pattern source | Hue rotation | Positive lightness lift | Opacity |
|---|---|---|---|
| Noise | Existing noise hue modifier | Existing shimmer | Organic alpha mask |
| Directional sine | Hue bands | Lightness bands | Smooth alpha bands |
| Smooth grid | Hue grid | Lightness grid | Grid alpha mask |

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
**Amount zero is identity for every mechanic.** No mechanic changes the source's
range or asks it to return a target-specific value.

The colorizer order is palette sampling, mechanic application, existing RGB
brightness envelope, then alpha multiplication by coverage and value opacity.
The result remains straight alpha. Final premultiplication stays in
[the scan shader](../../core/render/scan/shader.h), after the chain; opacity
must not multiply RGB early. Legacy noise hue remains an adapter into these
internals, retaining its specialized kernels and existing ordering.

## Direct and baked sphere evaluation

For the initial SphereDirection specialization, evaluation strategy is explicit
in the effect or workbench recipe. Normalization here applies only to sphere
directions; a future Position3D implementation must preserve radius.

* **Direct:** normalize the input and evaluate the prepared source per sample.
  This preserves narrow or high-frequency features subject to scan resolution
  and floating-point math, and incurs per-sample arithmetic.
* **Baked cube:** evaluate the source into the existing six-face layout, clamp
  and round to signed bytes with the existing quantization, then bilinearly
  sample it. The cube table becomes reusable scalar storage. Preserve old
  noise names as compatibility aliases or wrappers during extraction.

There is no automatic strategy switch based on frequency, load or source kind.
The same recipe must select the same path on native, WASM and firmware. The
legacy noise path remains baked. New analytic presets initially choose direct
evaluation; baked analytic recipes are opt-in and need an error qualification.

Cube storage cannot reproduce arbitrary frequencies or sharp boundaries.
Quantization contributes at most 1/254 scalar error at a bake point; spatial
interpolation error depends on the source and is additional. Each baked recipe
declares an admitted parameter range and measured approximation bounds. Reject
unqualified baked configurations rather than quietly blurring them or switching
paths. Reuse cube face conventions and duplicated edge evaluation so new
sources do not introduce seams.

## Ownership, cache invalidation and replay

Prepared state and tables are owned by the effect/runtime arena; views borrow
them through the draw and any retained interval. Direct endpoint state must
contain its complete immutable source inputs, including any generator handle
whose resource stays immutable through that interval. Own/copy direct generator
configuration, or borrow an explicitly immutable retained generator; a lifetime-valid
pointer to a subsequently mutated generator is insufficient. Never retain a pointer into
mutable live controls or overwrite the previous table before its final read.

The current noise cache keys scale and phase under a documented precondition
that its generator configuration and destination remain unchanged. Preserve
that optimized legacy contract. A general field cache additionally keys source
kind/recipe version, all source parameters, seed and generator configuration,
phase, and bake layout/version. A configuration generation may replace the
configuration portion only if every mutation increments it. Reset cache validity
on ownership changes and restore. Mechanic kind, amount and other mechanic
controls, and palette changes, do not invalidate the spatial table. Any
mechanic/color cache has a separate identity: palette changes invalidate the
hue-rotation color table; spatial pattern changes do not.

Snapshots store the source recipe, parameters, explicit evaluation strategy,
seed where applicable, current loop phase, and the independent mechanic kind
and parameters. Derived tables and pointers are
not serialized. Restoring rebuilds derived state deterministically before
scanning; invalid or unavailable recipes fail explicitly. Exact continuation of
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

Preserve the existing noise interpolation semantics. The colorizer reads the
previous and current spatial tables at the same column direction, using the
column's shared interpolated palette mapping and hue amount, samples each
endpoint's hue-rotation/palette resources, then interpolates endpoint RGB.
Alpha follows the existing temporal alpha path. It does **not** interpolate
the scalar first and hue-rotate once, nor evaluate independent complete old and
new color-control sets. Preserve existing endpoint selection, Q16 weights,
rounding and stationary-resource fast paths for legacy noise.

For a new pattern interval with unchanged source kind, strategy and mechanic,
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
Preserve legacy noise alpha fast paths where their alpha-preserving premise holds.
Existing specialized noise kernels retain their current float/rounding order;
algebraic rearrangement of those kernels is not assumed bit-identical.

A source kind, evaluation strategy, mechanic kind or custom recipe-version change is
a topology change, not a numeric parameter blend. For the initial pattern
implementation, hold the **current** source and mechanic selections and prepared
states/table through that interval, including its fraction-zero boundary.
Hold the mechanic's current typed parameters too: hue turns must never be
interpolated as opacity amount. Same-kind continuous controls interpolate normally.
Apply that mechanic to both palette endpoints and interpolate their modulated
colors after shared downstream controls, as above. The switch may therefore be
discontinuous; expose and test that policy rather than claiming cross-pattern
morphing. Ordinary unchanged-topology intervals retain exact endpoint behavior.
The fraction-zero result uses the previous palette with the held current
modulation, including its alpha effect; it need not equal the previous fully
modulated frame. Never reinterpret an old source's or mechanic's resources
under the new type.
Do not change legacy noise topology behavior as part of this feature.

## Compatibility and delivery

Use an independent modulation descriptor containing source selector/parameters,
evaluation strategy and mechanic selector/parameters. It is the canonical
generic representation, not a new non-hue operation hidden in a hue-mode enum.
Keep `HueMode::NONE = 0`, `NOISE = 1` and `PATH_LENGTH = 2`, including the
existing displacement alias. Decode old NOISE data as noise plus hue rotation;
no descriptor means existing behavior. Keep legacy field IDs such as `hue-noise-scale` and
`hue-noise-speed`, defaults, seed initialization, phase advancement, noise-loop
construction and optimized noise bake unchanged. PATH_LENGTH remains its
existing independent hue behavior; initially reject configurations requesting
it simultaneously with the generic slot. Do not silently stack mechanics or
reinterpret path length as a new pattern. Preserve decoding and
encoding of the existing color-clock-v1 snapshot; add versioned pattern state
instead of altering its meaning. Append new IDs and update activation gates,
parameter metadata and workbench controls together.

The existing modifiers must fall out of the generic implementation. Schematically:

```cpp
template<class PaletteSource>
using NoiseHueImplementation = PatternModulatedPalette<
    PaletteSource, BakedSphereField<Noise>, HueRotation>;
template<class PaletteSource>
using NoiseShimmerImplementation = PatternModulatedPalette<
    PaletteSource, BakedSphereField<Noise>, PositiveLightnessLift>;
```

NoiseHuePalette and NoiseShimmerPalette preserve their public names and bind/LUT
APIs as wrappers around these specializations. Their existing table arguments
bind the generic baked field view. Consolidate implementation, retaining
specialized fast paths, rather than maintaining duplicate independent engines.
Preset names need no migration. Existing one-dimensional NoiseWarpModifier and
RippleModifier remap the coordinate before palette lookup; they are not Color4
mechanics and need no forced merger. A scalar-domain pattern may later drive
the same hue/lightness/opacity mechanics through the generic protocol.

Source and mechanic selectors have independent stable identifiers; custom
policies are compiled in and selected through registered serializable recipes,
never persisted callback pointers.

Implement in stages: extract reusable scalar storage behind compatibility
entry points; add the generic source/mechanic protocol and all three mechanics;
then wire metadata, snapshots and workbench support. Exercise noise, sine,
interference and grid through that slot. No effect preset migration is required.

## Budgets and acceptance

One current signed-byte cube costs 3,456 bytes; retaining both endpoints costs
6,912 bytes. Optional float copies cost another 13,824 bytes per endpoint.
Reuse existing tables for migrated noise users. Direct policies allocate no
cube, and effects with no pattern allocate no pattern resources. Do not enable
float mirrors or instantiate every built-in for every effect by default.

Implementation must report incremental RAM1 code, RAM1 variables and FLASH
data against the baseline, passing the repository's current firmware gates.
For new direct presets, measure frame preparation and scan/column timing on
Teensy under the segmented driver; establish an admitted frequency/configuration
range that stays inside its timing budget. No fixed cycle estimate is claimed
by this proposal. Baking performs 3,456 source evaluations when its key changes,
and none when the complete key is unchanged. Scan-time allocation is forbidden.

Acceptance requires:

* Unchanged legacy noise table bytes and existing render goldens across seeds,
  positive/negative speeds, zero amount and temporal paths.
* Independent known-value tests for sine, interference and grid, output bounds,
  periodic phase, positive direction-scale invariance and invalid axes.
* Cross-axis tests using one source with all mechanics and one mechanic with
  different sources. Amount zero is identity throughout; scalar zero follows
  each mechanic's definition. Verify shimmer's existing lightness/alpha behavior.
* Compile-time domain checks rejecting mismatched coordinate types and implicit
  Position3D-to-SphereDirection conversion. Future volumetric fields must not
  normalize away radius; these checks require no additional renderer plumbing.
* Asymmetric orientation cases proving projection-frame versus camera-frame
  provenance, and cube edge/corner seam tests from both adjacent faces.
* Direct-versus-baked error measurements across each admitted baked range,
  separating quantization from spatial error; no blanket arbitrary-pattern
  equivalence claim.
* Temporal tests whose expected endpoint RGB differs from scalar-first blending,
  including changed phase/scale, stationary resources and topology holds.
  Verify opacity's modulated alpha at both endpoints and interior fractions,
  including equal base alphas but unequal masks; scan premultiplies exactly once.
* Cache tests for every key component, palette-only and mechanic-only changes,
  with no spatial rebake when only the mechanic changes; snapshot replay
  after cache reset, old snapshot compatibility, and unknown-recipe rejection.
* Native and WASM validation plus firmware memory/timing evidence for the
  implementation. This proposed-spec change itself needs documentation checks
  only.
