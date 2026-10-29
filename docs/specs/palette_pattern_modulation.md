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

Initial delivery supplies the sphere specialization; other renderer plumbing
and presets are deferred. The generic boundary admits vectors whose computed
float squared length is finite and at least `math::EPS_NORMALIZE_SQ` (currently
1e-12), matching [Vector normalization](../../core/math/3dmath.h). Explicitly
check finiteness before normalization; the existing lower-bound check alone
does not reject infinity. Apply this admission to directions and source axes.
Construct the unit SphereDirection tag only at this validated boundary.

Normalize once per sample and share that coordinate across retained endpoints.
Cube baking supplies unit directions; cube lookup may elide normalization using
scale-invariant ratios. The legacy raw cube path keeps its existing precondition
of finite direction with positive float squared length; the generic boundary
is deliberately narrower. Test scale invariance only while both vectors remain
admitted. Neither path interprets direction length as radius or normalizes
Position3D; `FieldSample` gains no unit-length requirement.

Prepare validates finite parameters and admitted axes.
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

struct SphereDirection; // Constructed by the validated domain adapter.

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
mechanic protocol as Opacity. Nonzero hue/color LUT fusion is a qualified
approximation of this reference, not bitwise equivalence. New generic accelerated
routes must still honor the exact identity rules below.
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
For the generic reference and new generic accelerated routes, **amount zero is
exact Color4 identity for every mechanic**. A zero resolved hue rotation,
including zero after periodic wrapping, returns the exact base palette sample;
bypass the hue/color LUT. Positive lightness lift also returns its exact input
when lift is zero. These identities apply before downstream brightness/coverage.
No mechanic changes the source range or asks for a target-specific value.

Legacy active LUT routes may resample palette color even at zero hue shift,
because their palette axis is quantized. Compatibility adapters preserve that
existing resampling, rounding and output as a legacy numeric profile within
the common generic engine, rather than claiming exact identity or creating a
second engine. New explicit descriptors use the new generic identity contract.

Zero amount may skip field evaluation and preparing unused mechanic/color
resources, while source clocks continue. Keep any still-needed retained endpoint
resources valid. Configuration validation remains mandatory even when skipped.

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

Snapshots store the enabled state, source recipe/parameters, evaluation strategy,
seed where applicable, generic loop phase, independent mechanic kind/parameters
and retained dormant legacy phase. Derived tables and pointers are not serialized.
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

Preserve the existing noise interpolation semantics. The colorizer reads the
previous and current spatial tables at the same column direction, using the
column's shared interpolated palette mapping and hue amount, samples each
endpoint's hue-rotation/palette resources, then interpolates endpoint RGB.
Alpha follows the existing temporal alpha path. It does **not** interpolate
the scalar first and hue-rotate once, nor evaluate independent complete old and
new color-control sets. Preserve existing endpoint selection, Q16 weights,
rounding and stationary-resource fast paths for legacy noise.

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
Preserve legacy noise alpha fast paths where their alpha-preserving premise holds.
Existing specialized noise kernels retain their current float/rounding order;
algebraic rearrangement of those kernels is not assumed bit-identical.

A source kind/custom recipe version, evaluation strategy, mechanic kind,
enabled state, effective descriptor ownership/presence or numeric-profile change
is a topology change, even when other selectors match. Hold the **entire current
effective modulation**: ownership, numeric profile, source parameters and
clock-selected state/resources, mechanic typed parameters/state and enabled
state, including fraction zero. This applies to both descriptor insertion and
removal. Continuous mechanic controls interpolate only when the entire topology
is unchanged, even if the mechanic kind alone matches. For example, sine-to-grid with opacity amount
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
Intervals wholly within legacy ownership retain existing noise topology behavior;
transitions across descriptor ownership follow the current-effective-owner policy above.

## Compatibility and delivery

Use an independent modulation descriptor with source selector/parameters,
evaluation strategy, mechanic selector/parameters and enabled flag (default true).
Keep `HueMode::NONE = 0`, `NOISE = 1`, `PATH_LENGTH = 2` and its displacement alias.
Resolve slot ownership before preparing resources:

| Explicit descriptor | Legacy hue mode | Effective behavior |
|---|---|---|
| Absent | NONE | Existing no-modulation behavior. |
| Absent | NOISE | Derive noise plus hue rotation with the legacy numeric profile. |
| Absent | PATH_LENGTH | Existing independent path-length hue behavior. |
| Present | NONE or NOISE | Descriptor solely owns the slot; legacy noise controls are dormant. |
| Present | PATH_LENGTH | Reject the configuration. |

A disabled explicit descriptor is identity; it never falls back to legacy NOISE.
Validate recipe and parameters regardless of enabled state or amount. Dormant
legacy controls cannot stack, drive generic resources or advance their clock.
Keep legacy IDs, defaults, seed setup, noise loop and optimized bake unchanged.
Encode old descriptor-absent data in its unchanged legacy shape. Encode an
explicit descriptor once; any preserved legacy fields are explicitly inactive,
not a second effective representation. Preserve color-clock-v1 meaning and add
versioned generic state. Append new IDs and update gates, metadata and controls.

An explicit descriptor's generic phase starts at zero on creation, stays
normalized to [0, 1), and is retained across source, strategy and mechanic changes.
With that descriptor present, advance once per engine update at its source speed even while
amount is zero or enabled is false. Source changes use the new speed thereafter,
without resetting phase. Disabled modulation may skip sampling/bakes; reenable
uses the current phase and restores valid resources before use. Restore a
snapshot's phase rather than reinitializing it.

Freeze the legacy hue-noise phase while a descriptor is present. Removing it
resumes the retained legacy phase if the old mode is NOISE; never advance both
clocks. Palette-mapping oscillation clocks remain separate and unaffected.

The existing modifiers must fall out of the generic implementation. Schematically:

```cpp
template<class PaletteSource>
using NoiseHueImplementation = PatternModulatedPalette<
    PaletteSource, BakedSphereField<Noise>, HueRotation<LegacyNoiseNumerics>>;
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

* Unchanged legacy noise table bytes and existing render goldens across seeds,
  positive/negative speeds, zero amount/shift and temporal paths, including legacy
  active-LUT resampling. Separately assert exact generic Color4 identity at zero
  amount and resolved hue rotation; qualified nonzero approximations need bounds.
* Independent known-value tests for sine, interference and grid, output bounds,
  periodic phase, direction-scale invariance within admission, rejected below-floor
  or nonfinite squared lengths, and invalid axes.
* Cross-axis tests using one source with all mechanics and one mechanic with
  different sources. Generic zero amount and mechanic-specific scalar-zero rules
  hold; legacy numeric profiles preserve existing output. Verify shimmer behavior.
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
  Test source changes with same-kind amount changes, mechanic/strategy changes,
  enabled transitions and their fraction-zero current-state holds.
  Insert and remove descriptors with matching Noise/Baked/Hue selectors but
  different phases; include scalar zero to prove the current numeric profile
  governs fraction zero, preserving either legacy resampling or generic identity.
  Verify opacity's modulated alpha at both endpoints and interior fractions,
  including equal base alphas but unequal masks; scan premultiplies exactly once.
* Cache tests for every key component, palette-only and mechanic-only changes,
  with no spatial rebake when only the mechanic changes; snapshot replay
  after cache reset, old snapshot compatibility, and unknown-recipe rejection.
* Mixed descriptor/legacy decode, precedence, disabling/reenabling and removal
  tests; no fallback or dual clock advance, phase continuity through source
  switches, snapshot phase restoration and validation of disabled recipes.
* Native and WASM validation plus firmware memory/timing evidence for the
  implementation. This proposed-spec change itself needs documentation checks
  only.
