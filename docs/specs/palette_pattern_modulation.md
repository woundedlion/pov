# Spatial palette pattern modulation

**Status: PROPOSED.** This document specifies a future engine change. It does
not describe implemented pattern controls.

## Problem and scope

The spatial palette noise modifier already has the coordinates needed for
analytic patterns. It samples a scalar field over sphere directions and uses
that scalar to rotate palette hue. Generalize the field source so waves,
interference and custom C++ functions can drive the same operation.

The abstraction is **spatial scalar pattern plus color target**. Noise is one
source. Analytic evaluation is one way to produce the scalar, rather than a
different kind of palette. Keep ordinary palette lookup, `get(float)`, unchanged.

The first implementation adds pattern-driven hue rotation to composed pullback
effects and the shader workbench. It also exposes the scalar sampler so the
existing shimmer wrapper can consume it with its separate lightness policy.
It does not add a general modifier graph, user-entered expressions, projected
coordinate domains, palette-coordinate displacement, opacity modulation or
arbitrary ordering of multiple color targets. Existing raymarch and other
wrapper consumers remain compatible; migrating every effect is unnecessary.

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

The direction's frame is defined by
[the carrier contract](../../core/render/pullback/contract.h): a field sampled
from PLANE retains the projection-frame direction; a field sampled directly
from SPHERE retains the camera-frame direction. Plane endomorphisms do not
rewrite the retained sphere point. It is not universally a world-space point.
The new modifier inherits these semantics without adding another transform.

[Shimmer](../../core/color/noise_shimmer_palette.h) already shares the spatial
table but consumes only positive noise lobes to lift OKLab lightness toward
white, preserving alpha. This is a different target from signed hue rotation.

## Scalar and domain contracts

A pattern takes a **unit sphere direction** and immutable prepared state and
returns a finite scalar in [-1, 1]. Zero means no hue shift; positive and
negative values rotate in opposite directions. Pattern output is independent
of the palette coordinate, hue amount, palette contents and brightness gain.

The public sampling adapter accepts the same finite, nonzero directions as the
existing cube sampler. Direct evaluation normalizes once before calling the
pattern. Cube baking supplies normalized directions; cube lookup keeps its
ratio-based implementation and needs no normalization. Neither path interprets
input length as radius. Tests allow floating-point tolerance when checking
scale invariance. No new unit-length requirement is imposed on `FieldSample`.

Prepare validates finite parameters, normalizes axes and rejects zero axes.
Custom sources must satisfy the output contract; debug assertions diagnose
nonfinite output, and the shared scalar boundary clamps finite numerical
overshoot. Clamping is not a way to conceal an incorrectly normalized custom
pattern. Sampling performs no allocation, random-state advancement, mutation,
I/O or clock reads. All time and seed inputs are explicit prepared values.
Custom sample cost must be bounded; unbounded iteration is outside the protocol.

## Source API and initial patterns

Use a small policy protocol, consistent with the existing prepare/sample split
in [pullback sources](../../core/render/pullback/source.h). The following types
are schematic, not existing declarations:

```cpp
struct PatternClock {
  float phase_turns;
};

template<class Pattern> struct PreparedPattern {
  typename Pattern::Prepared state;
};

struct DirectionalSine {
  struct Params { math::Vector axis; float frequency; float offset_turns; };
  struct Prepared { math::Vector axis; float frequency; float phase_turns; };
  static Prepared prepare(const Params&, PatternClock);
  static float sample(const math::Vector& unit_direction, const Prepared&);
};

template<class Pattern> struct DirectSphereField;
template<class Pattern> struct BakedSphereField;
template<class Source, class Field> struct PatternHuePalette;
template<class Source, class Field> struct PatternShimmerPalette;
```

Provide three initial source policies:

| Source | Definition and controls |
|---|---|
| Noise | Adapter around the current configured noise generator, spatial scale and sphere-loop construction. Seed and generator configuration belong to source identity. |
| Directional sine | `sin(2*pi*(frequency*dot(axis,direction) + phase + offset))`; unit axis, nonnegative frequency and offset in turns. Frequency is cycles per unit axial coordinate, so the sphere spans twice that many cycles along the axis. |
| Two-wave interference | Half the sum of two directional sines, with separate axes, frequencies and offsets. Both use the same loop phase; the fixed half gain bounds the result without dynamic normalization. |

Phase is in turns and speed is turns per engine update, matching current
palette noise timing. Advance it once per update, wrap to [0, 1), and preserve
negative-speed behavior. Do not substitute wall-clock seconds. Both analytic
sources return to the same field after one phase turn. Source preparation
hoists validation and frame-constant work outside scanning.
Admission establishes upper frequency and argument bounds compatible with the
chosen trigonometric implementation, including the project's fast-trig limits.
Fold phase before evaluation and keep the resulting argument within that admitted
range. Qualify those bounds during implementation rather than inventing a
universal cutoff. Arbitrary C++ patterns do not imply unbounded trigonometric
arguments or detail free of spatial aliasing.

An effect may supply any custom C++ policy satisfying this protocol. Compile-time
instantiation provides extensibility; this is not an expression interpreter or
an unrestricted runtime callback registry. Built-ins may have a small tagged
selector, resolved during preparation into specialized sampling paths. Custom
sources need a stable recipe identifier and serializable parameter schema to
participate in workbench replay; a function pointer is not a snapshot recipe.

## Applying the scalar

Hue rotation resolves `hue_turns = amount * scalar` and uses the existing
hue-rotation machinery, including periodic hue lookup and gamut mapping.
Amount remains signed and independent of the source. Zero amount bypasses the
color adjustment. Coverage, palette alpha, field-based opacity and the existing
brightness envelope retain their current roles and ordering.
Legacy NOISE is a compatibility adapter into these same field/target internals,
with specialized kernels retained; it is not a second independent engine.

The shimmer adapter resolves `lift = amount * max(0, scalar)` for amount in
[0, 1], then uses the existing lightness-lift operation. It is not a brightness
multiplier and does not reinterpret negative lobes as darkening. Source reuse
does not require putting hue and shimmer into one runtime target enum. Share
prepared resources only when their full source identity and time agree; avoid
a second sample when a caller explicitly applies both to one point. A general
composition-order policy is deferred.

## Direct and baked evaluation

Evaluation strategy is explicit in the effect or workbench recipe:

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
on ownership changes and restore. Hue amount, palette changes and shimmer amount
do not invalidate the spatial table. Palette changes invalidate the separate
hue-rotation color table; spatial pattern changes do not.

Snapshots store the source recipe, parameters, explicit evaluation strategy,
seed where applicable, and current loop phase. Derived tables and pointers are
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

For a new pattern interval with unchanged source kind, strategy and target,
retain the two prepared source states or two baked tables. Evaluate each at
the column direction, apply the shared column amount and mapping to its color
resources, then use the existing endpoint RGB interpolation. Do not additionally
interpolate source phase/parameters inside those endpoint evaluations. This
avoids double interpolation and aligns direct and baked resource semantics.

A source kind, evaluation strategy, target or custom recipe-version change is
a topology change, not a numeric parameter blend. For the initial pattern
implementation, hold the **current** source selection and source state/table
through that interval, including its fraction-zero boundary. Blend palette
resources with shared controls as usual. The source change may therefore be
discontinuous; expose and test that policy rather than claiming cross-pattern
morphing. Ordinary unchanged-topology intervals retain exact endpoint behavior.
Alpha retains its palette/coverage endpoint behavior even when the source is
held. Never reinterpret an old source's resources under the new source type.
Do not change legacy noise topology behavior as part of this feature.

## Compatibility and delivery

Keep `HueMode::NONE = 0`, `NOISE = 1` and `PATH_LENGTH = 2`, including the
existing displacement alias. Add a distinct PATTERN mode without reinterpreting
stored values: the new mode distinguishes pattern recipes while preserving the
meaning of persisted ordinals 1 and 2. Keep legacy field IDs such as `hue-noise-scale` and
`hue-noise-speed`, defaults, seed initialization, phase advancement, noise-loop
construction and optimized noise bake unchanged. An absent pattern descriptor
means the existing behavior, not a new default pattern. Preserve decoding and
encoding of the existing color-clock-v1 snapshot; add versioned pattern state
instead of altering its meaning. Append new IDs and update activation gates,
parameter metadata and workbench controls together.

Implement in stages: extract reusable scalar storage behind compatibility
entry points; add source policies and hue integration; then wire metadata,
snapshots and workbench support. Introduce shimmer's generic adapter without
rewriting its existing noise-facing API. No effect preset migration is required
to enable the feature.

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
* Independent known-value tests for sine and interference, output bounds,
  periodic phase, positive direction-scale invariance and invalid axes.
* Asymmetric orientation cases proving projection-frame versus camera-frame
  provenance, and cube edge/corner seam tests from both adjacent faces.
* Direct-versus-baked error measurements across each admitted baked range,
  separating quantization from spatial error; no blanket arbitrary-pattern
  equivalence claim.
* Temporal tests whose expected endpoint RGB differs from scalar-first blending,
  including changed phase/scale, stationary resources and topology holds.
* Cache tests for every key component and palette-only changes; snapshot replay
  after cache reset, old snapshot compatibility, and unknown-recipe rejection.
* Native and WASM validation plus firmware memory/timing evidence for the
  implementation. This proposed-spec change itself needs documentation checks
  only.
