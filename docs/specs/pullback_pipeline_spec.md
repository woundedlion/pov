# First-class pullback pipeline and stage library

**Status: LANDED, revision 6 (2026-08-14).** The composition core, standard
carriers, provider concepts, stage combinators, and concrete operator catalog
specified here ship in `core/render/pullback.h`, which landed in `13186d7c`.
Its consumers are `effects/ShaderBall.h`, `effects/CurlLattice.h`, and
`effects/fixed/FixedLookRuntime.h` — and through the latter,
`effects/AlienOcean.h`, `effects/FacetGrid.h`, and the four
`effects/fixed/*Looks.h` banks. The
verification artifacts (`tests/test_pullback.h`, `tests/pullback_manifest_check.cpp`,
`tests/data/pullback/`, `tools/pullback_capture.py`) ship with it. Section 17
is the exception: everything listed there is a design record only, with no
implementation in the tree. This revision supersedes revision 5's deferred
framework plan.
The architecture review correctly rejected shared pixel/block drivers and the
GS/BZ/Voronoi migrations, but it coupled that rejection to the unrelated
question of whether pullback composition and concrete pullback operators belong
in the engine. This revision separates those decisions.

The proposed change promotes ShaderBall's typed pipeline, standard carriers,
and reusable concrete operator/stage families into `core/render/pullback.h`.
ShaderBall remains the first production consumer, but consumer count is not an
admission rule for deliberate engine vocabulary. This follows the precedent of
`core/render/filter.h`: a core rendering facility may contain both composition
plumbing and concrete visual operations even when some operations currently
have one or zero shipping consumers.

This revision does **not** introduce acquisition contexts, `Program`,
`PixelDriver`, `BlockDriver`, reduction policies, or finish policies. It does
not migrate BZReactionDiffusion, GSReactionDiffusion, or Voronoi. Their
traversals and shading loops remain unchanged. The independently useful A/B
capture and deterministic WASM verification work from revision 5 remains
endorsed as verification infrastructure and is retained here only to the
extent required to prove this migration.

The shipped behavior remains specified by
[inverse_sampling_pipeline_spec.md](inverse_sampling_pipeline_spec.md) (the
"v1 spec") and [shaderball_spec.md](shaderball_spec.md). This document amends
their location and abstraction-boundary decisions as described in Section 15.

## 1. Decision and motivation

ShaderBall currently owns two different kinds of code in one effect header:

1. effect policy: presets, editable slots, parameter interpolation, clocks,
   resource ownership, frame preparation, transition choreography, topology
   admission, the compiled-program manifest, and runtime program selection;
2. rendering vocabulary: typed pullback carriers, stage-chain validation,
   lenses, projections, surface maps, planar warps, source fields, material
   shaping, and colorization.

The first category is properly effect-owned. The second category is an engine
expressive surface. Keeping it private has three costs:

- a future effect cannot compose a pullback without depending conceptually on
  ShaderBall or recreating its framework;
- reusable looks such as Bonne projection, mirror tiling, wave shear, primitive
  lattice sampling, iso-contour shaping, and projection-edge coverage are not
  discoverable beside the engine's forward filter vocabulary;
- the simulator-only dynamic backend and the compiled program set can drift
  because both orchestrate effect-private kernels independently.

The earlier "second concrete consumer" threshold is therefore retired for
intentional render vocabulary. Reuse count remains evidence about an
abstraction, but it is not a prerequisite when all of the following are true:

- the operation has a stable mathematical meaning independent of an effect;
- its input, output, and prepared-state requirements can be documented without
  naming an effect;
- it can be instantiated and tested with a synthetic frame;
- unused template inventory emits no code or data;
- moving it preserves the shipping hot-path and ownership contracts.

This is the **core-admission test** used by this specification. A helper that
fails it stays in ShaderBall even if it is called by a core stage.

## 2. Goals

The design shall:

- provide `Pullback::Pipeline<Binding, Stages...>` as a first-class typed value
  pipeline in `core/render/pullback.h`;
- define a standard six-role pullback vocabulary and standard carriers for the
  sphere-to-source-to-color path;
- provide concrete, composable core policies for the lenses, projections,
  surface maps, planar warps, source fields, material operations, and color
  operations currently implemented by ShaderBall;
- keep mutable animation and resource ownership outside the hot pipeline;
- let a consumer map its immutable frame into core operations through narrow,
  compile-time state providers, without adopting a universal engine frame;
- let compiled programs and runtime/reference dispatch call the same core
  operator kernels;
- preserve ShaderBall's program IDs, topology keys, manifest rows,
  `ShadeFunction` ABI, selection behavior, stage instrumentation, output, and
  device performance;
- remain header-only, allocation-free, and free of virtual dispatch, RTTI, and
  per-pixel type erasure;
- instantiate only explicitly named compiled programs; the core catalog shall
  not create a Cartesian product of looks;
- make the public API and diagnostics sufficient for a new pullback effect to
  compose stages without copying ShaderBall internals.

## 3. Non-goals

This design does not:

- generalize scan traversal or replace any `Scan::Shader` entry point;
- migrate GSReactionDiffusion, BZReactionDiffusion, Voronoi, a forward
  renderer, or the ray-marched volume path;
- provide candidate acquisition, shared-stencil, block-coherence, reduction,
  or per-pixel finish abstractions;
- move ShaderBall presets, slot enums, `Config`, `Params`, `FrameState`,
  `TopologyKey`, `InversePipelineId`, manifest rows, resource preparation, or
  transition state into core;
- make every runtime slot combination a compiled template instantiation;
- require core operators to own mutable resources or animation state;
- merge `Pullback::Pipeline` with `Filter::Pipeline`; their cardinality,
  direction, state, and terminal contracts remain different;
- change a projection formula, warp order, approximation, floating-point
  operation order, palette result, or alpha convention;
- require a visual operator to be split into the smallest mathematically
  possible node when the shipping implementation needs a fused stage boundary;
- introduce a generic render graph, stored callback chain, virtual stage base,
  or runtime stage registry.

## 4. Location, namespace, and dependency rules

The public facility lives in `core/render/pullback.h` in namespace
`Pullback`, alongside `Scan`, `Filter`, and `SDF`.

Public groups are:

```text
Pullback::Pipeline                 typed six-role coordinator
Pullback::StageKind                standard semantic role enumeration
Pullback::CodeEmission             placement metadata
Pullback::ProjectionSample         sphere-to-plane carrier
Pullback::WarpResult               complete planar-warp carrier
Pullback::SurfaceResult            one sphere-space map result
Pullback::WarpStepResult           one planar-warp result
Pullback::SourceInput              projection + warp carrier
Pullback::MaterialInput            projection + warp + field carrier
Pullback::MaterialSample           terminal material carrier
Pullback::Stage::*                 six stage combinator families
Pullback::Surface::*               sphere-space map policies
Pullback::Lens::*                  lens policies
Pullback::Projection::*            projection policies
Pullback::Warp::*                  planar-warp policies
Pullback::Source::*                scalar-source policies
Pullback::Weight::*                signal-weight policies
Pullback::Transfer::*              value-transfer policies
Pullback::Coverage::*              coverage policies
Pullback::Color::*                 colorization policies and kernels
```

`pullback.h` may include headers from `core/math`, `core/color`, and the
minimal engine concept/profiling headers it needs. It shall not include an
`effects/` header, refer to `ShaderBall`, or require the effect registry.
No `core/` header may include `ShaderBall.h` as a consequence of this work.

Existing pure mathematical kernels remain in their natural owners:

- lenses remain in `core/math/lenses.h`;
- projection kernels and `projections::ProjectionKernelResult` remain in
  `core/math/projections.h`;
- shared noise and stereographic helpers remain in their existing core math
  headers;
- palette/gamut operations remain in `core/color`.

`pullback.h` supplies typed policies and orchestration around those kernels.
An effect-private pure kernel that passes the Section 1 admission test moves to
the corresponding core math/color owner or, when it is meaningful only as a
pullback operator, into `pullback.h`. The migration shall not duplicate a
kernel in core and ShaderBall after the final phase.

## 5. Standard data model

### 5.1 Stage roles

Core standardizes this evaluation order:

```text
OUTER_CAMERA -> SURFACE_PROJECT -> PLANAR_WARP
             -> SOURCE -> MATERIAL -> COLOR
```

`StageKind` has exactly those six enumerators in that order. They describe
pullback evaluation, not authored forward order. Version 1 still requires one
stage of each kind. Variation within a role is expressed by sub-policies and
fused stage combinators rather than by changing the role sequence.

| Kind | Exact input | Exact output |
|---|---|---|
| `OUTER_CAMERA` | `Vector` | `Vector` |
| `SURFACE_PROJECT` | `Vector` | `ProjectionSample` |
| `PLANAR_WARP` | `ProjectionSample` | `SourceInput` |
| `SOURCE` | `SourceInput` | `MaterialInput` |
| `MATERIAL` | `MaterialInput` | `MaterialSample` |
| `COLOR` | `MaterialSample` | `Color4` |

The fixed roles are intentional. They capture the stable semantic interfaces
that make concrete operators interoperable. A future renderer needing a
different cardinality or semantic sequence may add a sibling pipeline type; it
must not weaken this pipeline's diagnostics into an untyped arbitrary chain.

### 5.2 Public carriers

The following field sets are normative. Names may be adjusted only before
Phase A lands; field order and types may not change during migration without
refreshing the layout and disassembly evidence.

```cpp
struct ProjectionSample {
  Complex coords;
  uint8_t region_id;
  uint8_t component_id;
  uint8_t boundary_flags;
  float fade_edge_distance;
  float value_weight;
  uint8_t flags;
  uint8_t traits = 0;
  uint8_t edge_class = 0;
  float domain_coverage = 1.0f;
  Vector sphere = Vector();
  float surface_path_length = 0.0f;
};

struct SurfaceResult {
  Vector sphere;
  float path_length;
};

struct WarpStepResult {
  Complex coords;
  Complex delta;
  float path_length;
};

struct WarpResult {
  Complex coords;
  float path_length;
};

struct SourceInput {
  ProjectionSample projected;
  WarpResult warped;
};

struct MaterialInput {
  ProjectionSample projected;
  WarpResult warped;
  float field;
};

struct MaterialSample {
  float value;
  float coverage;
  Vector sphere;
  float path_length;
};
```

`ProjectionSample` is a generic topology-bearing projection result:

- `coords` are coordinates in the source plane/chart;
- `region_id`, `component_id`, `boundary_flags`, `traits`, `edge_class`, and
  `flags` preserve projection topology and parity;
- `fade_edge_distance` is finite, nonnegative distance to a fade-eligible cut
  or singular boundary;
- `value_weight` and `domain_coverage` are finite values in `[0, 1]`;
- `sphere` is the projection-local sphere direction that produced the sample;
- `surface_path_length` is accumulated by sphere-space maps before projection.

`WarpResult` carries only what a downstream stage consumes: the warped
coordinates and `path_length`, the sum of the step path lengths. No net
displacement vector is accumulated across the warp chain. Each policy still
reports a per-step `WarpStepResult::delta`, which is the value its own
`path_length` is derived from. `MaterialSample::path_length` is the total
surface plus
planar path length. It replaces ShaderBall's misleading member name
`warp_displacement`; this is a source-level rename only and shall preserve
layout and arithmetic.

The projection boundary bit definitions currently private to ShaderBall move
beside `ProjectionSample` as a scoped or prefixed public bit vocabulary. Their
numeric values remain unchanged. Projection-specific flag meanings remain
documented by the projection policy that produces them.

The defaults are load-bearing compatibility with the current
`ProjectedLookup` constructor: partial aggregate initialization must retain
unit domain coverage rather than value-initialize it to zero. Phase A records
`sizeof`, `alignof`, `std::is_trivially_copyable_v`, and every
field offset for the old and new carriers at both production/test
instantiations. The new records shall remain aggregate-friendly, trivially
copyable, and free of ownership.

### 5.3 No universal frame type

Core does not define a monolithic `Pullback::FrameState`. Consumers prepare
different parameters and resources, and forcing them into a common record
would either expose effect policy or add hot-path copies.

Instead, a pipeline has a `Binding` type:

```cpp
struct ExampleBinding {
  using FrameState = ExampleEffect::FrameState;
  using Instrumentation = Pullback::NoInstrumentation;
};
```

Concrete operators take narrower **state-provider** types. Each provider names
the same `FrameState` and exposes only the data required by that operation.
Providers are empty compile-time adapters with `always_inline` static accessors.
They neither own nor copy state.

This is the principal decoupling boundary: core owns algorithms and carriers;
the effect owns frame layout and maps it into those algorithms.

## 6. Pipeline and stage contract

### 6.1 Pipeline type

The public coordinator is:

```cpp
template <typename BindingT, typename... Stages>
struct Pipeline {
  using Binding = BindingT;
  using FrameState = typename BindingT::FrameState;

  __attribute__((always_inline))
  static Color4 evaluate(const Vector &view, const FrameState &frame);

  HS_FLASH_MEMBER
  static Color4 shade(const Vector &view, const FrameState &frame);
};
```

`evaluate` recursively invokes the six stages and is the entry point used by a
consumer wrapper that needs custom placement. `shade` preserves the ordinary
two-argument ABI and calls `evaluate`. ShaderBall may retain derived/named
wrappers where its accepted manifest requires `FASTRUN`, `noinline`, or
alignment attributes; those wrappers call `evaluate` and do not reimplement
any stage.

The coordinator and recursive carrier fold are `always_inline` and shall emit
no symbols, frames, or stored stage objects. `shade` is the only default
address-taken wrapper.

### 6.2 Stage contract

Every top-level stage provides:

```cpp
using Binding = /* pipeline binding */;
using FrameState = typename Binding::FrameState;
using Input = /* exact public carrier */;
using Output = /* exact public carrier */;

static constexpr Pullback::StageKind KIND;
static constexpr Pullback::CodeEmission EMISSION;
static constexpr bool APPROXIMATE;
static constexpr bool TERMINAL;
static constexpr bool NON_FLOATING_FIELDS_EXACT;
static constexpr Pullback::ApproximationOracleId ORACLE;
static constexpr std::array<Pullback::ApproximationMetric, N> METRICS;

static Output run(const Input &, const FrameState &);
```

Stages are empty policy types. They allocate nothing, retain no frame
reference, have no mutable static or function-local state, and read no runtime
selector for a choice already fixed in their type. Continuous values and
deliberately runtime color controls remain immutable frame reads.

`CodeEmission`, `ApproximationDomain`, `ApproximationAggregation`,
`ApproximationMetric`, and `ApproximationOracleId` move from ShaderBall into
`Pullback`. Initial oracle IDs are `NONE`, `PEIRCE_FAST_SQUARE`, and
`HUE_ROTATION_AND_NOISE_LUTS`; future core approximate operators extend the
enumeration with their oracle and tests in the same change.

An exact stage declares `ORACLE == NONE` and an empty metric list. An
approximate stage declares a non-`NONE` oracle, a non-empty metric list
including a final-framebuffer metric, and
`NON_FLOATING_FIELDS_EXACT == true`. Approximation approval belongs to the
operator/stage that calls the approximation, not to the consumer.

### 6.3 Validation

Validation retains the shipped two-level design. Detection-only Level 1 checks
arity and member presence/types without indexing a malformed pack. Level 2 is
formed only after Level 1 succeeds and checks:

- exactly six stages;
- each `Stage::Binding` is exactly the pipeline `Binding`;
- kinds occur in the Section 5.1 order;
- the first input is `Vector` and final output is `Color4`;
- each `run(const Input&, const FrameState&)` returns exactly `Output`;
- adjacent output/input carrier types match exactly;
- only the color stage is terminal;
- every stage is empty;
- approximation metadata is well formed;
- optional `Binding::ExtraValidation<Stages...>::value` is true.

Named `static_assert` diagnostics cover: wrong arity, missing/mistyped stage
contract, binding mismatch, wrong order, wrong return type, carrier mismatch,
non-empty policy, terminal misuse, malformed approximation metadata, and
consumer validation failure.

The validation booleans remain publicly testable without instantiating a
failing pipeline. `tests/test_pullback.h` contains positive and negative
constexpr predicate tests for every category. This specification does not
require a generated negative-translation-unit harness.

### 6.4 Topology matching

Topology admission is not a core-stage concern. A core stage describes what it
does; an effect decides which configuration selects it.

ShaderBall supplies empty derived wrappers or aliases that add:

```cpp
static constexpr bool implements(const TopologyKey &);
```

`Pipeline` conditionally exposes an `implements(key)` fold when every stage
provides a compatible predicate. The member is instantiated only when used.
`TopologyKey` therefore remains private to ShaderBall.

A thin matching wrapper may inherit `run`, carriers, and metadata from a core
stage, but it shall add no data and no alternate rendering code. Core may
provide a generic decorator for this pattern only if it remains independent of
the key type; it is not required by this revision.

ShaderBall's cross-stage rule that unconditional Peirce edge distance requires
edge-fade coverage moves to `Binding::ExtraValidation`. Coverage enums and
topology rules do not enter the pipeline coordinator.

## 7. State providers and instrumentation

### 7.1 Provider contract

A concrete operator is parameterized by a provider, not by an effect:

```cpp
struct OuterWarpState {
  using Binding = ShaderBallBinding;
  using FrameState = typename Binding::FrameState;

  static const auto &params(const FrameState &);
  static const auto &prepared(const FrameState &);
  static float phase(const FrameState &);
  static const FastNoiseLite &noise(const FrameState &);
};
```

Only the accessors needed by the selected operator are required. For example,
`Warp::MirrorTile` does not require `noise`, and `Lens::Glitch` requires no
provider at all. Accessor return requirements are structural and documented at
the operator declaration: a wave-shear parameter view must expose
`strength` and `frequency`; an affine prepared view must expose the fields its
formula reads. Providers may return existing consumer records by const
reference. Core shall not require construction of a per-pixel view object.

Every provider:

- is empty and trivially constructible;
- names its owning pipeline `Binding` and derives `FrameState` from that
  binding;
- names `FrameState` exactly;
- exposes only static, `always_inline`, const-frame accessors;
- returns const references/pointers or scalar values with lifetimes valid for
  the draw;
- performs no validation, allocation, mutation, or runtime dispatch;
- is checked by a provider-specific C++20 concept and a named diagnostic when
  its operator is instantiated.

Provider concepts are deliberately local to each operator. There is no giant
`PullbackBinding` concept requiring resources an effect does not use.

The initial provider surface is normative at the category level:

| Provider category | Accessors available to policies in that category |
|---|---|
| orientation | prepared inverse `conjugate(frame)` |
| surface map | `params(frame)`, `prepared(frame)`, `phase(frame)`, optional `noise(frame)`, and `path_length_required(frame)` |
| projection | prepared frame `conjugate(frame)` plus scalar `pole_fade`, `central_meridian`, `coordinate_scale`, `standard_parallel`, `layout_scroll`, and `edge_distance_required` accessors as required by the selected map |
| planar warp slot | `params(frame)`, `prepared(frame)`, `phase(frame)`, and optional `noise(frame)`; basis, envelope, integrator, and polar mode are template facts in a compiled policy |
| source | `params(frame)`, `prepared(frame)`, and optional noise resource/time accessors; the selected source policy determines the required subset |
| material | value/coverage scalar accessors (`iso_level`, `iso_width`, band values, cutout values, `edge_width`) required by the selected policies |
| color | immutable color parameters/clocks, generated palette binding, prepared hue-rotation LUT, prepared hue-noise LUT, and deliberately runtime mapping/brightness/hue mode values |

An operator's declaration narrows this table with a `requires` expression that
names every field/member it reads and no unrelated member. For example,
`WaveShear` requires only `strength`, `frequency`, prepared rotation sine/cosine,
and phase; `MirrorTile` requires only its cell/offset prepared values. These
requirements are part of the public doxygen contract. Adding a new hot-path
read therefore changes the provider concept and its tests in the same commit.

For compiled policies, a provider's `Binding` must exactly equal the enclosing
stage's `Binding`; provider concepts expose this as a testable boolean before a
`static_assert`. This is what turns accidental use of a provider from another
effect or frame layout into a named binding diagnostic rather than a deep
substitution error.

### 7.2 Instrumentation

Moving code into core shall not erase ShaderBall's stage buckets or bake
ShaderBall profile fields into core.

`Binding::Instrumentation` supplies an optional zero-state hook policy. The
required shape is:

```cpp
struct NoInstrumentation {
  struct Token {};
  static Token mark();
  template <Pullback::ProfileEvent> static void span(Token);
};
```

`ProfileEvent` covers the existing generic boundaries: `LENS`,
`SURFACE_NOISE`, `PROJECTION`, `PLANAR_WARP`, `MIRROR_TILE`, `SOURCE`,
`MATERIAL`, and `COLOR`. `NoInstrumentation` compiles to no statements.
ShaderBall supplies a hook policy that maps these events to the existing
`HS_SB_STAGE_MARK`/`HS_SB_STAGE_SPAN` counters. Instrumented builds shall
preserve the current boundary and nesting order. Shipping builds shall show no
hook symbol or instruction.

## 8. Concrete core operator and stage catalog

The catalog is public template inventory. Defining a policy in the header does
not instantiate it. Only a pipeline alias or a direct dynamic/reference call
emits code.

Every operator exposes a direct pure entry point (`apply`, `project`, or
`sample`, as appropriate). Top-level `Stage::*` combinators call those entry
points. ShaderBall's dynamic/reference backend calls the same entry points from
its runtime switches. There shall be one implementation of each formula after
migration.

Every sub-policy also exposes the approximation metadata block from Section
6.2, minus `KIND` and `TERMINAL`. A top-level combinator forwards the metadata
of its one approximate child, or declares itself exact when all children are
exact. The initial catalog permits at most one approximate child in a
top-level stage; composing two distinct approximation oracles fails under the
named diagnostic `pullback stage: multiple approximation oracles require an
explicit combined oracle`. This prevents a stage from silently dropping one
child's error contract. A future approximation-set design may relax the rule
with its own specification and tests.

Each combinator exports aliases naming its constituent policies so consumers
and tests do not pattern-match template arguments: surface/project exports
`PreLensSurface`, `LensPolicy`, `PostLensSurface`, and `ProjectionPolicy`;
planar warp exports its policy tuple and `policy_at<I>`; source exports
`SourcePolicy`; material exports `WeightPolicy`, `TransferPolicy`, and
`CoveragePolicy`; color exports `ColorPolicy`. `Pipeline` similarly exports
`Binding`, `FrameState`, `stage_at<I>`, and `STAGE_COUNT`. Static operator facts
needed by consumer validation, such as unconditional edge distance, are
forwarded by the owning top-level stage.

Top-level stage combinators are inline by default. Non-default placement uses
one public wrapper:

```cpp
namespace Pullback::Stage {
template <CodeEmission EmissionV, typename StageImplementationT>
struct Placed;
}
```

`EmissionV` is one of `CodeEmission::INLINE_ONLY`,
`CodeEmission::OUT_OF_LINE_FLASH`, or `CodeEmission::OUT_OF_LINE_ITCM`.
`Placed` forwards the complete stage contract, sets `EMISSION = EmissionV`, and
defines the public `run` through an enum-specialized entry: `always_inline` for
`INLINE_ONLY`, `HS_FLASH_MEMBER` for `OUT_OF_LINE_FLASH`, and
`FASTRUN HS_NOINLINE_NOCLONE` for `OUT_OF_LINE_ITCM`. The wrapped implementation
is inlined into that entry. Metadata that disagrees with actual placement is
therefore not representable. ShaderBall uses `Placed` only where the baseline
stage is not inline, and supplies the exact baseline enum during migration.

### 8.1 Outer camera

```cpp
namespace Pullback::Stage {
template <typename BindingT, typename OrientationState>
struct OuterCamera;
}
```

`OrientationState` supplies `FrameState` and
`conjugate(const FrameState&) -> const Quaternion&` (a by-value return is
permitted only when the consumer already stores it by value and disassembly
shows no copy penalty). The stage maps `Vector -> Vector` with `rotate(view,
conjugate(frame))`.

The initial catalog contains this prepared-inverse-orientation camera. Camera
animation, wandering, and conjugate preparation remain consumer responsibilities.

### 8.2 Surface/lens/projection

The fused stage is:

```cpp
namespace Pullback::Stage {
template <typename BindingT,
          typename PreLensSurface,
          typename LensPolicy,
          typename PostLensSurface,
          typename ProjectionPolicy>
struct SurfaceProject;
}
```

Its exact evaluation order is:

```text
pre  = PreLensSurface::apply(input, frame)
lens = LensPolicy::apply(pre.sphere, frame)
post = PostLensSurface::apply(lens, frame)
local = rotate(post.sphere, ProjectionPolicy::frame_conjugate(frame))
out = ProjectionPolicy::project(local, frame)
out.sphere = local
out.surface_path_length = pre.path_length + post.path_length
```

An identity surface returns `{input, 0}`. A lens returns a `Vector` and adds no
path length. A surface map returns `SurfaceResult`. A projection returns a
complete `ProjectionSample` except for `sphere` and
`surface_path_length`, which the fused stage sets as shown. This preserves the
single `SURFACE_PROJECT` carrier boundary and permits compile-time composition
inside it.

Required surface policies:

- `Surface::Identity`;
- direct sphere-space noise;
- curl sphere-space noise for every currently implemented integrator, including
  the shipping simplex/Euler specialization.

Noise basis and surface integrator are compiled topology facts. Public policy
shapes are `Surface::DirectNoise<State, NoiseBasisV>` and
`Surface::CurlNoise<State, NoiseBasisV, SurfaceIntegrator>`, where
`NoiseBasisV` is the core `::NoiseBasis` non-type value and
`SurfaceIntegrator` is one of `Surface::Euler`, `Surface::Midpoint`, or
`Surface::Midpoint2x`. Direct kernels accept basis/integrator as explicit
runtime values for the dynamic dispatcher. A compiled policy shall not read
either selector from its provider.

Noise surface providers expose `path_length_required(frame)`. When false, the
policy returns an exact zero path without performing the square root or
integration bookkeeping used only by path-sensitive hue. This preserves the
current cross-stage cost optimization without making the surface policy inspect
ShaderBall color enums. ShaderBall's provider computes the boolean from its
prepared hue state and runtime hue mode.

Required lens policies:

- identity;
- glitch;
- twist;
- ordinary kaleidoscope;
- Mobius;
- tetrahedral, octahedral, and dodecahedral kaleidoscopes;
- triangular, square, pentagonal, hexagonal, and octagonal prism
  kaleidoscopes.

`SurfaceLens::TANGENT_NOISE` is a serialization tombstone, not a rendering
policy. Valid configurations reject it and the current lens dispatcher has no
formula. The Section 8.7 census classifies it as `SENTINEL`/import-only,
parallel to `WarpStageKind::LEGACY_STEREO_NOISE`. It may enter the catalog only
after a formula and finite-domain contract are separately specified and
admitted.

Policies that are already pure functions in `core/math/lenses.h` are thin
empty wrappers. Parameterized policies use a narrow provider.

Required projection policies:

- stereographic;
- folded sinusoidal;
- folded/front/back gnomonic variants;
- north/south Bonne variants;
- Peirce quincuncial layouts, including the approximate fast-square policy;
- vertical/horizontal Airocean;
- equirectangular.

Compile-time layout/hemisphere choices are template arguments when the compiled
program fixes them. Continuous central meridian, scale, scroll, fade, and
prepared frame rotation come through providers. Projection policies convert
`projections::ProjectionKernelResult` into `ProjectionSample` in one shared
helper; scaling of coordinates and edge distance occurs in the same arithmetic
order as today.

The shared Peirce and Airocean direct kernels take an explicit
`edge_distance_required` boolean. A compiled policy whose topology proves the
answer supplies a boolean template argument and emits no provider read; the
dynamic policy obtains it from `ProjectionState::edge_distance_required(frame)`.
ShaderBall's dynamic provider preserves the current predicate over edge-fade
material coverage and both warp envelopes. This dependency is prepared or
adapted by the consumer; core projection code does not inspect ShaderBall
coverage or warp enums.

The Peirce fast-square policy owns its oracle metadata. Whether edge distance
is unconditional is exposed as a static policy fact for ShaderBall's
`ExtraValidation`; it is not a special case in the coordinator.

### 8.3 Planar warps

The stage is a variadic lookup-order sequence:

```cpp
namespace Pullback::Stage {
template <typename BindingT, typename... WarpPolicies>
struct PlanarWarp;
}
```

It maps `ProjectionSample -> SourceInput`. Starting from
`projected.coords`, it invokes each policy in template order. For each
`WarpStepResult`, it passes the output coordinates to the next policy and adds
the reported `path_length` to the running total. The empty pack is exact
identity.

This generalizes ShaderBall's two fixed positions without imposing two slots on
other consumers. ShaderBall still instantiates exactly its authored outer then
inner lookup order; its authoring model remains a two-slot program.

Required policies:

- identity;
- affine frame;
- wave shear;
- vortex;
- vector noise;
- curl flow for each implemented integrator;
- mirror tile;
- polar chart, linear and logarithmic.

`WarpStageKind::LEGACY_STEREO_NOISE` is not a policy. It is a serialization
tombstone rejected by `valid_config`; its current switch arm is unreachable and
has no admitted rendering formula. The Section 8.7 census classifies it as
`SENTINEL`/import-only. A future implementation would require a separately
specified formula and admission before it could enter the core catalog.

Compile-time choices use core tag types, not consumer authoring enums. The
initial tags are `Warp::FlatEnvelope`, `Warp::ProjectionWeightEnvelope`,
`Warp::EdgeFadeEnvelope`; `Warp::Euler1`, `Warp::Midpoint2`,
`Warp::Midpoint4`; and `Warp::LinearPolar`, `Warp::LogarithmicPolar`.
ShaderBall maps its serialized enums to these types when declaring a compiled
program or entering a dynamic switch arm.

Polar harmonic is an independent discrete topology fact, not part of the
continuous warp parameter view. The compiled policy is
`Warp::PolarChart<State, PolarMode, uint8_t Harmonic>` and requires
`1 <= Harmonic <= Pullback::Warp::MAX_POLAR_HARMONIC`, where the initial public
maximum is the current admitted value 16. Its shared direct kernel takes mode
and harmonic explicitly. A compiled stage supplies both from its type; the
dynamic dispatcher supplies the validated runtime values from ShaderBall's
`WarpStageSpec`. Neither path makes the core kernel read that private record.

Envelope policies (`FLAT`, projection weight, edge fade) are separate
compile-time policies where the compiled topology fixes them. They receive the
original `ProjectionSample`, not progressively warped metadata. A runtime
dynamic dispatcher may select the same policies by switch.

Closed-form policies use one shared helper to derive the step `delta` and its
`path_length` from the input and output coordinates. Integrating policies
report the net displacement over their sub-steps as `delta` while summing the
individual sub-step lengths into `path_length`, so the two stay distinct. Both
compute `path_length` only when the frame requests it. Exact-zero strength
bypasses all noise, trigonometry, and integration exactly as today.

### 8.4 Sources

The stage is:

```cpp
namespace Pullback::Stage {
template <typename BindingT, typename SourcePolicy>
struct Source;
}
```

It maps `SourceInput -> MaterialInput`, preserving the complete projection and
warp records and adding one signed scalar field. `SourcePolicy::sample` receives
the `SourceInput`, not only its final `Complex`, so spherical fields may read
`projected.sphere` without a separate stage shape.

Required source policies:

- twin wave;
- rings;
- spiral;
- grid/coupled grid;
- projected noise contour;
- spherical noise contour;
- primitive lattice.

Noise source shapes are
`Source::ProjectedNoise<State, NoiseBasisV>` and
`Source::SphericalNoise<State, NoiseBasisV>`, again using core
`::NoiseBasis` as a non-type topology value. Their shared direct kernels accept
an explicit runtime basis for dynamic dispatch. The operator census asserts
that every surface/source noise basis and integrator represented in
`TopologyKey` is either encoded in the compiled type or is being exercised by
the deliberately runtime dynamic arm.

Coordinate conditioning is an explicit source sub-policy. Twin wave, rings,
spiral, and grid/coupled-grid consume the existing stereographically
conditioned coordinate. Projected noise contour, spherical noise contour, and
primitive lattice do not. The compiled types encode the appropriate
conditioning; a source that consumes raw or spherical coordinates does not pay
for it. Direct dynamic/core comparisons cover all seven source enumerators.
Source phase and prepared trigonometric state come through providers. Source
policies own no clocks.

### 8.5 Material

The material stage is:

```cpp
namespace Pullback::Stage {
template <typename BindingT,
          typename WeightPolicy,
          typename TransferPolicy,
          typename CoveragePolicyT>
struct Material;
}
```

It maps `MaterialInput -> MaterialSample`. Evaluation order is normative:

1. weight the signed field (`NONE` or projection `value_weight`);
2. map to and clamp the unit interval;
3. apply the selected transfer;
4. calculate a policy factor from the original projection metadata and/or the
   transferred value, then multiply it by
   `projected.domain_coverage` (`OPAQUE` has factor 1);
5. emit `{value, coverage, projected.sphere,
   projected.surface_path_length + warped.path_length}`.

Required weight policies: none and projection weight.

Required transfer policies: linear, ridge, iso-contour, and smooth bands.

Required coverage policies: opaque, projection weight,
projection-weight-squared, value cutout, and edge fade. Edge width zero remains
a hard edge; positive width uses the current smooth ramp. Policies shall not
silently clamp invalid projection metadata that projection validation is
required to reject.

The current fused compiled combinations may be expressed as aliases over this
template. The compiler must see the fixed policies; no runtime switch is added
to a compiled material stage.

### 8.6 Color

The terminal stage is:

```cpp
namespace Pullback::Stage {
template <typename BindingT, typename ColorPolicyT>
struct Color;
}
```

It maps `MaterialSample -> Color4`, is the only terminal stage, and returns
straight alpha. `Scan::Shader` remains responsible for final premultiplication.

Required public color kernels/policies cover the current generated-palette
path:

- cup, bell, ascending/linear, and descending/reverse palette mappings;
- phase offset and oscillation;
- palette lookup;
- none, sphere-noise, and total-path-length hue shifts;
- brightness envelope and value-dependent opacity;
- the prepared hue-rotation and hue-noise LUT approximations with their
  current oracle metadata.

Core owns the runtime selector types consumed by these kernels:

```cpp
namespace Pullback::Color {
enum class PaletteMapping : uint8_t {
  CUP = 0, BELL = 1, LINEAR = 2, REVERSE = 3
};
enum class BrightnessEnvelope : uint8_t {
  NONE = 0, CUP = 1, BELL = 2, ASCENDING = 3, DESCENDING = 4
};
enum class HueMode : uint8_t {
  NONE = 0, NOISE = 1, PATH_LENGTH = 2
};
}
```

The values and semantics are stable core API, not aliases of ShaderBall's
serialized enum types. ShaderBall keeps its authoring enums and provides
constexpr conversions guarded by a `static_assert` for every enumerator/value
pair. Its color provider returns the core enum types, never a raw integer.
Because the admitted values are intentionally identical, the conversion is a
zero-cost underlying-value conversion; disassembly must show no per-pixel
conversion switch.

The shipping `GeneratedPalette` color policy may retain deliberately runtime
mapping, brightness, and hue selectors because ShaderBall does not compile
those facets into its topology key. Those switches call the public core
kernels and are documented as runtime policy, not accidental dispatch. A
future consumer may instantiate fixed sub-policies without changing the
terminal carrier.

Palette generation, LUT preparation, mutable provider state, and resource
lifetime remain outside core stages. A color provider exposes const prepared
bindings and scalars only.

#### 8.6.1 Prepared hue LUT contract

Core owns non-owning views and constants for the two prepared LUTs:

```cpp
namespace Pullback::Color {
struct HueRotationLutView {
  static constexpr int VALUE_STEPS = 64;
  static constexpr int HUE_STEPS = 16;
  static constexpr size_t SIZE = VALUE_STEPS * HUE_STEPS;
  const Pixel *data;
  bool active;
};

struct HueNoiseLutView {
  static constexpr int FACE_COUNT = 6;
  static constexpr int FACE_STEPS = 24;
  static constexpr int FACE_SIZE = FACE_STEPS * FACE_STEPS;
  static constexpr size_t SIZE = FACE_COUNT * FACE_SIZE;
  const int8_t *data;
  bool active;
};
}
```

The views own no storage. Providers return them by const reference from the
immutable frame, or by value when the compiler produces the same two-scalar
loads; the latter requires disassembly evidence. Effect-owned mutable spans are
passed only to cold preparation kernels before the frame view is published.

The rotation table is row-major `[value][hue]` in `Pixel` elements. Preparation
samples value row `i` at `ONE_BELOW_UNIT * i / 63` and hue column `j` at
`j / 16`, using the existing gamut-aware hue-rotation base and kernel. Sampling
linearly interpolates adjacent value rows and wrapping adjacent hue columns
with `frac_to_q16` weights and `Pixel::lerp16`, in the current arithmetic order.

The noise table is six consecutive row-major `24 x 24` signed-byte faces in
`+X, -X, +Y, -Y, +Z, -Z` order. Face directions before normalization are:

```text
+X ( 1, v,  u)    -X (-1, v, -u)
+Y ( u, 1,  v)    -Y ( u,-1, -v)
+Z ( u, v,  1)    -Z (-u, v, -1)
```

Preparation samples `u,v` inclusively over `[-1,1]`, clamps noise to `[-1,1]`,
multiplies by 127, rounds halves away from zero, and stores `int8_t`. Sampling
selects the dominant cube axis with tie precedence X, then Y, then Z, applies
the corresponding face coordinates above, bilinearly interpolates four signed
bytes, and multiplies by `1/127`.

`Pullback::Color` provides stateless
`prepare_hue_rotation_lut`, `prepare_hue_noise_lut`,
`sample_hue_rotation_lut`, and `sample_hue_noise_lut` entry points implementing
this contract. Preparation accepts fixed-extent mutable spans and the required
const palette/noise inputs; sampling accepts the views. ShaderBall retains
buffer ownership, arena budgeting, preparation scheduling, and `active`
decisions, but it does not retain any indexing, interpolation, face-mapping, or
quantization formula.

### 8.7 Catalog completeness and exclusions

Phase B builds a census from the current ShaderBall discrete operator enums and
classifies every enumerator as one of:

- `CORE_POLICY`: implemented by a named public policy above;
- `ALIAS`: another spelling/configuration of a named core policy;
- `SENTINEL`: invalid/count/import-only value with no rendering semantics;
- `EFFECT_ONLY`: rejected by the Section 1 admission test, with a written
  reason.

No currently reachable rendering enumerator may disappear from the census.
An `EFFECT_ONLY` classification requires explicit architecture-review approval;
"only ShaderBall uses it" is not a valid reason. This prevents the migration
from promoting only the convenient shipping subset while leaving the rest of
the expressive catalog private.

## 9. ShaderBall integration

### 9.1 What moves to core

ShaderBall replaces its private definitions with core types or thin aliases:

- `InverseStageKind` -> `Pullback::StageKind`;
- emission and approximation metadata types -> `Pullback` equivalents;
- `ProjectedLookup`, `SurfaceNoiseResult`, `PlanarWarpStageResult`,
  `PlanarWarpResult`, `SourceInput`, `MaterialInput`, and `MaterialSample` ->
  aliases of the Section 5 carriers (with the `path_length` member rename
  propagated);
- `InversePipeline` -> alias of `Pullback::Pipeline<ShaderBallBinding, ...>`;
- private stage bodies -> core stage aliases plus effect-side topology matchers;
- reusable pure operator bodies -> the Section 8 core catalog.

### 9.2 What remains in ShaderBall

The following remain private effect code:

- discrete authoring enums and serialized/imported slot values;
- `Config`, `Params`, interpolation, schema validation, and GUI field layout;
- `FrameState` and all mutable owners from which it is prepared;
- state-provider adapters mapping `FrameState` to core operators;
- clocks, walks, prepared transforms, resource union/preparation, and lifetime
  checks;
- `TopologyKey`, `InversePipelineId`, `ProgramDescriptor`, manifest rows,
  continuous preconditions, and resource-readiness predicates;
- preset choreography, topology transitions, `PreparedEndpoint`, `FrameShader`,
  and draw orchestration;
- profiling-hook mapping and test-only ownership generation checks.

This boundary is normative. Moving a listed effect-policy type to core requires
a separate specification explaining its independent semantics.

### 9.3 Assembly example

Illustrative assembly follows. Public `Pullback` family names and roles are
normative; private provider, matcher, and program-alias names may differ:

```cpp
using Surface = ShaderBallMatchedStage<
    Pullback::Stage::SurfaceProject<
        ShaderBallBinding,
        Pullback::Surface::Identity,
        Pullback::Lens::Kaleidoscope,
        Pullback::Surface::Identity,
        Pullback::Projection::BonneNorth<BonneProjectionState>>,
    BonneKaleidoscopeMatcher>;

using Warp = ShaderBallMatchedStage<
    Pullback::Stage::PlanarWarp<
        ShaderBallBinding,
        Pullback::Warp::MirrorTile<OuterWarpState>>,
    OuterMirrorMatcher>;

using Source = ShaderBallMatchedStage<
    Pullback::Stage::Source<
        ShaderBallBinding,
        Pullback::Source::PrimitiveLattice<SourceStateProvider>>,
    PrimitiveLatticeMatcher>;

using Material = ShaderBallMatchedStage<
    Pullback::Stage::Material<
        ShaderBallBinding,
        Pullback::Weight::Projection,
        Pullback::Transfer::Linear,
        Pullback::Coverage::EdgeFade<ValueStateProvider>>,
    ProjectionLinearEdgeFadeMatcher>;

using Outer = ShaderBallMatchedStage<
    Pullback::Stage::OuterCamera<ShaderBallBinding, OuterCameraState>,
    MatchEveryTopology>;

using Color = ShaderBallMatchedStage<
    Pullback::Stage::Color<
        ShaderBallBinding,
        Pullback::Color::GeneratedPalette<ColorStateProvider>>,
    MatchEveryTopology>;

using Program = Pullback::Pipeline<
    ShaderBallBinding, Outer, Surface, Warp, Source, Material, Color>;
```

Every ShaderBall pipeline stage is decorated, including stages whose matcher
always returns true. Consequently every assembled program satisfies the
conditional `Pipeline::implements(key)` contract and the unchanged manifest
`static_assert`; core stages themselves remain topology-agnostic.

The manifest still binds `&Program::shade` or a placement-specific wrapper to
the existing `Color4 (*)(const Vector&, const FrameState&)` ABI. It keeps the
same 11 rows, IDs, keys, and topology/precondition/resource checks.

### 9.4 Dynamic/reference backend

The simulator-only dynamic backend remains a runtime dispatcher because its
purpose is experimentation outside the closed compiled roster. Its switches
shall call the same core operator entry points used by typed policies.

For example, the dynamic wave-shear case calls the public wave-shear kernel
with the selected provider data; it does not retain an effect-private copy.
The reference backend may compose operators manually to preserve independent
orchestration, but it may not duplicate formulas. Exact-reference functions
for approximate kernels remain independent as required by the oracle contract.

Dynamic-versus-compiled equivalence tests therefore check composition and
selection differences, while operator unit tests check formulas.

## 10. Ownership, lifetime, and mutation

- The consumer's `FrameState` is immutable for the duration of a draw.
- Providers borrow only state reachable from that frame.
- Mutable `FastNoiseLite`, palette cyclers, animation objects, generated
  palettes, LUT storage, and arenas remain consumer-owned.
- Frame resource pointers are const bindings whose owners outlive the draw.
- No core stage calls a resource setter, advances a clock, steps an animation,
  prepares a transform, or allocates.
- No provider returns a reference to a temporary. Scalar returns are permitted;
  aggregate parameter/prepared records are returned by const reference.
- The coordinator and policies contain no objects, so a pipeline has no
  lifetime independent of the frame.
- Debug/test owner-generation stamping remains in ShaderBall's frame accessors;
  extraction shall not bypass it.
- Transition rendering remains sequential: prepare and consume one endpoint
  before shared backing storage is overwritten.

The existing ShaderBall stack, persistent arena, RAM2, and effect-heap budgets
remain unchanged. Public carriers add no allocation or hidden ownership.

## 11. Code generation and instantiation contract

### 11.1 Header-only inventory

`pullback.h` is header-only. A policy that is never instantiated emits no text,
rodata, vtable, RTTI, constructor, or registration record. Tests include an
empty-framework device build and compare its ELF with the baseline.

### 11.2 Inlining and placement

- coordinator recursion, state-provider accessors, trivial policy wrappers,
  and stages marked `INLINE_ONLY` are `always_inline`;
- existing large flash leaves retain `HS_FLASH_MEMBER`/equivalent placement;
- an existing measured ITCM leaf retains its placement unless a new device
  acceptance explicitly changes it;
- address-taken program wrappers retain their current individual placement,
  `noinline`, and alignment attributes;
- moving a function to core does not authorize blanket O3 or forced inlining;
- no runtime discrete dispatch is introduced into a compiled stage;
- no carrier move may introduce sret traffic or an unexpected stack frame.

Template policy wrappers should forward explicit values into shared
non-template/policy-independent leaves where that preserves codegen. Binding
types shall not cause one clone of a large mathematical leaf per consumer.
The ELF audit checks `.isra`, `.constprop`, and IPA clones as well as source
names.

### 11.3 Program-set budget

ShaderBall continues to instantiate only the closed manifest. Adding public
inventory does not add a manifest row. The final ELF shall contain:

- the same address-taken compiled-program wrapper count as the baseline;
- no private `InversePipeline` coordinator symbol;
- no emitted `evaluate`/carrier-fold/provider-accessor symbols;
- no second copy of a moved operator formula;
- no simulator-only dynamic dispatcher in the Teensy image;
- no change to the selectable topology census.

Per-`W,H` duplicate accounting from the v1 spec remains in force. Mangled names
may change; byte and placement changes require attribution.

## 12. Correctness and test requirements

### 12.1 Core unit tests

`tests/test_pullback.h` is created in Phase A and wired into the ordinary
native test module/CMake/CI assertion accounting in the same change.

#### 12.1.1 Phase A composition tests

The Phase A surface covers:

- every pipeline validation predicate, positive and negative;
- exact stage invocation order and one invocation per stage;
- carrier forwarding and exact terminal return;
- custom binding/frame types proving there is no ShaderBall dependency;
- `NoInstrumentation` producing no observable calls and a counting hook
  receiving the exact event order;
- public carrier compatibility defaults, finite/range invariants, layout, and
  partial-initialization assertions.

#### 12.1.2 Phase B provider and operator tests

Starting in Phase B, the same module additionally covers:

- provider concept success and one malformed case per accessor category;
- identity surface, lens, warp sequence, weight, transfer, coverage, and color
  policies;
- a two-warp sequence proving progressive coordinates, original projection
  metadata, net-delta addition, and path-length addition;
- direct-operator versus top-level-stage equality;
- direct dynamic/core equality for every catalog operator and selector;
- every approximate operator against its independent exact oracle and
  checked-in manifest.

Concrete operator tests use deterministic inputs and directly pin formula
results or compare with the pre-migration ShaderBall helper during the
migration phase. Phase A gates only Section 12.1.1; the complete Section 12.1
suite becomes mandatory from Phase B onward.

### 12.2 Catalog and ShaderBall tests

The existing `tests/test_shaderball.h` suite is updated, not weakened. It adds:

- the Section 8.7 operator-enum census;
- every manifest row assembled solely from core top-level stages plus empty
  ShaderBall matching wrappers;
- all state providers empty and bound to the exact `FrameState`;
- provider accessors selecting the expected frame member/resource for outer and
  inner slots;
- dynamic/reference dispatch and compiled policies calling identical core
  operator entry points;
- all 11 program keys, selection IDs, continuous preconditions, resource
  readiness, rejection, and transition endpoints unchanged;
- stage-boundary comparison for every public carrier;
- seam, singularity, region, component, parity, and edge metadata coverage;
- zero/nonzero warp strength and exact-zero bypasses;
- zero-width/nonzero-width edge fade, zero/partial alpha, and hue modes;
- both roster resolutions and every manifest entry in framebuffer sweeps.

The test-only legacy stage path may remain behind a migration-only define until
Phase C acceptance. It is deleted before landing Phase D. No shipping build
option restores it afterward.

### 12.3 Cross-commit framebuffer equivalence

Phase P creates the previously aspirational manifests under
`tests/data/pullback/`: one versioned `programs.json` entry per compiled
program and separate oracle manifests for Peirce fast-square and the combined
hue-rotation/noise LUT approximation. Each program entry records the stable ID
and topology key, corpus generator and version, parameter default/endpoint/
interior cases, sample count, seam/singularity/boundary probes, every carrier
tolerance (absence means bit identity), final-framebuffer tolerance, and the
base SHA/toolchain used to measure it. Each oracle manifest additionally
records its exact callable, natural-domain metrics, non-floating exactness,
measured baseline results, accepted limits, and limit provenance.

`tools/generate_pullback_manifest_header.py` validates the JSON schema and
generates a build-directory C++ header consumed by native tests; the
cross-checkout harness consumes the same JSON directly. CMake depends on the
JSON and generator, and CI rejects stale or schema-invalid data. There is no
second hand-maintained C++ manifest.

The Phase P harness builds isolated base and candidate checkouts with the same
pinned compiler/toolchain and ccache disabled. It captures every manifest entry
at both roster resolutions for:

- native Debug/reference-oracle configuration;
- release WASM with the shipping optimization/LTO flags;
- parameter defaults and the endpoint/interior cases named by each Phase P
  program manifest;
- seam/singularity probes and transition endpoints.

Native output is bit-identical. A release-WASM mismatch is rebuilt in the
strict-FP diagnostic configuration described by revision 5
(`-ffast-math`/LTO removed, contraction disabled, NaN/Inf behavior retained).
If strict FP also differs, the migration is incorrect. A release-only licensed
reassociation may be accepted only with a disassembly attribution and the
revision 5 bound: absolute per-channel delta at most 4 in 16-bit linear units
and at most 0.1% of pixels differing in any frame, further bounded by any
tighter stage oracle.

The harness records base/head SHAs, toolchain pins, configurations, capture
manifests, and hashes. Its self-tests are wired through the repository's tools
test convention.

### 12.4 Approximation oracles

Phase P creates and baselines the Peirce fast-square and prepared hue LUT oracle
manifests described in Section 12.3. Their exact callables remain independent
of the optimized functions. Moving metadata or a call site does not reset or
widen a limit. Natural-domain and final-framebuffer metrics must match the
baseline manifest exactly.

Any additional approximate policy promoted by the Section 8.7 census must
bring, in the same phase, an exact callable, deterministic corpus, limits with
provenance, non-floating exactness assertion, and final-framebuffer metric.

## 13. Device and performance acceptance

This migration promises no visual or performance change. It is accepted only
if the final core-backed ShaderBall is non-regressing under the shipped v1
Section 14 capture protocol as explicitly amended by Sections 13.3 and 15.

### 13.1 Cheap gates for every phase

- native Debug and release tests;
- release WASM smoke at both resolutions;
- full-roster `teensy_gate` size/headroom pass;
- local size-trail query recorded in phase notes;
- expected ITCM delta of zero; any nonzero delta blocks until attributed;
- flash delta beyond +/-64 bytes attributed by symbol/section;
- carrier layout and `-fstack-usage` comparison;
- `nm`, `readelf`, and `objdump` checks from Section 11.

### 13.2 Shade-path disassembly

Each production migration phase archives before/after disassembly for all
address-taken wrappers and moved out-of-line leaves. Reordered addresses and
mangled names are expected. For inline stage arithmetic, the instruction
sequence, calls, branch structure, and load offsets shall match. A difference
requires source-level attribution and framebuffer evidence; an unexplained
difference blocks the phase.

### 13.3 Matched device protocol

Phase C runs the v1 Section 14 matched capture protocol, as amended here, on
one pinned Teensy 4.0 with the real segmented 288x144 driver, 600 MHz DWT
cycles, 16-frame windows, 110-second fast choreography, all authored
IDs/transitions, and the last-to-first wrap.

Control and candidate are separately linked images from the same migration
revision, selected by a branch-only compile-time old/core switch so each image
contains exactly one implementation. Captures are run in `A->B`, `B->A`,
`A->B` order with one unrecorded complete wrap after every flash. Each stream
is arm-stamped and validated before comparison. The switch and legacy arm are
removed before landing.

Phase P adds profile-only telemetry before this protocol is usable. Each image
emits exactly one boot record:

```text
Pullback arm: LEGACY|CORE|LANDED sha=<short-sha>
```

The branch build supplies `LEGACY` or `CORE`; a normal build defaults to
`LANDED`. ShaderBall also emits an event whenever the tuple changes:

```text
Pullback program: preset=<i>/<N> pipeline=<InversePipelineId|NONE> \
endpoint=<steady|from|to>
```

The event is profile-only and deduplicated by the last emitted tuple. It is
generated from the selected/prepared endpoint, not reconstructed from the
current preset after the fact. Thus both frames of a two-frame through-clear
transition are observable even though timing counters aggregate in 16-frame
windows. Presence/identity of each transition half is a coverage assertion;
the surrounding window remains the timing sample.

`tools/parse_profile.py` gains `--expected-pullback-arm` and
`--shaderball-program-manifest` options for `validate`. The protocol command is:

```text
python tools/parse_profile.py <capture.log> validate --scope frame \
  --expected-pullback-arm <LEGACY|CORE|LANDED> \
  --shaderball-program-manifest <generated-manifest.json>
```

Validation rejects a missing/duplicate/wrong arm stamp, SHA mismatch with the
capture metadata, an unknown program ID, a preset/program mismatch, a missing
steady program, either missing transition half, an unexpected dynamic `NONE`,
or a missing last-to-first wrap. Parser self-tests cover every rejection and a
complete synthetic cycle. These records supplement the existing preset,
monotonic-frame, counter-tree, and root/wall checks; they do not replace them.

The normative root-cycle/wall-sum agreement threshold is **5 ppm**. This
supersedes both v1's unattainable 0.7 ppm text and the current parser's overly
permissive 100 ppm check. The value is grounded in the accepted shipping
captures (3.5 and 3.7 ppm) with margin for integer/serial measurement
quantization. Phase P changes the parser and its boundary tests to 5 ppm and
records a fresh unmodified-base capture that passes. A base capture above 5 ppm
blocks migration and requires an evidence-backed spec amendment; it is not
silently grandfathered.

Acceptance uses the v1 capture shape with a behavior-neutral, one-sided
non-regression rule. For each metric, the noise band is
`max(max(control) - min(control), min(0.3 ms, 0.5 * control median))` across the
three control captures. For ShaderBall's own render and clean-shader metrics,
candidate median and maximum shall each be at or below the corresponding
control value plus the band. For every other roster ID, candidate maximum
shall be at or below control maximum plus the band. A faster candidate passes;
there is no material-gain requirement for a code-ownership refactor.

Additionally:

- every program ID and transition half appears as expected;
- the full cycle has zero new spill, cadence, or deadline failure;
- the full-roster memory and section gates pass;
- a final separately linked landed image receives a confirmation capture.

A same-image dual-arm build may be used to attribute flash-layout effects but
does not replace separately linked acceptance. Any layout exception follows v1
Section 14 and is recorded with wrapper/leaf addresses and disassembly.

## 14. Migration plan

Phases land in order `P -> A -> B -> C -> D`. Each phase is independently
reviewable and leaves the tree buildable. Work may be prepared in isolated
branches, but integration is serial so carrier/API changes have one owner.

### Phase P - verification baseline

- Create `tests/data/pullback/programs.json`, the Peirce and combined-hue
  oracle manifests, their JSON-schema validator/header generator, CMake
  dependency, and CI stale/schema checks from Sections 12.3-12.4. Populate
  measured results and limits from the unmodified base implementation before
  any core migration.
- Record the base SHA, native/WASM/Teensy toolchain pins, current manifest and
  selectable-topology census, carrier layouts, wrapper/leaf ELF census,
  `-fstack-usage`, full-roster size/headroom, and accepted shipping device
  profile.
- Implement or finish the isolated two-checkout framebuffer harness and strict-
  FP diagnostic path from Section 12.3.
- Add all-program, both-resolution capture entry points if the current harness
  lacks them.
- Add the profile-only pullback arm/program/transition telemetry, generated
  manifest, `parse_profile.py` validation options, capture-metadata SHA check,
  5 ppm root/wall threshold, and parser self-tests specified in Section 13.3.
  Wire the branch-only arm and short-SHA build defines without adding either
  old/core switch to master.
- Self-test toolchain mismatch refusal, per-checkout build isolation, manifest
  completeness, and deterministic replay.

Phase P is additive verification infrastructure. It changes no render path.

### Phase A - carriers and composition core

- Create `core/render/pullback.h` with metadata enums, standard carriers,
  `Pipeline`, recursive evaluation, validation predicates, and
  `NoInstrumentation`.
- Alias ShaderBall's carrier and metadata names to core, propagating the
  `MaterialSample::path_length` rename.
- Define `ShaderBallBinding`, its extra validation hook, and the instrumentation
  adapter.
- Re-express private `InversePipeline` as a core alias while retaining the
  existing private concrete stage bodies temporarily. Each temporary stage
  gains the two type aliases required by the core contract:
  `using Binding = ShaderBallBinding` and
  `using FrameState = typename Binding::FrameState`; its `run` body and
  metadata remain unchanged.
- Add and wire `tests/test_pullback.h`; re-point ShaderBall white-box validator
  access through a narrow shim.

Gates: Sections 12.1.1, 12.3, 13.1, and shade-path disassembly. Manifest rows
and wrapper count remain unchanged.

### Phase B - concrete catalog

- Add the Section 8 stage combinators, provider concepts, direct operator entry
  points, and required concrete policy inventory.
- Implement core counterparts of reusable effect-private pure helpers in their
  correct core math/color or pullback owner. The legacy helpers remain
  temporarily unchanged because production ShaderBall still uses them in this
  phase; temporary duplication is deleted in Phases C-D.
- Build the Section 8.7 enum census and tests.
- Add synthetic-frame tests for every family and deterministic formula tests
  for every operator; add missing approximation oracles before admitting an
  approximate policy.
- Keep production ShaderBall on its temporary private stage bodies during this
  phase. Direct test/reference calls may compare the new core operators with
  the old helpers.

Because the catalog is uninstantiated in production, the empty-framework ELF
comparison must show no device text/data change beyond explicitly shared
helpers already used elsewhere.

Gates: the complete Section 12.1 suite, the Phase P manifests/oracles, Section
13.1, and the empty-framework ELF comparison above.

### Phase C - ShaderBall stage migration

- Add narrow state providers for every selected outer-camera, surface, lens,
  projection, warp-slot, source, material, and color policy.
- Replace all 11 private compiled stage bodies with core stage aliases plus
  empty topology-matching wrappers.
- Route the simulator dynamic/reference dispatcher through the same core direct
  operator entry points.
- Preserve program wrappers, manifest rows, IDs, keys, ABI, placement,
  instrumentation, and draw site.
- Retain a branch-only legacy/core switch solely for cross-commit and device
  acceptance, with exactly one arm per shipping-style image.
- Run all Sections 12 and 13 gates, including the full matched device protocol.

Phase C does not land until every program and authored transition passes. A
policy that cannot preserve semantics or device acceptance blocks Phase C; it
is not silently left private, because mixed duplicate ownership would defeat
the catalog contract. The exception requires revising this specification with
the precise `EFFECT_ONLY` rationale.

### Phase D - cleanup and closure

- Delete legacy private coordinator, carriers, concrete stage bodies, moved
  helpers, branch-only switch, and migration-only tests.
- Prove no duplicate formulas and no legacy renderer code in either Teensy or
  WASM release targets. Prove that no simulator-only dynamic dispatcher is
  present in the Teensy image. The release-WASM dynamic backend remains
  available and contains only switches/orchestration over core operator
  kernels.
- Run the final landed-image confirmation capture and full authored cycle.
- Update README architecture text, the v1 amendment, ShaderBall spec
  abstraction boundary, API documentation, and profile/census records.
- Record the final core catalog and any approved `EFFECT_ONLY` entries.

Rollback before Phase D selects the last accepted pre-migration image. After
Phase D, rollback is the ordinary commit revert; there is no permanent runtime
or build option retaining the old private implementation.

## 15. Amendments to existing specifications

When Phase A lands, the v1 spec is provisionally amended:

- Sections 8.2 and 17 no longer require private nested coordinator/validation
  types or a second consumer before promotion;
- the coordinator, standard carriers, emission metadata, and approximation
  metadata live in `core/render/pullback.h`;
- `InversePipeline` is a ShaderBall alias of `Pullback::Pipeline`;
- the stage signature remains the two-argument
  `run(input, const FrameState&)`; revision 5's unshipped `Context` parameter is
  discarded;
- validator rules are pinned by constexpr predicate tests rather than the
  aspirational negative-TU harness described in v1;
- all program-set, state, selection, and size obligations remain in force;
- for this ownership-only migration, v1 Section 14's capture mechanics remain
  in force but its comparison rule is superseded by Section 13.3: the floored
  noise band and one-sided median/maximum non-regression rule apply, and no
  material-gain threshold is required;
- v1 Section 14's root-cycle/wall-sum threshold is superseded by the 5 ppm
  requirement in Section 13.3, which Phase P enforces consistently in the
  parser, parser tests, and new baseline captures.

When Phase D lands, `shaderball_spec.md` Section 0.10's abstraction threshold
and ownership table are amended:

- concrete pullback stages and operators that pass the Section 1 core-admission
  test are engine vocabulary even with one current production consumer;
- projected-plane warps, source functions, projection slots, material shaping,
  and color kernels are no longer categorically effect-local;
- ShaderBall still owns animation, authoring, frame preparation, resources,
  topology selection, and the fixed two-slot authoring model;
- the typed operator family's "prepared params + immutable input -> pure
  result" convention is implemented by core policies plus consumer state
  providers.

Historical migration records are not rewritten. A short status note points to
this specification as the superseding architectural decision.

## 16. Required deliverables

Implementation is incomplete without:

- `core/render/pullback.h` containing the composition core, carriers, provider
  concepts, stage combinators, and concrete catalog;
- a complete Section 8.7 operator census;
- ShaderBall aliases/providers/matchers with no rendering formula duplicated;
- unchanged manifest and selectable-topology censuses;
- dynamic/reference reuse of core operator entry points;
- validator, provider, carrier, operator, stage, catalog, selection, seam,
  transition, and framebuffer tests;
- checked-in per-program equivalence manifests and approximation manifests for
  every approximate core policy, all consumed from the same JSON artifacts by
  tests and harnesses;
- native and WASM cross-commit capture evidence at both resolutions;
- carrier-layout, stack, ELF, section, clone, and wrapper-count reports;
- full-roster size/headroom gates;
- matched Teensy device captures and landed-image confirmation;
- v1 spec, ShaderBall spec, README, and API documentation updates;
- deletion of migration switches and private duplicate implementations.

## 17. Explicitly deferred work and revival triggers

The following remain deferred:

- a shared acquisition/reduction/finish `Program` abstraction;
- `PixelDriver` and `BlockDriver`;
- GS, BZ, and Voronoi migrations;
- a core runtime compiled-program manifest/selector;
- arbitrary stage graphs or variable role sequences;
- a spawned-entity `PlaneTransformer` lifecycle abstraction.

They may be proposed independently when a concrete need appears:

- a second renderer needs the same coarse acquisition or block driver;
- two effects share loop orchestration that existing `Scan` entry points do not
  already provide;
- a second effect needs runtime selection among multiple compiled pullback
  topologies;
- independently spawned planar transforms need shared prepare/reclaim
  lifecycle;
- a real consumer cannot be expressed by the six-role carrier contract.

Such a proposal builds on the core operators here. It shall not be justified
merely by classifying unrelated inverse lookups under one vocabulary.
