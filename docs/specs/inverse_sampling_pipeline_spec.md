# Variadic inverse-sampling pipeline

**Status: SHIPPING.** This specification defines the sole shipping renderer for
ShaderBall-style pullback rendering. The equivalence, size, and device-timing
gates below have passed and the monolithic runtime-dispatch renderer is gone
from the device and WASM release images: it survives only behind
`HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES` as the host-test oracle every
authored preset is compared against.

## 1. Purpose

ShaderBall answers one question for every visible sphere sample:

> Which authored source sample contributes to this view direction, and what
> final straight-alpha color does it produce?

The current implementation has a generic runtime-dispatched path and several
hand-written static paths. The static paths are substantially faster, but each
duplicates orchestration across lens, projection, warp, source, material, and
color operations. Adding another specialization requires copying that
orchestration and manually preserving metadata semantics.

The inverse-sampling pipeline replaces both forms with one reusable
compile-time abstraction. The shipping program contains the closed set of
pipelines reachable from the authored roster and supported GUI topology edits.
It contains no generic shading fallback.

## 2. Goals

The design shall:

- express the pullback as typed stages with compile-time ordering checks;
- inline stage boundaries for selected authored topologies;
- preserve exact projection metadata, coverage, deformation, and alpha;
- preserve the current single-sample surface/lens projection behavior;
- resolve every discrete selector represented by a program policy once per
  prepared frame, never per pixel; continuous-value branches already present
  in the selected kernel remain permitted and do not become template
  dimensions;
- exclude the runtime-dispatch renderer and its inactive switch arms from the
  final release link;
- keep the former generic behavior as a host-only reference oracle;
- expose stage traits for code placement, profiling, and specialization
  accounting;
- prevent an unrestricted Cartesian product of template instantiations;
- remain compatible with `Scan::Shader::draw` and its segmented driver;
- require measured frame-time benefit under the full-roster ITCM ceiling.

## 3. Non-goals

This design does not:

- turn ShaderBall into a forward renderer;
- replace `Scan::Shader` or change pixel traversal;
- reuse `Filter::Pipeline` as the implementation type;
- make every GUI combination a template instantiation;
- accept a discrete GUI topology that has no compiled pipeline;
- move continuous parameters into template arguments;
- approximate a stage merely because it is statically selected;
- allocate per pixel, use virtual dispatch, or require RTTI;
- change the authored ordering described by the ShaderBall specification;
- promise that a static pipeline is faster without on-device evidence;
- introduce lens mixing or multi-branch shading in the first version.

## 4. Relationship to `Filter::Pipeline`

`Filter::Pipeline<W, H, Stages...>` is a forward command pipeline. A stage
receives a plotted sample, may emit zero or more downstream samples through a
callback, and eventually reaches a canvas sink. Its traits describe spatial
domain, history, terminal behavior, and segmented-rendering requirements.

The inverse pipeline is a pullback value pipeline. It starts with one view
direction, transforms typed intermediate values toward an authored source,
and returns one `Color4`. It has no canvas sink, history flush, plot callback,
or world-to-screen domain transition.

The designs may share small generic facilities when their contracts are truly
identical:

- compile-time type-list lookup;
- dependent-false diagnostics;
- trait-fold helpers;
- stage uniqueness checks.

They shall not share a base class or stage vocabulary. Forward filters and
inverse stages have different cardinality, lifetime, and terminal contracts.

## 5. Authored order and evaluation order

The authored spatial field order remains:

```text
Source field -> Planar Warp -> Projection -> Surface -> Outer Camera
```

The pullback evaluates the inverse relationship from a view sample:

```text
OuterCamera -> SurfaceProject -> PlanarWarp -> Source -> Material -> Color
```

The template parameter order is the pullback evaluation order. Stage names
must describe pullback behavior and must not imply that authored order changed.
Material shaping and colorization consume the sampled field during shading;
they are not invertible nodes in the authored spatial chain.

## 6. Data model

The first version standardizes the existing ShaderBall carriers rather than
inventing narrower per-specialization structs:

| Stage boundary | Value |
|---|---|
| scan to outer camera | `Vector` view direction |
| outer camera to surface/project | `Vector` outer-local direction |
| surface/project to planar warp | `ProjectedLookup` |
| planar warp to source | `PlanarWarpResult` plus `ProjectedLookup` |
| source to material | signed `float` field plus lookup metadata |
| material to color | `MaterialSample` |
| terminal | straight-alpha `Color4` |

Two aggregate carriers make dependencies explicit:

```cpp
struct SourceInput {
  ProjectedLookup projected;
  PlanarWarpResult warped;
};

struct MaterialInput {
  ProjectedLookup projected;
  PlanarWarpResult warped;
  float field;
};
```

Implementations may pass these values in registers or decompose them during
inlining. Their conceptual fields are part of the contract even if the
compiler removes the aggregate.

### 6.1 Unary cardinality

Version 1 carries exactly one value across every boundary. Current production
shading applies the selected lens once and projects the result once;
`surface_lens.mix`, `join_projected`, and `blend_outputs` do not participate in
that path. Existing GUI compatibility keeps accepting the mix parameter, but
program selection shall not read it or introduce new mixing.

Multi-branch lens blending is a possible later semantic feature, not pipeline
infrastructure hidden in this proposal. It requires separately approved
reference semantics for surface-noise placement on each branch, ordered
endpoint weights, projected-versus-color join rules, and independent
framebuffer tests before the pipeline can model it.

## 7. Stage contracts

Stages are policy types used only through static functions. All mutable and
prepared data lives in the frame object. `InversePipeline` contains no policy
subobjects and is never instantiated, so policy types add no per-effect or
per-frame storage. `std::is_empty_v` may reject accidental data members, but
the ban on mutable static and function-local state remains a review rule.

### 7.1 Common requirements

Every stage shall provide:

```cpp
static constexpr InverseStageKind KIND;
static constexpr CodeEmission EMISSION;
static constexpr bool APPROXIMATE;
static constexpr bool TERMINAL;
using Input = /* exact carrier type */;
using Output = /* exact carrier type */;
```

`KIND` is compared with the expected kind at its semantic position; validation
does not depend on enum ordinal values. `EMISSION` is one of `INLINE_ONLY`,
`OUT_OF_LINE_FLASH`, or `OUT_OF_LINE_ITCM`. It describes emission intent, not
a section owned by a fully inlined stage. `TERMINAL` is false except for the
color stage. `APPROXIMATE == true` requires a non-`NONE`
`ApproximationOracleId` and a constexpr list of named error metrics. Each
metric declares its comparison domain, unit, aggregation (`MAXIMUM`, `MEAN`,
or another explicitly named statistic), and numeric limit. The metadata also
declares whether every non-floating output field is exact and names a final
framebuffer metric. Exact stages declare the `NONE` oracle and an empty metric
list.

Stages shall not retain references to a `FrameState`, allocate storage, mutate
the frame, or read slot tags already fixed by their template type. Continuous
parameters and prepared resources remain frame reads.

### 7.2 Outer camera

```cpp
struct OuterCameraStage {
  static Vector run(const Vector &view, const FrameState &frame);
};
```

The result must match `outer_camera_lookup`. Identity may be a distinct stage
whose `run` returns the input; it is not represented by a runtime branch.

### 7.3 Surface and projection

```cpp
struct SurfaceProjectStage {
  static ProjectedLookup run(const Vector &outer_local,
                             const FrameState &frame);
};
```

This stage owns:

- surface noise placement before or after the lens;
- the selected surface lens;
- projection-frame rotation;
- the selected projection kernel;
- exact region, component, boundary, edge, trait, sphere, weight, and domain
  metadata.

Fusing surface and projection is deliberate. Surface-noise placement and
projection-frame rotation share the selected lensed direction, and keeping the
boundary fused prevents policy combinations that reorder them.

### 7.4 Planar warp

```cpp
struct PlanarWarpStage {
  static SourceInput run(const ProjectedLookup &projected,
                         const FrameState &frame);
};
```

It applies authored warp stages in pullback order, accumulates `net_delta`,
`deformation`, and `path_length`, and preserves the original projected
metadata. A no-warp policy emits a zero-delta result without runtime slot
inspection.

### 7.5 Source

```cpp
struct SourceStage {
  static MaterialInput run(const SourceInput &input,
                           const FrameState &frame);
};
```

It owns source-coordinate conditioning and one selected source function. A
source policy must not read `frame.slots.function`. Noise basis or integrator
choices that alter control flow may be template parameters only when the
program set explicitly instantiates them. A large leaf may remain out of line,
but its discrete algorithm is still chosen by its policy type rather than a
runtime selector.

### 7.6 Material

```cpp
struct MaterialStage {
  static MaterialSample run(const MaterialInput &input,
                            const FrameState &frame);
};
```

It owns signal weighting, value transfer, and coverage. Projection metadata
must remain available until this stage completes. Policies shall preserve the
zero-width edge behavior and out-of-domain coverage of `shape_material`.

### 7.7 Color

```cpp
struct ColorStage {
  static Color4 run(const MaterialSample &material,
                    const FrameState &frame);
};
```

It returns straight-alpha `Color4` under the existing renderer contract. It may
use prepared palette or hue fields owned by the frame's resource state.
`Scan::Shader` performs premultiplication at the final canvas write.

## 8. Pipeline type

The coordinator is genuinely variadic:

```cpp
template <typename... Stages>
struct InversePipeline {
  HS_FLASH_MEMBER static Color4 shade(const Vector &view,
                                      const FrameState &frame);
};
```

Version 1 accepts exactly six stages. The pack form makes carrier validation,
trait folds, instrumentation, and future fused-stage experiments generic while
the validator still requires the six semantic kinds in this exact order:

```text
OUTER_CAMERA, SURFACE_PROJECT, PLANAR_WARP, SOURCE, MATERIAL, COLOR
```

An internal `StageAt<I, Stages...>` selects each type. A recursive carrier fold
invokes `StageAt<I>::run`, requires its exact return type to equal `Output`, and
requires `Output` to equal the next stage's `Input`. `shade` supplies `Vector`
and requires the final output to be exactly `Color4`.

`shade` performs:

```cpp
return run_stage<0, Stages...>(view, frame);
```

`shade` is address-taken and therefore explicitly `HS_FLASH_MEMBER`. Every
recursive coordinator and carrier-fold helper is `always_inline` and must not
produce a symbol. Stage bodies marked `INLINE_ONLY` are also `always_inline`
and execute in the wrapper. Large kernels marked `OUT_OF_LINE_FLASH` remain
`HS_FLASH_MEMBER` leaves. A measured small leaf may be marked
`OUT_OF_LINE_ITCM`, but the wrapper itself is not promoted. Thus "fusion"
means removal of dispatch and intermediate orchestration around selected
stage bodies, not forced inlining of every projection, noise, or color kernel.

### 8.1 Compile-time validation

Validation is a two-level process. Level 1 uses detection-only traits to check
pack length and the presence and types of required members without forming any
`StageAt<I>` outside the pack. Only when Level 1 succeeds may Level 2 form
positions, expressions, and adjacent carrier pairs. Named `static_assert`s
shall then fail with a targeted diagnostic unless:

- the pack contains exactly six stages;
- every stage has the required trait names and static `run` function;
- each stage `KIND` equals its expected semantic kind;
- each `run(Input, const FrameState&)` returns exactly `Output`;
- every adjacent `Output` and `Input` type is identical;
- `SurfaceProject` returns one `ProjectedLookup`;
- only the color stage declares terminal output;
- approximate stages expose a named oracle identifier and a well-formed metric
  list containing a final-framebuffer metric.

Position-specific concepts are documentation and overload aids; the named
asserts are the stable diagnostic surface. Compile-fail tests separately cover
wrong arity, missing traits, wrong order, wrong return type, carrier mismatch,
terminal misuse, and malformed approximation metadata on native Clang and the
shipping Teensy GCC. They pin the named category, not compiler substitution
text.

### 8.2 Location and access

Version 1 originally defined the coordinator, policies, and carrier aliases as
private nested types of `ShaderBall<W, H>`. The pullback-catalog amendment moves
the reusable coordinator, standard carriers, stage metadata, instrumentation
hooks, and operator policies to `core/render/pullback.h`. ShaderBall adapts its
private immutable `FrameState` through an empty binding and narrow state
providers; its program set, topology keys, resource checks, and continuous
preconditions remain effect-owned. Code examples retain their original short
names for readability.

The program set is consequently instantiated per `W, H` specialization. ELF and
WASM accounting must check for duplicate bodies across the production and
test resolutions. The coordinator is parameterized by a consumer binding, so
promotion does not require a universal engine frame or a second shipping
consumer.

## 9. Closed program set and selection

Only explicitly listed program entries are instantiated. The constexpr
manifest contains exactly one entry for every canonical topology reachable
from the shipping preset roster or an enabled discrete GUI edit. It contains
no speculative Cartesian product and no fallback entry. CI compares the
manifest count and IDs with the generated selectable-topology census and the
ELF census of address-taken pipeline wrappers. A new selectable topology must
add one measured manifest entry and one optimization-ledger row.

```cpp
using DodecahedralNoiseGridPipeline = InversePipeline<
    Outer::Animated,
    SurfaceProject::DirectNoiseAfter<
        SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
        Projection::STEREOGRAPHIC>,
    PlanarWarp::MirrorOuter,
    Source::Grid,
    Material::LinearEdgeFade,
    Color::Liquid>;

using DodecahedralNoiseGrid = ProgramEntry<
    InversePipelineId::DODECAHEDRAL_NOISE_GRID,
    DodecahedralNoiseGridPipeline>;

using ShaderBallPrograms = InverseProgramSet<
    DodecahedralNoiseGrid,
    /* every other selectable topology */>;
```

Program types and IDs use semantic names rather than preset numbers. Preset
numbers are unstable and must not appear in type names or selection
conditions. `ProgramEntry<ID, Pipeline>` owns the entry identity, canonical
`TopologyKey`, continuous precondition, resource-readiness predicate, and
non-null `&Pipeline::shade`. Two entries may not use the same ID or pipeline
type.

`make_topology_key(const Config&)` constructs a value-initialized key and
copies every active discrete discriminator owned by a stage:

- function, projection, projection-frame policy, and surface lens;
- signal-weight, value-transfer, coverage, and colorizer policies;
- Peirce layout, Airocean layout, Bonne hemisphere, and gnomonic hemisphere;
- surface-noise kind, placement, basis, and curl integrator;
- source noise basis;
- both warp stages' kind, basis, envelope, polar mode, curl integrator, and
  polar harmonic.

Noise seeds and resource IDs select prepared data but not code paths; policies
read the bound resources and do not key on them. Continuous parameters remain
frame reads. A continuous precondition is present only when a pipeline omits
or changes a branch from the reference kernel. Otherwise the stage reads the
continuous value and preserves the reference branch exactly.

Key equality is memberwise; byte hashing and `memcmp` are forbidden because
padding is not semantic. Canonicalization replaces inactive discriminators
with their enum-zero defaults: non-selected projection layouts and hemisphere
policies; surface-noise placement, basis, and integrator when surface noise is
`NONE`; and, for each warp stage, every discriminator its selected kind does
not read — basis unless the kind samples noise, envelope unless the kind
scales its amplitude by one, integrator unless the kind is curl flow, and
polar mode and harmonic unless the kind is the polar chart. Source noise basis is canonicalized when the selected
source does not use noise. The canonicalization table is exhaustive over
`Slots`, source-noise policy, surface-noise policy, and both `WarpStageSpec`s;
adding a discriminator fails a census test until the table classifies it.

`InverseProgramSet::find(key, config)` runs only when an authored preset is
validated or a discrete GUI edit is proposed. It returns an
`InversePipelineId` only for an exact canonical-key match whose continuous
precondition passes. It has no default result. Authored presets store the
validated ID beside their `Config`; continuous GUI edits retain it, while a
discrete edit is committed only after `find` accepts the candidate. Rejected
edits leave the previous selected configuration active and report
`UNSUPPORTED_TOPOLOGY`.

The GUI derives enabled discrete choices from
`InverseProgramSet::available(current_key, edited_field)` so it normally never
offers an unsupported combination. Continuous controls remain unchanged.
Supporting free Cartesian composition would require compiling that Cartesian
program set and is outside the space budget; topology controls are therefore
curated by the manifest.

The effect carries a selected configuration explicitly:

```cpp
using ShadeFunction = Color4 (*)(const Vector &, const FrameState &);

struct ProgramDescriptor {
  InversePipelineId id;
  TopologyKey key;
  ShadeFunction shade;
  bool (*continuous_parameters_supported)(const Config &);
  bool (*resources_ready)(const FrameState &);
};

struct SelectedConfig {
  Config config;
  InversePipelineId pipeline;
};

struct PreparedEndpoint {
  FrameState frame;
  ShadeFunction shade;
  InversePipelineId pipeline;
  float alpha;
};
```

`prepare_endpoint` verifies that the stored ID still matches the candidate's
canonical key and continuous precondition before preparing the frame.
Successful `prepare_resource_union` remains a caller-side prerequisite. After
preparation, the entry's resource predicate validates every binding or table
the pipeline can dereference, including liquid LUT readiness. Any topology,
precondition, or resource failure rejects endpoint preparation; it never
invokes another renderer.

`FrameShader` receives a non-null `ShadeFunction` and calls it unconditionally
for each sample. This indirect call selects one already-fused compiled
pipeline; it is not a runtime-dispatch shading path. Dispatching once around
`Scan::Shader::draw` would instantiate the complete scan loop for every
pipeline type, so version 1 preserves the single shared scan instantiation and
measures the per-sample call cost.

Debug and test builds attach owner-generation tokens to prepared backing
stores. `FrameShader::operator()` validates the token before dereferencing the
frame. Version 1 retains sequential transition lifetime: prepare one endpoint,
draw it to completion, and only then overwrite shared prepared storage for the
other endpoint. Two simultaneously live prepared endpoints require distinct
backing stores and are outside this contract.

### 9.1 Effect-side usage

Preset declarations bind authored values to a compiled program ID:

```cpp
static constexpr SelectedConfig PRESETS[] = {
    {make_dodecahedral_noise_grid_config(),
     InversePipelineId::DODECAHEDRAL_NOISE_GRID},
    /* remaining authored presets */
};
```

The effect prepares and draws one endpoint through the program set:

```cpp
HS_COLD_MEMBER bool prepare_endpoint(const SelectedConfig &selected,
                                     const LookRuntime &look, float alpha,
                                     PreparedEndpoint &out) const {
  const ProgramDescriptor *program =
      ShaderBallPrograms::get(selected.pipeline);
  if (program == nullptr ||
      program->key != make_topology_key(selected.config) ||
      !program->continuous_parameters_supported(selected.config))
    return false;

  out.frame = prepare_frame(selected.config, look);
  if (!program->resources_ready(out.frame))
    return false;

  out.shade = program->shade;
  out.pipeline = program->id;
  out.alpha = alpha;
  return true;
}

HS_FLASH_MEMBER void draw_endpoint(Canvas &canvas,
                                   PreparedEndpoint &prepared) const {
  FrameShader shader{&prepared.frame, prepared.alpha, prepared.shade,
                     prepared.pipeline};
  Scan::Shader::draw<W, H>(canvas, shader);
}
```

The steady-state draw site is consequently small:

```cpp
PreparedEndpoint prepared;
if (!prepare_endpoint(active_config, runtime, 1.0f, prepared)) {
  report_prepare_failure();
  return;
}
draw_endpoint(canvas, prepared);
```

`FrameShader::operator()` contains no fallback:

```cpp
HS_FLASH_MEMBER Color4 operator()(const Vector &view) const {
  Color4 color = shade_function(view, *frame);
  color.alpha *= alpha;
  return color;
}
```

Transition code repeats `prepare_endpoint` then `draw_endpoint` for each
visible half in the established through-clear order. GUI code submits a
candidate `Config` to `ShaderBallPrograms::find`, with enabled discrete values
enumerated by `available`; on failure it preserves the current
`SelectedConfig`. Startup validates the compiled preset table and uses the
first authored preset if persisted state names an unavailable topology. No
effect call site refers to an individual stage after declaring the program
set.

## 10. State, ownership, and memory

- `FrameState` is immutable for the duration of a draw.
- Resource pointers in the frame borrow effect-owned persistent storage whose
  lifetime exceeds the draw.
- Large prepared tables live in effect persistent state, not `FrameState` or
  the call stack.
- The coordinator stores no policy subobjects; policy types are required to be
  empty as an additional accidental-member check.
- The selector and program set allocate nothing.
- No pipeline may increase the device effect heap, persistent arena, stack, or
  RAM2 use without updating the corresponding budget test and report.
- Device records include `sizeof` and alignment for every carrier. The ARM
  build retains `-fstack-usage` files and wrapper disassembly; emitted fused
  boundaries may not introduce hidden aggregate-return storage, unexpected SP
  frames, or other sret traffic. Device stack watermarking bounds the real ARM
  call chain; the native stack-paint test remains a separate regression gate.
- Transition rendering prepares and consumes one endpoint frame before any
  shared prepared backing store is overwritten. Supporting two live frames
  requires distinct backing stores and a new lifetime test.
- Debug and test builds stamp prepared backing stores and borrowing frames with
  an owner-generation token. Preparing endpoint B invalidates endpoint A's
  token. Tests prepare and shade A, prepare B, prove A is rejected after the
  overwrite, and then exercise the exact through-clear prepare/draw sequence.

## 11. Code placement and instantiation budget

Template fusion is constrained by the full Phantasm image, not the
single-effect profile image.

Resource attribution uses four same-commit, same-toolchain link classes: the
former runtime renderer as a migration baseline, an empty pipeline framework,
a separate one-entry build for every manifest entry, and the complete shipping
program set. Each build retains its linker map, `readelf`, `nm`, and
wrapper/leaf disassembly. Bytes are classified by VMA and section rather than
symbol type. The one-entry-minus-empty delta is that entry's marginal cost;
the full-minus-empty delta is the program-set cost; and replacement savings
are the final shipping image versus the migration baseline. Shared leaves and
read-only data are charged once to a named shared pool. A diagnostic
`-fno-ipa-icf` build, or an equivalent address-alias report, discloses
identical-code folding.

Each program entry shall report:

- unique and shared flash bytes;
- unique and shared ITCM bytes;
- runtime-dispatch renderer bytes retained in the final image, which must be
  zero;
- frame-time delta for every topology using it;
- whether identical-code folding or compiler cloning changed the accounting;
- wrapper and out-of-line leaf VMAs and sizes.

Each address-taken pipeline wrapper is `HS_FLASH_MEMBER`. Large projection,
noise, color-preparation, and other statically selected kernels may remain
out-of-line flash leaves. ITCM placement is permitted only for a measured leaf whose full-roster
build passes `teensy_gate`'s derived RAM1-bank calculation and explicit
headroom ratchets. The gate records ITCM bytes consumed, ITCM bytes free before
the next FlexRAM bank boundary, rounded bank allocation, DTCM variables, DTCM
local free bytes, stack headroom, RAM2 allocator free bytes, persistent-arena
free bytes, and flash bytes. The 196,608-byte bank seen in the current image is
not invariant when DTCM variables change. Unless a separately reviewed
resource reallocation changes the ratchets, the final roster retains at least
648 bytes of ITCM boundary headroom, 12,864 bytes of DTCM local free space,
4,224 bytes of RAM2 allocator free space, and a 12,288-byte stack-headroom
floor. `CodeEmission` must agree with emitted wrapper/leaf attributes and ELF
VMAs; an inline-only stage has no independent byte or placement result.

ELF validation fails on unexpected emitted `run_stage`, carrier-fold,
`.isra`, or `.constprop` bodies and reports their VMA and size. It also fails
when the full program-set wrapper census differs from the checked manifest.

The final ELF must contain no former top-level runtime-dispatch `shade`, no
fallback selector, and no switch helper reachable only from that renderer.
Shared mathematical kernels used by typed stages are not considered retained
renderer code.

The implementation budget is exactly the deduplicated set of selectable
shipping topologies, not every legal slot combination. A new entry is accepted
only when a new topology is intentionally exposed, the full-roster resource
gate passes, and its flash/ITCM cost is recorded in the optimization ledger.

## 12. Instrumentation

Every stage position maps to one existing ShaderBall DWT bucket:

| Pipeline position | Bucket |
|---|---|
| outer camera | `outer_camera` if added; otherwise root residual |
| surface/project | `lens`, `surface_noise`, and `projection` sub-buckets |
| planar warp | `planar_warp` and optional detail buckets |
| source | `source` |
| material | `material` |
| color | `color` |

Instrumentation is compiled out of shipping images. An instrumented pipeline
must preserve the reference renderer's stage boundaries so ratios remain
comparable. Absolute timing from per-pixel instrumented images is diagnostic
only.

The selector shall expose the selected `InversePipelineId` to the profile
harness. A cycle report must prove that every authored topology selected its
declared program. GUI tests must prove that an unsupported discrete edit is
rejected without changing the active program.

## 13. Correctness requirements

### 13.1 Reference-oracle equivalence

Each program entry provides a checked-in deterministic equivalence manifest
naming its corpus generator and version, exact parameter cases, sample count,
seam and singularity probes, and a numeric tolerance for every floating
carrier member. Absence of an explicit tolerance means bitwise equality. The
same manifest compares the static pipeline with the host-only reference at every
relevant stage boundary and at final output. Changing a corpus, sample count,
or tolerance requires review with the entry; a wrapper receives no additional
tolerance merely because it changes inlining.

Tests shall cover:

- all program-set topology keys;
- the exact endpoint and interior continuous-parameter cases named by each
  entry's equivalence manifest;
- projection cuts, singularities, regions, components, traits, and edge
  classes;
- the accepted but shading-inert lens-mix parameter;
- surface noise before and after the lens;
- zero and nonzero warp strengths;
- zero-width and nonzero-width edge coverage;
- zero-alpha and partial-alpha color output;
- transition endpoint preparation and unsupported-topology rejection.

The square Peirce surface/project policy is approximate because
`peirce_projection_fast_square` has nonzero error against the exact Peirce
kernel. It must declare the existing exact Peirce oracle, coordinate and edge
budgets, and a final-framebuffer budget even when the host reference and the
pipeline call the same fast kernel. Reference equivalence alone is not an
approximation oracle.

Approximation classification follows the implementation called by the policy,
not only the difference introduced by the wrapper. Every liquid-color policy
that may sample `PreparedLiquidHue` is approximate and references the direct
hue-rotation oracle, including maximum-channel, aggregate-channel, and
final-framebuffer metrics. Sharing the same approximation with the host
reference does not make the stage exact.

### 13.2 Approximate stages

An approximate stage requires:

- an exact reference callable independent of the optimized implementation;
- a checked-in deterministic oracle manifest naming its corpus generator,
  sample count, seams, and domain-boundary cases;
- its complete named metric list, including natural-output and framebuffer
  domains;
- an explicit exactness declaration for non-floating metadata;
- an identifier linking the program entry to its oracle test.

Approximation approval belongs to the stage, not `InversePipeline`.

Selector predicate tests are table-driven. Every continuous precondition is
tested at its accepted value, at the adjacent representable value on each
side, at domain endpoints, and with NaN and infinities where the carrier can
hold them. Tests assert both the selected ID and final equality with the host
reference. Unsupported GUI topologies must leave the previous selected
configuration active and must never call an unrelated pipeline wrapper.

### 13.3 Cardinality regression

Selector and framebuffer tests shall prove that every version-1 pipeline
evaluates one surface/project sample per view direction. Dormant join helpers
must not become reachable as an incidental consequence of the refactor.

## 14. Performance acceptance

Timing uses the real segmented 288x144 driver, 600 MHz DWT cycles, 16-frame
windows, the 110-second fast-cycle choreography, two-frame through-clear
transitions, and a 1,400-revolution epoch. Candidate and control use the same
source revision, compiler, flags, roster manifest, prepared inputs, and phase.
A pair is two consecutive complete captures, one per arm. Execute three pairs
in order `A -> B`, `B -> A`, `A -> B`, yielding exactly three validated
110-second captures per arm. Run one complete unrecorded roster wrap after
each flash or reset. An invalid capture invalidates its pair, which is repeated
in the same order. When both implementations fit in one image, select the arm
once per prepared frame and use that image for attribution. Final acceptance
repeats the protocol with separately linked shipping images because final
flash layout is itself part of the result.

`parse_profile.py validate` must report monotonic frames, no epoch reset,
complete root and stage counter trees, and root-cycle/wall agreement within
0.7 ppm. A generated roster manifest is the capture oracle: every authored ID,
both visible halves of every transition, and the last-to-first wrap must appear
with the expected compiled pipeline ID. Reports retain per-ID median,
p95, maximum render and clean-shader time, spill count, and transition
ownership. Observations are paired by pair index. For each per-ID metric, the
noise band is `max(control) - min(control)` across the three control captures.
Candidate gain is the median of the three control observations minus the
median of the three candidate observations. A specialization is materially
positive only when its clean-shader gain is at least the greater of 0.50 ms or
twice that noise band. Every per-ID candidate maximum must remain at or below
the control maximum plus its noise band, and the final linked image must meet
the 62.5 ms display deadline with zero spills.

A program entry is accepted only when all of the following hold:

1. Host-reference equivalence and candidate fixed-preset captures validate.
2. Candidate render and clean-shader distributions pass the matched-capture
   non-regression and material-gain rules above.
3. The full authored shipping cycle and last-to-first transition have no new
   spill, cadence, or non-target-topology regression after final linking.
4. The full Phantasm size gate passes with required ITCM, RAM1, RAM2, stack,
   and persistent-arena headroom.
5. ELF inspection confirms expected section placement and no unexpected IPA
   clones or duplicate inactive arms.
6. Native stage, pipeline, seam, stack, arena, and case-census tests pass.
7. Release WASM smoke passes at both roster resolutions.

Global O3 is a diagnostic comparison only. It cannot override a regression in
the shipping selective-O3 image.

## 15. Migration plan

### Phase A: framework and reference harness

- Move the monolithic renderer behind a host-test-only reference interface.
- Add stage concepts, validated variadic `InversePipeline`, canonical topology
  key, program entry, and closed program set.
- Implement one identity/no-op pipeline and compare it with the host reference.
- Add compile-fail tests for stage order and terminal contracts.

### Phase B: complete the selectable program set

- Re-express every deduplicated shipping and GUI-selectable topology as a
  program entry, beginning with the existing static specializations.
- Keep the old static wrappers only in migration builds until their program
  replacements pass stage-boundary, framebuffer, and matched-device gates.
- Generate the selectable-topology census and require an exact manifest match.

### Phase C: effect integration

- Store `InversePipelineId` in every authored `SelectedConfig`.
- Route preset, transition, and GUI candidate preparation through the program
  set as shown in Section 9.1.
- Add selection and rejection tests for every authored preset and supported or
  unsupported GUI topology.
- Emit pipeline IDs in profile markers.

### Phase D: remove the runtime renderer — executed

- Remove the top-level runtime-dispatch renderer and obsolete hand-written
  wrappers from all release targets.
- Keep the reference renderer in host tests only.
- Prove zero retained renderer bytes by the final ELF/WASM symbol census, then
  run the full acceptance protocol.

## 16. Required deliverables

An implementation is incomplete without:

- the framework and stage policies;
- the closed program set and exact selectable-topology census;
- the effect integration API in Section 9.1;
- compile-time contract tests;
- stage-boundary and final-color equivalence tests;
- selection, rejection, preset-census, and transition tests;
- stack, arena, full-roster size, native, and WASM gates;
- shipping fixed and full-cycle device captures;
- ELF symbol/section accounting;
- an accepted/rejected ledger entry for every new program entry;
- ELF and WASM proof that the runtime-dispatch renderer is absent.

Rollback republishes the last accepted pipeline-only image. There is no
shipping build option that restores the runtime-dispatch renderer.

## 17. Implementation decisions

The pullback-catalog amendment places its type-list, validation helpers,
standard carriers, and policies in `core/render/pullback.h`. It remains
independent of `filter.h`; consumers bind their own frame type and
instrumentation adapter. Public policies are admitted by reusable semantics,
not by the number of current shipping consumers.

The native CMake suite runs generated negative translation units through a
scripted `try_compile`, requires compilation to fail, and matches the stable
diagnostic category. The same corpus is compiled by the shipping Teensy GCC.

A same-image migration/reference selector is an attribution aid, not an
acceptance dependency and never ships. Use it only when the compiler retains
both arms and the report records their addresses. The separately linked
matched shipping protocol in Section 14 remains mandatory and decides
acceptance.
