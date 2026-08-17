# Shader workbench and fixed-pipeline effects

**Status: LANDED, revision 9 (2026-08-17), except §10.2 and §10.3 (code
generation and generated-file layout), which are specified but deferred — see
the status note on §10.** The migration shipped in
`69d4751c`. This document defines the authoring and product architecture that
follows the extraction of the pullback pipeline and operator catalog. The
underlying formulas and rendering behavior remain specified by
[pullback_pipeline_spec.md](pullback_pipeline_spec.md),
[inverse_sampling_pipeline_spec.md](inverse_sampling_pipeline_spec.md), or
[shaderball_spec.md](shaderball_spec.md); Section 14 records the completed
migration.

## 1. Decision

Holosphere shall separate pullback exploration from shipped effect identity:

```text
Pullback stage and operator library
                |
                v
Shader: dynamic, snappable authoring workbench
                |
                | freeze/export
                v
Fixed-pipeline effect: fixed structure, morphable presets
```

`Shader` is the user-facing name of the development tool. Its implementation
name shall be unambiguous in code, for example `ShaderWorkbench` or
`PullbackShader`. It is not a shipping effect in the normal firmware roster.

A fixed-pipeline effect is an ordinary first-class `Effect` whose concrete
`Pullback::Pipeline`, compact parameter type, preset bank, and transition/state
policy are declared together. Presets update that effect's parameters; they do
not select stages. Structural changes are authored by snapping in Shader and
are exported as a different effect, not disguised as a preset morph.

There is no `Family` runtime abstraction, `family_id`, or second identity layer.
The effect's stable `effect_id` is its only shipping identity. The term
"family" is reserved for informal product discussion and has no schema or code
meaning.

`ShaderBall` remains only as the dynamic implementation behind the
simulator-only `Shader` workbench and as a legacy WASM alias. It is excluded
from firmware, where the retained looks are first-class fixed-pipeline effects.

## 2. Motivation

The original ShaderBall abstraction claimed that substantially different
looks were points in one morphable parameter space. Its compiled form no
longer has that property. At the 2026-08-15 census, 14 presets select 13 exact
compiled pipelines; only the two sinusoidal curl presets share a pipeline
identifier.
Even a shared compiled pipeline is not sufficient for a morph: ShaderBall's
stable-topology rule additionally requires identical discrete slots and
resource-defining fields.

Consequently ShaderBall currently contains a second effect registry inside an
effect:

- `TopologyKey` describes an implementation structure;
- `InversePipelineId` and `ProgramDescriptor` form a closed program manifest;
- runtime selection finds the compiled program for a configuration;
- unsupported combinations exposed by the broad configuration schema are
  rejected in release builds;
- topology-changing preset transitions fade through clear and never
  structurally interpolate.

This is useful workbench behavior but a misleading shipping-effect boundary.
The proposed split makes both sides honest: Shader explores structures and
snaps between them; each exported effect has a fixed structure and genuinely
morphs its presets.

## 3. Terminology

- **Operator**: one reusable mathematical operation from the `Pullback`
  catalog, such as a lens, projection, planar warp, source, transfer, or
  coverage policy.
- **Stage graph**: the ordered operators and any explicit branch/join nodes
  that define one pullback computation.
- **Effect semantic descriptor**: authoring/build-time canonical metadata for
  every immutable fact that determines graph evaluation or the path between
  two parameter values: stage graph, parameter and interpolation schema,
  clocks, preparation, resources, approximation policy, serialization, and
  state handoff. It does not contain presets or published transition edges and
  does not create a runtime object or identity separate from the effect.
- **Preset bank manifest**: the mutable, transactionally generated collection
  of preset IDs and values, directed transition edges, per-edge easing and
  duration, absent-edge behavior, dwell/order choreography, and generated array
  order for one effect. It schedules descriptor-defined paths but cannot
  define a new interpolation path.
- **Continuous parameters**: finite numeric values that a fixed graph can
  interpolate without changing operator identity, carrier shape, resource
  topology, or control schema.
- **Preset**: one immutable `preset_id`, display metadata, and a concrete
  continuous-parameter value belonging to exactly one effect.
- **Endpoint identity**: one stable `effect_id` for any roster effect, with
  optional descriptor, preset-bank, preset, and custom-parameter data when the
  endpoint is a fixed-pipeline effect.
- **Snap**: an atomic structural replacement with no claim of visual
  continuity.
- **Freeze/export**: validation and deterministic conversion of a Shader
  document into a fixed-pipeline effect or a preset for an existing effect.

## 4. Goals

The design shall:

- provide a dynamic authoring tool capable of composing the public pullback
  operator catalog without adding every composition to firmware;
- make discrete edits immediate and explicit rather than pretending they are
  lerpable parameters;
- support ordered collections such as multiple sequential lenses and planar
  warps in the authoring model;
- retain one serializable, reimportable source document for every exported
  effect and preset;
- deterministically decide whether an export belongs to an existing effect or
  requires a new effect structure;
- generate concrete typed pipelines without runtime discrete dispatch in the
  per-pixel shipping path;
- let each compiled effect expose only the controls meaningful to its fixed
  structure;
- preserve one-endpoint-at-a-time fade-through-clear behavior for changes
  between structurally different effects;
- prove that each compiled export matches its dynamic preview;
- avoid a universal frame type, universal parameter superset, or Cartesian
  product of template instantiations;
- retain the current device timing, ITCM, flash, RAM, stack, and arena gates.

## 5. Non-goals

The first version does not:

- ship Shader's dynamic dispatcher in normal firmware;
- make arbitrary runtime graphs available on the device;
- require every experimental Shader document to become a first-class effect;
- define a general-purpose GPU shading language;
- make structural snapping visually continuous;
- promise that two arbitrary lenses, projections, sources, or branches can be
  interpolated coherently;
- require a fully variadic or graph-shaped `Pullback::Pipeline` before effect
  extraction begins;
- generate hand-tuned code placement or accept a performance regression merely
  because an effect was generated;
- place presets with different effect semantic descriptors in one effect;
- use generated C++ as the only editable source of an authored pattern.

## 6. Fixed-pipeline effect identity and morph contract

A fixed-pipeline effect declares its concrete pipeline stages, compact
parameters, preset bank, and transition/state policy in its effect spec. The
canonical effect semantic descriptor is the authoring and generator form of
those declarations. Exact descriptor equality decides whether Shader may add
a preset to an existing effect; similar appearance or reuse of a shade function
does not. Semantic facts include:

- operator identities, order, and branch/join shape;
- projection layout and hemisphere policies;
- lens, surface-map, warp, source, material, coverage, and color policy types;
- noise bases and numerical integrators selected as compiled policies;
- resource kinds, counts, binding roles, alias partitions, seeds, and stable
  resource identities whenever they cannot transition continuously;
- discrete parameter facts that affect the generated algorithm, carrier, or
  safe interpolation path;
- every parameter's storage representation, binding, semantic unit, admitted
  domain, and interpolation trait;
- versioned path policies that map one edge progress value to each parameter
  or parameter group, including any ordered staggering;
- clock identities, domains, update equations, wrap rules, sharing, pause/reset
  behavior, parameter bindings, and the transition-progress clock/unit;
- frame-preparation and provider policy identities;
- persistent-state and cross-effect handoff policy identities;
- serialization schema and field identities;
- exact or versioned approximation policies and their oracle identities.

Names, descriptions, preset values, choreography, and performance placement
are not part of semantic equality. Whether an otherwise available continuous
parameter is shown in a particular authoring UI is effect-owned presentation
metadata and does not create a new effect.

Every document field shall be classified exhaustively as one of:

- **semantic**: participates in the effect semantic descriptor;
- **preset**: a value interpreted by a semantic parameter binding;
- **effect metadata**: effect-owned display and control-exposure policy;
- **product metadata**: gallery grouping and scheduling policy that does not
  affect the effect's morph contract;
- **study metadata**: authoring-only suggestions that cannot change an
  existing effect during `Add preset`.

Seeds, resource identities and aliasing, integer counts, enum-like values,
constrained coefficients, and clock topology are semantic unless their
operator declares and tests an explicit continuous transition. Palette mapping
is the first such enum-like exception: it is preset state in the shared color
stage, not pipeline identity. A schema census shall fail when a field lacks a
classification.

Every transitionable parameter declares one versioned interpolation trait.
Version 1 scalar storage is binary32. Palette mapping uses enum8 storage with
the closed values `CUP`, `BELL`, `LINEAR`, and `REVERSE`; other integer and enum
fields remain semantic.
For exact endpoints, every reference function returns stored `a` when `t <= 0`
and stored `b` when `t >= 1`. For interior `t`:

- `LINEAR`: `a + (b - a) * t`;
- `LOG_POSITIVE`: `exp(lerp(log(a), log(b), t))`, with both endpoints finite
  and strictly positive;
- `SHORTEST_PERIODIC(period)`: choose the unique delta in
  `[-period/2, period/2)`; an exact half-period tie therefore chooses the
  negative direction, then wrap `a + delta * t` into `[0, period)`;
- `NORMALIZED_LINEAR(group)`: normalize `(1-t)*a + t*b` over the complete named
  vector group and reject a path whose unnormalized norm reaches the catalog
  epsilon. Antipodal endpoints are rejected unless another versioned trait
  defines their path;
- `MIXED_ENUM`: retain both enum endpoints and `t` in the prepared carrier. For
  palette mapping, the color stage evaluates both mapping coordinates and
  linearly blends those coordinates before palette sampling. Equal endpoints
  collapse to the single exact mapping.

Inputs and final results use the common core reference implementation; final
storage conversion is IEEE-754 round-to-nearest, ties-to-even with no implicit
saturation. Safety validation runs over that quantized path. A catalog operator
may add a versioned trait only with a pure reference function, degeneracy
rules, and path validator. Neither Shader nor generated code may substitute an
implementation-local lerp. Easing and duration live in the preset bank
manifest and do not change the geometric path.

The descriptor contains a closed set of versioned path policies. Version 1
provides `PARALLEL` and `STAGGERED_ORDERED(groups)`. A policy maps canonical
edge progress to the progress supplied to every interpolation group; this map,
including group order and slice boundaries, is part of effect semantic
identity. A preset-bank edge references one declared path-policy ID. It cannot
carry an ad hoc stagger flag or other path-shaping data.

Version 1 bank easing IDs are `LINEAR` and `EASE_IN_OUT_SIN`, implemented by
shared pure core reference functions. They map `[0,1]` monotonically into
`[0,1]`, return exact zero and one at the endpoints, and cannot overshoot.
Duration is a positive integer count of the descriptor-declared transition
ticks; version 1 within-effect transitions use one tick per successfully completed
`draw_frame()`. Additional easing or progress-clock semantics require a
versioned descriptor/catalog amendment.

For duration `D >= 1`, a newly started edge renders evaluations `k = 0..D`.
Before evaluation `k`, raw progress is exactly `r = k / D` in the shared
binary32 reference operation, eased progress is `u = easing(r)`, each group
receives `g = path_policy(group, u)`, and its stored value is
`interpolation_trait(a, b, g)`. Thus easing is applied before path scheduling.
The effect then advances its descriptor-defined clocks exactly once using
those evaluated values, prepares the frame, and renders it. Only a successful
completed draw increments `k`; pause leaves it unchanged. Evaluation zero is
the exact source endpoint, evaluation `D` is the exact destination endpoint,
and the edge becomes steady only after evaluation `D` completes. `D == 1`
therefore renders source and destination on two successive completed frames
with no interior evaluation.

Every advertised preset transition within an effect shall satisfy the effect's
continuous-path validator. If the UI or effect API permits arbitrary preset
selection with morphing, every ordered preset pair shall be admitted. An effect
may instead publish an explicit transition graph. Version 1 assigns each
transition origin an absent-edge fallback of `SNAP` or `REJECT`; it never
silently invokes the parameter lerp. A same-effect fade is deferred until it
has an explicit envelope, swap point, clock/state, interruption, and duration
contract.

The set of presets, published directed edges, edge references to declared path
policies, easing, duration, absent-edge fallback, dwell/order choreography, and
generated order belongs to the preset bank manifest and does not participate
in effect semantic equality. The descriptor defines all paths an edge may
select; the manifest defines which edges exist and schedules them.

`Add preset` validates both directed paths between the candidate and every
existing preset for an all-to-all effect. For an explicit graph, the export
request names every new directed edge and no edge is inferred. The generated
artifact contains the complete transition graph and assertions or host checks
for every published edge. Verification samples endpoints and adversarial
interior progress values through the canonical interpolant with live clock and
resource state.

The fixed-pipeline effect boundary is deliberately stricter than visual
similarity. Two looks that appear related but require different stage graphs
are separate effects or separate non-morphing variants in a product group.

## 7. Shader workbench

### 7.1 Availability and ownership

Shader is enabled in the WASM/simulator authoring build and may be enabled in a
dedicated diagnostic build. It shall be excluded from normal device rosters
and shall not pull a runtime stage registry or dynamic switch arms into a
shipping image.

Shader owns mutable authoring state, graph editing, control presentation,
diagnostics, and document import/export. Mathematical kernels remain owned by
the core pullback catalog. The dynamic evaluator and compiled policies shall
call the same kernels rather than maintain effect-private formula copies.

### 7.2 Graph model

The authoring model describes semantics rather than the current ShaderBall
slot layout. At minimum it shall represent:

```text
view -> outer transforms -> sphere-space operations -> projection
     -> planar operations -> source -> material -> color
```

Collections are ordered. In particular, sphere-space operations may contain
zero or more surface maps and lenses, and planar operations may contain zero
or more warps. An operation has a stable operator identifier, compile-time
policy facts, continuous values, resource bindings, and an optional exposed
control definition.

Branching is explicit. A branch node names its child paths and a join node
defines the domain and rule by which they recombine. A graph containing an
unsupported branch may be saved, but it cannot be exported as a shipping
effect until the compiled catalog can represent and validate that join.

The version 1 document grammar is a rooted, acyclic, typed, structured
series-parallel graph with one view input and one color output. Ports are
ordered and named by the operator catalog. A branch owns ordered child regions
and one paired join; shared subgraphs and arbitrary fan-out are rejected. A
child node's canonical identity is its branch-anchor identity plus ordered port
path, and the paired join derives its identity from that anchor rather than
from one of its incoming paths. The grammar also rejects cycles, duplicate
edges, multiple producers for one input, unbalanced branches/joins, unreachable
nodes, and carrier mismatches. Source labels and node/edge declaration order
are non-semantic; canonical node and edge tables are sorted by derived
identity. Symmetric unordered joins are unsupported until a future grammar
defines their canonical labeling.

Phase A evaluates and exports only the linear six-role subgrammar. Known
structured branch regions outside that subset round-trip losslessly with a
`VALID_BUT_UNSUPPORTED` diagnostic and cannot preview or export until their
typed dynamic evaluator and compiled composite land.

Documents are untrusted input. The schema shall bound document bytes, node and
edge counts, nesting depth, string length, parameter count, and requested
resources. Operator and resource kinds are closed catalog identifiers. File
references resolve beneath fixed authoring roots; user strings are never
interpolated as paths, C++ identifiers, or C++ tokens. Unknown semantic fields
round-trip losslessly but block preview and export. Only a versioned extension
namespace declared non-semantic may be ignored safely.

### 7.3 Snapping and editing

A structural edit builds one immutable candidate snapshot containing the
descriptor, parameter values, control schema, prepared bindings, and edit
generation:

1. construct and validate the candidate descriptor;
2. resolve or prepare all required resources;
3. discard the candidate if a newer edit generation now exists;
4. at a frame boundary, atomically publish the entire snapshot only when a
   complete next frame can render;
5. retain the previous snapshot until every render reader has released it.

No intermediate partially rebound graph or mismatched control schema may be
visible. Failed, cancelled, or stale preparation leaves the active snapshot
unchanged. Structural edits snap by default. Shader may offer a preview fade,
but that fade is presentation and does not assert structural morphability.

Continuous edits update values within the active descriptor and may be
smoothed or animated. Changing a field that participates in the descriptor is
a structural edit even when its storage type is numeric or enum-like.

### 7.4 Diagnostics

Shader shall explain, before export:

- malformed carrier or stage ordering;
- missing resources or invalid resource aliases;
- invalid parameter ranges and unsafe interpolation paths;
- an operator unavailable in the compiled catalog;
- a graph unavailable to the current compiled pipeline grammar;
- estimated resource counts and known target restrictions;
- whether the descriptor matches an existing fixed-pipeline effect.

Diagnostics have stable phase and node/port-path reason codes. Schema-invalid,
semantically invalid, valid-but-unsupported, target-unavailable, ambiguous
effect match, and operational export failures are distinct results; none mutates an
authoritative source or generated artifact.

### 7.5 Product integration and document lifecycle

Shader is a dedicated workbench route/build mode, not a normal effect-card
entry. Opening a legacy custom ShaderBall configuration may route into Shader
through the legacy importer. The workbench assigns a document identity, marks
uncommitted edits dirty, autosaves recoverable drafts locally, and requires an
explicit discard or successful save before resolution/effect navigation can
replace them. Reload restores the last committed document plus any newer draft
as a recoverable choice. Unsupported semantic nodes remain visible and
round-trippable but disable preview/export with stable diagnostics.

The release-exclusion gate is defined at the build-feature and link-symbol
level: firmware rosters cannot register Shader, and the final ELF may not
contain its parser, graph evaluator, authoring registry, or generator support.

## 8. Multiple lenses and graph semantics

"Two lenses" has three distinct meanings and shall not be represented by one
ambiguous capacity-two slot.

### 8.1 Sequential composition

For lenses `A` then `B`, evaluation is `B(A(input))`. The compiled catalog
shall provide an allocation-free composition policy equivalent to:

```cpp
Pullback::Lens::Sequence<A, B, ...>
```

This is compatible with the existing fused surface/project stage and is the
minimum required core extension for exporting a sequential multi-lens graph.
Instrumentation may report the sequence as one lens bucket unless a diagnostic
build requests child attribution.

### 8.2 Lens-result interpolation

A blend between `A(input)` and `B(input)` is an explicit blend operator with a
declared normalization and singularity policy. It evaluates both branches.
The exporter shall not infer this operation from two adjacent lenses.
Kaleidoscopic and discontinuous maps may fail visual or mathematical morph
admission even when both endpoints are individually valid.

### 8.3 Independent projected branches

Branches that apply different lenses or projections and sample independently
before joining form a graph, not a linear stage list. Their join must specify
whether it combines directions, projected coordinates, material samples, or
final premultiplied colors. A generic variadic linear pipeline does not by
itself provide this semantic.

The first implementation may reject independent projected branches on export.
It shall preserve them in the authoring document so support can be added later
without changing their intended meaning.

## 9. Authoring document and canonical descriptor

### 9.1 Source of truth

The editable source is a versioned data document, conventionally:

```text
patterns/<effect-or-study>.shader.json
```

A document contains one effect semantic descriptor and one or more concrete
presets. An unpromoted study normally begins with one preset; after promotion,
the effect document is the source of truth for the complete generated preset
bank. Presets with another descriptor cannot be appended to that document.

The exact serialization may change before implementation, but the document
shall contain:

- document and operator-catalog schema versions;
- immutable `effect_id` and `preset_id` values where assigned;
- graph nodes, edges, order, and branch/join declarations;
- stable operator identifiers and compiled policy facts;
- parameter classifications, storage representations, interpolation traits,
  values, units, bounds, defaults, and exposure suggestions;
- clock, preparation, persistent-state, and handoff policies;
- resource keys, seeds, and binding roles;
- one versioned preset bank containing ordered preset records (`preset_id`,
  values, and display metadata), directed edge records (source/destination IDs,
  descriptor-declared path-policy ID, versioned easing ID, positive duration
  and unit), per-transition-origin absent-edge fallback, dwell/order
  choreography, and explicit generated order;
- suggested effect metadata for an unpromoted study, or authoritative effect
  metadata for a promoted document;
- deterministic initial clock/state values used for verification;
- optional authoring-only notes and preview settings.

Generated C++ is derived output and shall carry a do-not-edit marker. Exported
effects must remain reimportable from the data document; parsing generated C++
is not a supported round trip.

### 9.2 Canonicalization

Canonicalization shall:

- resolve operator aliases to stable catalog identifiers;
- apply schema defaults explicitly;
- remove inactive fields that cannot affect the selected operator;
- normalize policy values and resource roles;
- preserve semantically significant operation order;
- assign deterministic node and parameter identifiers;
- reject unknown fields unless a versioned forward-compatibility rule retains
  them losslessly;
- produce byte-stable canonical output for equal descriptors.

Semantic equality is exact comparison of canonical descriptors. A digest may
be stored for lookup and diagnostics, but a hash match is not the authority for
equality. The digest algorithm and descriptor schema version must accompany
the digest.

Continuous preset values, display names, descriptions, and choreography are
excluded when calculating the effect semantic descriptor.

Canonical input is UTF-8 without a byte-order mark. Duplicate object keys,
non-finite numbers, values outside their declared storage representation, and
unknown semantic fields are errors. Object keys are sorted by Unicode code
point after NFC normalization. Arrays retain order only when their schema field
marks it semantic; node and edge declaration arrays are sorted by derived
canonical identity. Numeric values are
converted to their declared fixed-width integer or IEEE-754 binary32 storage
once, before preview, matching, or generation. Binary32 signed zero is
normalized to positive zero unless the operator declares its sign semantic;
generated C++ uses an exact hexadecimal floating literal. The descriptor
format version defines remaining escaping, line-ending, and digest-byte rules.

Document syntax, canonical descriptor format, operator semantic ABI, generator
output, and serialization schema have separate versions. Older documents are
interpreted with their original defaults, migrated to one current canonical
descriptor, and only then compared. A formula or default with observable
semantic change receives a new operator revision or an explicit tested
migration; a spelling-stable operator identifier alone does not imply semantic
compatibility.

A checked effect registry maps each fixed-pipeline effect's immutable
`effect_id` to one source path and maps each descriptor digest plus exact
descriptor to exactly one `effect_id`. Duplicate descriptors under different
effect IDs, duplicate IDs, and ambiguous aliases are validation errors
independent of filesystem order. File moves and display-name changes do not
change identity. Preset IDs are immutable and unique within an effect;
generated ordering is by explicit source order, while persistence, URLs,
captures, and tooling use IDs rather than array positions.

Factories, URLs, saved state, captures, profiles, product groups, and restore
tokens use the general registry's `effect_id`. Fixed-pipeline metadata is
attached to that ordinary registry row; there is no family registry or alias
layer.

Once a preset record is published, its ID and parameter value are immutable.
Changing the authored value creates a new preset ID. A removed preset becomes a
tombstone retaining its last value and legacy aliases for restore/URL migration
but is excluded from new selection and choreography. Bank revisions may add
presets or change edges, scheduling, and display metadata; they may not silently
retarget a persisted preset ID. Custom effect states persist `effect_id`, no
`preset_id`, the descriptor version/digest, and their serialized parameters.

Committed bank manifests and generated restore metadata are retained for every
snapshot version the compatibility policy still supports. A token that pins an
older bank digest restores against that exact retained manifest or a versioned
migration proven to reproduce its preset, pending transition, choreography
position, and next-frame state. `EXACT_STATE` is valid only for those two
outcomes. When support for a bank revision is intentionally retired, restore
returns a stable `UNSUPPORTED_BANK_REVISION` result and may offer an explicit
reset/custom-state import; it may not silently reinterpret the token under the
current bank.

Each registry row also declares promotion evidence keys of `(effect_id,
descriptor_digest, preset_bank_digest, capability_profile_id)`. Factory, URL,
restore, gallery, and roster lookup
distinguish `UNKNOWN_ID` from `KNOWN_UNAVAILABLE`. An unavailable effect remains
first-class metadata but exposes no creator for that profile; restore retains
the requested stable identity and offers an explicit compatible target or
fallback rather than silently selecting another effect.

## 10. Export and generated artifacts

**Implementation status: only §10.1 classification ships.** `classifyExport` in
`scripts/shader_workbench.mjs`, reachable as `shader_workbench_cli.mjs classify`,
performs the non-mutating classification below. The generation, staging, and
publication machinery of §10.2 and §10.3 is unimplemented and deferred: there is
no `generated/` tree, no `export` or `generate` command, and no registry-manifest
transaction. Every shipping fixed-pipeline effect header is hand-written and
hand-owned; its `DESCRIPTOR_DIGEST` and `PRESET_BANK_DIGEST` are maintained by
hand against the authoring document. §10.2 and §10.3 state the intended design,
not current behavior.

### 10.1 Classification and outcomes

Export first performs a non-mutating semantic classification:

1. **Add preset candidate**: the canonical semantic descriptor equals exactly
   one registered fixed-pipeline effect's descriptor.
2. **Create effect candidate**: the descriptor is compilable but matches no
   registered effect.
3. **Rejected**: schema-invalid, semantically invalid, ambiguous, unsupported
   by the compiled grammar, or unavailable for the requested target capability
   profile.

An add-preset candidate becomes valid only after the candidate value and every
new published directed transition pass range and path admission. A create-
effect candidate remains a draft until identity, effect metadata,
choreography, build, verification, and target gates pass. Parse errors, user
cancellation, stale-source conflicts, generator/build failures, and review
rejection are operational outcomes distinct from semantic classification.

Every export request pins a versioned capability profile containing target,
resolution, build features, operator catalog, enabled approximations, and
resource ceilings. Structural compilability is target-independent; promotion
approval is target-specific and records the profiles it passed.

The exporter shall never choose a merely similar effect, discard an operator,
or replace a branch with an approximation without an explicit versioned
operator and approximation oracle.

### 10.2 Generated effect content (deferred)

For a new effect, deterministic generation shall provide ordinary C++ effect
declarations rather than a family runtime abstraction. At minimum it provides:

- the concrete `Pullback` pipeline or composite policy type;
- the compact continuous `Params` definition or schema adapter;
- immutable preset initializers;
- narrow providers mapping the effect frame into core operators;
- compile-time resource requirements and binding validation;
- parameter registration metadata;
- serialization field identifiers and schema version;
- authoring/build artifacts containing the canonical descriptor and digest
  used for future matching; release firmware need not retain the descriptor;
- dynamic-versus-compiled verification manifest entries;
- static assertions that every generated preset belongs to the effect and is
  in range.

The generated structure should remain as direct as the hand-written pipeline
declarations it replaces, for example:

```cpp
struct CurlLatticeSpec {
  using Pipeline = Pullback::Pipeline<
      Binding, OuterCamera, SinusoidalProjection, CurlWarp,
      LatticeSource, Material, Color>;

  struct Params {
    float curl_strength;
    float scale;
    float phase;
  };

  static constexpr Preset PRESETS[] = { /* generated values */ };
};
```

`Pipeline` is resolved at compile time. Presets update `Params`; no preset
stores or selects a pipeline identifier.

The promoted effect document is the sole authority for its preset bank,
including edge and dwell/order choreography. Generated runtime tables are
purely derived from it. The small hand-owned effect wrapper retains the visual
name and description when product policy does not generate them, playlist
membership, unusual lifecycle behavior, and reviewed code-placement decisions.
It shall not shadow preset-bank fields. A separate product playlist overlay may
assign effect-level airtime or grouping, but it cannot change within-effect
paths or preset order without naming a different bank revision.

For an existing effect, the effect document is authoritative for parameter
presentation, defaults, preset-bank choreography, and control exposure. An
imported study contributes values by canonical parameter ID and explicit transition
edges; conflicting study metadata is diagnosed and cannot silently alter the
effect API. A new preset requires a stable ID and a complete choreography
entry, whether explicitly authored or accepted from an effect default.

Export is transactional. Add-preset transactions lock by `effect_id` and
verify the caller's expected source digest. Create-effect transactions hold the
registry-manifest lock across final classification and publication, then repeat
descriptor, ID, alias, and source-path uniqueness checks under that lock; two
equal descriptors cannot publish under different new IDs. Both paths write a
complete draft bundle beneath a staging root and run generation, build, static
validation, and host equivalence against that bundle. They refuse to overwrite
an unexpected hand-owned path. After all required checks and human-owned fields
are complete, publication atomically advances the registry manifest to the new
bundle, using compare-and-swap against its expected digest. Failure preserves
the previous committed bundle and retains the draft with diagnostics.
Generation shall be idempotent and independent of filesystem enumeration
order, locale, address layout, and wall-clock time. Concurrent and
crash-injection tests cover equal-descriptor creation, lost updates, and
partial-publication cases.

Because adding a preset or edge changes the bank digest, publication also
invalidates all prior promotion evidence. The transaction must re-run every
capability profile it intends to keep enabled and atomically publish the new
evidence, explicitly remove profiles not reapproved, or refuse publication.
Host-only generation success cannot preserve stale device availability.

### 10.3 Generated-file layout (deferred)

No `generated/` tree exists; the layout below describes the deferred design. The
implementation may refine locations, but ownership shall remain visible, for
example:

```text
patterns/CurlLattice.shader.json       editable authoring source
generated/CurlLattice.generated.h      deterministic pipeline/preset output
effects/CurlLattice.h                  hand-owned effect identity and policy
tests/data/pullback/CurlLattice/...    exported verification captures
```

## 11. Fixed-pipeline effect contract

Each extracted effect directly declares or includes its concrete pipeline,
parameter type, per-field interpolants and path validators, frame preparation,
clocks, resource traits, parameter schema, state-handoff policy, and presets.
Generated Shader preview and compiled effects call the same reference
interpolants.

No generic effect shell is required for extraction. After at least two concrete
effects exist, genuinely repeated lifecycle, preset-selection, scan,
serialization, or transition code may move into small shared helpers. Such
helpers do not introduce a family identity, runtime descriptor, universal
configuration, or hot-path type erasure.

Each effect keeps a compact immutable frame and exposes only the narrow
providers required by its operators. A shipping fixed-pipeline effect:

- calls its fixed compiled shade function directly;
- stores no active pipeline identifier;
- performs no topology-key search;
- exposes no unsupported discrete stage combinations;
- allocates only resources named by its fixed graph;
- morphs only along admitted continuous paths;
- selects presets by immutable ID and treats generated array position as an
  implementation detail;
- contains no dynamic fallback when the normal release features are selected.

## 12. Pipeline evolution

Effect extraction and pullback-pipeline generalization are separate projects.
The current six-role `Pullback::Pipeline` can express all current compiled
ShaderBall programs and therefore does not block their promotion.

Sequential multi-lens export should first use `Lens::Sequence` inside the
existing surface/project stage. If real authored effects subsequently require
independent sphere-space stages, core may evolve toward an adjacent-carrier-
checked variadic chain. Such a change must retain semantic role validation,
terminal validation, approximation metadata, instrumentation, and code-
placement behavior.

A linear variadic chain shall not claim to solve branches. Branch/join support
requires an explicit typed composite or render-graph contract with bounded
evaluation count and named join semantics.

No pipeline-generalization step is accepted on the promise of simpler compile
time alone. It must reduce source complexity or enable a retained effect while
passing compile-time, code-size, disassembly, and device timing measurements.

## 13. Transitions and effect lifecycle

### 13.1 Within a fixed-pipeline effect

Published morph edges interpolate the effect's continuous values and
persistent rates, never stage identities or accumulator positions. Exact
endpoint values must be produced at transition boundaries. The effect
validator rejects an unsafe continuous path before it begins. Selection of a
preset with no directed morph edge follows the preset bank manifest's explicit
`SNAP` or `REJECT` policy and never falls through to interpolation.

### 13.2 Between effects

Structural changes are controller transactions. They never render two effects
in one frame and never require two live `Effect` instances, but they do require
an explicit output-envelope, presentation-fence, handoff, and failure contract.
The base effect/controller API shall provide a target-independent transition
envelope applied to the complete published output; an effect-private shader
alpha is not that contract.

The logical states are:

```text
STEADY_OUT -> FADING_OUT -> CLEAR_PRESENTED -> CONSTRUCTING
           -> PREPARING_FIRST_FRAME -> PUBLISHING_HIDDEN_FRAME
           -> HIDDEN_FRAME_PRESENTED -> COMMIT_READY
           -> FADING_IN -> STEADY_IN

CONSTRUCTING|PREPARING_FIRST_FRAME|PUBLISHING_HIDDEN_FRAME
  -- recoverable failure --> RESTORING_OUT
RESTORING_OUT -> RESTORE_FRAME_READY -> FADING_BACK -> STEADY_OUT
RESTORING_OUT -- failure --> CLEAR_FAILSAFE
CLEAR_FAILSAFE -- new request --> CONSTRUCTING
```

Before fade-out, the controller validates destination identity, factory,
resolution, static budgets, serialized preset, and handoff compatibility
without constructing a second effect. At `CLEAR_PRESENTED`, it waits for the
target's display/DMA presentation fence before releasing outgoing resources.
It then destroys the outgoing effect, constructs the incoming effect, imports
the admitted handoff, and prepares its first complete frame while output
remains clear. `PUBLISHING_HIDDEN_FRAME` publishes that incoming frame with an
exact zero output envelope, so the display remains clear, and waits for its
presentation/DMA acknowledgement. Only then does
`HIDDEN_FRAME_PRESENTED -> COMMIT_READY` occur. Reaching `COMMIT_READY` means
all fallible work and the incoming presentation fence have succeeded; its
transition to `FADING_IN` commits controller identity, URL, and persisted
selection and is infallible under controller invariants. Failure before
teardown cancels the transaction; failure after teardown remains
deterministically clear and either
reconstructs the outgoing effect from its restore token or follows the
target-specific fail-safe policy.

Every registry effect has an endpoint identity and declares one restore
capability: `EXACT_STATE`, `RESET_ONLY`, or `NONE`. A fixed-pipeline-effect
token may contain effect/preset IDs, descriptor and bank digests, serialized parameters
and pending state, and handoff/persistent state. A generic token always contains
the stable effect ID, RNG seed identity and visit position, animation pause,
transition origin, and any effect-owned serialization supported by its restore
capability. `RESET_ONLY` reconstructs documented initial state; `NONE` requires
the target's clear fail-safe after teardown and must be disclosed by preflight.
Public selection and URL remain the outgoing identity until `COMMIT_READY`
commits the destination; rollback therefore exposes no temporary destination
identity. A restored outgoing effect must prepare and
present one complete fenced frame before `FADING_BACK`. Failure to reconstruct
or prepare that frame enters `CLEAR_FAILSAFE`, whose retry, diagnostic, and
fail-stop behavior is target-defined. A retry is a new request, which the
controller accepts from `CLEAR_FAILSAFE` and resumes at `CONSTRUCTING`: teardown
has already run, output is already clear, and no outgoing effect remains to
restore, so that request's restore capability is `NONE`.

Phase D introduces status-bearing controller adapters for destination
preflight, effect construction/init, handoff import, and first-frame
preparation; the existing `void init()` and `void draw_frame()` APIs alone do
not satisfy this transaction. Recoverable statuses include unavailable target,
invalid restore or handoff, resource-capacity rejection, and first-frame
preparation rejection. Contract violations, memory corruption, and existing
`HS_CHECK` invariants remain fail-stop. WASM reconstructs the captured outgoing
endpoint according to its advertised capability after a recoverable
post-teardown failure; Phantasm remains clear and
follows its synchronization specification's retry or fail-stop decision. No
target may report a committed switch before first-frame success.

The handoff is a versioned semantic state, not a memory-layout copy. It may
carry compatible clocks, orientations and walk state, palette sequence state,
choreography position, and other catalog-defined persistent values. Each
effect declares which fields it exports and accepts. A missing or incompatible
field follows an explicit reset policy. Migration acceptance compares complete
legacy ShaderBall transition sequences when continuity is claimed; a fresh
effect initialization is not considered equivalent.

Manual, automatic, synchronized, restore, and authoring transition origins are
distinct. The policy for each origin specifies snap versus fade, pause
behavior, request replacement, and whether a queued request may supersede the
destination before or after `CLEAR_PRESENTED`. `SYNCHRONIZED` preserves the
base `PresetChangeOrigin::SYNCHRONIZED` semantics unless a versioned controller
policy explicitly maps it elsewhere. Migrated manual ShaderBall selection
retains its current snap/reset behavior unless a separately approved product
change says otherwise.

Phantasm currently performs an epoch-synchronized teardown followed by a fixed
`K`-revolution black construction window and publishes incoming frame zero at
the commit boundary. Implementing fades there requires an amendment to
`phantasm_frame_sync_spec.md` defining envelope timing in revolutions,
repeat-loss and late-join behavior, presentation acknowledgement, interruption,
and commit deadlines. Until that amendment lands, promoted effects remain
off-roster comparison effects and Phantasm retains its existing handoff.

The single-device Holosphere `POVDisplay` adapter shall implement the same
envelope, clear-presentation acknowledgement, status-bearing construction,
teardown, restore, and `CLEAR_FAILSAFE` contracts around its sequential
`show<E>()` lifecycle. If that adapter is not implemented in Phase D, capability
profiles for that target mark promoted effects `KNOWN_UNAVAILABLE` and no
Holosphere roster may name them.

Shader's development UI may snap without this transition. Preview fades are
optional and must not affect the exported morph contract.

## 14. Migration plan

### Phase A - document and export foundation

- Introduce stable effect and preset identities, display metadata, and
  legacy aliases in the registry, URL, save-state, capture, profile, and roster
  protocols before any ShaderBall identity is removed.
- Define the versioned Shader document schema, field census, linear six-role
  grammar, canonical semantic descriptor, and transactional draft export.
- Adapt the existing dynamic/reference pullback backend to evaluate linear
  documents through core kernels. Generic branch code generation is deferred.
- Generalize capture and profiling selection from the literal `ShaderBall`
  name and numeric preset indices to stable effect/preset IDs.
- Convert current ShaderBall presets into source documents without changing
  shipping code.
- Produce a census mapping all 14 current presets, and any later additions
  discovered at Phase A start, to canonical descriptors and retained/retired
  decisions.
- Capture a fresh current-source baseline covering every retained preset and
  every transition edge. Resolve or explicitly waive any already-red preset
  before defining "non-regression."

### Phase B - pilot effect

- Promote the two sinusoidal curl presets, which already share one compiled
  pipeline, into the first off-roster fixed-pipeline comparison effect.
- Prove continuous path admission, serialization, export regeneration, dynamic
  equivalence, resource bounds, and device timing.
- Keep ShaderBall's copies during the comparison period.
- Run full-roster flash, ITCM, RAM1, RAM2, stack, heap, and arena gates even
  while both copies exist. If the duplicate image cannot link safely, keep the
  pilot in a dedicated comparison build rather than weakening a device limit.

### Phase C - second effect and shared helpers

- Promote a structurally different current program as a second off-roster
  comparison effect, including a singleton preset if necessary.
- Factor only lifecycle and preparation behavior genuinely shared by the first
  two effects into small helpers. Keep their pipeline and parameter
  declarations concrete.
- Compare the complete full-roster section map and device timing against the
  unsplit ShaderBall implementation to detect scan, binding, wrapper, and
  template duplication early.

### Phase D - controller transition

- Move one-endpoint fade-through-clear scheduling to the effect/show boundary.
- Define and implement the controller state machine, output-envelope API,
  presentation fence, restore token, versioned handoff, and failure policy.
- Amend the Phantasm synchronization specification before enabling device
  fades; verify effect-local resources are reclaimed only after clear is
  presented.
- Exercise manual, automatic, synchronized, restore, and authoring origins;
  all ordinary-effect and fixed-pipeline-effect endpoint combinations with
  every restore capability; pause and request replacement in every state;
  init/restore failure; lost/repeated sync symbols; and exact transition
  endpoints.

### Phase E - complete promotion

- Promote each retained canonical descriptor as a concrete effect and author
  additional morphable presets where the visual identity warrants them.
- Assign user-facing visual names rather than mechanically exposing full
  pipeline type names.
- Provide a versioned legacy importer that runs before effect-name validation.
  It consumes the legacy effect alias, schema version, accepted and requested
  configurations, pending field IDs, optional preset index, runtime clocks, and
  persistent transition state. Classification is total over every valid legacy
  snapshot:

    - an exact curated value maps to `effect_id` plus `preset_id` and translated
      handoff state;
    - a fixed-pipeline effect descriptor with non-preset values maps to
      `effect_id`, no preset ID, serialized custom parameters, and translated
      handoff state;
    - an accepted configuration with a structurally different pending request,
      or any valid configuration with no compiled effect, opens in Shader with
      accepted, requested, pending, and runtime state preserved separately and
      an explicit notice.

  Invalid legacy snapshots follow their historical rejection/default policy
  and record a diagnostic. URL and persisted identity rewrite only after the
  destination's first successful frame. Retain a `ShaderBall` tombstone alias
  while any supported legacy snapshot version remains.
- Replace positional RNG seeding with a stable seed identity before roster
  expansion, or explicitly approve and capture the resulting changes to every
  later effect.
- Define a product group for the promoted effects. Its gallery discovery,
  child ordering, device duration weights, reachable presets, and total airtime
  are explicit; replacing one 120-second ShaderBall row with many independent
  120-second rows is not the default migration.
- Remove ShaderBall from normal rosters after every retained look and external
  identifier has an intentional destination.
- Retain Shader as the authoring tool and the dynamic backend as its evaluator
  and test oracle.

### Phase F - vocabulary growth

- Add `Lens::Sequence` when the first retained multi-lens document requires
  export.
- Add a variadic carrier chain or explicit branch composite only in response to
  a concrete retained graph and its acceptance tests.

## 15. Verification and acceptance

### 15.1 Authoring and generation

Tests shall cover:

- document parse, validation, canonicalization, and byte-stable regeneration;
- canonical equality under source-label, node-declaration, edge-declaration,
  object-key, and valid topological-order permutations;
- round trips that preserve all semantic and unknown retained fields;
- rejection of duplicate keys, cycles, alias loops, path escapes, oversized
  documents/graphs/resources, raw-token injection, and unknown semantic fields;
- structural equality and inequality for every descriptor field;
- registry uniqueness for effect IDs, preset IDs, descriptors, and aliases;
- preset-bank schema, digest, edge references, absent-edge origin policies,
  easing/duration units, choreography, and generated-order round trips;
- exact v1 transition sequences for `D == 1`, ordinary durations, pause, and
  both endpoint frames;
- retained-bank exact restore, tested bank migration, tombstone lookup, and
  unsupported-bank rejection;
- add-preset versus create-effect classification;
- unsupported-export diagnostics for unknown operators and branches;
- deterministic identifiers independent of document key order;
- sequential lens ordering and the distinction between sequence, blend, and
  branch nodes.

### 15.2 Dynamic-versus-compiled equivalence

Every export adds a deterministic comparison case. The dynamic Shader graph
and generated effect shall receive identical view samples, prepared values,
resources, clocks, seeds, and palettes. Exact stages compare bitwise where
declared by the operator catalog; approved approximations use the exact
descriptor-selected approximation revision, registered oracle, metric version,
and thresholds. The exporter cannot widen a catalog tolerance.

Framebuffer capture shall include every exported preset, transition endpoints,
singularities, seams, projection cuts, lens boundaries, and resource-dependent
paths. Every published directed morph edge is sampled in both persistent clock
evolution and at adversarial interior `t` values, including periodic boundaries
and log-domain extrema. A generated effect is not accepted merely because its
types compile.

Each verification manifest pins stable effect and preset IDs; descriptor and
preset-bank digests; operator, generator, corpus, and oracle versions; target
and feature profile;
resolution and sample ordering; toolchain and floating-point configuration;
seeds, initial clocks, palettes, and resources; exactness classifications; and
all thresholds. Regeneration that changes an approved baseline fails pending
explicit review rather than silently updating it.

### 15.3 Release exclusion

Release ELF inspection shall prove that Shader's graph evaluator, authoring
schema machinery, operator registry, and inactive runtime switch arms are not
linked into normal firmware. Each fixed-pipeline effect must contain only its explicitly
named pipeline and required resources.

### 15.4 Resource and performance gates

The split is accepted only if:

- every promoted preset remains within its device frame-time gate;
- full-roster ITCM, flash, RAM1, RAM2, stack, effect heap, and persistent arena
  remain within their existing budgets;
- no effect allocates the former ShaderBall resource superset without need;
- aggregate wrapper/template duplication is measured and reviewed;
- every enabled capability profile has promotion evidence for the exact
  descriptor and preset-bank digests being published;
- every migration phase that retains duplicate ShaderBall/extracted-effect code passes a
  full-roster link/section gate or remains confined to a diagnostic build;
- selectively placed hot shader bodies retain their measured placement and
  timing unless a new measurement approves a change;
- the controller transition never renders two heavy endpoints in one frame.
- the incoming hidden frame is published at zero envelope and acknowledged
  before destination identity commits.

Conceptual simplification is expected, but no compile-time or binary-size
improvement is assumed without measurement.

## 16. Naming and promotion policy

Pipeline types may retain precise structural names such as
`GnomonicDodecahedralGridWaveMirrorPipeline`. First-class effects should use
short visual identities that remain meaningful as continuous presets are
added. The exporter may suggest a name but a human approves it before roster
entry.

The registry separates immutable machine identity from presentation and C++
type identity. `effect_id` keys save state, URLs, captures, profiles, seeding,
and aliases. Display name and description may change without breaking those
surfaces. `preset_id` is likewise independent of generated array position and
display label.

Not every Shader study is promoted. Promotion requires:

- a retained visual identity;
- a compilable and validated canonical descriptor;
- at least one curated preset;
- a plausible continuous authoring space, even if the first revision has one
  preset;
- dynamic equivalence and target resource evidence;
- an intentional position in galleries and device rosters.

An exact pipeline with one preset may begin as a singleton effect. It shall not
be merged with a structurally different effect merely to avoid a singleton.
Conversely, implementation similarity alone does not require two separately
named user-facing effects when a higher-level product grouping can present
them without claiming a morph.

## 17. Rejected alternatives

### Keep ShaderBall as the permanent shipping container

Rejected because it preserves a broad but mostly unavailable configuration
surface, an internal program registry, and topology changes masquerading as
preset transitions.

### Create one effect per current preset

Rejected as the governing rule. The correct unit is one fixed-pipeline effect;
multiple presets with the same descriptor belong directly to that effect.

### Group by visual similarity and fade internally

Rejected as the effect definition. Product navigation may group related
effects, but an effect advertising morphable presets must not hide structural
replacement behind its preset API.

### Use one universal effect frame and parameter superset

Rejected because it recreates ShaderBall's inactive fields, broad resource
ownership, and runtime selection in a different type.

### Ship the dynamic graph evaluator

Rejected for normal firmware because it adds inactive dispatch arms, weakens
the closed compiled roster, and competes with the hot-path and memory budgets.

### Require a generic render graph before extracting any effect

Rejected because the current six-role pipeline already represents all current
compiled programs. Graph generalization should be justified by a retained
authoring case, especially for branch semantics.

## 18. Specification relationship

[pullback_pipeline_spec.md](pullback_pipeline_spec.md) remains authoritative
for the ownership and contracts of the core pullback catalog. This document
does not move effect lifecycle, mutable resources, or universal frame state
into core. `Shader` is the dynamic authoring consumer, and fixed-pipeline
effects are the shipping consumers.

[shaderball_spec.md](shaderball_spec.md) and
[inverse_sampling_pipeline_spec.md](inverse_sampling_pipeline_spec.md) now
serve as historical compatibility and migration records for the retired
shipping container. Each promoted effect has a canonical generated descriptor
plus hand-owned effect policy defining its behavior.

## 19. Completion record

The landed migration satisfies these criteria:

- Shader saves, loads, snaps, validates, and previews structural documents;
- canonical effect semantics include graph, interpolation, clocks,
  preparation, resources, serialization, approximation, and handoff policies;
- stable effect/preset IDs and legacy aliases replace names and positions on
  persisted and tooling surfaces;
- export classification deterministically separates an added preset from a new
  effect, or rejects it, without changing the committed bundle; the staging,
  generation, and publication transaction of §10.2 and §10.3 remains deferred;
- fourteen structurally different fixed-pipeline effects ship through the
  common contracts without dynamic per-pixel dispatch;
- a sequential two-lens graph can be represented and, once retained, exported
  through an explicit sequence policy;
- dynamic-versus-compiled captures are generated and pass for every promoted
  effect, from its hand-owned header rather than from a generated bundle;
- inter-effect fade-through-clear has a fenced controller transaction and an
  approved Phantasm synchronization amendment;
- all retained ShaderBall looks have migration destinations and compatibility
  mappings, including full configuration and runtime-state imports;
- ShaderBall is removed from normal shipping rosters;
- device timing and full-roster memory gates pass with the final effect set.
