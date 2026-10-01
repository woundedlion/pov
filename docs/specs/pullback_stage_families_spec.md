# Pullback stage families: arbitrary chains over ranked carriers

**Status: §§1–6 and §8 LANDED; §7 PARTIAL.** The static ranked pipeline and
its migration (the §6 contract cut-over) ship, as does the preview
interpreter (§8: `core/render/pullback/interpreter.h`, `workbench/shader/chain_host.h`,
and the `setShaderChain` binding). Promotion/verification (§7) is landed per
sub-section: the operator authority (§7.1), the field-table half of §7.2, and
the manifest half of §7.4 ship; the allocator, binding table, promotion pin,
and roster-derived acceptance registry are design, and each sub-section
carries its own banner. Location, provider, instrumentation, and ownership
contracts are collected in §§10–12. The document is three systems with distinct invariants and
failure modes, layered in
order: **the static ranked pipeline and its migration (§§1–6)**, which
is self-contained; **promotion and verification (§7)**, a resource
allocation and CI layer over it; and **the preview interpreter's engine
contract (§8)**, a runtime mirror of it. Each later part depends only on
the parts before it. The authoring tool is separate:
[shader_workbench_chain_spec.md](shader_workbench_chain_spec.md).

## 1. Historical problem

The ranked pipeline replaced the fixed six-slot model. The original design
is available in Git history; §§2–6 describe the shipped model.

## 2. Families and the structural model

A **family** is the geometric domain a value lives in as it flows from
view vector to pixel color — the pullback analogue of World → Screen →
Pixel:

| Rank | Family   | Domain                                    | Endomorphisms (today's vocabulary)      |
|------|----------|-------------------------------------------|-----------------------------------------|
| 0    | `SPHERE` | unit view directions                      | camera rotation, lenses, surface noise  |
| 1    | `PLANE`  | projected complex coords + provenance     | planar warps                            |
| 2    | `FIELD`  | unit field value + coverage               | transfer, coverage shaping              |
| 3    | `COLOR`  | straight-alpha color                      | (new territory: grading, tone ops)      |

**The structural model is one mechanism: typed carriers.** Each family
has exactly one carrier (§3), and the carrier set is **closed with a
single authority**: one ordered type list,
`CarrierList = TypeList<SphereSample, PlaneSample, FieldSample, Color4>`,
in rank order. Everything else derives from it — `CanonicalCarrier<T>`
(membership: `T` appears in the list), the constexpr `family_of<T>`
rank (`T`'s list index — **total**: a nonmember yields a sentinel
rank rather than a formation error, which §5's staged diagnostics
rely on), and the interpreter's maximum slot size and
alignment (§8's `static_assert` folds over the same list). Nothing is
stored twice, so membership, rank, and the runtime ABI cannot drift
when carriers evolve; and neither the list nor its derived traits is a
customization point a consumer could specialize — a new carrier
arrives only by spec revision (§3's evolution contract) extending the
list. Membership is validated for every leaf's Input and Output (§5).
Closure is load-bearing, not tidiness: the interpreter sizes
its evaluation slots from this set (§8), so an open universe would let
a stage pass static validation while violating the runtime ABI. Every stage names an Input and an Output
carrier, and a chain is structurally legal iff

- adjacent carriers agree (`Stage[i]::Output == Stage[i+1]::Input`),
- every stage satisfies `family_of<Input> <= family_of<Output>`,
- the first Input is `SphereSample` and the last Output is `Color4`.

There is no other structural metadata — though structure is one concern
among several: callable contracts, provider validation, and the semantic
invariants of §3 remain distinct checks layered on top. Monotonicity is
what carrier adjacency alone cannot express — without it, a consumer-authored
down-crossing stage (`FieldSample → PlaneSample`) could re-enter PLANE and
sample again — and together the rules give the once-per-boundary law: each
family boundary is crossed at most once, and the chain can never return to
a lower family. Anything that "looks back" (a warp modulated by a field
value) is a policy input within a stage, not a rank decrease.

A stage whose ranks are equal is an **endomorphism**; one whose Output
rank exceeds its Input rank is a **crossing**. The distinction is
vocabulary, not machinery — no contract member or validation row depends
on it. Rank skips are legal, exactly as a World filter may be followed
directly by a Pixel filter: a 3D-noise source sampled on the sphere is a
`SPHERE→FIELD` crossing that never projects; `Pullback::RayStage<Renderer>` in `core/render/pullback/ray.h` ships the
`SPHERE→COLOR` crossing.

Family *names* survive as the rank's vocabulary — in diagnostics ("pullback
pipeline: a stage may not decrease its family rank — families are ordered
Sphere, Plane, Field, Color", keeping the filter pipeline's prose-assert style)
and in the workbench's band rendering. The mechanism is the rank function.

COLOR carries **straight alpha**: `GeneratedPalette` scales `alpha` by
coverage and opacity but never premultiplies the channels — final
premultiplication stays in `Scan::Shader`, after the chain. COLOR
endomorphisms are therefore provenance-free `Color4 → Color4` operators —
they may read `FrameState`, but anything that needs `value`, `coverage`,
`sphere`, or `path_length` (PATH_LENGTH hue, brightness envelopes) lives
in the Colorize policy, because those inputs are consumed at the
crossing. The family ships empty — its cost is one rank entry; the
crossing it terminates already exists.

## 3. One canonical carrier per family

> **Historical migration terminology:** References below to `ProjectionSample`, `WarpResult`, `SourceInput`, `MaterialInput`,
> `MaterialSample`, and `SurfaceProject` describe the pre-cut-over API.
> The migration is complete; the shipped carriers and combinators are in
> `core/render/pullback/contract.h` and `core/render/pullback/stage.h`.

Free chaining within a family requires every stage in that family to
speak one type. The slot-specific carriers collapse to four, each
following the same grammar — working state, provenance record, trace
accumulator:

```cpp
struct SphereSample {         // rank 0 — must stay all-float: as a 4-float
  Vector dir;                 // homogeneous aggregate it is returned in
  float path_length;          // s0-s3 across out-of-line boundaries
};

struct ProjectionProvenance { // exactly what a projection computes:
  uint8_t region_id;          // regions, fades, weights — nothing more
  uint8_t component_id;
  uint8_t boundary_flags;
  float fade_edge_distance;
  float value_weight;
  uint8_t flags;
  uint8_t traits = 0;
  uint8_t edge_class = 0;
  float domain_coverage = 1.0f; // defaulted: projections routinely omit
};                              // trailing fields (stereographic,
                                // from_kernel) and rely on full coverage

struct ProjectionResult {     // the projection policy protocol
  Complex coords;
  ProjectionProvenance provenance;
};

struct PlaneSample {          // rank 1 (~44 B, parity with today's boundary)
  Complex coords;             // THE planar coordinate; warps advance it
  ProjectionProvenance provenance;
  Vector sphere;              // sample point; combinator-written (§4)
  float path_length;
};

struct FieldSample {          // rank 2 (~24 B)
  float value;                // field value; value in [0,1] is a carrier
                              // invariant, established by the crossing
  float coverage;             // accumulated coverage
  Vector sphere;              // sample's sphere point, read by Colorize
  float path_length;
};

// rank 3: Color4, straight alpha (endomorphisms are Color4 -> Color4)
```

`SphereSample` carries the working `Vector` and accumulated path length;
`SurfaceResult` remains the surface policy protocol. Historical carrier mapping:
`ProjectionSample`+`WarpResult`+`SourceInput` → `PlaneSample`;
`MaterialInput`+`MaterialSample` → `FieldSample`. `ProjectionSample`
is **replaced** by `ProjectionResult` — a protocol with no ignored
state: no `sphere` field for a policy to fail to set (today policies
leave it defaulted and the fused stage overwrites it; the new protocol
cannot express the mistake), and no `surface_path_length` (its job
moves to the carrier accumulator). `Stage::Project` assembles the
carrier: the result's `coords` become the working coordinate, its
`provenance` embeds unchanged, and `sphere` is combinator state written
from the pre-projection point. The carrier holds one planar coordinate,
not two: the only consumer of the embedded copy today is the warp
chain's seed (`Stage::Warp::run`), which is precisely the split the crossing
now performs. `SurfaceResult` and `WarpStepResult` remain the surface and warp
policy protocols, respectively.

The `PLANE→FIELD` crossing **consumes** projection provenance rather
than carrying it: the weight policy eats `value_weight` into the value,
the crossing's coverage policy eats provenance into `coverage` alongside
`domain_coverage` — edge fade included, since `EdgeFade` reads only
`fade_edge_distance` and its width parameter, never the field value —
and the planar topology fields are PLANE-family vocabulary no downstream
policy touches. Exactly one upstream fact outlives the crossing:
`sphere` — the sample point the `Project` combinator recorded in the
carrier (a sphere-domain crossing writes `dir` directly, its true
sample point, not a fabrication), Colorize's hue-noise and palette
input. FIELD carries no sentinels and no fabricated provenance.
The shipped value-dependent coverage policy is `ValueCutout`, whose
contract is that it reads the **current FIELD `value`** and nothing
else — it remains a FIELD stage consuming only carrier state, and a
chain may legally place it before, between, or after transfers (each
placement reads a different value and means something different).
Migration of today's shipped behavior places it after the migrated
`Transfer` stage, because today's cutout reads the post-transfer
value — an ordering fact of the migration tables, not of the stage
contract. The carriers are **curated shading
contexts, not a closed algebraic ideal**: `sphere` and `path_length`
ride in FIELD because the shipped colorizer consumes them, and a future
consumer of a fact the crossing discards (a `region_id`-driven
colorizer, say) widens the canonical carrier — a deliberate, global,
spec-level change, never a per-effect patch. Carrier evolution being
explicit and versioned is the contract; carrier immutability is not.
Widening carries two standing obligations: every canonical carrier
stays trivially destructible, and the interpreter's slot size and
alignment derive from the carrier set by `static_assert` (§8), so a
widening recompiles the slots rather than silently overflowing them.

The scan still hands the pipeline a `Vector`: `evaluate` seeds
`SphereSample{view, 0.0f}` at entry — the ENTRY rule checks the first
stage against `SphereSample`, and the adapter is the pipeline's, not a
stage.

Every stage is fully inlined into the scan today, so carriers are
compile-time fiction after SROA; only `Stage::Placed` (§4) passes them
across a real call boundary. A placed sphere-only run passes 16-B
register-returned `SphereSample`s — narrower than today's fused-surface
spill — and a placed run ending at `Project` returns a ~44-B
`PlaneSample` by memory, parity with today. Placement boundaries are
chosen per effect, and §6's gates prove each one.

Path length becomes a single accumulator in the carrier, replacing
`surface_path_length` plus a separate warp sum joined at the material
stage. The `path_length_required` provider gate and the policies'
constant-zero returns survive unchanged, so for effects that never track
path the accumulator constant-folds to zero through the inlined chain.

Structural legality is what the typed-carrier mechanism guarantees;
semantic invariants are not — they are **obligations**, each established
at a named point, preserved by every shipped combinator, pinned by
tests, and assumed by any custom stage. The type system cannot promise
them: nothing prevents a consumer-authored endomorphism from emitting
NaN. The rule is generic, stated once: **every stage may assume its
Input carrier's invariants and must establish its Output carrier's
invariants**; endomorphisms additionally preserve the monotonic and
immutable fields named below. The middle column records only the
*shipped* stage that first establishes each carrier — a
consumer-authored crossing owes the same output obligation (a custom
`SphereSample → Color4` sky stage establishes the `Color4` row itself).
The normative table:

| Carrier | First established by (shipped) | Endomorphisms must preserve |
|---|---|---|
| `SphereSample` | entry adapter | `dir` unit-length; `path_length` finite, ≥ 0, non-decreasing |
| `PlaneSample`  | `Project`     | `provenance` and `sphere` immutable — endomorphisms advance `coords` and `path_length` only; provenance values in range (`value_weight`, `domain_coverage` ∈ [0, 1]; `fade_edge_distance` finite), a projection-policy obligation; `path_length` finite, ≥ 0, non-decreasing |
| `FieldSample`  | `Sample`      | `value ∈ [0, 1]`; `coverage ∈ [0, 1]`, non-increasing; `sphere` immutable; `path_length` finite, ≥ 0, non-decreasing |
| `Color4`       | `Colorize`    | straight alpha; channels and alpha stay in [0, 1] |

`FieldSample.coverage ∈ [0, 1]` is established without a clamp:
`ProjectionCoverage` policies are obliged to return [0, 1], and
`value_weight`/`domain_coverage` are in range by the projection
obligation above, so `Sample`'s coverage product is in range by
construction.

### Numeric note

Path length re-associates: `((pre + post) + w1) + w2` replaces
`(pre + post) + (w1 + w2)`, differing only when two or more warp stages
contribute nonzero path. Device builds compile with `-ffast-math`, which
already licenses reassociation, so the capture-manifest and golden re-bake
at §6 step 2 is required.

## 4. Stage vocabulary

> **Historical migration terminology:** References below to `ProjectionSample`, `WarpResult`, `SourceInput`, `MaterialInput`,
> `MaterialSample`, and `SurfaceProject` describe the pre-cut-over API.
> The migration is complete; the shipped carriers and combinators are in
> `core/render/pullback/contract.h` and `core/render/pullback/stage.h`.

The policy layer — `Surface::*`, `Lens::*`, `Projection::*`, `Warp::*`,
`Source::*`, `Weight::*`, `Transfer::*`, `ProjectionCoverage::*`,
`ValueCoverage::*`, `Color::*` — keeps its callable logic. The mechanical
signature edits — this list is
normative and intended to be exhaustive; if the cut-over finds another,
the spec is what changes: projection policies retype their return
`ProjectionSample` → `ProjectionResult` (they never set the removed
fields); planar source policies retype `SourceInput` → `PlaneSample` (field
path `input.warped.coords` → `input.coords`), while spherical source policies
consume `SphereSample`; warp and
weight policies retype their
`const ProjectionSample &` parameter to `const ProjectionProvenance &`;
color policies retype `MaterialSample` → `FieldSample` (`sample.sphere`
keeps its spelling); the coverage vocabulary moves to the crossing's
`ProjectionCoverage` signatures and `ValueCutout` drops the parameter it
ignores (§ the `Sample` bullet below); and the projection *helper*
free functions (`projection.h`: `stereographic`, `folded_sinusoidal`,
`equirectangular`, `gnomonic`, …) retype their `ProjectionSample`
return to `ProjectionResult`, which rewrites each flat positional
aggregate return into the nested `{coords, {…}}` provenance form —
initializer rewrites, not just signature retypes, and the largest
single block of mechanical edits in the cut-over. No logic changes
anywhere in this list.

The stage combinators in `stage.h` become one family-typed wrapper per
policy role. Combinators are verbs; the policy namespace each wraps is
its noun:

```
SPHERE  Stage::Rotate<OrientationProvider>    SphereSample -> SphereSample
        Stage::Displace<SurfacePolicy>        (wraps Surface::*)
        Stage::Lens<LensPolicy>               (wraps Lens::*)
crossing Stage::Project<ProjectionPolicy>     SphereSample -> PlaneSample
         Stage::SampleSphere<SourcePolicy>    SphereSample -> FieldSample
         RayStage<Renderer>                  SphereSample -> Color4
PLANE   Stage::Warp<WarpPolicy>               PlaneSample -> PlaneSample
crossing Stage::Sample<SourcePolicy,
                       WeightPolicy = Weight::Projection,
                       CoveragePolicy = ProjectionCoverage::Weight>
                                              PlaneSample -> FieldSample
FIELD   Stage::Transfer<TransferPolicy>       FieldSample -> FieldSample
        Stage::ApplyCoverage<ValueCutout<...>>     FieldSample -> FieldSample  (value-dependent only)
crossing Stage::Colorize<ColorPolicy>         FieldSample -> Color4
COLOR   (none yet; the family exists so they can)
```

Each wrapper's transformation is **normative** — the carrier types alone
cannot distinguish add-from-replace or preserved-from-dropped state, and
capture equivalence hangs on exactly those details:

```
Rotate    out = {rotate(in.dir, Provider::conjugate(frame)), in.path_length}
Displace  r = Policy::apply(in.dir, frame[, prepared])
          out = {r.sphere, in.path_length + r.path_length}   // adds, never replaces
Lens      out = {Policy::apply(in.dir, frame), in.path_length}
Project   local = rotate(in.dir, Policy::frame_conjugate(frame))
          r = Policy::project(local, frame)          // ProjectionResult
          out = {r.coords, r.provenance, /*sphere=*/local, in.path_length}
Warp      step = Policy::apply(in.coords, in.provenance, frame[, prepared])
          out = {step.coords, in.provenance, in.sphere,      // provenance, sphere unchanged
                 in.path_length + step.path_length}
Sample    normative pseudocode in the bullet below
SampleSphere out = {clamp_unit((Policy::sample(in) + 1) / 2), 1,
                    in.dir, in.path_length}
Transfer  out = in;  out.value = Policy::apply(in.value, frame)
ApplyCoverage out = in;  out.coverage = in.coverage * Policy::apply(in.value, frame)
Colorize  Policy::apply(in, frame) -> Color4
```

The projection protocol has no `sphere` field (§3): the sample point is
combinator state, written by `Project` from `local` — a policy cannot
even accidentally control it, which strengthens today's convention
where policies leave the field defaulted and the fused stage overwrites
it after projecting. These transformations are implemented **once, as shared carrier
kernels** — free functions over the canonical carriers — called by both
the template combinators and the interpreter's erased adapters (§8), so
the semantics cannot fork between the two execution paths.

- `Stage::Project` is the tail of today's `SurfaceProject::run` (rotate by
  `frame_conjugate`, project, fill provenance). The pre/post-lens surface
  pairing, and the Emscripten flattened path, dissolve: a chain lists
  `Displace` before or after `Lens`, or not at all. (The flattened path
  is computation-identical to `[Displace, Project]` with no lens stage;
  retiring it changes nothing numerically.)
- **`Stage::Sample` establishes the FieldSample invariant.** The crossing
  samples the source, applies the weight policy to the raw signed field,
  ramps, and seeds the carrier — normatively:

  ```cpp
  const auto source_span = Instrumentation::mark();
  const float raw = SourcePolicy::sample(input, frame /*, prepared */);
  Instrumentation::template span<ProfileEvent::SOURCE>(source_span);
  const auto material_span = Instrumentation::mark();
  const float weighted = WeightPolicy::apply(raw, input.provenance, frame);
  const FieldSample out{
      .value = Detail::clamp_unit((weighted + 1.0f) * 0.5f),
      .coverage = CoveragePolicy::apply(input.provenance, frame)
                * input.provenance.domain_coverage,
      .sphere = input.sphere,
      .path_length = input.path_length,
  };
  Instrumentation::template span<ProfileEvent::MATERIAL>(material_span);
  return out;
  ```

  `Stage::SampleSphere` establishes the same invariant without projection.
  Its source consumes `SphereSample`; the signed field is ramped directly,
  coverage is `1`, and the input direction and path length are preserved.

  This is the distillation today's Material combinator performs
  (including the `domain_coverage` seed), relocated to the boundary
  where plane provenance becomes field state. The crossing's coverage
  policy is one slot from a dedicated vocabulary —
  `ProjectionCoverage::{None, Weight, WeightSquared, EdgeFade<Provider>}`,
  signature `apply(const ProjectionProvenance &, const FrameState &) →
  float` — preserving today's **mutual exclusivity**: coverage modes
  never stack, so a migrated edge-fade chain does not silently acquire a
  projection-weight factor. The two consumers name their modes
  differently, so the normative migration is two tables. ComposedEffect
  `ProjectionCoverageMode`: `NONE` → `None`, `WEIGHT` → `Weight`,
  `WEIGHT_SQUARED` → `WeightSquared`, `EDGE_FADE` → `EdgeFade`. ShaderWorkbench
  `CoveragePolicy`: `OPAQUE` → `None`, `PROJECTION_WEIGHT` → `Weight`,
  `PROJECTION_WEIGHT_SQUARED` → `WeightSquared`, `EDGE_FADE` →
  `EdgeFade`, `VALUE_CUTOUT` → `None` at the crossing plus a
  `Stage::ApplyCoverage<ValueCutout>` FIELD stage. Its policy's
  signature is `apply(float value, const FrameState &) → float`
  and it multiplies into the accumulated `coverage`. Because the boundary is crossed
  exactly once (§2), a double ramp or an unramped value reaching Colorize
  is unrepresentable in the shipped vocabulary: no stage emits a signed
  value into the chain and there is no second ramp to apply. (`value ∈
  [0, 1]` is a semantic invariant of the built-in stages, not a numeric
  refinement type — a consumer-authored FIELD endomorphism must preserve
  it.) Weighting is a crossing policy rather than a stage because the raw
  signed field exists only inside the crossing; `Weight::Projection` is
  the default and `Weight::None` the alternative, matching the dynamic
  backend's existing signal-weight slot. A sphere-domain source crossing
  does the same with neutral weight. `RayStage` is a complete one-stage chain;
  the field path uses `SampleSphere -> Colorize`, and the plane path requires
  `Project -> Sample -> Colorize`. `Rotate` is optional in either path.
- `PlanarWarp`'s variadic policy list becomes N consecutive `Warp` stages.
- `Identity` policies stop appearing in pipelines — absence of a stage
  *is* the identity.
- The five shipped crossings are not a closed set: any combinator whose
  Output rank exceeds its Input rank is admitted by the same rules, which
  also admits new crossings without a schema change. `RayStage` supplies the
  shipped `SPHERE→COLOR` crossing.
- **`Stage::Placed<EMISSION, Stages...>`** is the single grouping
  construct: a contiguous run of stages that is itself a stage (Input =
  first's Input, Output = last's Output, internal adjacency and
  monotonicity validated, Prepared a nested tuple), emitted as one call
  unit under the given `CodeEmission`. The single-stage form is today's
  `Placed`; `INLINE_ONLY` is a transparent forwarder, which is what lets a
  derivation layer make placement a computed value rather than a
  structural branch (§6). Nesting a sub-chain is sound here where the
  filter pipeline forbids it (`is_pipeline`) because a pullback run is a
  pure value function, while a filter pipeline owns the canvas sink.
  ComposedEffect's prescribed displacement placement preserves today's
  single flash-call boundary; without multi-stage placement, lens and
  projection code would return to the inlined hot scan loop of every
  displacing effect. Placement is always written by the author — never
  inferred — because placements are deliberate, ITCM-ledger-scored
  decisions. The pipeline therefore has **two views**: the flattened
  *semantic leaf list* — where carrier adjacency, monotonicity, binding,
  approximation metadata, and consumer predicates live, and which has no
  notion of emission — and the *structural placement tree*, where
  emission lives exclusively: `Placed` nodes carry `CodeEmission`, and a
  bare stage in the pipeline list is an implicit inline placement node.
  `Placed` contributes nothing to the semantic view, so adding or
  removing a placement wrapper can never change whether a chain
  validates or what a predicate matches; placement affects exactly
  prepared-state layout and code emission. A `Placed` whose post-filter
  stage list is empty collapses to `void` itself, so a fully-conditional
  group vanishes exactly like its members. The two views come from
  **one normalization**, run once at pipeline assembly:
  `descriptors + Binding → bound execution tree → semantic leaf list` —
  the bind step of the contract paragraph below is part of it, and the
  tree is the authority: normalization builds the bound execution tree
  first and projects the leaf list from it. Ownership is fixed per view,
  **with separate indexings** — the views describe different tree shapes
  and share no indices. The leaf list owns every validation fold, the
  consumer predicates, and the public introspection surface:
  `STAGE_COUNT` and `stage_at<I>` count and index *leaves*. The
  execution tree owns `prepare` and `run` through its own internal
  indexing (`node_at<I>`, node count): `evaluate` is an index recursion
  over execution *nodes*, `PreparedTuple` is indexed by node and nests
  per placement node so an out-of-line call receives one contiguous
  prepared sub-tuple, and a `Placed` node runs its children by the same
  recursion over its own nested pack and tuple. Today's coupling of
  `stage_at<Index>` / `STAGE_COUNT` / `std::get<Index>(prepared)` in one
  recursion therefore splits: the recursion keeps its shape but runs on
  node indices; the leaf indices are a read-only public surface that
  never drives execution. Per-leaf contract facts (`RUN_RETURNS`,
  `PREPARES`) are validated on bound leaves; the tree only composes
  them. No API may mix the views: anything semantic reads leaves,
  anything executional reads nodes.

The shipped `Stage::Contract` declares no `KIND` or `TERMINAL`;
`Input`/`Output`/`Prepared`/approximation metadata remain per stage,
while `EMISSION` leaves the stage contract for the placement tree
(above). The binding machinery is rebuilt, not retained: the policy
surface is mixed — provider-parameterized policies export `using
Binding` (warps, sources, provider-bound lenses), while parameterless
policies (`Transfer::Ridge`, `Lens::Glitch`, `Weight::Projection`) and
provider-templated ones without the alias (`GeneratedPalette`) do not —
so combinators are **unbound descriptors** declaring carriers and
`Policies`, and pipeline assembly normalizes recursively, tree-first:
(1) remove `void` entries while preserving placement nodes (an empty
`Placed` collapses to `void`); (2) bind every leaf *in place in that
tree* via an internal `Descriptor::Bind<Binding>`, producing bound
stages with concrete `FrameState`, `Instrumentation`, and `Prepared`;
(3) derive the semantic leaf list as a projection of the bound tree.
The execution tree is the one authoritative representation; the leaf
list is a view of it, so placement grouping is never flattened away and
then reconstructed. Binding agreement is checked by a **non-asserting
predicate** — `descriptor_bindable<Descriptor, Binding>()` — evaluated *before* the
bind step instantiates anything, so the pipeline's named `BINDINGS`
assertion (§5) is what reports a foreign binding, with today's
message, rather than a template-formation abort inside `Bind`
preempting the named surface; `Bind`'s internal `PROVIDER_VALID`
asserts remain as a backstop, unreachable through pipeline assembly.
Execution uses only bound types. This is
what keeps binding services available inside combinators: `mark`/`span`
resolve through the bound stage's binding exactly as
`BindingT::Instrumentation` does today. `StageDescriptor`,
`descriptor_bindable()` and `Detail::DescriptorPrepares` validate descriptors
and bound leaves; the pipeline reports `CONTRACTS`, `BINDINGS` and `PREPARES`.

**The descriptor contract is public — and it is an execution contract
only.** Consumer-authored combinators — including the rank-skipping
crossings §2 admits — implement the same contract the shipped ones do;
extension edits no pipeline machinery. A descriptor provides:
`Input`/`Output` (drawn from the **closed set** of canonical carriers,
§2 — a type with a homemade `family_of` specialization is not a
carrier); `Policies` (a tuple, possibly empty); and
`template <typename Binding> Bind`, yielding the bound stage with
`Prepared` (trivially destructible),
`static Prepared prepare(const FrameState &)`,
`static Output run(const Input &, const FrameState &, const Prepared &)`
(always-inline), and the approximation metadata (`APPROXIMATE`,
`ORACLE`, `METRICS`, `NON_FLOATING_FIELDS_EXACT` — defaulted by
`ApproximationDefaults`, aggregated by `CombinedApproximation` over
`Policies`);
provider-identity asserts (`PROVIDER_VALID`) fire inside `Bind`.
Promotion metadata is deliberately **not** part of this contract: an
operator's provider requirements live in the promotion catalog (§7),
and an operator without that metadata simply is not promotable —
hand-written static stages carry no tool-facing declarations. The forwarding mechanism is
normative: the author writes `run` (and optionally `prepare`) as
**binding-templated statics on the descriptor**, and the helper base's
`Bind<Binding>` forwards into them — that is how an unbound descriptor
can define execution that needs the binding-dependent `FrameState` and
`Instrumentation`. A complete consumer-authored combinator, doubling as
the rank-skip example:

```cpp
struct SkyGradient
    : Pullback::Stage::Contract<SkyGradient, SphereSample, Color4> {
  using Policies = std::tuple<>;   // ApproximationDefaults apply

  template <typename Binding>
  static Color4 run(const SphereSample &in,
                    const typename Binding::FrameState &frame,
                    const NoPrepared &) {
    const auto start = Binding::Instrumentation::mark();
    const Color4 out = shade_sky(in.dir, frame);
    Binding::Instrumentation::template span<ProfileEvent::COLOR>(start);
    return out;
  }
};
```

Nothing else is required: `Contract` derives carriers and defaults,
`Bind` forwards, and the pipeline validates it like any shipped stage.
The shipped combinators are implemented against this exact contract;
that they need nothing more is the proof it suffices.

Instrumentation lives inside combinator and policy bodies (the
MirrorTile pattern), pinned to `NoInstrumentation` by composed effects —
but the responsibility split moves span boundaries, so event ownership
is normative, chosen to keep today's report buckets meaningful:

| Stage | Event |
|---|---|
| `Rotate` | none (uninstrumented today) |
| `Displace` | `SURFACE_NOISE` |
| `Lens` | `LENS` |
| `Project` | `PROJECTION` |
| `Warp` | `PLANAR_WARP` (MirrorTile keeps its policy-internal `MIRROR_TILE`) |
| `SampleSphere` | `SOURCE` |
| `Sample` | `SOURCE` around the source-policy call; `MATERIAL` around weight + ramp + projected coverage |
| `Transfer`, `ApplyCoverage` | `MATERIAL` |
| `Colorize` | `COLOR` |

`MATERIAL` thus becomes the sum of up to three spans (Sample's tail,
Transfer, ApplyCoverage) covering exactly the work today's single Material
span covers, so per-bucket cycle reports keep their meaning across the
migration; splitting the fused surface stage likewise recovers today's
`LENS`/`SURFACE_NOISE`/`PROJECTION` granularity without touching the
pipeline recursion.

For consumer validation hooks, every combinator exports its policies
uniformly (`using Policies = std::tuple<...>` — a single-element tuple
for most, `{SourcePolicy, WeightPolicy, CoveragePolicy}` for `Sample`,
so every policy participates in provider validation, approximation
aggregation, and predicates — plus descriptive aliases for each), and
the pipeline exposes trait folds over its flattened leaf list —
`any_stage<Predicate>`, `stage_matching<Predicate>` — with predicates
matching over `Policies`; `Placed` is invisible to them by the
transparency rule above. ShaderWorkbench's `ExtraValidation` uses its own hand-written fold for the
edge-fade/projection compatibility check; it does not use these trait folds.

## 5. Validation

The staged structure (earlier levels gating later ones) and the
named-boolean surface stay, with one more level than today. The order
is explicit: descriptor contract shape, then `CANONICAL`, then the
rank rows (`MONOTONE`, `CARRIERS`, `ENTRY`, `EXIT`), then the bound
callable checks. Ordering is load-bearing at the second step: the rank
rows evaluate only after `CANONICAL` passes, and `family_of<T>` is
total (§2's sentinel rank), so a malformed descriptor cannot detonate
template formation inside a rank fold before the named `CANONICAL`
assertion reports — the diagnostic leads with the actual cause even
under eager instantiation. The rows become §2's rules:

| Today          | Relaxed                                                        |
|----------------|----------------------------------------------------------------|
| `ARITY` (== 6) | `NONEMPTY` (gates the other rows, as `ARITY` does today)       |
| `ORDER` (fixed kind table) | `MONOTONE`: per stage, `family_of<Input> <= family_of<Output>` |
| `CARRIERS` (slot-typed)    | `CARRIERS` (adjacent Output/Input equality), plus `ENTRY` (first Input is `SphereSample`) and `EXIT` (last Output is `Color4`) |
| `TERMINALS`    | retired — subsumed by `EXIT` + monotonicity; a mid-chain `Color4` producer followed by COLOR endomorphisms is now a feature |
| `EMPTY_DESCRIPTORS`, `RUN_RETURNS`, `EXTRA_VALIDATION` | unchanged, evaluated over the flattened leaf stages |
| `CONTRACTS`, `BINDINGS`, `PREPARES` | reformulated over the descriptor contract (§4): every leaf is bound by the pipeline's bind step; a descriptor whose policy names a foreign binding fails there |
| — | `DESCRIPTOR_IDENTITY`: bound stages preserve their descriptor identity |
| `APPROXIMATIONS` | unchanged per leaf stage (`CombinedApproximation`'s ≤1 oracle within one stage's policy list) |
| —              | `CANONICAL`: every leaf's Input and Output satisfies `CanonicalCarrier<T>` (§2) — membership in the closed `CarrierList`, the same list `family_of` ranks over |

All folds run over the flattened leaf list, so placement wrappers cannot
perturb them.

Cross-stage approximation correctness is **not provable from stage
metadata**: errors compose across domains (projected-coordinate error
perturbs the palette's input, so an approximate projection and an
approximate color stage still interact at the framebuffer), and
conversely two same-domain approximations can be jointly fine. The
authority is pipeline-level acceptance — capture-manifest thresholds
measured over the pipeline's actual output, after every stage including
any COLOR endomorphisms. Its representation, the roster-derived
registries, and the CI gates covering both static and interpreted chains
are specified in §7.4; per-leaf metadata remains natural-domain
bookkeeping here.

### What deliberately stays rigid

- **Pipeline mechanics.** `PreparedTuple`, `prepare`/`prepare_into`/
  `shade`/`shade_prepared`, and the index-recursive always-inline
  `run_stage` are already length-generic (the recursion mentions no
  count; repeated `NoPrepared` tuple elements are distinct subobjects
  and cost roughly a byte each plus padding — a few bytes per chain,
  not zero; `PREPARED_BLOB_BYTES` stays as the backstop, and
  uniquely-indexed empty prepared types are the named remedy if the
  bytes ever matter). They keep their exact *shape* — always-inline index
  recursion with `std::get<Index>` over a pack-indexed tuple — while the
  indices become execution-node indices per §4's two-view split; the
  stage-dispatch inline cliff is real, and no dispatch loop or IIFE may
  replace the recursion.
- **Snapshot layout.** Each composed effect's `ParamsFor<Spec>` is a concrete,
  fixed layout derived from its declared provider instances. Its C++ snapshot
  version changes whenever that layout changes; the effect and preset identities
  in promoted documents remain stable.

### Scope boundary

The ranked pipeline has no fixed capacity for parameterized stages. A composed
Spec declares `template <typename B> using Pipeline = Pullback::Pipeline<B, ...>`;
`ComposedEffect<W, H, Derived, Spec>` derives present-only parameters, clocks,
and noise resources from that declaration. Firmware arena and prepared-state
budgets still bound each concrete instantiation. Promotion's catalog requirements,
parameter binding and topology acceptance remain specified in ?7.

## 6. Pipeline assembly and consumer migration

The authoring surface is one spelling. The relaxed pipeline keeps the
`Pipeline` name — there is no separate chain type — and a concrete
pipeline is its binding, once, followed by its stages:

```cpp
using GnomonicGridMirror = Pullback::Pipeline<ShaderBinding,
    Stage::Rotate<CameraProvider>,
    Stage::Project<Projection::Gnomonic<...>>,
    Stage::Warp<Warp::MirrorTile<...>>,
    Stage::Sample<Source::Grid<...>>,   // weight + projected coverage by default
    Stage::Colorize<Color::GeneratedPalette<...>>>;  // no Transfer/Coverage:
                                                     // absence is identity
```

Two rules keep it that flat:

- **Binding appears exactly once.** `Pipeline`'s first parameter binds
  the whole list; the stages are descriptors (§4) — provider-bound ones
  must agree with it, parameterless ones are bound by instantiation.
  Nothing is threaded per stage.
- **No conditional vocabulary.** `Pipeline` and `Placed` ignore `void`
  entries — the entire hook a derivation layer needs — and validation
  runs on the post-filter leaf list so diagnostics name the stages an
  author actually sees. Concrete pipelines never contain a conditional.

The only conditional assembly in the codebase is ComposedEffect's
derivation, which already computes per-family policies with
`conditional_t` and now yields `void` where a family is absent, with
placement selected by the author's `SurfacePlacement` template argument.
A composed Spec declares only its actual sphere stages, including displacement
placement and code emission:

```cpp
template <typename B>
using Pipeline = Pullback::Pipeline<B,
    Stage::Rotate<OuterCameraProvider<B>>,
    Stage::Placed<CodeEmission::OUT_OF_LINE_FLASH,
        Stage::Lens<Lens::HexagonalPrismKaleidoscope>,
        Stage::Displace<Surface::DirectNoise<
            SurfaceProvider<B, DirectSurfaceParams, true>,
            math::NoiseBasis::SIMPLEX>>,
        Stage::Project<Projection::Stereographic<ProjectionProvider<B>>>>,
    Stage::Sample<Source::Grid<SourceProvider<B, GridSourceParams>>>,
    Stage::Colorize<Color::GeneratedPalette<ColorProvider<B, HUE, BRIGHTNESS>>>>;
```

`Stage::Placed` remains the explicit boundary for moving a contiguous run out of
ITCM. Resource discovery descends through that run's leaves; it does not flatten
away its emission boundary.

### Landing plan

**No transitional compatibility surface** — no legacy carrier typedefs,
no validation-boolean aliases, no six-slot sugar, no `KIND` shim. The
landings form a progression of complete states; anything a landing breaks
migrates inside that landing. The contract cut-over is one wide landing,
wide because the old surface has real readers: the tests read today's
validation booleans and positional slot typedefs, ShaderWorkbench's dynamic
backend aggregate-initializes the old carriers field-by-field, and its
test-pinned projection-join facility lerps `surface_path_length` between
branches. Enumerating those readers is what makes the landing plannable.

1. **Prelude (independent, final).** Delete the never-instantiated
   `GnomonicDodecahedralGridVectorMirrorPipeline::shade` — its
   two-argument `run_stage<0>(view, frame)` call predates the
   prepared-state parameter and surfaces the moment the rewrite
   instantiates it. Pure dead-code removal, standalone.
2. **The contract cut-over, one landing.** Families, the four carriers,
   the family-typed combinators, variadic `Placed`, and the Binding-once
   void-filtering `Pipeline` with its new validation replace the old
   contract outright, and every reader migrates in the same landing:
   - `composed_effect.h` → the relaxed `Pipeline` with the prescribed
     placement;
   - ShaderWorkbench template pipelines → chains (dropping `Identity` slots);
     `ExtraValidation` retains its hand-written compatibility fold;
   - ShaderWorkbench's dynamic backend → the canonical carriers. The
     projection-join facility (`join_projected` /
     `projection_join_compatible`) is gone: the lens-blend transition it
     was reserved for no longer exists, and no render path ever called
     it. When a transition feature does land, the chain-level
     representation is fixed now: compatible projections may join inside
     a curated compound `SPHERE→PLANE` operator that lerps `coords`,
     selects the discrete topology and edge fields (`region_id`,
     `component_id`, `boundary_flags`, `fade_edge_distance`, `flags`,
     `traits`, `edge_class`) from the nearer branch (`mix < 0.5` →
     direct), **recomputes** `value_weight` from the blended coordinates
     via pole attenuation, lerps `domain_coverage` and path length, and
     normalized-lerps `sphere` (`nlerp_unit`). Strict projections
     (Bonne, Peirce, Airocean) refuse the plane join and require
     complete-output blending, which a linear chain cannot express and
     which is therefore an **evaluator-level two-pass** — shade both
     branches, blend the resulting `Color4`s — outside the chain program
     entirely. Blending a single plane sample is not equivalent
     (everything downstream is nonlinear), so no `PlaneSample`
     interpolation may be substituted. `PREPARED_BLOB_BYTES` re-verified
     against the larger prepared tuples;
   - `test_pullback.h` / `test_shader_chain.h` validation and
     positional-typedef reads;
   - the deletions land here too: the six-slot combinators, the old
     carriers including `ProjectionSample` itself (replaced by
     `ProjectionResult`), the Emscripten flattened path,
     now-unreferenced `Identity` policies, and `Lens::Sequence` —
     consecutive `Lens` stages are the composition mechanism, so it
     duplicates the pipeline's own sequencing and adds no expressive
     power.
   Gates: capture manifests with the path-length re-bake (§3), clean
   builds of both phantasm sides, teensy size trail, a map-file grep for
   standalone stage symbols (the inline-cliff regression tell), and an
   on-device A/B profile **on a displacing effect specifically** — that
   is where placement equivalence is at risk. The codegen goal:
   identical for chains that replicate today's placement *and* are
   untouched by the path-length reassociation (§3) — single-warp or
   untracked-path effects; for chains the reassociation touches, every
   difference is attributed and accepted through the capture and device
   evidence, not assumed away.
3. **Companion layers, sequenced after, each landing final** — the
   interpreter (§8), promotion/verification tooling (§7), and the
   [workbench](shader_workbench_chain_spec.md); they need the landed C++
   vocabulary but add no C++ compatibility surface.

## 7. Promotion and verification

**Status: PARTIAL; each sub-section carries its own banner.** The operator
authority and the field tables it reads ship; the allocator, the binding
table, the promotion pin, and the roster-derived acceptance registry with its
`any_approximate` and `APPROXIMATION_DOMAINS_DISJOINT` folds are design.

A layer over the static model, with its own invariants and failure
modes: allocation can fail where chain validation succeeds, and its
gates run in CI rather than at compile time.

### 7.1 One operator authority

**Status: PARTIAL.** `OperatorDescriptor` and `make_operator_descriptor()`
(`core/render/pullback/operator_model.h`), the `OPERATOR_TABLE` built from
those records (`core/render/pullback/operator_table.h`), and the generated
catalog (`core/render/pullback/catalog_export.h`) golden-pinned by
`tests/test_shader_chain.h`. The shipped factory derives one schema per model
with the topology enum8s as ordinary fields rather than instantiating
per-variant recipes; the provider-requirement function of topology and the
binding-table conformance test wait on the §7.2 allocator.

A promotable operator is visible to three systems — promotion (its
provider requirements), the interpreter (its runtime ABI entry), and
the authoring tool (its catalog entry). These are **views of one C++
record**, an `OperatorDescriptor`: a **type-level stage recipe** —
`template <typename Binding, typename Assignment, typename Topology>
Stage` — rather than one concrete stage type, because both remaining
degrees of freedom are chosen *after* the catalog is written:
allocation selects providers (composed policies encode their instance key in the
type — `WarpProvider<B, Outer>` versus `Inner` — and `Assignment` is
what picks it), and promotion pins the document's bank-invariant
topology values (§7.3) as `Topology`, the provider-free static
configuration that becomes template arguments (weight mode, coverage
mode, noise basis). Only the **carrier pair and the provider-free
logical policy identity** are invariant across instantiations — the
concrete C++ policy types are not, since provider selection is part of
the type — so catalog views read that metadata from the record itself,
never from an instantiation. Alongside the recipe: the parameter
schema — **topology-invariant by construction**: the union of every
variant's fields plus the topology enum8s themselves. Registration
happens at chain compile, before any value (topology values included)
arrives through the value channel (§8), so the schema cannot depend
on them; a field the current topology value deactivates (edge-fade
width under `Weight` coverage) keeps **full schema citizenship** —
registered, stored, defaulted, interpolated, validated,
preset-required, digest-bearing — and is merely unread by the active
variant's `prepare`/`run`, exactly the status mode-gated fields
already have in shipped composed effects. Also in the record: the
optional provider-requirement
function of topology (§7.2; absent = not promotable), the runtime
adapter callbacks with their
block sizes, and worst-case approximation metadata (§7.4). Authority
is **derivation, not colocation** — fields that merely sit in the same
record can still disagree with each other — and derivation needs a
source that actually carries the runtime facts, which the stage
recipe alone does not: §4's descriptor contract is execution-only
(carriers, policies, `Prepared`, `prepare`, `run`) and knows nothing
of parameter blocks, instance state, or lifecycles. The authoring
unit is therefore an explicit **operator model** owning four members:
the shared semantic kernel (§4's free functions over the carriers);
the static stage recipe (the §4 descriptor, for template pipelines —
the model wraps it, the contract itself stays execution-only); the
runtime types — the `FIELDS`-bearing parameter family and the
instance-state type with its lifecycle (`init`/`migrate`/`destroy`/
`advance`), hand-authored exactly where an operator has state and
defaulted to zero-size state with no-op lifecycle where it does not;
and the erased adapters over the kernel. The `OperatorDescriptor` is
produced by a **factory over the model**, deriving every derivable
field from the member that owns it: carrier ids from the recipe's
`Input`/`Output`, block sizes and alignments from `sizeof`/`alignof`
of the model's parameter, prepared, and instance-state types
(worst-cased across topology variants by the same fold §7.2 names),
the runtime callbacks as the model's adapters and lifecycle, the
parameter schema from the family `FIELDS` tables and their stable ids
(§7.2), and approximation metadata as the fold over the variants'
`CombinedApproximation`s. Hand-authored facts are only what no type
carries — names, documentation, the requirement declarations — plus
the state logic itself, which the model localizes rather than
derives: "one record" eliminates duplicated *facts*, and is honest
that behavior is written, once, in the model. The factory
instantiates **every topology variant** at
compile time — conformance is instantiation, so a variant that
violates a policy contract fails the build rather than the first
preset that selects it — and a test asserts the §7.2 binding table
and the interpreter's registered parameter schema describe the same
field set. The
interpreter's operator table is built from these records; the tool
catalog is *generated* from them, with the golden pin guarding only
the C++-to-tool generation step — which, once everything inside the
record is derived, is all that is left to guard; the promotion
emitter reads them directly. Integrating a new operator is authoring
one model — which is what makes §9's "one catalog entry" claim true
rather than aspirational.

### 7.2 Resource allocation

**Status: PARTIAL.** Named per-instance composed storage, stable machine
`Field::id`, field validation/interpolation and canonical typed registration ship.
Catalog `requirements(Topology)`, promotion's parameter-binding table and the
promotion emitter remain design.

A composed Spec declares its ranked pipeline explicitly. Parameter-consuming
providers carry a compile-time string `ResourceKey`, their parameter family and
resource kind. `WarpProvider<B, "ripple", WaveShearParams>` and
`LensProvider<B, "lens-b">`, for example, address separate instances. Repeating a
family with different keys creates independent parameters, clocks and resources;
three clocked warps and two Mobius lenses are ordinary valid compositions.

`ParamsFor<Spec>` discovers provider requirements through each stage's policies,
including compound policies and the leaves of `Stage::Placed`. Parameterless
policies contribute no storage. It constructs one typed block per key. Reusing
one key is explicit sharing and requires the same family and resource kind;
conflicting declarations fail compilation. `get<"key">()` addresses one block;
`visit` drives validation, interpolation, registration, clock advancement and
frame preparation. Each noise-consuming instance owns its own persistent noise
field. Each affine warp owns its own rotation accumulator.

The standard source, projection and color keys carry the shared composed
lifecycle's source/camera/palette roles. All authored effects address parameters
through `get<"key">()`; `Family<"key">` and `HAS<"key">` describe their blocks.
Absent instances have no parameter members or placeholder families. Standard
instance names retain their storage and
registration order; additional keys follow declaration order within each resource
kind. New keys and repeated families qualify display controls as
`<key>.<display-name>` to keep targets distinct. Ordinary controls use the
family's canonical `ParamSpec` through shared typed registration. The exact
registered parameter count and qualified-name bytes enter the persistent arena
budget; no fixed slider capacity or warp/lens slot capacity remains.

Promotion remains an explicit resource-allocation problem with multiplicity.
Each promotable operator's catalog exposes `requirements(Topology)`: the
parameter families and state its recipe's policies actually read at the pinned
topology. Compound policies combine their children's requirements. An operator
without this metadata is not promotable. The interpreter budgets worst-case
state and prepared footprints across topology variants; a promoted concrete
pipeline budgets the variants and instances it actually declares.

The emitter assigns a deterministic key to each document instance in chain order
and writes providers already specialized with that key and family. Distinct
instances receive distinct keys by default; identical families are not collapsed.
The assignment owns a parameter-binding table mapping each `(instance_id,
field_id)` to the concrete typed parameter block and registered control. It is
injective except for an explicitly shared requirement. Shared requirements must
agree in defaults, every preset and transition scheduling; disagreement refuses
promotion. Requirements with disjoint fields may choose a composed family, but
independent instances never require such a merge merely because they have the
same resource kind. `Bind` and the execution stage contract receive no assignment.

Conformance compares the emitted binding field set with the interpreter's
registered schema. Topology, approximation and budget checks remain promotion
conditions; the number of clocked warps or Mobius lenses is not a rejection rule.

### 7.3 Topology parameters

**Status: PARTIAL.** The interpreter half ships: topology parameters are
enum8 runtime switches (§8). The promotion pin is design.

**Topology parameters** — the catalog-flagged enum8 class (weight mode,
coverage mode, noise basis) — are ordinary runtime switches in the
interpreter, free to vary between presets; but promotion pins them as
template arguments in one typedef, so a promotable document must hold
them invariant across its preset bank. A bank that varies them stays
interpreter-only or splits into one document per topology.

### 7.4 Approximation acceptance

**Status: PARTIAL.** The manifest half ships:
`tools/generate_pullback_manifest_header.py` validates the manifests under
`tests/data/pullback/` and generates the native header,
`tests/pullback_manifest_check.cpp` (`unit_pullback_manifest`) checks their
identity properties. No gate currently pairs manifest programs with a live
program table.
The registry derived from the composed-effect roster, the `any_approximate`
fold, the interpreted-chain approximation aggregate, and
`APPROXIMATION_DOMAINS_DISJOINT` are design.

The acceptance representation is a capture-manifest entry — capture key
plus final framebuffer thresholds — measured by the oracle harness over
the pipeline's actual output, after every stage including any COLOR
endomorphisms, so the bound is a property of the program, not of the
Colorize leaf. A leaf's `FRAMEBUFFER` metric declares its intended
contribution to that bound; `has_final_framebuffer_metric` checks only
that the declaration exists, and the manifest's thresholds are the
binding ones. Enforcement is a named CI gate, not a compile-time assert,
and its scope is a registry **derived from the authoritative rosters**,
not maintained beside them — a hand-kept table could only prove that
registered pipelines have manifests, never that every shipping pipeline
is registered. The effect roster that already gates device builds
enumerates every composed effect, each naming its `RenderPipeline`, and
ShaderWorkbench's program table enumerates its studies and dynamic programs;
the capture registry is generated from those two declarations, so a
shipping pipeline absent from the registry is unrepresentable rather
than checklist-caught. The completeness check walks the derived
registry, reads each pipeline's `any_approximate` fold, and fails CI
for any approximate entry absent from the manifest. Interpreted chains
have no compile-time fold, so they take a **second path**: program
compilation aggregates approximation metadata from the operator-table
entries it resolves. That aggregate stays conservative under later
preset changes because operator-table approximation metadata is
**topology-invariant by construction** — each entry declares the worst
case across its topology-parameter variants, and topology values apply
after compilation (§8), so no preset can make a compiled program more
approximate than its aggregate claims. The enumerable shippable set is
the shipped document catalog — each cataloged document is compiled in CI and its
runtime aggregate checked against the manifest. Ad-hoc workbench chains
are previews, not shipping artifacts, and are explicitly outside the
gate, as is a pipeline outside both rosters — a private, unshipped
experiment. Those boundaries are stated, not implied. Per-leaf metadata
remains natural-domain bookkeeping; the informational
`APPROXIMATION_DOMAINS_DISJOINT` fold (error domains = `METRICS` domains
excluding `FRAMEBUFFER`) stays available for consumers to assert via
`EXTRA_VALIDATION` as a lint against unreviewed oracle stacking — the
shipping `PEIRCE_FAST_SQUARE` + `GeneratedPalette` pairing satisfies
it — but it is a diagnostic, not a correctness gate.

## 8. Preview interpreter: engine contract

Template-instantiating arbitrary chains at runtime is impossible, so the
workbench's dynamic preview becomes a **stage-program interpreter** — the
generalization of the per-sample switch dispatch ShaderWorkbench's dynamic
backend already does, walking an array instead of fixed slots. Only its
engine contract lives here; routing and editing are the tool spec's
concern.

- **The engine trusts nothing across the boundary.** The wire payload
  of `setShaderChain` is an ordered list of `{instance, operator}`
  and nothing else — no offsets, no family tags, nothing layout-shaped,
  and **no parameter values**: values flow through the existing
  bulk document channel (`setShaderChainParameters`, keyed `instance.field`,
  validated against the operator's schema) *after* compilation, per the
  apply order below, so the transaction boundary is the chain compile
  alone and values always apply to a committed program. `instance`
  is the document's chain-entry label: operator ids alone cannot
  distinguish `warp1` from `warp2`, and the engine needs the instance
  identity to register per-instance parameter definitions before values
  apply. Structural variants (weight mode, coverage mode, noise basis)
  are enum8 parameters arriving through that same value channel, so a
  chain entry carries no field the engine drops. **Three identities
  are distinct, related by projection**: the *descriptor digest* covers
  the descriptor — chain, parameter schemas and defaults with display
  units excluded, serialization fields, path policies — and not the
  preset bank, which digests separately as `preset_bank_digest`, so a
  preset edit leaves the descriptor digest, and the parity toggle it
  arms, untouched; the *program-shape identity* is the ordered
  `{instance, operator}` list, exactly what `setShaderChain`
  consumes — many documents share one program shape and differ only in
  the values they then apply; the *instance-state identity* is a
  single entry's `(instance_id, operator_id)` pair, the migration key
  below. The digest refines the shape; neither collapses into the
  other, and the tool's bypass toggle is an **ephemeral program-shape
  override** — it compiles a shape omitting one entry while the
  document and its digest are untouched, which is exactly why the
  companion spec keeps bypass session-only and never serialized.
  `setShaderChain` *compiles* the shape: it
  resolves each operator against the engine's own operator table (the
  C++ ground truth the catalog is pinned to), rejects duplicate or
  malformed `instance` labels (they own parameter namespaces), validates
  existence, carrier adjacency, entry/exit, and the arena, length and
  parameter-count budgets — whose authoritative limits (a single arena's
  capacity, the chain-length cap, the chain-wide parameter cap, alongside
  the per-op block sizes already cataloged) are
  **exported with the catalog**, so the editor can account for aligned
  totals before offering an edit; the transactional double-buffering
  below is engine-internal and never inflates or halves this figure —
  and only then lays out the program. Per-chain monotonicity
  needs no separate walk: each operator's
  `family_of(in) <= family_of(out)` is proven once at operator-table
  construction and pinned by a test, and adjacency over monotone
  operators yields a monotone chain. The program is an arena-backed
  array of ops whose param-block and prepared-state offsets are
  computed internally and never cross the API. Family is not stored at
  all; it derives from the operator's carrier pair. Compilation is
  **transactional**: any failure is a structured refusal (code +
  offending `entry_index`, −1 for the chain as a whole) that leaves the
  previous program, the registered
  parameter definitions, the parameter generation, and all live
  instance state unchanged. The tool's catalog validation is editor UX,
  not the trust boundary.
- **Each operator-table entry is a runtime operator ABI**: parameter-block
  size and alignment, prepared-state size and alignment, instance-state
  size and alignment (compilation arena-allocates, migrates, and
  destroys those blocks, so their layout is table data exactly like
  the other two — declared as the worst case across topology variants,
  matching the eager per-variant construction below), the callbacks
  below. Resource dependencies are implemented by the callbacks rather than
  declared as descriptor metadata:

  ```cpp
  construct_params(void *params);
  param_address(void *params, uint16_t schema_index) -> void *;
  validate(const void *params) -> const char *;
  // InstanceId carries the (instance_id, operator_id) pair identity.
  init(void *dst, InstanceId);                   // infallible; construct owned resources
  migrate(void *dst, const void *src, InstanceId) -> Status;  // src untouched
  destroy(void *state);
  advance(void *state, const uint8_t *params);   // per frame, steps clocks
  prepare(const FrameContext &, const uint8_t *params,
          const void *state, uint8_t *prepared);
  run(const void *in, void *out, const FrameContext &,
      const uint8_t *params, const uint8_t *prepared);
  ```

  **Instance state** is what parameter and prepared blocks cannot
  cover: persistent per-instance accumulators and owned resources —
  phase clocks, initialized noise generators — the runtime analogue of
  ShaderWorkbench's bounded arrays of pre-initialized noise resources, and
  it is reachable from execution: `advance` steps clocks once per
  frame, then `prepare` reads the updated state to derive the frame's
  prepared block (`run` needs only `prepared`). Instance state never
  caches parameter-derived values — anything parameter-derived is
  `prepare`'s job, recomputed per frame — which is what makes `init`
  independent of the apply-values-after-compile ordering: `init`
  constructs owned resources (noise seeded from a stable hash of the
  `InstanceId`) and zeroes accumulators, nothing more. `init` is
  **infallible by contract** (nothrow, no `Status`): it
  placement-constructs deterministic resources into storage whose
  capacity the compile transaction already proved from the ABI's
  instance-state size, and it allocates nothing itself — so the
  transactional guarantee needs a failure path only from `migrate`,
  where duplicating live accumulated state is genuinely fallible. An
  operator whose fresh construction could fail has no home in this
  ABI; admitting one is a spec revision, not a throw out of `init`.
  And `init`
  constructs resources **eagerly for every topology variant** the
  operator declares, which is what the worst-case budget actually pays
  for: a runtime topology switch selects among already-constructed
  resources and never constructs or destroys, so the state layout is
  topology-invariant by construction (ShaderWorkbench's pre-initialized
  noise arrays are the precedent). State identity
  is the `(instance_id, operator_id)` pair: a structural edit
  `migrate`s the state of instances whose pair survives (an unchanged
  warp keeps its phase and warmed noise across an edit elsewhere),
  while the same label with a different operator gets a fresh `init`,
  never a migration; removed instances get `destroy`. `migrate`'s
  postconditions are explicit: on success, `dst` is a **transactional
  clone** — owned resources duplicated, never moved and never shared
  (`src` is const, and a move would mutate it; the commit's destruction
  of the losing arena is what releases the old resources), so `dst`
  remains valid and destroyable after `src` is destroyed; on failure,
  `dst` is left **unconstructed** — the callee
  rolls back any partial construction before returning its `Status`,
  and the arena treats a failed `dst` as never constructed (it is not
  `destroy`ed). `migrate` must also
  leave `src` untouched until commit: compilation builds a candidate
  state arena (`init` for new pairs, `migrate` for surviving ones — a
  failing `Status` aborts the compile transactionally, and on abort
  every candidate state that compile already constructed, `init`s and
  completed `migrate`s alike, is destroyed before the refusal returns:
  the candidate arena tears down as a unit), swaps only on success, and
  `destroy`s the losing side. That transaction implies **two live
  storage sets** — old and candidate op arrays, param, prepared, and
  state blocks, and owned resources coexist until commit — so the
  memory model is explicit: **two fixed arenas with an active index**,
  each sized to the full exported budget. Compilation lays the
  candidate out in the inactive arena while the active one stays
  untouched and running; success flips the index and tears down the
  loser; refusal tears down the candidate only. The steady-state
  transaction therefore never allocates and cannot fail for memory —
  the budget check runs at validation, before layout, against **one
  arena's capacity**, which is the figure the catalog exports; the
  doubling is engine cost, invisible to the editor's accounting, and
  cheap where the interpreter actually runs (Emscripten and
  test-oracle builds only — never the device). Compilation budgets instance
  state worst-case for operators whose topology parameters can switch
  resource needs, so a preset change never demands an allocation
  mid-run.

  **`FrameContext`** carries the base projection orientation, three borrowed
  baked-palette pointers, and the hue-rotation and hue-noise LUT pointers.
  Parameters and instance state reach callbacks separately. Borrowed pointers
  remain valid within the frame that owns them.

  The callbacks are thin adapters over the **shared carrier kernels of
  §4** — the same free functions the template combinators call — with
  the adapter reading the op's parameter block where a static provider
  would read a named `FrameState` slot (the static provider wrappers
  are compile-time-bound to named instances and cannot serve an arbitrary
  third instance; §7.2's allocation limit is the same fact seen from
  the promotion side). ShaderWorkbench's dynamic backend is the existing
  precedent for this shape.
- **The erased carrier ABI is explicit**, because a homogeneous op array
  cannot invoke heterogeneously-typed callbacks unaided. Evaluation owns
  two carrier slots whose size and alignment are derived by
  `static_assert` from the closed carrier set (§2) — today that means
  sized for `PlaneSample` — and ping-pongs between them. Every
  canonical carrier is trivially destructible (§3's widening
  obligation), so placement-constructing into a slot reuses its storage
  and ends the previous carrier's lifetime with no destructor call;
  each adapter placement-constructs its output into `out`, and the
  reinterpretation of `in` is sound because compile-time adjacency
  validation guarantees the previous op's output carrier is exactly
  this op's input carrier — the validated carrier pair in the operator
  table is what licenses the erasure. Sharing the canonical carrier
  structs prevents **layout** drift only; semantic parity with the
  template path is engineered, not assumed — the shared kernels above
  are the mechanism (Project's sphere assignment, Sample's
  weight/ramp/coverage sequence, and Coverage's accumulation exist
  once), and a **static-versus-erased parity test covers every
  operator × topology variant** — variant coverage is not optional,
  because a default-variant pass would miss exactly what varies: enum
  dispatch, variant-gated parameters, per-variant resource selection,
  and prepared-state differences — with parameter values at
  representative boundaries (defaults and range endpoints) in each
  variant. Per-op prepared state is arena-sized at compile
  (structural edit) time, like the param blocks — the static-blob +
  `static_assert` pattern that guards template pipelines cannot cover
  unbounded chains.
- `setShaderChain` replaces the thirteen structural enum parameters. It
  synchronously rebuilds the registered parameter definitions and bumps
  the parameter generation before returning — an async rebuild would let
  preset values apply against a stale definition snapshot. The
  registered definitions are each operator's **full union schema**
  (§7.1), which is what makes registration computable from the wire
  payload at all: topology values arrive only later, through the same
  value channel as everything else, and flipping one switches which
  fields the variant *reads*, never which fields *exist*. Apply order
  is fixed: `setShaderChain` → `setShaderChainParameters` →
  `syncEffectGui` → `invalidate`.
- Device exclusion is named: the interpreter lives behind the workbench
  build-flag pattern (`HS_ENABLE_*` defaulting to 0 outside
  Emscripten/test-oracle builds, `#error` under `ARDUINO`), its effect is
  excluded from `HS_PHANTASM_EFFECT_LIST`. No release-ELF symbol inspection
  gate currently verifies interpreter exclusion.
- Promotion of a document into a composed effect follows §7: the
  typed instance allocation and budgets must succeed and topology parameters must be
  bank-invariant. The interpreter is never shipped to the device.

## 9. What this unlocks

At the pipeline layer (hand-written effects, ShaderWorkbench studies, the
workbench interpreter), **immediately, with the shipped vocabulary
recombined**: lens sandwiches and double Mobius; displacement at any
depth relative to lenses; any warp count and order; transfer and
value-cutout chains; projection-free sphere-sampled sources use the
`SPHERE→FIELD` `Stage::SampleSphere` crossing, and ray shaders use the shipped
`SPHERE→COLOR` `RayStage` crossing. **Admitted by the rules but
awaiting a new combinator**: COLOR grading stages (the family ships empty).
The rules make these one-combinator additions instead of schema changes;
they are not day-one capabilities. At the ComposedEffect layer: any
chain whose provider allocation succeeds (§7) — every shipping effect,
plus recombinations with independently parameterized warp and lens instances.
Topology, catalog conformance and memory budgets govern promotion acceptance.

Combining two sources is a policy-level concern, not a chain-shape gap:
monotonicity correctly forbids a second `Sample` crossing, and source
policies already receive the full carrier, so product or max combinator
policies over two sources compose fields with zero pipeline change at
the template layer; none ships. The interpreter would ship such a
compound source as its own curated operator entry — a scalar
source-expression ABI (nested operators, each with its own parameter and
prepared blocks) is deliberately out of this revision, and joins true
Fork/Join branching as a reserved future extension, taken up only if the
curated set grows unwieldy.

The integration surface for any new operator is one `OperatorDescriptor`
record (§7.1), from which the promotion, interpreter, and tool views are
all generated — no schema change, no new banks.

## 10. Location, namespace, and dependency rules

The public facility lives under `core/render/pullback/`, in namespace
`Pullback`, alongside `Scan`, `Filter`, and `SDF`. `core/render/pullback.h`
is the umbrella over the composition core only: the carrier contract, the
field tables, the surface, lens, projection, warp, source, material, and
color policy families, the ray stage (`core/render/pullback/ray.h`), and the
stage combinators. The chain interpreter
(`core/render/pullback/interpreter.h`, `core/render/pullback/operator_model.h`,
`core/render/pullback/operator_table.h`, `core/render/pullback/operators.h`
with its per-family `operators_*.h` headers, and
`core/render/pullback/catalog_export.h`), the composed-effect base
(`core/render/pullback/composed_effect.h`), and the shared runtime seeds
(`core/render/pullback/runtime_seeds.h`) are not reachable from the umbrella;
their consumers include them directly.

Public groups are:

```text
Pullback::Pipeline                 typed ranked-chain coordinator
Pullback::CodeEmission             placement metadata
Pullback::SphereSample             rank-0 sphere carrier
Pullback::PlaneSample              rank-1 planar carrier
Pullback::SurfaceResult            one sphere-space map result
Pullback::WarpStepResult           one planar-warp result
Pullback::FieldSample              rank-2 scalar carrier
Color4                            rank-3 color carrier
Pullback::Field / Fields           field-table records and their curve, interpolation, and validity helpers
Pullback::Stage::*                 ranked stage combinators
Pullback::RayStage                 fused ray-query and appearance stage
Pullback::Kernel                   shared carrier kernels the combinators and erased adapters call
Pullback::Surface::*               sphere-space map policies
Pullback::Lens::*                  lens policies
Pullback::Projection::*            projection policies
Pullback::Warp::*                  planar-warp policies
Pullback::Source::*                scalar-source policies
Pullback::Weight::*                signal-weight policies
Pullback::Transfer::*              value-transfer policies
Pullback::ProjectionCoverage::*    projection-coverage policies
Pullback::ValueCoverage::*         value-coverage policies (`ValueCutout`)
Pullback::Color::*                 colorization policies and kernels
Pullback::Interp                   chain interpreter: operator model and table, `Op::*` operators, catalog export
```

Carrier declarations live in `core/render/pullback/contract.h`; `Color4` is
defined in `core/color/pixel.h`. See §3.

`pullback.h` may include headers from `core/math`, `core/color`,
`core/animation` (`Animation::RippleParams` is defined in
`core/animation/params.h`), and the minimal engine concept/profiling headers
it needs. `core/render/pullback/contract.h` includes color and math headers; the engine dependency
arrives through color headers and `core/animation/transformer.h`. It shall not include an
`effects/` or `workbench/` header, refer to `ShaderWorkbench`, or require the
effect registry.

Existing pure mathematical kernels remain in their natural owners:

- lenses remain in `core/math/lenses.h`;
- projection kernels and `projections::ProjectionKernelResult` remain in
  `core/math/projections.h`;
- shared noise and stereographic helpers remain in their existing core math
  headers;
- palette/gamut operations remain in `core/color`.

`pullback.h` supplies typed policies and orchestration around those kernels.
Core and consumer adapters share the same mathematical kernels.


## 11. No universal frame type

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
Providers are empty compile-time adapters with inline static accessors
(`always_inline` in ComposedEffect).
They neither own nor copy state.

This is the principal decoupling boundary: core owns algorithms and carriers;
the effect owns frame layout and maps it into those algorithms.

### 11.1 Provider contract

A concrete operator is parameterized by a provider, not by an effect:

```cpp
struct OuterWarpState {
  using Binding = ShaderWorkbenchBinding;
  using FrameState = typename Binding::FrameState;

  static const auto &params(const FrameState &);
  static auto prepare(const FrameState &);
  static float phase(const FrameState &);
  static const FastNoiseLite &noise(const FrameState &);
  static bool path_length_required(const FrameState &);
};
```

Only the accessors needed by the selected operator are required. For example,
`Warp::MirrorTile` does not require `noise`, and `Lens::Glitch` requires no
provider at all. Accessor return requirements are structural and documented at
the operator declaration: a wave-shear parameter view must expose
`strength` and `frequency`; an affine prepared record must expose the fields its
formula reads. `params(frame)` returns an existing consumer record by const
reference; `prepare(frame)` returns the prepared record by value, and the
policy names that return type `Prepared`. Core shall not require construction
of a per-pixel view object.

Every provider:

- is empty and trivially constructible;
- names its owning pipeline `Binding` and derives `FrameState` from that
  binding;
- names `FrameState` exactly;
- exposes only inline static const-frame accessors (`always_inline` in ComposedEffect);
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
| surface map | `prepare(frame)` and `path_length_required(frame)`, plus the subset of `params(frame)`, `phase(frame)`, `noise(frame)`, `scale(frame)`, and `strength(frame)` the selected map reads |
| projection | prepared frame `conjugate(frame)` plus scalar `singularity_fade`, `central_meridian`, `coordinate_scale`, `standard_parallel`, and `layout_scroll` accessors as required by the selected map; edge-distance demand is the `EdgeDistanceRequired` template argument of the Peirce and Airocean policies, not a provider read |
| planar warp slot | `params(frame)` and `prepare(frame)`, plus `path_length_required(frame)`, `phase(frame)`, and `noise(frame)` where the selected policy reads them; basis, envelope, integrator, and polar mode are template facts in a compiled policy |
| source | `params(frame)`, `prepare(frame)`, and optional noise resource/time accessors; the selected source policy determines the required subset |
| material | value/coverage scalar accessors (`iso_level`, `iso_width`, band values, cutout values, `edge_width`) required by the selected policies |
| color | immutable color parameters/clocks, generated palette binding, prepared hue-rotation LUT, prepared hue-noise LUT, and deliberately runtime mapping/brightness/hue mode values |

An operator's declaration narrows this table with a `requires` expression that
names every field/member it reads and no unrelated member. For example,
`WaveShear` requires `strength`, `frequency`, `phase(frame)`, prepared rotation
sine/cosine, `path_length_required(frame)`, and `edge_width` only under the
edge-fade envelope; `MirrorTile` requires `cell_x`, `cell_y`, prepared rotation
sine/cosine and mirror transform, and `path_length_required(frame)`. These
requirements are part of the public doxygen contract. Adding a new hot-path
read therefore changes the provider concept and its tests in the same commit.

For compiled policies, a provider's `Binding` must exactly equal the enclosing
stage's `Binding`; provider concepts expose this as a testable boolean before a
`static_assert`. This is what turns accidental use of a provider from another
effect or frame layout into a named binding diagnostic rather than a deep
substitution error.

### 11.2 Instrumentation

Moving code into core shall not erase ShaderWorkbench's stage buckets or bake
ShaderWorkbench profile fields into core.

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
`MATERIAL`, and `COLOR`. `MIRROR_TILE` is nested inside `PLANAR_WARP`;
its cycles are a subset and must not be added again when totaling stage time.
`NoInstrumentation` compiles to no statements.
ShaderWorkbench currently supplies a no-op hook policy. The event-to-counter
mapping described by this design was not implemented.

## 12. Ownership, lifetime, and mutation

- The consumer's `FrameState` is immutable for the duration of a draw.
- Providers borrow only state reachable from that frame.
- Mutable `FastNoiseLite`, palette cyclers, animation objects, generated
  palettes, LUT storage, and arenas remain consumer-owned.
- Frame resource pointers are const bindings whose owners outlive the draw.
- No core stage calls a resource setter, advances a clock, steps an animation,
  prepares a transform, or allocates.
- No provider returns a reference to a temporary. Parameter records are returned
  by const reference; prepared records and LUT views by value.
- The coordinator and policies contain no objects, so a pipeline has no
  lifetime independent of the frame.
- Transition rendering remains sequential: prepare and consume one endpoint
  before shared backing storage is overwritten.

The existing ShaderWorkbench stack, persistent arena, RAM2, and effect-heap budgets
remain unchanged. Public carriers add no allocation or hidden ownership.

## 9. Executable snapshots and shader host retirement

**IMPLEMENTED.** `ShaderChain` is the sole simulator authoring host. The slot
host, fixed-slot parameter configuration, admission fold and full-configuration
WASM channel have been removed. Only current effect identities and document formats are accepted.
The program, instance clocks, noise seeds, walk state and all generated palette
cycles share the typed transactional contract in
[chain_snapshot_spec.md](chain_snapshot_spec.md). Capture producers execute
versioned chain fixtures and compare their frame bytes and exact-kernel oracle
metrics against the frozen before corpus at both supported resolutions.
