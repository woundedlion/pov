# Pullback stage families: arbitrary chains over ranked carriers

**Status: §§1–6 and §8 LANDED; §7 PROPOSED.** The static ranked pipeline and
its migration (the §6 contract cut-over) ship, as does the preview
interpreter (§8: `core/render/pullback/interpreter.h`, `workbench/ShaderChain.h`,
and the `setShaderChain` binding). Promotion/verification (§7) is not yet
implemented. Where this spec and
[pullback_pipeline_spec.md](pullback_pipeline_spec.md) disagree about the
stage model (stage kinds, slot count, carrier records), this spec describes
what ships. The document is three systems with distinct invariants and
failure modes, layered in
order: **the static ranked pipeline and its migration (§§1–6)**, which
is self-contained; **promotion and verification (§7)**, a resource
allocation and CI layer over it; and **the preview interpreter's engine
contract (§8)**, a runtime mirror of it. Each later part depends only on
the parts before it. The authoring tool is separate:
[shader_workbench_chain_spec.md](shader_workbench_chain_spec.md).

## 1. The problem

`Pullback::Pipeline` validates a hard-wired shape: exactly six stages, one
each of `OUTER_CAMERA, SURFACE_PROJECT, PLANAR_WARP, SOURCE, MATERIAL,
COLOR`, in that order, with slot-specific carrier types between them
(`contract.h` `ORDER`/`ARITY`/`CARRIERS`). All flexibility lives *inside* a
slot as policy parameters: `SurfaceProject` owns four policy slots
(pre-lens surface, lens, post-lens surface, projection), `PlanarWarp` a
variadic policy list, `Material` a fixed weight/transfer/coverage triple.

Consequences:

- No chain the slots didn't anticipate: two lenses with a rotation between
  them, displacement on both sides of the lens, a second transfer, any
  color processing after the palette.
- Unused capability is expressed as `Identity` policies filling mandatory
  slots, which is why `SurfaceProject::run` carries an Emscripten-only
  flattened path for the identity-lens, identity-post-surface case.
- The carrier zoo (`SurfaceResult`, `WarpResult`, `SourceInput`,
  `MaterialInput`, `MaterialSample`) exists to serve slot boundaries, not
  domains, so a new stage shape means a new carrier and new validation rows.
- The shader workbench mirrors the rigidity: thirteen fixed option banks,
  most set to "None", and a document schema whose graph is pinned to "the
  linear six-role graph" despite reserving 256-node limits.

The filter pipeline solved the same problem differently: stages carry a
`domain_rank` (World = 0, Screen = 1, Pixel = 2), the chain must be
non-decreasing in rank, and *within* a domain any stage sequence is legal.
That is the model to mirror.

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
`SPHERE→FIELD` crossing that never projects; a `SPHERE→COLOR` sky shader
is expressible the day someone writes that combinator.

Family *names* survive as the rank's vocabulary — in diagnostics ("a
Plane stage may not follow a Field stage — families are ordered Sphere,
Plane, Field, Color", keeping the filter pipeline's prose-assert style)
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
  uint8_t traits;
  uint8_t edge_class;
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

Mapping from today: `Vector`+`SurfaceResult` → `SphereSample`;
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
chain's seed (`stage.h:221`), which is precisely the split the crossing
now performs. `WarpStepResult` survives as the warp policy protocol
type.

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
The one genuinely deferred coverage policy is `ValueCutout`, whose
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
  `CoverageKind`: `PROJECTION` → `Weight`, `PROJECTION_SQUARED` →
  `WeightSquared`, `EDGE_FADE` → `EdgeFade`. ShaderWorkbench
  `CoveragePolicy`: `OPAQUE` → `None`, `PROJECTION_WEIGHT` → `Weight`,
  `PROJECTION_WEIGHT_SQUARED` → `WeightSquared`, `EDGE_FADE` →
  `EdgeFade`, `VALUE_CUTOUT` → `None` at the crossing plus a
  `Stage::ApplyCoverage<ValueCutout>` FIELD stage. That deferred stage's
  signature is `apply(float value, const FrameState &) → float`
  (dropping the `ProjectionSample` parameter today's version ignores),
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
  does the same with neutral weight. The minimal legal chain over the
  shipped crossings is `Project → Sample → Colorize` — `Rotate` is the
  camera, present in every real effect but not required by ENTRY, which
  `Project`'s `SphereSample` input already satisfies. (A future
  `SPHERE→COLOR` combinator would make a one-stage chain legal; the
  rules, not this list, define the minimum.)
- `PlanarWarp`'s variadic policy list becomes N consecutive `Warp` stages.
- `Identity` policies stop appearing in pipelines — absence of a stage
  *is* the identity.
- Today's three crossings are not a closed set: any combinator whose
  Output rank exceeds its Input rank is admitted by the same rules, which
  is how a `SPHERE→FIELD` noise source or a `SPHERE→COLOR` sky stage
  arrives — as a new combinator, not a schema change.
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

Contract changes: `Detail::StageContract` drops `KIND` and `TERMINAL`;
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
predicate** — `Bindable<Descriptor, Binding>` — evaluated *before* the
bind step instantiates anything, so the pipeline's named `BINDINGS`
assertion (§5) is what reports a foreign binding, with today's
message, rather than a template-formation abort inside `Bind`
preempting the named surface; `Bind`'s internal `PROVIDER_VALID`
asserts remain as a backstop, unreachable through pipeline assembly.
Execution uses only bound types. This is
what keeps binding services available inside combinators: `mark`/`span`
resolve through the bound stage's binding exactly as
`BindingT::Instrumentation` does today. `HasStageContract`/
`HasTypedStageContract`, `BINDINGS`, and `PREPARES` are reformulated
over descriptors and bound leaves — replaced, not carried over.

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
| `Sample` | `SOURCE` around the source-policy call; `MATERIAL` around weight + ramp + projected coverage |
| `Transfer`, `Coverage` | `MATERIAL` |
| `Colorize` | `COLOR` |

`MATERIAL` thus becomes the sum of up to three spans (Sample's tail,
Transfer, Coverage) covering exactly the work today's single Material
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
transparency rule above. That is the mechanism by which ShaderWorkbench's
`ExtraValidation` ("`EDGE_DISTANCE_UNCONDITIONAL` on the projection
requires an edge-fade coverage") is re-expressed once positional slot
typedefs are gone.

## 5. Validation

The staged structure (earlier levels gating later ones) and the
named-boolean surface stay, with one more level than today. The order
is explicit: descriptor contract shape, then `CANONICAL`, then the
rank rows (`MONOTONE`, adjacency, `ENTRY`, `EXIT`), then the bound
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
| `CARRIERS` (slot-typed)    | adjacency fold, plus `ENTRY` (first Input is `SphereSample`) and `EXIT` (last Output is `Color4`) |
| `TERMINALS`    | retired — subsumed by `EXIT` + monotonicity; a mid-chain `Color4` producer followed by COLOR endomorphisms is now a feature |
| `EMPTY_DESCRIPTORS`, `RUN_RETURNS`, `EXTRA_VALIDATION` | unchanged, evaluated over the flattened leaf stages |
| `CONTRACTS`, `BINDINGS`, `PREPARES` | reformulated over the descriptor contract (§4): every leaf is bound by the pipeline's bind step; a descriptor whose policy names a foreign binding fails there |
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
- **`Pullback::Params` and its eight parameter groups** — source,
  projection, outer_warp, inner_warp, surface, lens, value, color
  ("groups", not "families": they are unrelated to the carrier
  families of §2). Preset schemas, slider registration, and
  interpolation need a fixed, nameable layout. The relaxation is a
  *pipeline-shape* concern; the parameter-schema layer assembles a chain
  instead of six slots. Do not try to make `Params` variadic.

### Scope boundary

Promotion — turning an authored chain into a ComposedEffect — is a
resource-allocation problem over `Params`'s finite provider slots,
specified in §7. The boundary that matters here: chains whose allocation
fails are the territory of the workbench interpreter (§8) and
hand-written effects; the pipeline model itself imposes no such limit.

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
placement a computed value. The sphere run keeps both displacement
slots — today's derivation places curl displacement *before* the lens
and `DirectNoise` *after* it (`DIRECT_SURFACE`), and capture equivalence
requires preserving that:

```cpp
using PreDisplaceStage = std::conditional_t<
    HAS_SURFACE_NOISE && !DIRECT_SURFACE, Stage::Displace<CurlPolicy>, void>;
using PostDisplaceStage = std::conditional_t<
    DIRECT_SURFACE, Stage::Displace<DirectNoisePolicy>, void>;

using SphereRun = Stage::Placed<
    HAS_SURFACE_NOISE ? CodeEmission::OUT_OF_LINE_FLASH
                      : CodeEmission::INLINE_ONLY,
    PreDisplaceStage, LensStage, PostDisplaceStage, ProjectStage>;
```

Those aliases live in `composed_effect.h` beside the `*PolicyFor`
metafunctions they extend.

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
     `ExtraValidation` → the §4 trait folds;
   - ShaderWorkbench's dynamic backend → the canonical carriers. The
     projection-join facility (`join_projected` /
     `projection_join_compatible`) has no render-path caller today — it
     is a test-pinned facility reserved for the lens-blend transition —
     and it retypes to a `PlaneSample` join utility preserving its exact
     field-by-field algorithm: `coords` lerped; the discrete topology
     and edge fields (`region_id`, `component_id`, `boundary_flags`,
     `fade_edge_distance`, `flags`, `traits`, `edge_class`) selected
     from the nearer branch (`mix < 0.5` → direct); `value_weight`
     **recomputed** from the blended coordinates via pole attenuation;
     `domain_coverage` and path length lerped; `sphere`
     normalized-lerped (`nlerp_unit`). Strict projections (Bonne,
     Peirce, Airocean) still refuse the plane join. When the transition
     feature lands, the chain-level representation is fixed now:
     compatible projections may use the join inside a curated compound
     `SPHERE→PLANE` operator, while strict projections require
     complete-output blending, which a linear chain cannot express and
     which is therefore an **evaluator-level two-pass** — shade both
     branches, blend the resulting `Color4`s — outside the chain program
     entirely. Blending a single plane sample is not equivalent
     (everything downstream is nonlinear), so no `PlaneSample`
     interpolation may be substituted. `PREPARED_BLOB_BYTES` re-verified
     against the larger prepared tuples;
   - `test_pullback.h` / `test_shader_workbench.h` validation and
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

A layer over the static model, with its own invariants and failure
modes: allocation can fail where chain validation succeeds, and its
gates run in CI rather than at compile time.

### 7.1 One operator authority

A promotable operator is visible to three systems — promotion (its
provider requirements), the interpreter (its runtime ABI entry), and
the authoring tool (its catalog entry). These are **views of one C++
record**, an `OperatorDescriptor`: a **type-level stage recipe** —
`template <typename Binding, typename Assignment, typename Topology>
Stage` — rather than one concrete stage type, because both remaining
degrees of freedom are chosen *after* the catalog is written:
allocation selects providers (today's policies encode their slot in the
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

Promotion is an explicit **resource-allocation problem**: a
capacity-aware assignment of stage instances to provider slots. Not a
stage count, and not set inclusion either — a set loses multiplicity,
and three independently parameterized warps must not collapse into a
requirement two slots could satisfy. Requirements are **promotion-catalog
metadata, not stage-contract members** (§4), and they are a **function
of topology**: the stage recipe varies with `Topology` (§7.1), so its
resource needs vary with it — `Sample` requires edge-fade width state
only when its coverage mode is `EdgeFade`. One unconditional list per
operator would either over-reject (demand the EdgeFade slot for a
`Weight`-coverage document) or under-allocate. Each promotable
operator's catalog entry therefore exposes two views:
`requirements(Topology)`, the parameter-schema and resource instances
the recipe's policies read at those topology values, which promotion
evaluates at the document's pinned topology (§7.3) — recursively
composed for compound policies (a compound
source like `Source::Multiply<A, B>` concatenates its children's) —
and the worst-case resource and state footprint across variants, which
is what the interpreter budgets (the same worst-case-across-variants
rule its eager per-variant construction and topology-invariant
approximation metadata already follow, §7.4/§8). An operator without
requirement metadata simply is not promotable, which is the safe
default. The slots `Params` and `FrameState` back have fixed
capacities: clocked warp × 2 (`outer_phase`/`inner_phase`,
`outer_rotation` on the first), Mobius lens family × 1, orientation × 1,
surface family × 1, source family × 1, projection family × 1, value
family × 1, color family × 1. The value slot has a synthesis rule rather
than a count: promotion constructs `ValueT` as the family covering the
*union* of the value fields the chain's requirements name, so
`IsoContour` + `EdgeFade` promotes (their field sets are disjoint) while
two independently parameterized `IsoContour` stages do not (one family
cannot carry two `iso_level` sets). The mechanism is a reusable
composition, not an undefined merge: `Value::Combine<Families...>`
inherits each required family, so field access by name resolves through
the base and `ValueProvider` is unchanged. The field-table half rests
on one language fact, not a redesign of `fields.h`: `Field<Owner>`
holds an owner-typed member pointer in a homogeneous array, and a
`float Base::*` converts implicitly — including in constant
expressions — to `float Combined::*`, so `Combine`'s `FIELDS` is a
constexpr concatenation that rebuilds each base entry as a
`Field<Combined>`: a genuine homogeneous
`std::array<Field<Combined>, N>` that `Fields::interpolate`, `valid`,
and slider registration iterate unchanged. Disjointness cannot hang on
the display `name` — it is nullable and cosmetic — so `Field` gains a
**stable machine `id`**: non-null, unique within its family, and the
same token the catalog parameter schema exposes as the `field` segment
of `<label>.<field>` ids (§7.1's schema generation reads it from the
table — one authority for both spellings). The compile-time assert is
id-disjointness across the combined table; repeating one family
(`Combine<IsoContour, IsoContour>`) is ill-formed before any assert
runs, as duplicate direct base classes — which *is* the
double-`IsoContour` rejection. Two instances share one slot only by
explicit declaration of a shared requirement; distinct instances are the
default. Slot assignment alone is type-correct but not value-correct,
so the `Assignment` also owns a **parameter-binding table**: a mapping
from every `(instance_id, field_id)` the chain's schemas declare to
the concrete `Params` member and control registration the emitter
writes. The mapping is **injective** — one document parameter, one
storage slot — except where a shared requirement declares an explicit
**alias**, and aliased parameters must agree everywhere the document
speaks: identical defaults, equal values in every preset, and
identical transition scheduling, checked at promotion with
disagreement a refusal. Nothing upstream enforces this — the workbench
gives every instance an independent parameter namespace, so two
"shared" instances can freely diverge in the document, and a promoted
effect with one storage slot would render that divergence differently
than the interpreter previews it. The binding table is also one side
of §7.1's conformance test: promotion bindings and the interpreter's
registered parameter schema must describe the same field set. Requirements declare their **compatible slot sets**, because slots are
not interchangeable: only the outer warp slot carries `outer_rotation`,
so a rotation-consuming warp (`AffineFrame`) is outer-only while a
plain clocked warp accepts either — today an affine warp bound to the
inner provider silently receives zero rotation, which is exactly the
mis-assignment chain-order first-fit would produce for
`WaveShear → Affine`. Allocation computes a **deterministic matching**
over compatible slots, and determinism is by canonical construction,
not a tie-break phrase: slots are ordered as this section declares
their capacities (outer warp before inner, and so on), instances in
chain order, and the assignment is the **lexicographically minimal
complete matching** — each instance in chain order takes the earliest
compatible slot that still leaves a complete matching for the
remaining instances. (Chain order alone does not choose uniquely among
multiple complete matchings, and the choice is semantic: the
assignment determines the generated provider types and parameter
mappings.) The instance
labels record the mapping, and promotion succeeds iff a complete
matching exists — `WaveShear → Affine` therefore emits Affine → outer,
WaveShear → inner. The assignment never enters the
pipeline type system: the promotion emitter resolves each requirement to
its slot and writes descriptors *already specialized* with that slot's
provider (`WarpProvider<B, Outer>`, …) into the generated
`Pipeline<...>` typedef — instantiations of the operator's stage recipe
(§7.1) with the allocation's assignment and the document's pinned
topology values. `Bind` receives no assignment, and hand-written
pipelines have always been written in exactly this already-specialized
form. Operators with empty requirements consume
nothing: a second parameterless lens (`Glitch`, the kaleidoscopes —
stacked as consecutive `Lens` stages) promotes today, as does a lens
sandwich whose
only parameterized member is the single Mobius. A chain whose
assignment fails (third clocked warp, second Mobius family) is
interpreter-only until a per-instance provider design exists (indexed
`WarpProvider<B, Slot>`, a phase array in `FrameState`, per-instance
registration — a named follow-on, not part of this spec).

### 7.3 Topology parameters

**Topology parameters** — the catalog-flagged enum8 class (weight mode,
coverage mode, noise basis) — are ordinary runtime switches in the
interpreter, free to vary between presets; but promotion pins them as
template arguments in one typedef, so a promotable document must hold
them invariant across its preset bank. A bank that varies them stays
interpreter-only or splits into one document per topology.

### 7.4 Approximation acceptance

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
  of `setShaderChain` is an ordered list of `{instance_id, operator_id}`
  and nothing else — no offsets, no family tags, nothing layout-shaped,
  and **no parameter values**: values flow through the existing
  per-parameter channel (`setParameter`, keyed `instance_id.field`,
  validated against the operator's schema) *after* compilation, per the
  apply order below, so the transaction boundary is the chain compile
  alone and values always apply to a committed program. `instance_id`
  is the document's chain-entry label: operator ids alone cannot
  distinguish `warp1` from `warp2`, and the engine needs the instance
  identity to register per-instance parameter definitions before values
  apply. Structural variants (weight mode, coverage mode, noise basis)
  are enum8 parameters arriving through that same value channel, so a
  chain entry carries no field the engine drops. **Three identities
  are distinct, related by projection**: the *document digest* covers
  the whole descriptor — chain, parameter schemas and defaults,
  presets, path policies; the *program-shape identity* is the ordered
  `{instance_id, operator_id}` list, exactly what `setShaderChain`
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
  malformed `instance_id`s (they own parameter namespaces), validates
  existence, carrier adjacency, entry/exit, and the arena/length
  budget — whose authoritative limits (a single arena's capacity, the
  chain-length cap, alongside the per-op block sizes already
  cataloged) are
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
  offending instance) that leaves the previous program, the registered
  parameter definitions, the parameter generation, and all live
  instance state unchanged. The tool's catalog validation is editor UX,
  not the trust boundary.
- **Each operator-table entry is a runtime operator ABI**: parameter-block
  size and alignment, prepared-state size and alignment, instance-state
  size and alignment (compilation arena-allocates, migrates, and
  destroys those blocks, so their layout is table data exactly like
  the other two — declared as the worst case across topology variants,
  matching the eager per-variant construction below), the callbacks
  below, and the shared frame
  resources (noise fields, palette, LUTs) the operator reads:

  ```cpp
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

  **`FrameContext`** is the interpreter's per-frame snapshot: the
  interpolated parameter view, the spatial transforms, and borrowed
  shared resources (the palette bake, LUTs, noise fields). Its pointers
  alias engine-owned storage and are valid only within the
  `draw_frame()` that built it — the same lifetime contract as
  ComposedEffect's `FrameState`.

  The callbacks are thin adapters over the **shared carrier kernels of
  §4** — the same free functions the template combinators call — with
  the adapter reading the op's parameter block where a static provider
  would read a named `FrameState` slot (the static provider wrappers
  are compile-time-bound to fixed slots and cannot serve an arbitrary
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
  is fixed: `setShaderChain` → apply preset values by id →
  `syncEffectGui` → `invalidate`.
- Device exclusion is named: the interpreter lives behind the workbench
  build-flag pattern (`HS_ENABLE_*` defaulting to 0 outside
  Emscripten/test-oracle builds, `#error` under `ARDUINO`), its effect is
  excluded from `HS_PHANTASM_EFFECT_LIST`, and the release-ELF symbol
  inspection extends to the new entry points.
- Promotion of a document into a composed effect follows §7: the
  capacity assignment must succeed and topology parameters must be
  bank-invariant. The interpreter is never shipped to the device.

## 9. What this unlocks

At the pipeline layer (hand-written effects, ShaderWorkbench studies, the
workbench interpreter), **immediately, with the shipped vocabulary
recombined**: lens sandwiches and double Mobius; displacement at any
depth relative to lenses; any warp count and order; transfer and
value-cutout chains; projection-free sphere-sampled sources use the
`SPHERE→FIELD` `Stage::SampleSphere` crossing. **Admitted by the rules but
awaiting a new combinator**: sky shaders
(`SPHERE→COLOR`), and COLOR grading stages (the family ships empty).
The rules make these one-combinator additions instead of schema changes;
they are not day-one capabilities. At the ComposedEffect layer: any
chain whose provider allocation succeeds (§7) — every shipping effect,
plus
recombinations including parameterless-lens stacks; only chains needing
a slot that does not exist stay interpreter-only.

Combining two sources is a policy-level concern, not a chain-shape gap:
monotonicity correctly forbids a second `Sample` crossing, and source
policies already receive the full carrier, so `Source::Multiply<A, B>` /
`Source::Max<A, B>` combinator policies compose fields with zero
pipeline change at the template layer. The interpreter ships compound
sources as their own curated operator entries — a scalar
source-expression ABI (nested operators, each with its own parameter and
prepared blocks) is deliberately out of this revision, and joins true
Fork/Join branching as a reserved future extension, taken up only if the
curated set grows unwieldy.

The integration surface for any new operator is one `OperatorDescriptor`
record (§7.1), from which the promotion, interpreter, and tool views are
all generated — no schema change, no new banks.
