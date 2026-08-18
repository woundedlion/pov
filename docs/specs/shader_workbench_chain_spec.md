# Shader workbench chains: document schema v2 and the chain editor

**Status: PROPOSED.** The tool half of
[pullback_stage_families_spec.md](pullback_stage_families_spec.md): the
document schema, migration, and editor for authoring the chains that spec
defines. Lands after that spec's C++ cut-over; nothing here is
implemented. The workbench code lives in the daydream repo
([shader/shader_workbench.mjs](https://github.com/woundedlion/daydream/blob/master/shader/shader_workbench.mjs),
[tools/shader_documents.js](https://github.com/woundedlion/daydream/blob/master/tools/shader_documents.js)).

## 1. Document schema v2

**Encoding: an ordered `chain` array, not nodes + edges.** v1's
nodes/edges graph only validates because roles are unique; its
canonicalizer sorts nodes by a role index, strips labels, and rewrites
edges as role pairs — all three moves collapse under duplicate operators
(two warps in either order canonicalize identically, a digest collision
between semantically distinct documents, while a shuffled identical
document digests differently). A guaranteed-linear chain says so:

- `descriptor.chain` is an ordered array of `{label, operator}` —
  exactly the engine wire format's `{instance_id, operator_id}`, so
  the chain projection of a document *is* its compiled program shape,
  and the two cannot structurally drift; the full document digest
  refines that shape with parameter schemas, presets, and path
  policies (the main spec §8 names the three identities — digest,
  program shape, instance state). A chain entry carries nothing the
  engine would drop.
  v1's per-node `policy` and `resources` fields do not survive:
  structural variation is expressed either as a distinct operator id or
  as an enum8 parameter declared in the operator's catalog schema and
  carried as an ordinary parameter value (the interpreter treats these
  as runtime switches, as the dynamic backend already does; promotion
  maps them to template arguments), and an operator's resource needs
  live in the engine's operator table, not the document. The v1
  importer consumes the old `policy` objects into operator selection
  plus enum8 values. Array order is both the semantic and the
  **canonical** order; there is no `edges` field and no role sort.
  Graph generality for a hypothetical DAG future is not paid for
  silently — if Fork/Join ever lands it arrives as a versioned schema
  change.
- **Labels are canonical and digest-bearing**: they are the instance
  identity parameters bind to, so renaming an instance is a descriptor
  change.
- **Operator catalog v2**: each operator declares its carrier pair
  (input/output carrier, from which family rank and all editor legality
  derive) and its parameter schema. The catalog is the single source the
  validator and the editor read, and it is generated from — or
  golden-pinned against — the C++ combinator typedefs, so the tool
  cannot drift from the engine's ground truth. Compound sources
  (`Source::Multiply<A, B>`) appear as their own curated catalog
  entries; operator nesting is reserved for a future schema revision.
- **Bindings become real.** v1's `parameter.binding` is validated as a
  bare identifier and cross-checked against nothing (shipped documents
  contradict any label convention in both directions). v2: a parameter id
  is `<label>.<field>`; the label segment must resolve to a chain entry
  and the field to that operator's catalog schema, else
  `UNBOUND_PARAMETER`. The existing `ID_PATTERN` admits these ids; v2
  fixes `.` as exclusively the label/field separator (`-` only within a
  segment).
- **Diagnostics keep v1's discipline** — phase/code/path on every
  failure. New semantic codes, at minimum: `ENTRY_FAMILY`, `EXIT_FAMILY`,
  `FAMILY_ORDER`, `UNKNOWN_OPERATOR`, `DUPLICATE_LABEL`,
  `UNBOUND_PARAMETER`. The import path surfaces the full diagnostic
  list, not the first line — with a valid-by-construction editor, import
  is where diagnostics are actually consumed.
- Digest-stability tests include a duplicate-operator chain (reordered
  document → same digest; swapped stage order → different digest).

## 2. v1 → v2 identity

Expansion rewrites every descriptor digest (roles → stages, identities
dropped) and every parameter id (label namespacing), which cascades into
preset values, staggered path-policy groups (v1's `STAGGERED_ORDERED`
transition scheduling), and `serialization.fields`. Requirements:

- Expansion is the **single code path**: the re-exported v2 catalog
  documents are by definition its output, gated by a test that pins them
  byte-identical after canonicalization.
- A frozen v1-digest → v2-digest migration table serves the registry and
  the workbench's fixed-vs-dynamic preview routing, both of which match
  by digest today.
- Deterministic label assignment for expanded instances (v1 slot order:
  `warp1`, `warp2`, …) and a complete parameter-id rewrite map.
- v1 documents that expand to the same chain (distinct only by
  identity-policy spelling) collide by design — expansion canonicalizes
  them — and the registry migration merges their entries deliberately
  rather than reporting ambiguity.

**Engine control names.** The hand-written alias table
(`engineParameterNames`) cannot be deleted outright: fixed-preview
routing writes into compiled composed effects whose slider names come
from `Params` registration. Convergence rule: newly promoted effects
register label-derived control names, so document ids and engine names
coincide by construction; the alias table survives only for effects
promoted before this spec, and shrinks as they are re-registered.

## 3. The chain editor

The thirteen fixed folder banks are replaced by direct manipulation of
the chain: a vertical rail of stage cards grouped into four family bands
with crossings drawn as band boundaries, the default scratch chain being
`Rotate → Project → Sample → Colorize` (the minimal legal chain over the
shipped crossings is `Project → Sample → Colorize`; the default adds the
camera), and
selecting a card showing that stage's parameters labeled by instance
("Wave Shear · warp2"). An operator's parameter list is its full union
schema (main spec §7.1): fields the current topology value deactivates
(edge-fade width while coverage is Weight) may render dimmed, but they
stay in the document, every preset, and the digest — deactivation
changes what the engine reads, never what the document carries.

**One edit rule generates every operation.** A structural edit replaces a
contiguous span of the chain with a catalog sequence whose endpoint
carriers match the span's neighbors:

- *Insertion* replaces an empty span (the palette at any point lists
  exactly the operators that fit — legality is computed from catalog
  carrier pairs **and the engine's exported budgets**, accounting
  aligned block totals against arena capacity and chain length, never
  checked after the fact, so the editor has no invalid-document state
  and no over-budget refusal at apply).
- *Removal* replaces with the empty sequence — automatically legal for an
  endomorphism, automatically illegal across a crossing, which is why
  crossings are removed by replacement (e.g. Project + Sample → a
  sphere-domain source) rather than deletion.
- *m-for-n replacement* is the general case, and makes skip-collapse and
  skip-expansion symmetric: without it, plain insertion could never add a
  crossing back, and the rail could author chains it cannot edit out of.
- *Bypass* is a session-only skip of one endomorphism — an ephemeral
  **program-shape override** in the main spec §8's identity terms: the
  engine compiles a shape omitting the entry while the document and
  its digest are untouched, which is why it is never serialized and
  never part of the canonical descriptor — A/B toggling cannot change
  document identity. Crossings get no toggle.

**Structural edits reconcile the whole document atomically.** The
validator requires every preset to carry every parameter and staggered
path policies to schedule every group, so an unreconciled edit
invalidates the entire preset bank. On insert: backfill every preset with
the operator's catalog defaults. On remove: drop the instance's values
from presets, staggered groups, and `serialization.fields`, and drop
transition edges that become degenerate. On relabel: rewrite ids in
place. The editor never produces a document that fails its own
validator, and each structural edit plus its reconciliation is one atomic
undoable step — document-level snapshot history suffices (documents are
≤1 MB) and is required, because removal-with-replacement is otherwise a
one-misclick preset-bank data-loss machine.

**Discoverability and access.** A browsable catalog panel (family-grouped,
illegal entries visible but disabled with the reason) accompanies
contextual insertion — the old banks were ugly but advertised the whole
vocabulary, and contextual palettes alone would hide it. Keyboard and AT
parity is required: cards support move-up/move-down, the insertion
palette is a listbox, the selected card carries `aria-current`, and focus
restores to the edited card after the parameter-GUI rebuild. The existing
nav module's build-once observer assumption does not survive a mutable
rail; it is replaced, not adapted. A read-only breadcrumb of the chain
may sit above the canvas; as a second interactive navigation surface it
is deferred until long chains demonstrate the need.

**Engine boundary.** Applying a document follows the interpreter's engine
contract (main spec §8): `setShaderChain` → apply preset values by id →
`syncEffectGui` → `invalidate`, with the synchronous
definition-rebuild-and-generation-bump guarantee making the value
application safe. Presets and transitions are otherwise unchanged: the
machinery is parameter-id-driven and survives — given the reconciliation
rules above, which are what keep id-driven machinery coherent under id
churn.
