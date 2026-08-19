# Shader workbench chains: document schema v2 and the chain editor

**Status: §§1–4 LANDED 2026-08-19.** The tool half of
[pullback_stage_families_spec.md](pullback_stage_families_spec.md): the
document schema, migration, and editor for authoring the chains that spec
defines. All four sections ship in the daydream repo
([shader/shader_workbench.mjs](https://github.com/woundedlion/daydream/blob/master/shader/shader_workbench.mjs),
[tools/chain_document_store.js](https://github.com/woundedlion/daydream/blob/master/tools/chain_document_store.js),
[tools/chain_strip.js](https://github.com/woundedlion/daydream/blob/master/tools/chain_strip.js),
[tools/chain_library.js](https://github.com/woundedlion/daydream/blob/master/tools/chain_library.js),
[tools/chain_dock.js](https://github.com/woundedlion/daydream/blob/master/tools/chain_dock.js)).
§4's pipeline-strip workbench is the editor's *surface*; every editing
semantic in §3 carries forward beneath it.

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

**§4 supersedes this section's layout** (the vertical rail, the side
catalog panel, the sidebar document controls). Everything else — the one
edit rule, legality-before-gesture, atomic reconciliation, bypass, the
engine boundary — is the semantic contract §4 builds on and does not
restate.

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

## 4. The workbench surface: pipeline strip, live canvas, stage library

**Status: LANDED 2026-08-19.** A view-layer replacement of §3's rail. The document
store, the schema, digesting, migration, and the apply path are untouched;
what changes is where the chain, the vocabulary, and the render live on
the screen, and which gestures name the store's spans. §3's rail put the
program in a sidebar and gave the render the leftover space — backwards
for a tool whose entire feedback loop is *watching the render while
changing the program*. The redesign inverts it.

### 4.1 Layout

Three stacked regions plus one dock, replacing the sidebar:

- **Toolbar** (top edge, one slim row): document source picker, preset
  picker, Open…/Save/Save As, the descriptor digest (abbreviated,
  click-to-copy), and the document status output. The engine memory/
  compute stats stay in this row.
- **Pipeline strip** (below the toolbar): the loaded chain, left to
  right in execution order, grouped into the four carrier domains.
- **Canvas** (center): the real engine rendering the authored chain.
  The largest region by construction — it owns all vertical space the
  strip, library, and toolbar don't need. Perf/segment overlays keep
  their corners.
- **Stage library** (bottom): the whole operator catalog, grouped into
  the same four domains in the same order, drag sources for the strip.
- **Parameter dock** (right, collapsible): the selected stage instance's
  sliders. Docked beside the canvas rather than over it, because tuning a
  slider is exactly when the render must stay visible.

Each carrier domain has a fixed hue used by both its strip band and its
library group — the color is the visual statement of "stages from this
group go up there". On narrow viewports the strip and library scroll
horizontally and independently; the dock overlays the canvas's right edge;
the library collapses to domain tabs showing one group at a time.

### 4.2 The pipeline strip

One band per carrier — Sphere, Plane, Field, Color — always all four
rendered, left to right in catalog carrier order, so the strip reads as
the pipeline even when a band is empty.

- **Endomorphism stages** render as chips inside their carrier's band, in
  chain order. A band holds any number of chips in sequence — the strip
  is the chain, not a slot-per-domain form.
- **Crossings** render as **socket chips** spanning the boundary between
  their input and output bands, labeled with the carrier pair. The three
  boundaries are always occupied in a valid document (a chain enters on
  sphere and exits on color), which the strip makes structural: sockets
  are the fixed joints of the pipeline, bands are the variable runs
  between them.
- **Chip anatomy**: operator display name, instance label ("Wave Shear ·
  warp2"), and on the chip: a **✕ remove** icon and a **bypass** toggle
  (endomorphisms only), a **swap** icon (crossings only). The selected
  chip is outlined and carries `aria-current`; a bypassed chip renders
  dimmed.
- **Remove**: ✕ commits `replaceSpan(i, 1, [])` — legal by construction
  for an endomorphism, which is why only endomorphisms carry it.
  Crossings are removed by replacement (§3), so sockets carry swap
  instead: the replacement palette for that span, same-pair operators
  enabled, everything else present-but-disabled with the reason.
- **Insertion**: gaps between chips are the store's chain indices. Each
  band carries one persistent **+** affordance opening the insertion
  palette for that band's context gap (after the band's selected chip,
  else the band's last gap); during a drag, every legal gap materializes
  as a highlighted target. Both routes go through `legalInsertions` —
  the strip has no legality of its own.
- **Reorder**: chips drag within their band (a crossing doesn't reorder,
  as in §3's drag rule); Alt+Arrow is the keyboard equivalent. The move
  commits as the label-preserving m-for-m span replacement, so parameter
  values survive reorder.
- The strip rebuilds whole after every committed edit with keyboard focus
  restored to the edited chip — §3's render discipline, rotated.

### 4.3 The stage library

The catalog panel's discoverability requirement survives: the whole
vocabulary stays visible, grouped by the carrier each operator consumes,
with crossings showing their pair. What changes is the geometry and the
drop model:

- **Drag a library chip onto the strip to insert it.** On drag start the
  strip highlights every gap the store accepts the operator at — carrier
  legality *and* the exported budgets (arena bytes, op count, param
  count), computed before the gesture, so there is no over-budget refusal
  at drop.
- **Bands are coarse drop targets.** Dropping on a highlighted gap
  commits exactly there; dropping anywhere else on a band commits at the
  band's nearest legal gap. A band with no legal gap for the dragged
  operator shows a refusing state while hovered, with the reason in the
  shared status region. Coarse drops are what make the strip feel like
  "drag a stage onto a domain" instead of "hit a 6-pixel seam".
- **Click** an enabled entry to insert without a drag, at the context
  gap (after the selection, else the first legal gap) — §3's behavior.
- Entries that currently fit nowhere render disabled with the computed
  reason, never hidden. A text filter narrows all groups by name/id;
  filtering never hides an entry's disabled state.

### 4.4 Parameters

Committing an insert **selects the new instance and opens the parameter
dock on it** — adding a stage and immediately hearing its sliders is the
core authoring loop. Selecting any chip retargets the dock.

- Binary32 fields render as sliders with the catalog schema's min/max/
  curve; topology enum8s render as dropdowns (live structural switches on
  the chain path). Controls are labeled by instance.
- The union-schema discipline (§3) applies: fields the current topology
  values deactivate render dimmed, never dropped.
- **A dock edit is a document edit**: it writes the active preset's value
  for that `<label>.<field>` id, with the engine write as the side
  effect (per-control coalesced for undo, in the same history as
  structural edits). The document stays the single source of truth; Save
  never has to "capture" state that only the engine knows.

### 4.5 Use case: authoring a new effect

The scratch document opens as the default chain (`Rotate → Project →
Sample → Colorize`, §3) with one preset of catalog defaults — a valid,
rendering document from the first frame. Authoring is then: drag stages
in, tune the dock, iterate against the live render, name presets.

**Export** is the workflow's terminal step and has one format: **Save
produces the canonical v2 serialization** (`exportShaderDocumentJson`) —
download or re-save to the opened file handle. The toolbar digest is the
document's identity across the export boundary. Because the editor is
valid-by-construction (§3), any exported document validates cleanly and
sits within the engine's exported budgets; there is no "export failed"
state. The exported document is also the promotion input: main spec §7
consumes exactly this file to produce a compiled composed effect, and
shipping it as a pattern is the ordinary registry install of the saved
file. The workbench itself neither promotes nor installs — it guarantees
the artifact both consumers trust.

### 4.6 Use case: modifying an existing effect

Loading must populate the configurator **exactly as if the effect had
been authored from scratch in this session** — same strip, same dock,
same edit affordances, no separate read-only mode:

- The toolbar source picker lists the registry's chain documents; Open…
  imports a file. A v1 document expands through the single expansion path
  (§2) on import; import failures surface the full diagnostic list (§1).
- Authoring always routes the preview through the interpreter
  (`setShaderChain`) so the loaded chain is live-editable. When the
  loaded descriptor digest matches a promoted fixed effect (via the
  migration/registry table), the toolbar offers a **parity toggle** to
  the compiled build for A/B verification; the first descriptor-changing
  edit breaks the match and disables the toggle — bypass, being a
  program-shape override that never touches the digest (§3), does not.
- **Save As** writes a new `document_id` and leaves the loaded original
  untouched; plain Save re-exports over the opened handle. Either way the
  output is the same canonical serialization as §4.5 — modifying an
  effect and creating one converge on the same export.

### 4.7 Keyboard and AT

Full parity with the pointer gestures, rotated to the horizontal:

- The strip is a `listbox` with `aria-orientation="horizontal"`;
  Left/Right roam chips (roving tabindex), Alt+Left/Right move the
  focused endomorphism, Enter/Space selects, Delete removes an
  endomorphism or opens a crossing's swap palette, Insert opens the
  insertion palette at the following gap. Band `+` buttons and palettes
  keep §3's listbox behavior; focus restores to the edited chip after
  every rebuild.
- Library entries are buttons; disabled entries carry `aria-disabled`
  plus a described-by reason. Domain groups are labeled `group`s; the
  filter is a labeled searchbox.
- One shared live status region announces every refusal — strip, library,
  and dock report through it.

### 4.8 Module boundaries

The rail view in `chain_editor.js` and `chain_catalog_panel.js` are
replaced by a strip module and a library module implementing the same
store-facing contract (the `ChainStore` surface, the drag controller
owned by the strip, gap hit-testing generalized with band hit-testing for
coarse drops). `deactivatedParamNames`, the document store, `chain_apply`,
and the schema module are untouched. `shader.html` reflows to the §4.1
regions; the workbench stylesheet is rewritten rather than adapted — the
rail's vertical assumptions are load-bearing throughout it. The existing
editor tests port to the strip contract; new coverage: coarse-drop
snapping to the nearest legal gap, socket swap, ✕ legality (absent on
crossings), dock edits writing the active preset and participating in
undo, and the parity toggle disarming on the first descriptor edit.
