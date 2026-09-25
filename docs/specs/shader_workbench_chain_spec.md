# Shader workbench chains: document schema v2 and the chain editor

**Status: §§1–4 LANDED 2026-08-19, except §4.3.** The tool half of
[pullback_stage_families_spec.md](pullback_stage_families_spec.md): the
document schema, migration, and editor for authoring the chains that spec
defines. The schema, validation, canonical identity, and v1 expansion ship
here as `scripts/shader_workbench.mjs` (daydream's `shader/shader_workbench.mjs`
is the copy its engine-bundle installer writes); the document store and the
editor ship in the daydream repo
([tools/chain_document_store.js](https://github.com/woundedlion/daydream/blob/master/tools/chain_document_store.js),
[tools/chain_strip.js](https://github.com/woundedlion/daydream/blob/master/tools/chain_strip.js)).
§4's pipeline-strip workbench is the editor's *surface*; every editing
semantic in §3 carries forward beneath it.

## 1. Document schema v2

**Encoding: an ordered `chain` array, not nodes + edges.** v1's
nodes/edges graph only validates because roles are unique; its
canonicalizer sorts nodes by a role index, strips labels, and rewrites
edges as role pairs — all three moves collapse under duplicate operators
(two warps in either order canonicalize identically, a digest collision
between semantically distinct documents, while a shuffled identical
document digests differently). A guaranteed single-path chain says so:

- `descriptor.chain` is an ordered array of `{label, operator}` —
  exactly the engine wire format's `{instance, operator}` (the two keys
  `setShaderChain` reads; a missing or non-string one is
  `MALFORMED_PAYLOAD`), so
  the chain projection of a document *is* its compiled program shape,
  and the two cannot structurally drift; the descriptor digest
  refines that shape with parameter schemas and defaults, serialization
  fields, and path policies, display units excluded, while the preset
  bank digests separately as `preset_bank_digest` (the main spec §8
  names the three identities — digest, program shape, instance state).
  A chain entry carries nothing the engine would drop.
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
  cannot drift from the engine's ground truth. A compound source, when
  one ships, appears as its own curated catalog entry; operator nesting
  is reserved for a future schema revision.
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

- Expansion is the **single code path** for loading a v1 document. The
  committed v2 pattern documents are engine-owned artifacts, not
  expansion output: each is pinned as its own canonical re-export, and
  the five identity-frame replacements deliberately differ from what
  expanding their v1 fixture yields.
- A v1-digest → v2-digest migration table maps each v1 fixture digest
  onto the digest of the committed document of the same name. It is
  recomputed by
  [scripts/generate-shader-v2-documents.mjs](https://github.com/woundedlion/daydream/blob/master/scripts/generate-shader-v2-documents.mjs),
  which writes only the table, and a completeness test fails when the
  committed table drifts from what a rerun writes. Preview routing does
  not read it: a loaded document matches a promoted fixed effect on its
  v2 digest directly.
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

The semantic contract beneath §4's surface: the one edit rule,
legality-before-gesture, atomic reconciliation, bypass, and the engine
boundary. §4 owns the layout and the gestures and does not restate these.

The chain is edited by direct manipulation over the three editable carrier
domains, with Color represented only by the terminal output crossing. The
default scratch chain is `Rotate → Project → Sample → Colorize` (the minimal
legal chain is `SampleSphere -> Colorize`; the default uses the plane path
and adds the camera). Every instance the editor creates — a scratch
chain's, and each one an edit inserts — declares its operator's full union
schema (main spec §7.1); a loaded document's instances keep the
declarations it carries, which may be narrower. Fields the current
topology value deactivates (edge-fade width while coverage is Weight) may
render dimmed, but a declared field stays in the document, every preset,
and the digest — deactivation changes what the engine reads, never what
the document carries.

**One edit rule generates every operation.** A structural edit replaces a
contiguous span of the chain with a catalog sequence whose endpoint
carriers match the span's neighbors:

- *Insertion* replaces an empty span (the palette at any point lists
  exactly the operators that fit — legality is computed from catalog
  carrier pairs **and the engine's exported budgets**, accounting the
  engine's running-cursor block layout against arena capacity and chain
  length, never checked after the fact, so the editor has no invalid-document state
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

**Engine boundary.** Applying a document follows the interpreter's engine
contract (main spec §8): `setShaderChain` → apply preset values by id →
`syncEffectGui` → `invalidate`, with the synchronous
definition-rebuild-and-generation-bump guarantee making the value
application safe. Presets and transitions are otherwise unchanged: the
machinery is parameter-id-driven and survives — given the reconciliation
rules above, which are what keep id-driven machinery coherent under id
churn.

## 4. The workbench surface: pipeline strip and live canvas

**Status: LANDED 2026-08-19.** Where the chain, the vocabulary, and the
render live on the screen, and which gestures name the store's spans; the
document store, the schema, digesting, migration, and the apply path are
§§1–3's. The render owns the space: the tool's entire feedback loop is
*watching the render while changing the program*. The stage library of
§4.3 is deferred; the catalog reaches the strip through the band insertion
palettes.

### 4.1 Layout

Three stacked regions:

- **Toolbar** (top edge, one slim row): document source picker, preset
  picker, Open…/Save/Save As, the descriptor digest (abbreviated,
  click-to-copy), and the document status output. The engine memory/
  compute stats stay in this row.
- **Pipeline strip** (below the toolbar): the loaded chain, left to
  right in execution order, grouped into the Sphere, Plane, and Field
  carrier domains and ending at the Color output.
- **Canvas** (center): the real engine rendering the authored chain.
  The largest region by construction — it owns all vertical space the
  strip and toolbar don't need. Perf/segment overlays keep their
  corners.

A stage's parameters live **on its own chip** (§4.4), not in a side
panel: the chain and its tuning are one surface. The engine's global
controls float over the canvas region alone — never over the toolbar
or the strip, which are the page's primary UI.

Each editable carrier domain has a fixed hue for its strip band. Bands
size to their contents instead of sharing the strip's width or height.
On narrow viewports the strip scrolls horizontally without exposing a native
scrollbar. Edge arrow buttons, a mouse wheel over the strip, and Left/Right on
the focused strip background move the viewport.

### 4.2 The pipeline strip

The strip has three editable carrier bands — Sphere, Plane, and Field —
rendered left to right in catalog carrier order. Color is the terminal
output type, not a fourth carrier band; the final `Field → Color` socket
names that output directly.

- **Endomorphism stages** render as chips inside their carrier's band, in
  chain order. A band holds any number of chips in sequence — the strip
  is the chain, not a slot-per-domain form.
- **Crossings** render as **socket chips** after their input bands, labeled
  by the function of the carrier the crossing produces — `Projection:` for
  plane, `Source:` for field, `Color:` for color — ahead of the replacement
  selector that names the operator; the carrier pair appears only in the
  chip's accessible name. The Sphere → Plane and Plane → Field sockets sit
  between bands; the Field → Color socket ends the strip. A chain enters on
  sphere and exits on color, which the strip makes structural: sockets are
  the joints of the pipeline, bands are the variable runs between them. A
  crossing may skip a band: the catalog sphere → field sources
  (`sample.spherical-rings.v3`, `sample.spherical-noise.v3`) stand in for
  both the projection and the sample crossing, the plane band is then
  absent with no gap to insert at, and that socket has a replacement
  selector offering the projection/sample pairs that re-open it.
- **Chip anatomy**: operator display name and a row of icon buttons in the
  chip header: a **◉ bypass**
  toggle, a pair of **← →** reorder buttons and a **× remove**
  button (endomorphisms only), or a **replacement selector** (crossings
  only). The instance label appears in the chip's accessible name and rename
  field, not its visible heading. The selected chip is outlined and carries `aria-current`; a
  bypassed chip renders dimmed. Every chip carries its stage's parameter
  controls inline, built from the document's parameter declarations over
  the active preset's values; a chip discloses them transiently while the
  pointer hovers it or keyboard focus is inside it, and pins them open
  while selected.
- **Remove**: × commits `replaceSpan(i, 1, [])` — legal by construction
  for an endomorphism, which is why only endomorphisms carry it.
  Crossings are removed by replacement (§3), so sockets carry a
  selector of the same-pair operators the store accepts for that span;
  Delete opens the same set as a palette.
- **Insertion**: gaps between chips are the store's chain indices. Each
  band carries one persistent **+** affordance opening the insertion
  palette at the gap after the band's last stage. Insert opens it at the
  gap after the focused chip. Both routes go through `legalInsertions` —
  the strip has no legality of its own. A chain with no transfer preserves
  the field value.
- **Reorder**: a chip's ← → buttons move it within its band (a crossing
  doesn't reorder); Alt+Arrow is the keyboard equivalent. The move
  commits as the label-preserving m-for-m span replacement, so parameter
  values survive reorder.
- The strip rebuilds whole after every committed edit with keyboard focus
  restored to the edited chip.

### 4.3 The stage library — deferred

**Not shipped.** No library panel, drop target, or drag controller exists in
the tool, and pointer gestures never start a chip drag. Insertion runs
through §4.2's band palettes.

### 4.4 Parameters

A stage's controls live **inline in its chip**. Parameterized chips start
collapsed; selecting a chip opens its controls, and selecting it again
collapses them. Each numeric value is also directly editable. Committing an
insert selects the new instance, whose controls are
on screen: adding a stage and immediately hearing its sliders is the core
authoring loop, and the chain never has to be read in one place and tuned
in another.

- Controls are built from the **document**, not the engine: the
  declaration's storage and domain give the control its shape, the active
  preset's value gives its position. Binary32 fields render as sliders
  over the declared domain; topology enum8s render as dropdowns (live
  structural switches on the chain path). A row is labeled by its field
  segment alone — the chip already names the instance.
- The declaration discipline (§3) applies: a declared field the current
  topology values deactivate renders dimmed, never dropped.
- An instance that declares no parameters gets no control group, and a
  chip carrying one is wide. Every parameter remains visible within the chip;
  parameter groups never get their own scrollbar. The pipeline viewport
  scrolls horizontally rather than crushing its neighbours.
- **A stage edit is a document edit**: it writes the active preset's
  value for that `<label>.<field>` id, with the engine write as the side
  effect (per-control coalesced for undo, in the same history as
  structural edits). The document stays the single source of truth; Save
  never has to "capture" state that only the engine knows.

### 4.5 Use case: authoring a new effect

The scratch document opens as the default chain (`Rotate → Project →
Sample → Colorize`, §3) with one preset of catalog defaults — a valid,
rendering document from the first frame. Authoring is then: insert
stages from the band palettes, tune them on their chips, iterate against
the live render, name presets.

**Export** is the workflow's terminal step and has one format: **Save
produces the canonical v2 serialization** (`exportShaderDocumentJson`),
downloaded under the document's filename. The toolbar digest is the
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
been authored from scratch in this session** — same strip, same inline
controls, same edit affordances, no separate read-only mode:

- The toolbar source picker lists the registry's chain documents; Open…
  imports a file. A v1 document expands through the single expansion path
  (§2) on import; import failures surface the full diagnostic list (§1).
- Authoring routes the preview through the interpreter (`setShaderChain`) so
  the loaded chain is live-editable, and applies the selected effect preset to
  the document controls and interpreter parameters. When the loaded descriptor
  digest matches a promoted fixed effect's registry entry,
  the toolbar offers a **parity toggle** to the compiled build for A/B
  verification; the first descriptor-changing edit breaks the match and
  disables the toggle — bypass, being a program-shape override that never
  touches the digest (§3), does not.
- **Save As** writes a new `document_id` and leaves the loaded original
  untouched; plain Save re-downloads the loaded document under its own
  filename. Either way the output is the same canonical serialization as
  §4.5 — modifying an effect and creating one converge on the same
  export.

### 4.7 Keyboard and AT

Full parity with the pointer gestures, rotated to the horizontal:

- The strip is a `toolbar` with `aria-orientation="horizontal"` whose
  chips are `group`s, so their inline parameter controls stay exposed;
  the selected chip carries `aria-current`. Left/Right roam chips
  (roving tabindex), Alt+Left/Right move the focused endomorphism,
  Enter/Space selects, Delete removes an endomorphism or opens a
  crossing's swap palette, Insert opens the insertion palette at the
  following gap. Band `+` buttons open the insertion palette, a listbox;
  focus restores to the edited chip after every rebuild. When the
  strip background itself has focus, Left/Right scroll the viewport rather
  than moving a chip.
- Palette entries are `role="option"` rows in a listbox. Operators invalid at
  the active gap are omitted rather than rendered disabled.
- One shared live status region announces every refusal.

### 4.8 Module boundaries

`chain_strip.js` is the view over the store-facing contract of
`chain_document_store.js`, whose `setPresetValue` is the preset-value write
§4.4 needs. `chain_presentation.js` supplies `deactivatedParameterIds`,
which reads document declarations, and shared presentation rules.
`shader_documents.js` owns the document toolbar, Save and parity toggle.
`chain_apply.js` holds the engine boundary and `shader_workbench.mjs` the schema. `shader.html` lays
out the §4.1 regions with `shader.css`. The strip suite
(`chain_strip.test.js`) covers band and socket layout, socket replacement,
× legality (absent on crossings), button and Alt+Arrow reorder, undo and
redo through the apply path, and the inline controls' read of the active
preset and edit callbacks.

**Pointer behaviour needs a real browser.** The unit suite renders into a
DOM fake with neither layout nor pointer capture, which cannot see an
absolutely placed palette landing far from the control that opened it,
nor a pointer travel swallowing the click that should have selected a
chip. A headless-Chrome probe drives the strip's gestures with a real
mouse (palette placement, chip selection, reorder, inline control
round-trip) and gates alongside the page smoke.
