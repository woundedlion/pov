# Shader documents

Shader studies use `*.shader.json` as their editable source. The authoring
backend validates the version 2 chain document against the operator catalog
(`scripts/engine_catalog.json`), derives a canonical effect semantic
descriptor, and computes its SHA-256 digest without including preset values or
display metadata.

Validate a document:

```text
node scripts/shader_workbench_cli.mjs check patterns/example.shader.json
```

Print the descriptor identity, the digest preimage: the canonical descriptor
with every parameter's `unit` removed, serialized with sorted keys. Its
SHA-256 is the `descriptor_digest` that `check` reports. Committed documents
carry `unit` on every parameter, so this output is not a diff base for a
document's `descriptor` section.

```text
node scripts/shader_workbench_cli.mjs descriptor patterns/example.shader.json
```

Classify it against an effect registry and capability profile:

```text
node scripts/shader_workbench_cli.mjs classify <document> <registry> <profile>
```

Version 2 documents encode the descriptor as an ordered chain of
`{label, operator}` entries validated against the operator catalog. Labels and
chain order are canonical and digest-bearing, structural variation is an
operator id or a topology enum8 parameter, and every parameter id binds
`<label>.<field>` against the catalog's operator schema. Only the current
document schema and current operator identities are accepted.

The descriptor is a semantic identity over the chain, not a complete
reproduction of the effect. It describes `chain`, `parameters`,
`path_policies` and `serialization`; noise seeds,
runtime clocks, prepared frames and approximation oracles live only in the C++
effect.

Scalar parameters use binary32 storage. Discrete choices — color mappings and
topology fields — use enum8 storage with
`MIXED_ENUM` interpolation, which carries both endpoints and a blend weight
through a transition. Preset records provide one
value for every parameter. Transition edges name a descriptor-owned path policy
and a bank-owned easing and positive duration.

## Generated and hand-authored sources

Twelve of the nineteen documents are generated. `node
scripts/generate_promoted_shader_documents.mjs` rewrites `alien_ocean`,
`grid_space`, `cosmic_eyeball`, `kaleidoscope_flowers`, `kaleidoscope_mandala`,
`alien_core`, `kaleidoscope_hex_bright`, `kaleidoscope_hex_soft`, `mobius_grid`, `kaleidoscope_pent_bright`,
`alien_brain` and `kaleidoscope_stained_glass` from the effect specs the script holds, so a
hand edit to those files is lost on the next run — change the spec instead. The specs use current chain IDs and named parameter values.

Seven documents are authored directly in this directory: `ash_cloud`,
`chromatic_lichen`, `mermaid_skin`, `example`, `kaleidoscope_hex_oil`,
`kaleidoscope_smooth` and `lattice_melt`. Edit these version 2 documents here.

After editing any document, run `node scripts/generate_composed_presets.mjs`
to refresh its effect header. Documents use canonical serialization: sorted
keys, two-space indentation and LF endings.
`node scripts/generate_promoted_shader_documents.mjs --check` verifies that
serialization and generated-document freshness.

The engine installs promoted pattern documents (every `source_documents` entry;
`example.shader.json` stays local) and `catalog.json` into
`daydream/generated/shader/patterns/`. The catalog is a manifest:
`source_documents` maps each current effect identity to its authoring document,
and `product_group` carries its gallery grouping. The shader compiler accepts
current documents only. `scripts/shader_workbench.test.mjs` gates membership:
every document backing an effect appears in `source_documents`.

A promoted document is also applied to its compiled effect by control name, one
writable parameter at a time, so every writable parameter id must resolve to a
control that effect registers. Chain labels are the vocabulary that resolution runs on — `camera`,
`lens`, `surface`, `project`, `warp1`/`warp2`, `sample`, `transfer`, `cutout`,
`colorize`. Generated and hand-authored documents use the same labels. `scripts/wasm_smoke.mjs` resolves writable promoted ids against
the running module's registered controls. Fixed topology fields and
`camera.spin-speed`, which `AshCloud` holds as a compile-time constant, are exempt.
Derived values are also exempt when the compiled controls used to derive them
match the document.

`lattice_melt.shader.json` is the editable source for the `LatticeMelt`
effect. Its two presets share one descriptor and vary only the
logarithmically interpolated sphere-noise scale (`LOG_POSITIVE`).

`chromatic_lichen.shader.json` is the editable source for the
`ChromaticLichen` effect. It combines a glitch lens, post-lens curl
displacement, and a low-frequency gnomonic grid.

`mermaid_skin.shader.json` is the editable source for the `MermaidSkin`
effect. It combines a folded grid, maximum palette chroma, and sphere-space
curl displacement.

`ash_cloud.shader.json` is the editable source for the `AshCloud` effect. Its
value cutout follows a primitive lattice displaced by sphere-space curl noise
and folded through a dodecahedral kaleidoscope.

`kaleidoscope_smooth.shader.json` describes the fixed stereographic dodecahedral-grid
pipeline used by `KaleidoscopeSmooth`. Its four presets differ only in source,
projection, warp, and color parameters.

`example.shader.json` carries no `effect_id`: it is the CLI's sample document
and backs no effect.

Every other document is the editable source of a `Pullback::ComposedEffect`
specialization. A document maps to its effect by `effect_id` == the effect's
`EFFECT_ID`. Each effect lives in its own header,
`effects/<ClassName>.h`:

| Document | Effect |
| --- | --- |
| `alien_ocean` | `AlienOcean` |
| `grid_space` | `GridSpace` |
| `cosmic_eyeball` | `CosmicEyeball` |
| `ash_cloud` | `AshCloud` |
| `lattice_melt` | `LatticeMelt` |
| `chromatic_lichen` | `ChromaticLichen` |
| `mermaid_skin` | `MermaidSkin` |
| `kaleidoscope_flowers` | `KaleidoscopeFlowers` |
| `kaleidoscope_smooth` | `KaleidoscopeSmooth` |
| `kaleidoscope_mandala` | `KaleidoscopeMandala` |
| `alien_core` | `AlienCore` |
| `kaleidoscope_hex_bright` | `KaleidoscopeHexBright` |
| `kaleidoscope_hex_soft` | `KaleidoscopeHexSoft` |
| `mobius_grid` | `MobiusGrid` |
| `kaleidoscope_pent_bright` | `KaleidoscopePentBright` |
| `kaleidoscope_hex_oil` | `KaleidoscopeHexOil` |
| `alien_brain` | `AlienBrain` |
| `kaleidoscope_stained_glass` | `KaleidoscopeStainedGlass` |
