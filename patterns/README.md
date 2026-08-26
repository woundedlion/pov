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

Print the canonical descriptor:

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
`<label>.<field>` against the catalog's operator schema. A `schema_version` 1
six-role graph document expands through `expandV1Document` before validation —
the single code path both generations share.

The descriptor is a semantic identity over the chain, not a complete
reproduction of the effect. Expansion keeps only `chain`, `parameters`,
`path_policies` and `serialization`; a v1 document's `resources`, `clocks`,
`preparation` and `approximation` sections do not survive it, so noise seeds,
runtime clocks, prepared frames and approximation oracles live only in the C++
effect.

Scalar parameters use binary32 storage. Discrete choices — color mappings and
the topology fields v1 expansion bakes out — use enum8 storage with
`MIXED_ENUM` interpolation, which carries both endpoints and a blend weight
through a transition. Preset records provide one
value for every parameter. Transition edges name a descriptor-owned path policy
and a bank-owned easing and positive duration.

## Generated and hand-authored sources

Twelve of the eighteen documents are generated. `node
scripts/generate_promoted_shader_documents.mjs` rewrites `alien_ocean`,
`grid_space`, `cosmic_eyeball`, `kaleidoscope_flowers`, `kaleidoscope_mandala`,
`alien_core`, `kaleidoscope_hex_bright`, `kaleidoscope_hex_soft`, `mobius_grid`, `kaleidoscope_pent_bright`,
`alien_brain` and `kaleidoscope_stained_glass` from the effect specs the script holds, so a
hand edit to those files is lost on the next run — change the spec instead. The
specs are written in the v1 six-role shape and the committed file is their
canonical v2 expansion.

The remaining six are hand-authored and no rerun writes them: `ash_cloud`,
`lattice_melt`, `chromatic_lichen`, `kaleidoscope_smooth` and
`kaleidoscope_hex_oil`, promoted from workbench snapshots, and
`example.shader.json`.

`shaderball_migration.json` is a manifest, not a shader document, and the CLI
above rejects it as one. `source_documents` maps each `effect_id` to the
document it is authored in, `product_group` carries the gallery grouping and
its device-second budget, and `destinations` maps every retired legacy
`ShaderBall` preset to the effect and preset that replaced it.
`scripts/shader_workbench.test.mjs` gates the tree against it: a document that
backs an effect must appear in `source_documents`.

A promoted document is also applied to its compiled effect by control name, one
parameter at a time, so every parameter id must resolve to a control that effect
registers. Chain labels are the vocabulary that resolution runs on — `camera`,
`lens`, `surface`, `project`, `warp1`/`warp2`, `sample`, `transfer`, `cutout`,
`colorize`, the labels the v1 expansion assigns — and a hand-authored document
takes the same ones. `scripts/wasm_smoke.mjs` resolves every promoted id against
the running module's registered controls; the topology fields a fixed build
bakes in are exempt, as is `camera.spin-speed`, which `AshCloud` holds as a
compile-time constant.

`lattice_melt.shader.json` is the editable source for the `LatticeMelt`
comparison effect. Its two presets share one descriptor and vary only the
linearly interpolated sphere-noise scale.

`chromatic_lichen.shader.json` is the editable source for the
`ChromaticLichen` effect. It combines a glitch lens, post-lens curl
displacement, and a low-frequency gnomonic grid.

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
