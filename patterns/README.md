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

Twelve of the sixteen documents are generated. `node
scripts/generate_promoted_shader_documents.mjs` rewrites `alien_ocean`,
`contour_lattice`, `cosmic_eyeball`, `equator_grid`, `facet_wave`,
`glitch_grid`, `hex_wave`, `kaleido_wave`, `mobius_grid`, `prism_lattice`,
`signal_weave` and `vector_facets` from the effect specs the script holds, so a
hand edit to those files is lost on the next run — change the spec instead. The
specs are written in the v1 six-role shape and the committed file is their
canonical v2 expansion.

The remaining four are hand-authored and no rerun writes them: `curl_lattice`,
`facet_grid` and `prism_spiral`, promoted from workbench snapshots, and
`example.shader.json`.

`shaderball_migration.json` is a manifest, not a shader document, and the CLI
above rejects it as one. `source_documents` maps each `effect_id` to the
document it is authored in, `product_group` carries the gallery grouping and
its device-second budget, and `destinations` maps every retired legacy
`ShaderBall` preset to the effect and preset that replaced it.
`scripts/shader_workbench.test.mjs` gates the tree against it: a document that
backs an effect must appear in `source_documents`.

`curl_lattice.shader.json` is the editable source for the `CurlLattice`
comparison effect. Its two presets share one descriptor and vary only the
linearly interpolated sphere-noise scale.

`facet_grid.shader.json` describes the fixed stereographic dodecahedral-grid
pipeline used by `FacetGrid`. Its four presets differ only in source,
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
| `contour_lattice` | `ContourLattice` |
| `cosmic_eyeball` | `CosmicEyeball` |
| `curl_lattice` | `CurlLattice` |
| `equator_grid` | `EquatorGrid` |
| `facet_grid` | `FacetGrid` |
| `facet_wave` | `FacetWave` |
| `glitch_grid` | `GlitchGrid` |
| `hex_wave` | `HexWave` |
| `kaleido_wave` | `KaleidoWave` |
| `mobius_grid` | `MobiusGrid` |
| `prism_lattice` | `PrismLattice` |
| `prism_spiral` | `PrismSpiral` |
| `signal_weave` | `SignalWeave` |
| `vector_facets` | `VectorFacets` |
