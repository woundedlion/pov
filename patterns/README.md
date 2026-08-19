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

Scalar parameters use binary32 storage. Discrete choices — color mappings and
the topology fields v1 expansion bakes out — use enum8 storage with
`MIXED_ENUM` interpolation, which carries both endpoints and a blend weight
through a transition. Preset records provide one
value for every parameter. Transition edges name a descriptor-owned path policy
and a bank-owned easing and positive duration.

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
