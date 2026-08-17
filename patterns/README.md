# Shader documents

Shader studies use `*.shader.json` as their editable source. The authoring
backend validates the version 1 linear pullback graph, derives a canonical
effect semantic descriptor, and computes its SHA-256 digest without including
preset values or display metadata.

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

Version 1 graphs contain the six compiled roles in order: `outer_camera`,
`surface_project`, `planar_warp`, `source`, `material`, and `color`. Source
labels and declaration order do not enter semantic identity. Stage policies,
parameter schemas, clocks, preparation, resources, serialization,
approximation, and handoff do.

Scalar parameters use binary32 storage. Discrete color mappings use enum8
storage with `MIXED_ENUM` interpolation, which carries both mapping endpoints
and a blend weight through the shared color stage. Preset records provide one
value for every parameter. Transition edges name a descriptor-owned path policy
and a bank-owned easing and positive duration.

`curl_lattice.shader.json` is the editable source for the off-roster
`CurlLattice` comparison effect. Its two presets share one descriptor and vary
only the linearly interpolated sphere-noise scale.

`facet_grid.shader.json` describes the fixed stereographic dodecahedral-grid
pipeline used by `FacetGrid`. Its four presets differ only in source,
projection, warp, and color parameters.

`example.shader.json` carries no `effect_id`: it is the CLI's sample document
and backs no effect.

Every other document is the editable source of a `FixedLook::Runtime`
specialization. A document maps to its effect by `effect_id` == the effect's
`EFFECT_ID`; directory placement predicts neither the base class nor the number
of effects in a header, so the table below is the map. `CurlLattice` and
`FacetGrid` are the two bespoke effects; `AlienOcean` shares the runtime with
the `effects/fixed/` group but keeps a single-effect file.

| Document | Effect | Header |
| --- | --- | --- |
| `alien_ocean` | `AlienOcean` | `effects/AlienOcean.h` |
| `curl_lattice` | `CurlLattice` | `effects/CurlLattice.h` |
| `facet_grid` | `FacetGrid` | `effects/FacetGrid.h` |
| `cosmic_eyeball`, `equator_grid`, `glitch_grid` | `CosmicEyeball`, `EquatorGrid`, `GlitchGrid` | `effects/fixed/GridMirrorLooks.h` |
| `facet_wave`, `signal_weave`, `vector_facets` | `FacetWave`, `SignalWeave`, `VectorFacets` | `effects/fixed/GridWarpLooks.h` |
| `contour_lattice`, `prism_lattice` | `ContourLattice`, `PrismLattice` | `effects/fixed/LatticeLooks.h` |
| `hex_wave`, `kaleido_wave`, `mobius_grid` | `HexWave`, `KaleidoWave`, `MobiusGrid` | `effects/fixed/TwinWaveMirrorLooks.h` |
