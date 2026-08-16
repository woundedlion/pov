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

Every continuous parameter uses binary32 storage and declares its unit, domain,
and interpolation trait. Preset records provide one value for every parameter.
Transition edges name a descriptor-owned path policy and a bank-owned easing
and positive duration.
