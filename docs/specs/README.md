# Design specs

Every spec carries its own status banner; the lines below repeat it. Where two
specs cover the same ground, the banner in the newer one says which half it
owns — this index records those relationships so a reader lands on the
authoritative document first.

## Rendering architecture

The shader-family spec:

| Document | Status and scope |
|---|---|
| [pullback_pipeline_spec.md](pullback_pipeline_spec.md) | LANDED. The composition core, standard carriers, provider concepts, stage combinators, and operator catalog that ship in `core/render/pullback.h`. Field-level source of truth for the carrier records the shader family names. §17 is a design record with no implementation. |
| [pullback_stage_families_spec.md](pullback_stage_families_spec.md) | LANDED §§1–6 (static ranked pipeline and its migration) and §8 (preview interpreter); §7 promotion/verification remains PROPOSED. Relaxes the six-slot stage-kind pipeline to arbitrary chains of stages over ranked carriers (Sphere → Plane → Field → Color), including the preview interpreter's engine contract. Where it and the pipeline spec disagree on the stage model, this spec describes what ships. |
| [shader_workbench_chain_spec.md](shader_workbench_chain_spec.md) | LANDED §§1–4. The tool half of the stage-families spec: shader document schema v2 (ordered chain array, digest-bearing labels), v1 migration and digest identity, the chain editor, and the pipeline-strip workbench surface. Ships in the daydream repo. |

The rasterizer and driver specs:

| Document | Status and scope |
|---|---|
| [segmented_stateful_effects_spec.md](segmented_stateful_effects_spec.md) | IMPLEMENTED. Cross-segment reach for history-reading effects: the compile-time filter traits, the `Effect::needs_full_frame()` query, and the two driver boundaries (`targets/wasm/engine_bindings.h` `setClip`, `hardware/pov_segmented.h` `clip_to_segment`) that honour it. |
| [congruence_class_lut_spec.md](congruence_class_lut_spec.md) | FACILITY ONLY. Congruence-class clustering and canonical distance LUTs (`core/mesh/mesh_classes.h`); landed and gate-green but wired to no effect. §11–§12 carry the measurements behind the deformation restriction and why IslamicStars was unwired. |

## Mesh morphing

| Document | Status and scope |
|---|---|
| [opchain_morph_spec.md](opchain_morph_spec.md) | LANDED. The shipped op-by-op build contracts: recipe model, `Animation::OpLeg` leg kinds, truncate edge cases, the smooth kis/needle bridge, the renderer constraints the design turns on, and the measured dead ends. Field-level source of truth is `effects/IslamicStars.h`, `core/animation/opleg.h`, `core/mesh/recipe.h` and `core/mesh/solids.h`. |
| [conway_morph_spec.md](conway_morph_spec.md) | SUPERSEDED by opchain_morph_spec.md for the design of record; retained for the §3 edge table and the §7 test plan the Conway-morph suites are numbered against. Symbol names and tuning constants in it no longer track the tree. |

## Phantasm hardware

| Document | Status and scope |
|---|---|
| [phantasm_pcb_spec.md](phantasm_pcb_spec.md) | SPECIFIED; the routed board is committed, with two recorded deviations — the §11.1 hand-solder lands and R-PWR-7's J1 keying are unmet by the shipped copper, each carrying a deviation block. Source of truth for the KiCad schematic and layout of the per-segment carrier board (`hardware/phantasm/`). One identical PCB ×4 is qualified; the N=8 profile is compile-tested only. |
| [phantasm_frame_sync_spec.md](phantasm_frame_sync_spec.md) | IMPLEMENTED, and describes the implementation. One-wire flywheel sync; protocol core `hardware/pov_sync.h`, device shell `hardware/pov_segmented.h`. |

Related indexes: [on-device effect profiles](../profiles/README.md) and the
ledgers under [`docs/ledgers/`](../ledgers).
