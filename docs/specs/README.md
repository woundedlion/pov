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
| [pullback_stage_families_spec.md](pullback_stage_families_spec.md) | LANDED §§1–6 (static ranked pipeline and its migration); §7 promotion/verification and §8 preview interpreter remain PROPOSED. Relaxes the six-slot stage-kind pipeline to arbitrary chains of stages over ranked carriers (Sphere → Plane → Field → Color), including the preview interpreter's engine contract. Where it and the pipeline spec disagree on the stage model, this spec describes what ships. |
| [shader_workbench_chain_spec.md](shader_workbench_chain_spec.md) | PROPOSED. The tool half of the stage-families spec: shader document schema v2 (ordered chain array, digest-bearing labels), v1 migration and digest identity, and the chain editor. Lands after the stage-families C++ cut-over. |

## Phantasm hardware

| Document | Status and scope |
|---|---|
| [phantasm_pcb_spec.md](phantasm_pcb_spec.md) | Source of truth for the KiCad schematic and layout of the per-segment carrier board (`hardware/phantasm/`). One identical PCB ×4 is qualified; the N=8 profile is compile-tested only. |
| [phantasm_frame_sync_spec.md](phantasm_frame_sync_spec.md) | IMPLEMENTED, and describes the implementation. One-wire flywheel sync; protocol core `hardware/pov_sync.h`, device shell `hardware/pov_segmented.h`. |

Related indexes: [on-device effect profiles](../profiles/README.md) and the
ledgers under [`docs/ledgers/`](../ledgers).
