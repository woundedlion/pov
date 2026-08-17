# Design specs

Every spec carries its own status banner; the lines below repeat it. Where two
specs cover the same ground, the banner in the newer one says which half it
owns — this index records those relationships so a reader lands on the
authoritative document first.

## Rendering architecture

The shader-family specs stack, narrowest scope last:

| Document | Status and scope |
|---|---|
| [pullback_pipeline_spec.md](pullback_pipeline_spec.md) | LANDED. The composition core, standard carriers, provider concepts, stage combinators, and operator catalog that ship in `core/render/pullback.h`. Field-level source of truth for the carrier records every other shader spec names. §17 is a design record with no implementation. |
| [inverse_sampling_pipeline_spec.md](inverse_sampling_pipeline_spec.md) | SHIPPING. The sole shipping renderer for ShaderBall-style pullback rendering: typed program selection, code-emission and ELF rules, and the device resource ratchets. The WASM images keep the dynamic backend deliberately. |
| [shaderball_spec.md](shaderball_spec.md) | LANDED. §0 is authoritative for ShaderBall's authored vocabulary, presets, and choreography, and for nothing beyond them — the renderer architecture is the inverse-sampling spec's. §§1–13 are the original merge record and are historical wherever they disagree with §0 or the code. |
| [shader_workbench_fixed_pipeline_effects_spec.md](shader_workbench_fixed_pipeline_effects_spec.md) | LANDED, except §10.2–§10.3 (code generation and generated-file layout), which are specified and deferred. Defines the authoring workbench and the fixed-pipeline product effects; formulas and rendering behavior stay specified by the three specs above. |

## Phantasm hardware

| Document | Status and scope |
|---|---|
| [phantasm_pcb_spec.md](phantasm_pcb_spec.md) | Source of truth for the KiCad schematic and layout of the per-segment carrier board (`hardware/phantasm/`). One identical PCB ×4 is qualified; the N=8 profile is compile-tested only. |
| [phantasm_frame_sync_spec.md](phantasm_frame_sync_spec.md) | IMPLEMENTED, and describes the implementation. One-wire flywheel sync; protocol core `hardware/pov_sync.h`, device shell `hardware/pov_segmented.h`. |

Related indexes: [on-device effect profiles](../profiles/README.md) and the
ledgers under [`docs/ledgers/`](../ledgers).
