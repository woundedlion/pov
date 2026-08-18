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

## Phantasm hardware

| Document | Status and scope |
|---|---|
| [phantasm_pcb_spec.md](phantasm_pcb_spec.md) | Source of truth for the KiCad schematic and layout of the per-segment carrier board (`hardware/phantasm/`). One identical PCB ×4 is qualified; the N=8 profile is compile-tested only. |
| [phantasm_frame_sync_spec.md](phantasm_frame_sync_spec.md) | IMPLEMENTED, and describes the implementation. One-wire flywheel sync; protocol core `hardware/pov_sync.h`, device shell `hardware/pov_segmented.h`. |

Related indexes: [on-device effect profiles](../profiles/README.md) and the
ledgers under [`docs/ledgers/`](../ledgers).
