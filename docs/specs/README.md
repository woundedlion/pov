# Design specs

Every spec carries its own status banner; the lines below repeat it.

## Phantasm hardware

| Document | Status and scope |
|---|---|
| [phantasm_pcb_spec.md](phantasm_pcb_spec.md) | SPECIFIED; the body describes committed rev 1.1 artifacts; the committed rev 1.2 project requires placement and routing. The routed board is committed, with two recorded deviations — the §11.1 hand-solder lands and R-PWR-7's J1 keying are unmet by the shipped copper, each carrying a deviation block. Source of truth for the KiCad schematic and layout of the per-segment carrier board (`hardware/phantasm/`). One identical PCB ×4 is qualified; the N=8 profile is compile-tested only. |
| [phantasm_frame_sync_spec.md](phantasm_frame_sync_spec.md) | IMPLEMENTED, with section 11.3 hardware measurement of the worst-case interrupt-mask window still open. One-wire flywheel sync; protocol core `hardware/pov_sync.h`, device shell `hardware/pov_segmented.h`. |

Related indexes: [on-device effect profiles](../profiles/README.md) and the
ledgers under [`docs/ledgers/`](../ledgers).
