# PHANTASM board projects

Use **[1.2/](1.2/)** for the current Quilter placement run.
Upload these three files together from that directory:

- [phantasm.kicad_pcb](1.2/phantasm.kicad_pcb): unplaced rev 1.2 board.
- [phantasm.kicad_sch](1.2/phantasm.kicad_sch): matching rev 1.2 schematic.
- [phantasm.kicad_pro](1.2/phantasm.kicad_pro): placement and design rules.

Keep the accompanying symbol and footprint libraries with the project when
opening it in KiCad 10. Every upload must use the board and schematic from the
same revision directory. Rev 1.2 requires placement, routing, and validation
before fabrication.

| Directory | Contents |
|---|---|
| `1.2/` | Current complete rev 1.2 project for Quilter |
| [1.1/](1.1/README.md) | Routed rev 1.1 project, libraries, and detailed technical reference |
| `gen/` | Generators and validation tools |

Regenerate the current project from the repository root, schematic first:

```sh
python hardware/phantasm/gen/board.py --force
python hardware/phantasm/gen/pcb.py --unplaced --force --force-teensy-library
```

These commands replace the generated rev 1.2 files in `1.2/`.
Before uploading a project that was opened in KiCad, restore its rule floors:

```sh
python hardware/phantasm/gen/heal_clearance.py hardware/phantasm/1.2/phantasm.kicad_pro
```

In Quilter, preserve the uploaded four-layer stackup and select **Preserve copper
on internal layers**. The board's schematic links identify related components.

The stackup is **Signal / GND / GND / Signal**. Both inner layers have GND
pours and the user names `GND` and `Ground`, matching
[Quilter's ground-layer naming rules](https://docs.quilter.ai/using-quilter/prepare-your-input-board-file).
Confirm both inner layers are classified as Ground in Quilter's stackup editor.
Its no-power-layer info is expected: this design uses two ground reference
planes, with the supply rails routed on the outer layers.

Electrical and mechanical requirements live in the
[PCB specification](../../docs/specs/phantasm_pcb_spec.md).
