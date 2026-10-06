# PHANTASM board projects

Use **[1.2/](1.2/)** for the current Quilter placement run.
The experimental **[1.3/](1.3/README.md)** project uses differential sync.
Rev 1.2 remains the generator default; select rev 1.3 explicitly.
Upload these three files together from that directory:

- [phantasm.kicad_pcb](1.2/phantasm.kicad_pcb): unplaced rev 1.2 board.
- [phantasm.kicad_sch](1.2/phantasm.kicad_sch): matching rev 1.2 schematic.
- [phantasm.kicad_pro](1.2/phantasm.kicad_pro): placement and design rules.

Keep the accompanying symbol and footprint libraries with the project when
opening it in KiCad 10. Every upload must use the board and schematic from the
same revision directory. Rev 1.2 requires placement, routing, and validation
before fabrication.

The upload board fixes **J1 power at the left edge** and **J2 LED, J3A SYNC IN,
J3B SYNC OUT in a column at the right end**, in that top-to-bottom order.
The column sits inboard of the right mounting holes, in the rev 1.1 connector
area. Their wire openings face left toward the hub. Keepouts reserve each
terminal's full plastic body plus 0.5 mm assembly clearance, preventing SMD
parts beneath the connector. Power's wire-access reservation extends off the
left edge. The outline remains **58.28 x 32 mm**, and all four mounting holes
retain the rev 1.1 coordinates. The Teensy shifts 3.5 mm left to clear the
terminal bodies; D_BUS is released for placement outside their keepouts.

J2/J3A/J3B in revisions 1.2 and 1.3 have no wire-access reservation.
The right-column connectors' wire openings sit about 1.5 mm from the Teensy PCB edge. Before fabrication,
verify harness insertion at the module's mounted height, including any ferrules.
If the harness cannot be inserted with the Teensy seated, land the J2/J3A/J3B
wires before seating it or use a socket tall enough for the harness to clear;
the placement and DRC gates do not check this.

Rev 1.2 omits the J4 debug header and SERIAL1_TX connection; Teensy pin 1 is
unconnected. R_MEN and MASTER_EN remain part of the sync control circuit.
Probe power directly at the Teensy; diagnostics use its USB connector.

Replace the uploaded files in Quilter with this complete revision pair and
confirm all four connectors appear inside the outline as pre-placed parts.
[Quilter preserves their positions and rotations](https://docs.quilter.ai/design-parameters/pre-placed-components).
Existing routed candidates do not acquire these constraints retroactively.

| Directory | Contents |
|---|---|
| `1.2/` | Current complete rev 1.2 project for Quilter |
| [1.3/](1.3/README.md) | Experimental differential-sync project; requires placement, routing, and firmware support |
| [1.1/](1.1/README.md) | Routed rev 1.1 project, libraries, and detailed technical reference |
| `gen/` | Generators and validation tools |

Regenerate the current project from the repository root, schematic first:

```sh
python hardware/phantasm/gen/board.py --revision 1.2 --force
python hardware/phantasm/gen/pcb.py --revision 1.2 --unplaced --force --force-teensy-library
```

These commands replace the generated rev 1.2 files in `1.2/`.
Omitting `--revision` also selects rev 1.2. To generate the experimental
project independently, use `--revision 1.3` on both commands; outputs go
to `1.3/`. Always select the same revision for schematic and PCB generation.
Running the PCB generator without `--unplaced` writes a separate
`phantasm-draft.kicad_pcb` and matching project for placement experiments.
Before uploading a project that was opened in KiCad, restore its rule floors:

```sh
python hardware/phantasm/gen/heal_clearance.py hardware/phantasm/1.2/phantasm.kicad_pro
```

Upload the `.kicad_pro` with the schematic and PCB so Quilter can import the
fabrication floors: 0.2 mm trace width and copper clearance, 0.6 mm via diameter,
0.3 mm drill, and 0.5 mm copper-to-edge clearance. The Default net class uses
0.3 mm tracks and 0.6/0.3 mm vias. Terminal footprint keepouts include 0.5 mm
body clearance. Quilter's global component-spacing setting is separate; review
its 0.5 mm value in the job setup rather than assuming it imports from KiCad.

The rev 1.2 schematic wires C_DEC1 directly to U_MCU VIN and C_DEC2 directly to U1 pin 14
to make their bypass assignments explicit. Both are 100 nF and preplaced on the
PCB. [Quilter prioritizes direct schematic wires when assigning bypass capacitors](https://docs.quilter.ai/placement-constraints/bypass-capacitors).
Start a fresh job to redetect the constraints; replacing files can retain an
existing incorrect assignment. Verify these two rows before submitting.

The rev 1.2 bypass table contains only those two decouplers. Remove C_IN
(bulk storage), C_LF (filter damping), and C_SYNC (signal filtering) from that
table. In Power Nets, keep 500 mA sizing for all five detected supply nets and
turn off **Attempt Power Pour** for each. These are Quilter job settings;
there is no documented KiCad property to import these exclusions or pour flags.
After correcting the job, duplicate it to reuse its design parameters and
physics constraints for subsequent runs.

In Quilter, preserve the uploaded four-layer stackup and select **Preserve copper
on internal layers**. Also add **GND plane In1.Cu** and **GND plane In2.Cu**
by name to the job's **Preserved Pours** table. Quilter deletes and regenerates
pours absent from that table; preserving an internal pour also requires preserving
the stackup. See [Quilter's preserved-pour setup](https://docs.quilter.ai/design-parameters/preserved-pours).
The board's schematic links identify related components.

The stackup is **Signal / GND / GND / Signal**. Both inner layers have GND
pours and the user names `GND` and `Ground`, matching
[Quilter's ground-layer naming rules](https://docs.quilter.ai/using-quilter/prepare-your-input-board-file).
Confirm both inner layers are classified as Ground in Quilter's stackup editor.
Its no-power-layer info is expected: this design uses two ground reference
planes, with the supply rails routed on the outer layers.

The generator emits 0.25 mm minimum pour thickness, 0.5 mm thermal gaps, and
0.5 mm thermal spokes. Check downloaded routed boards before selecting a
candidate; the candidate analyzer rejects zone settings below the fabrication
floors even if KiCad DRC reports no errors:

```sh
python hardware/phantasm/gen/heal_zones.py path/to/routed.kicad_pcb --check
```

If Quilter replaced the zone settings, repair a separate copy:

```sh
python hardware/phantasm/gen/heal_zones.py path/to/routed.kicad_pcb -o path/to/repaired.kicad_pcb
```

The repair raises undersized features to the generator defaults and removes
cached copper fills. It preserves tracks, vias, placement, and zone boundaries.
The output must differ from the source; replacing an existing output requires
`--force`. Keep the repaired board with its matching schematic, project rules,
and libraries. **Refill all zones in KiCad, save, then rerun DRC and fabrication
validation.** The settings check alone does not validate the refilled copper or
prove that ground pads remain connected.

Electrical and mechanical requirements live in the
[PCB specification](../../docs/specs/phantasm_pcb_spec.md).

The TBC05-02 and TBC05-03 terminal-block footprints in `1.2/phantasm.pretty/`
and TBC05-04 in `1.3/phantasm.pretty/` are hand-maintained source files.
The PCB generator reads and copies them; the 1.3 TBC05-02 and TBC05-03
footprints are generated copies of the 1.2 sources.
