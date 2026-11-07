# PHANTASM rev 1.2 board project

Rev 1.2 is the generator default. The [board index](../README.md) covers shared generation and [release validation](../README.md#validation-and-release).

## Quilter upload

Upload these three files together from that directory:

- [phantasm.kicad_pcb](phantasm.kicad_pcb): unplaced rev 1.2 board.
- [phantasm.kicad_sch](phantasm.kicad_sch): matching rev 1.2 schematic.
- [phantasm.kicad_pro](phantasm.kicad_pro): placement and design rules.

Keep the accompanying symbol and footprint libraries with the project when
opening it in KiCad 10. Every upload must use the board and schematic from the
same revision directory. Rev 1.2 requires placement, routing, and validation
before fabrication.

The rev 1.2 upload board fixes **J1 power at the left edge** and **J2 LED, J3A SYNC IN,
J3B SYNC OUT in a column at the right end**, in that top-to-bottom order.
The column sits inboard of the right mounting holes, in the rev 1.1 connector
area. Their wire openings face left toward the hub. Keepouts reserve each
terminal's full plastic body plus 0.5 mm assembly clearance, preventing SMD
parts beneath the connector. Power's wire-access reservation extends off the
left edge. The outline remains **58.28 x 32 mm**, and all four mounting holes
retain the rev 1.1 coordinates. The Teensy shifts 3.5 mm left to clear the
terminal bodies; D_BUS is released for placement outside their keepouts.

Rev 1.2 J2/J3A/J3B have no wire-access reservation. Their openings sit
about 1.5 mm from the Teensy PCB edge; J2.1 faces a header pin that remains
in its insertion path at any socket height. Trial-fit its harness before
fabrication. Placement and DRC do not qualify harness insertion. A connector whose pin
legend meets a mounting reserve has its complete legend on the back; rev 1.2
J2 and J3B use the back, while J3A retains its complete front legend.

Rev 1.2 omits the J4 debug header and SERIAL1_TX connection; Teensy pin 1 is
unconnected. R_MEN and MASTER_EN remain part of the sync control circuit.
Probe power directly at the Teensy; diagnostics use its USB connector.

Replace the uploaded files in Quilter with this complete revision pair and
confirm all four connectors appear inside the outline as pre-placed parts.
[Quilter preserves their positions and rotations](https://docs.quilter.ai/design-parameters/pre-placed-components).
Existing routed candidates do not acquire these constraints retroactively.

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

## Hardware and firmware

Rev 1.2 selects **ENIG** copper finish; the shipped rev 1.1 board uses
**Lead-Free HASL**. This is a revision-specific fabrication change.

### Terminal-block clearance

The generator uses **GCT TBC05-02-1-G-G** for J1 (board power) and
**TBC05-03-1-G-G** for J2 (LED signals), J3A and J3B (sync). J4 is absent from rev 1.2.
The [GCT mechanical drawing](https://www.farnell.com/cad/4513152.pdf) specifies
2.54 mm pitch, 1.3 mm PCB holes, and bodies **5.48 × 6.5 mm** (two positions)
or **8.02 × 6.5 mm** (three positions), 8.5 mm high.

Each terminal footprint reserves its body plus 0.5 mm on every side with a
front courtyard and a footprint-local component keepout. Traces, vias and
copper pours remain allowed beneath the plastic body. The placed draft packs
against the larger reservation; the unplaced output locks these connectors
at `TERMINAL_EDGE_PLACEMENTS`. Recheck
wire entry, screw access and rotor clearance after placement.

The committed rev 1.1 boards still carry pin-header footprints. Matching pin
pitch alone does not establish body or drill clearance for these terminals.
J1 remains pin 1 = +5 V, pin 2 = GND; the terminal block is not polarized.

### Sync circuit

Every master and follower uses the same PCB, populated parts, and pin map.
The ID straps select the role; `MASTER_EN` still gates both U1 channels C/D.

| Signal | Connections | Master | Follower |
|---|---|---|---|
| `FRAME_SYNC` | Teensy D3, R1/R2 divider, C_SYNC 220 pF | Receive-only, echo ignored | Receive with pad hysteresis |
| `SYNC_TX` | Teensy D4, U1 pin 9, R_TX 10 kΩ to GND | Sync pulse output | Output held LOW |
| `MASTER_EN` | Teensy D5, U1 pins 10/13, R_MEN pull-up | LOW after TX initialization | HIGH |

`R_TX` is populated on **every board**, including the master. It holds the AHCT
input LOW during reset while the GPIO is high impedance. The receive RC no
longer loads an AHCT input. `R_TX` uses the same 0603 10 kΩ catalog part as
`R_MEN` (C25804).

Rev 1.2 also closes rev 1.1's power-up feedback path. R_MEN pulls `/OE` to the
Teensy's 3V3 rail, which can rise after U1's 5 V supply; it cannot guarantee a
disabled driver throughout power-up. On rev 1.1 the briefly enabled receive-to-
transmit path can re-drive SYNC_BUS. Rev 1.2's separate SYNC_TX input and R_TX
hold its input LOW during that interval. This protection requires the revised
copper and firmware described below; the generator alone does not retrofit it.

**Firmware compatibility:** `HS_PHANTASM_BOARD_REV` comes from
`[phantasm_common] board_revision` in `platformio.ini`, which defaults to `11`.
Build rev 1.2 with `pio run -e phantasm_rev12`, or set `board_revision = 12`
for all Phantasm environments. The rev 1.2 pin map transmits on D4, configures
D3 as input with HYS, and initializes D4 LOW as an output on every board before
enabling the master. Followers keep D4 LOW; only the master emits pulses.
Regenerating the board does not change the firmware build flags or routed copper.

Generating rev 1.2 requires regenerating the schematic before either PCB draft.
The resulting PCB is unrouted and needs placement/routing validation before
fabrication; changing its revision label alone does not convert rev 1.1 copper.
