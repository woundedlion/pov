# PHANTASM rev 1.3 differential-sync prototype

Rev 1.3 is an experimental PCB project with an RS-485 electrical interface
for the existing count-coded sync protocol. It requires placement, routing,
electrical validation, and matching firmware before device use. Rev 1.2
remains the default generation target and retains its single-ended sync bus.

## Generate

From the repository root:

```sh
python hardware/phantasm/gen/board.py --revision 1.3 --force
python hardware/phantasm/gen/pcb.py --revision 1.3 --unplaced --force --force-teensy-library
```

Open or upload this directory's matching
[schematic](phantasm.kicad_sch), [board](phantasm.kicad_pcb), and
[project](phantasm.kicad_pro), with the accompanying libraries. These are
placement inputs, not routed fabrication artifacts. Keep downloaded candidates
and final routed files inside this revision's directory.

To regenerate rev 1.2, pass `--revision 1.2` to both commands, or omit the
selector. Generating one revision does not replace the other revision's files.

## Cable trunk

Use a specified 120-ohm twisted-pair cable with a separate reference-ground
conductor and, if fitted, a shield. Belden 8451 is nominally 45 ohms and is
not the cable specified for this interface.

Both sync connectors carry the same four nets:

| Pin | Signal | Wiring |
|---|---|---|
| 1 | A | First conductor of the twisted pair |
| 2 | B | Second conductor of the twisted pair |
| 3 | GND | Separate reference conductor; heavy wires carry LED returns |
| 4 | SHIELD | Shield/drain continuity |

Daisy-chain master to board 1 to board 2 to board 3. Each board's two
connectors join directly by net, with a short local transceiver branch.
Intermediate boards do not retransmit pulses. A powered-down board must
leave the bus usable; unplugging its cable interrupts downstream connectivity.
The cable may follow the rotor circumference, but ends at the last board:
do not connect the last board back to the master.

Enable the 120-ohm termination across A/B at the two physical ends only.
Leave intermediate termination jumpers open. Termination follows physical
cable position rather than the firmware segment ID. Carry shield continuity
through every board and close the shield-to-ground jumper only at the master.
The shield is separate from the reference-ground conductor.

The local transceiver branch has a matched 10 Ω pulse-proof resistor in
each conductor and a CDSOT23-SM712 protection array on the cable side.
The connectors and termination remain on the uninterrupted trunk side
of those resistors. Preserve pair symmetry and short protection return paths.

## Board layout

Rev 1.3 retains the 58.28 x 32 mm outline and mounting-hole centers.
The two four-position sync terminals need a different arrangement from
rev 1.2's three-connector column. Review the generated fixed placements,
terminal body keepouts, and hub-facing wire access before routing.
The fixed terminal origins are J2 at (36.6, 23.97), J3A at (48, 3.96),
and J3B at (48, 16.52), in millimeters. U1, R_D1, R_D2, and C_DEC2
are released for placement in the remaining area. The removed receive divider
and C_SYNC require no local filter traces.

For Quilter, set bypass assignments to C_DEC1 at U_MCU VIN, C_DEC2 at U1
pin 14, and C_DEC3 at U_SYNC pin 8. Their direct schematic wires identify
the parent supplies. C_BULK3 is bulk decoupling; place it near U_SYNC without
a bypass assignment. C_DEC2 and C_DEC3 need placement at their IC supply pins.

## Firmware contract

| Teensy pin | Function | Required behavior |
|---|---|---|
| D3 | FRAME_SYNC receive | Input; timestamp falling edges |
| D4 | SYNC_TX | Output; idle HIGH, active LOW pulses |
| D5 | SYNC driver enable | HIGH enables the driver; LOW disables it |

Initialize D5 LOW before configuring transmit or interrupts. Initialize D4
HIGH before enabling the master. Keep follower drivers disabled and receivers
enabled. The master drives the idle state continuously between pulses.
Arm receive interrupts after initialization; the physical leading edge is the
falling edge of the receiver output. The symbol count and pitch retain their
existing meanings.

The current firmware is not compatible with this pin and polarity contract.
The rev 1.2 MASTER_EN signal has the opposite enable polarity. Do not run a
legacy image on this prototype. The pulse-width redesign discussed separately
is not implemented by the PCB generator.

## Validation before fabrication and operation

Validate connector fit, wire access, component placement, and the mounting
envelope; perform routing, ground-plane refill, ERC, DRC, and PCB/schematic
net parity checks before manufacturing a prototype. Confirm the 3.3 V
regulator budget with the transceiver
driving both end terminations, including supply transients and decoupling.
Qualify protection and its clamp coordination against the transceiver limits
before treating the project as a fabrication design.
The fabrication exporter rejects rev 1.3 until its supplier catalog and
assembly metadata are qualified.

Measure A minus B and each bus conductor relative to local ground at the
nearest and farthest receivers. Exercise motor startup, speed changes, LED
load steps, reset, power loss at an intermediate board, and cable removal.
Check received edge counts, pulse width, voltage margin, common-mode limits,
and board-to-board timing skew. Integrated failsafe does not establish noise
margin on an undriven cable; assess whether one external bias network is needed.

The [PCB specification](../../../docs/specs/phantasm_pcb_spec.md#revision-13-experimental-differential-sync)
defines the electrical contract. Reference sources:

- [TI RS-485 design guide](https://www.ti.com/lit/pdf/slla272).
- [TI THVD1410 family datasheet](https://www.ti.com/lit/ds/symlink/thvd1450.pdf).
- [Belden 8451 specifications](https://www.belden.com/products/cable/audio-cable/analog-audio-cable/8451).
