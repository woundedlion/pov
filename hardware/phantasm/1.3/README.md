# PHANTASM rev 1.3 differential-sync design

Rev 1.3 uses RS-485 for the existing count-coded sync protocol. The electrical
review selects the components below for four boards with at most 1 m between
adjacent boards, entirely on the rotor (3 m total for three inter-board cables).
**This is a schematic and placement input, not a fabrication release.**
Placement, routing, power/harness qualification and matching firmware remain
required. Rev 1.2 remains the generator default.

## Selected components

The schematic carries Manufacturer, MPN and Datasheet fields for the differential
interface, decouplers and control resistors. These identify actual parts;
assembly-house stock, substitutions and CPL rotations still require verification.

| Reference | Manufacturer / orderable part | Selection |
|---|---|---|
| U_SYNC | Texas Instruments **THVD2410DR** | 3.3 V, SOIC-8, 500 kbps, slew-limited transceiver |
| D_SYNC | Bourns **CDSOT23-SM712** | RS-485 TVS array, SOT-23; pins 1/2 to A/B, pin 3 to GND |
| R_A, R_B | Vishay **CRCW060310R0FKEAHP** | 10 Ω, 1%, 0603, pulse-proof, 0.33 W at 70°C |
| R_TERM | Vishay **CRCW0805120RFKEAHP** | 120 Ω, 1%, 0805, pulse-proof, 0.5 W at 70°C |
| C_DEC1, C_DEC2, C_DEC3 | YAGEO **CC0603KRX7R9BB104** | 100 nF, 50 V, X7R, 10%, 0603 |
| C_BULK3 | YAGEO **CC0603KRX7R8BB105** | 1 µF, 25 V, X7R, 10%, 0603 |
| C_IN | Nichicon **UPW1H101MPD** | 100 µF, 50 V, 105°C, radial 8 × 11.5 mm, 3.5 mm lead pitch; retain mechanically |
| R_TX, R_MEN, R_DATA_PD, R_CLK_PD | YAGEO **RC0603FR-0710KL** | 10 kΩ, 1%, 0603; TX pull-up, DE and LED-input pull-downs |
| J3A, J3B | GCT **TBC05-04-1-G-G** | Four-position screw terminal; hand assembly and harness retention required |

Component references: [TI transceiver](https://www.ti.com/lit/ds/symlink/thvd2410.pdf),
[Bourns TVS](https://www.bourns.com/docs/product-datasheets/cdsot23-sm712.pdf),
[Vishay pulse resistors](https://www.vishay.com/docs/20043/crcwhpe3.pdf),
YAGEO [100 nF](https://www.yageogroup.com/download/specsheet/CC0603KRX7R9BB104),
[1 µF](https://www.yageogroup.com/download/specsheet/CC0603KRX7R8BB105), and
[10 kΩ](https://www.yageogroup.com/component-documentation/download/specsheet/RC0603FR-0710KL).
Resistor power ratings require temperature derating and suitable PCB heat spreading.

Inherited power/LED selections remain U1 **SN74AHCT125DR**, Q_REV **AO3401A**,
F1 **TLC-NSMD050**, FB **CBW321609U601T**, C_LF **CL31A226KAHNNNE**,
R_LF **RCA051R5JLF**, and R_D1/R_D2 **0805W8F330JT5E**. J1 is
**TBC05-02-1-G-G**, J2 is **TBC05-03-1-G-G**, and U_MCU is a **Teensy 4.0**.
The inherited supplier mappings are in [the fabrication tool](../gen/fab.py);
its rev 1.1 BOM (the only revision it exports) must not be used to assemble rev 1.3.
[Nichicon's C_IN selection](https://www.nichicon.com/en-us/part/upw1h101mpd/8471/)
matches the existing radial footprint; confirm its height against the rotor envelope.

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

Use a 120-ohm shielded twisted pair. **Belden 3106A, black jacket (color 010)**
includes a separate insulated reference-ground conductor and is the preferred
electrical selection. **Belden 3105A 010500** also has a suitable signal pair;
its reference must come from a separate insulated wire alongside the cable or
the existing common power-ground harness after the checks below. Belden 8451
is nominally 45 ohms and is not specified for this interface.

Both sync connectors carry the same four nets:

| Pin | Signal | Wiring |
|---|---|---|
| 1 | A | First conductor of the twisted pair |
| 2 | B | Second conductor of the twisted pair |
| 3 | GND | Reference conductor when fitted; heavy wires carry LED returns |
| 4 | SHIELD | Shield/drain continuity |

These are four electrical connections, not four insulated wires. A/B carry the
differential signal current. The reference bounds ground-potential differences;
the shield intercepts interference and is not a substitute for that reference.
With 3105A, use white/blue for A and blue/white for B. With 3106A, use
white/orange for A, orange/white for B and blue/white for reference GND.

To omit the separate reference wire with 3105A, all boards must retain a
low-impedance common power-ground connection, including while their positive
supplies are switched off. Leave terminal 3 empty; preserve shield continuity
on terminal 4 and bond it only at the master. Measure each A/B conductor
relative to each receiver's local ground under motor and LED disturbances.
Target -5 V to +10 V, providing 2 V headroom inside the populated TVS's
-7 V/+12 V working limits. This is an acceptance target, not a measured result.
If common ground can disconnect independently, this wiring is not qualified.
A separate reference wire also needs a fault-current review: parallel ground
paths can carry load current and must not replace the heavy returns.

[3105A](https://www.belden.com/products/cable/electronic-wire-cable/multi-pair-cable/3105a)
is 7.47 mm diameter with an 85.65 mm minimum bend radius;
[3106A](https://www.belden.com/products/cable/electronic-wire-cable/multi-pair-cable/3106a)
is 7.87 mm diameter with about a 91 mm minimum bend radius. Clamp the jacket
to the rotor structure and balance the harness. The terminal screws must not
carry the cable's centrifugal load. Electrical suitability does not establish
fit or rotational retention. The 010500 suffix specifies a black 500-foot reel,
not the required purchase quantity; obtain suitable cut lengths.

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

The local transceiver branch has a 10 Ω, 1% pulse-proof resistor in
each conductor and a CDSOT23-SM712 protection array on the cable side.
The connectors and termination remain on the uninterrupted trunk side
of those resistors. Preserve pair symmetry and short protection return paths.

## Electrical margins

THVD1410's ±18 V bus absolute limits do not cover the SM712's positive clamp
voltage. SM712 specifies +26 V/-14 V maximum clamp at 17 A; THVD2410's ±70 V
bus limits provide headroom. Trace inductance, pulse energy and protection-return
impedance still need qualification. The assembled board is **not** a ±70 V
continuous-fault interface: the TVS conducts earlier. Its usable ground-offset
envelope is constrained by the TVS's -7 V/+12 V standoff despite the IC's wider
common-mode range. These ratings come from the manufacturer references above.

The two 120 Ω terminations present 60 Ω nominal. A first-order DC budget uses
59.4 Ω minimum termination load, 20.2 Ω maximum total branch resistance,
and 1.5 V differential at the driver pins:

`V_AB = 1.5 × 59.4 / (59.4 + 20.2) = 1.12 V`.

This leaves about 0.92 V above the receiver's 0.2 V threshold magnitude before
cable resistance, leakage, temperature drift and disturbances. It is a resistive
estimate, not a simulated or measured worst-case guarantee. Require at least
1.0 V settled differential magnitude at every receiver in both driven states.

Master DE remains asserted between pulses, making idle actively driven. No
external bias is fitted: a weak network provides little guaranteed noise margin,
loads the active state and can inject current into a dead rail. Integrated
failsafe can take up to 18 µs to establish HIGH after loss of drive; it does not
guarantee rejection of arbitrary interference on an undriven cable. Firmware
must tolerate loss of master and reacquire from valid symbol bursts.

At 3.6 V, a conservative normal DC load estimate is
`3.6 / (59.4 + 19.8) + 0.0056 = 0.0511 A` for U_SYNC and the bus. Reserve
**60 mA external 3.3 V capacity**, plus other peripherals; this is below
[PJRC's 250 mA external allowance](https://www.pjrc.com/store/teensy40.html).
Budget **0.25 A for the complete logic card** until measured. Shorted-bus current
and regulator thermal behavior need separate testing. The inherited 0.15 A
input-drop calculation is insufficient for this master.

## Board layout

Rev 1.3 retains the 58.28 x 32 mm outline and mounting-hole centers.
The two four-position sync terminals need a different arrangement from
rev 1.2's three-connector column. Review the generated fixed placements,
terminal body keepouts, and hub-facing wire access before routing.
The fixed terminal origins are J2 at (36.6, 23.97), J3A at (48, 3.96),
and J3B at (48, 16.52), in millimeters. U1, R_D1, R_D2, and C_DEC2
are released for placement in the remaining area. The removed receive divider
and C_SYNC require no local filter traces.

Route J3A to J3B as an uninterrupted pair over a continuous ground plane.
Keep the connector-to-transceiver branch at most 20 mm and untwisted cable
ends at most 10 mm as layout targets. Put the TVS at the connector branch,
with a short, wide return and adjacent ground vias. Put R_A/R_B between that
protected trunk and the IC. Keep both paths symmetric; do not add a common-mode
choke or bus capacitors without measurements showing a need. Specify the
fabricator's stackup before selecting trace width/gap for 120 Ω differential
impedance; no universal geometry is implied by these files.

Keep U1's 33 Ω DATA/CLK resistors immediately at its outputs, and place
R_DATA_PD/R_CLK_PD at its inputs to define the reset state. Separate SPI traces
and LED/power-return paths from the sync branch. Place C_DEC3 within 3 mm
of U_SYNC pin 8 with a short ground return; place C_BULK3 nearby.

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

Wait at least 20 µs after valid supply and stable idle before arming receive
interrupts, then clear pending interrupt state. Use at least 5 µs LOW pulse
width as an initial hardware qualification target; confirm delivered pulse
width and single falling edges with the complete loaded harness. THVD2410's
specified differential edges are 240–600 ns under its datasheet test load;
the real cable and four protection arrays add capacitance.

The current firmware is not compatible with this pin and polarity contract.
The rev 1.2 MASTER_EN signal has the opposite enable polarity. Do not run a
legacy image on this prototype. The pulse-width redesign discussed separately
is not implemented by the PCB generator.

## Validation before fabrication and operation

Validate connector fit, wire access, component placement, and the mounting
envelope; perform routing, ground-plane refill, ERC, DRC, and PCB/schematic
net parity checks before manufacturing a prototype. Confirm the 3.3 V
regulator budget with the transceiver driving both end terminations, including
supply transients and decoupling. Require +3V3 to remain within 3.0–3.6 V and
+5V_LOGIC within U1's 4.5–5.5 V operating range at load and temperature. Recheck
F1's hot/post-trip resistance; a nominal 5 V source alone does not prove this.
The input retains a reverse-polarity FET and PTC but has **no qualified positive
overvoltage cutoff**. Verify slip-ring overshoot and establish upstream regulated
protection or redesign this input before electrical release. A generic 5 V TVS
whose clamp exceeds the downstream limits does not resolve that issue.

J1 remains an unkeyed screw terminal. A polarized, retained harness connection
and controlled terminal wiring are required: Q_REV cannot prevent a reversed
GND connection shorting through the sync, strip or USB grounds. Cut the Teensy
VIN/VUSB link before external-power operation with USB, per the assembly spec.
The fabrication exporter rejects rev 1.3 until its supplier catalog and
assembly metadata are qualified.

Measure A minus B and each bus conductor relative to local ground at the
nearest and farthest receivers. Exercise motor startup, speed changes, LED
load steps, reset, power loss at an intermediate board, and cable removal.
Check received edge counts, pulse width, voltage margin, common-mode limits,
and board-to-board timing skew. Record no extra or missing received pulses
during a sustained test with the real motors and LED loads. Repeat with each
follower unpowered, with the master reset/unpowered, and on reconnect. Test the
chosen reference-ground wiring explicitly. Power the prototype on a current-limited
bench supply before rotor operation. No EMC compliance or rotor qualification
is established by component ratings or schematic ERC.

The [PCB specification](../../../docs/specs/phantasm_pcb_spec.md#revision-13-experimental-differential-sync)
defines the electrical contract. Reference sources:

- [TI RS-485 design guide](https://www.ti.com/lit/pdf/slla272).
- [TI THVD2410 family datasheet](https://www.ti.com/lit/ds/symlink/thvd2410.pdf).
- [Belden 8451 specifications](https://www.belden.com/products/cable/audio-cable/analog-audio-cable/8451).
