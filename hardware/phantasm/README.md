# PHANTASM Segment Board — KiCad project

KiCad 10 schematic for the per-segment carrier board specified in
[../../docs/specs/phantasm_pcb_spec.md](../../docs/specs/phantasm_pcb_spec.md). One identical
PCB is built ×4 for the qualified configuration; a solder strap selects each
board's role (segment 0 = master/conductor, 1–3 = flywheel slaves). The default
firmware is **N = 4** (ID0/ID1 via `JP_ID0` / `JP_ID1`). The compile-tested
**N = 8** profile uses eight boards and reads `JP_ID2` (pin 23) for segments 4–7.
Its rotor mounting, balance, cabling, and swept envelope are not mechanically
qualified. The card is **logic-only (~0.15 A)** — each strip draws about 4.3 A
at N=4 or 2.2 A at N=8, injected **off-board** (spec §2.3), so nothing here
carries high current.

Contains the **schematic**, the corrected `C_IN`, `C_LF`, `C_DEC1/2`, and
`C_SYNC` placement/routing, and the completed Quilter routing. The main PCB is
the fabrication source of truth and has no unconnected pads.

## Files

| File | What |
|---|---|
| `phantasm.kicad_pro` | Project file |
| `phantasm.kicad_sch` | Schematic — all parts, values, footprints, full §10 connectivity. **Last written by KiCad, not by `gen/board.py`**, and not reproducible from it — see [Regenerating](#regenerating) |
| `phantasm.kicad_pcb` | Completed routed PCB with validated placement, control routing, planes, mounting, and service clearances |
| `quilter_incremental/` | Historical protected input snapshot used for the completed control-net routing — written by `gen/make_quilter_incremental.py`, which now refuses to run because the board is routed |
| `unplaced/phantasm_unplaced.kicad_pcb` | **4-layer** (SIG/GND/GND/SIG) outline + net-assigned footprints **staged below the board, unrouted** — for an autoplacer (Quilter). Stackup encoded in-file. A KiCad GUI re-save; `python gen/pcb.py --unplaced --force` regenerates it and discards that re-save |
| `unplaced/phantasm_unplaced.kicad_pro` | **Captured artifact — no generator writes it.** Quilter-facing DRC rules and net class for the unplaced board, wider than the routed project's. Pinned by `gen/constraints.py` (`UNPLACED_RULES`, `UNPLACED_DEFAULT_CLASS`) and `gen/tests/test_constraints.py`; restore from those if it is ever lost, because `gen/heal_clearance.py` only raises a value below a minimum and would leave the routed floors in place |
| `phantasm.kicad_sym` | Project symbol library: custom `Teensy4.0` + `+5V_RAW/+5V_LOGIC` power symbols |
| `phantasm.pretty/` | Project footprint library: generated `Teensy4.0` footprint (2×14 0.1″ THT) |
| `sym-lib-table` / `fp-lib-table` | Register the `phantasm` symbol / footprint libraries |
| `gen/` | The Python generators (schematic + PCB) — see [Regenerating](#regenerating) |

Open `phantasm.kicad_pro` in KiCad 10. Stock symbols/footprints come from the
standard KiCad libraries; the custom Teensy + power symbols and the Teensy
footprint come from the project `phantasm.kicad_sym` / `phantasm.pretty`.

## Validation

Every bullet but the ERC one maps to executed code. The gates that read the
committed board directly need no KiCad and run in CI
(`python -m unittest discover -s hardware/phantasm/gen/tests`,
`python gen/board_metadata.py --check`); the rest need a local KiCad 10
(`kicad-cli`) and run when the fab package is regenerated.

- **N=8 firmware:** `pio run -e phantasm8` compiles and links the optional
  eight-board profile; this is firmware validation, not rotor qualification.
- **ERC: 0 errors — verified once by hand.** No runner exists: `gen/` has no
  ERC step, `erc*.rpt` is gitignored, and neither CI nor the `just` recipes
  re-run it, so re-check it in KiCad (`kicad-cli sch erc`) after any schematic
  change. A warning-inclusive KiCad 10 run reports nine `lib_symbol_mismatch`
  notices for embedded copies of stock/custom symbols; the exported
  connectivity is verified separately below.
- **Netlist matches the electrical specification** — `gen/fab.py` holds the
  netlist it exports to the named-net table in `gen/check.py`, which also runs
  standalone against the committed schematic: every net in the table must match
  member-for-member, keyed on
  `(ref, pin)` so a connector or IC pinout permutation fails. A named net outside
  the table is reported as a `NOTE` and does not fail the gate. CI enforces the
  same table without KiCad: `gen/tests/test_check.py` applies it to the pad nets
  of the committed board, so copper that stops matching the spec fails on every
  push.
  The required connections are realized with the correct members
  (logic feed `J1 → F1 → Q_REV → FB → +5V_LOGIC`;
  series terminations `U1 out → R → J2`/bus; the pin-3 divider node ties Teensy D3,
  `U1` ch-C input, `R1`/`R2`/`C_SYNC`; ID0/ID1/ID2 straps; `MASTER_EN`; shield).
- **Copper-connectivity gate:** `gen/connectivity.py` unions the routed board's
  tracks, vias, pads and pour fills per net and rejects any net whose pads land
  in more than one island. Every other net gate in this list reads pad net
  *attributes*, which survive deleting the copper that realizes them; this one
  walks the copper. It needs no KiCad and runs in CI via
  `gen/tests/test_connectivity.py`; run it standalone with
  `python gen/connectivity.py`.
- **PCB geometry DRC: clean** (`kicad-cli pcb drc`): zero error-severity
  violations and zero unconnected pads.
- **Standard-cost via gate:** `gen/fab.py` rejects a routed board containing a
  via smaller than 0.45 mm or a drill smaller than 0.20 mm, and rejects any
  via pair with less than 0.15 mm of copper spacing (pad edge to pad edge).
- **Schematic parity gate:** `gen/fab.py` runs `kicad-cli pcb drc
  --schematic-parity` and rejects any board/schematic difference outside
  `KNOWN_PARITY_ITEMS` (mounting holes, jumper BOM attributes, the Teensy
  footprint field), so gerbers for stale copper cannot ship with a BOM
  exported from a newer schematic. Parity items are warning severity in
  KiCad, so this runs separately from the error-severity DRC gate.
- **Plot-origin gate:** the gerbers, Excellon drill, and CPL are all exported
  in absolute board coordinates, and `gen/fab.py` rejects a board carrying a
  non-zero drill/place origin (`aux_axis_origin`) — that would move only the
  origin-relative exports and place every assembled part off-board.
- **Zone-geometry gate:** `gen/fab.py` rejects a copper pour whose
  `min_thickness`, `thermal_gap`, or `thermal_bridge_width` is below the
  0.1016 mm (4 mil) minimum feature the fab resolves. KiCad DRC never flags
  these — thermal reliefs are same-net geometry — so a sub-process gap would
  export clean gerbers the fab fills as a solid pour, tying every
  through-hole GND pad to the full plane and starving the hand-soldered
  joints (R-ASM-6).
- **Shipped-land gate:** `gen/tests/test_pcb_lands.py` pins the pad geometry the
  routed board ships for every chip passive, so restoring a stock land over the
  widened sync-resistor pads fails in CI. It pins as-built geometry, not spec
  §11.1 geometry — see the lands note below. KiCad reports no parity difference for
  that edit — the widened pads keep the stock footprint id. The same file pins
  J1's footprint on both artifacts, which the assembly gate cannot see because
  `EXCLUDE_FP_SUBSTR` excludes every hand-soldered connector spelling.
- **Count-floor ratchets:** `gen/fab.py` rejects a board carrying fewer than
  100 vias or fewer than 2 copper pours — floors pinned to what the committed
  routed board holds (both counts are in the facts block below; the pours are
  the In1/In2 reference planes, and the mounting-hole keepout rule areas are
  counted separately because they pour no copper). Neither floor updates
  itself: after promoting a re-route, re-measure with `gen/board_metadata.py`
  and deliberately re-baseline the constants in `gen/fab.py`.
- **Assembly-metadata gate:** `gen/fab.py` rejects the assembly outputs
  unless the assembled references exactly match its LCSC assignment table,
  every rotation correction names an assembled part, and every centroid row
  is present, numeric, and on the top (assembly) side.
- **Layer-name gate:** `gen/tests/test_pcb_stack.py` rejects a declared layer
  name carrying whitespace on either board. KiCad builds each Gerber's filename
  from the layer name, so an Altium-style alias such as `Ground Layer 1` ships a
  space in the upload zip.
- **Board-revision gate:** `gen/tests/test_revision.py` reads the revision off
  the bottom silkscreen (`Phantasm Rev 1.1`) and requires the schematic title
  block, the routed board's title block and `gen/builder.py`'s `REVISION` to
  agree with it. The board title block is what KiCad writes into the Gerber X2
  `ProjectId` attribute; without it the gerbers ship `rev?`.
- **Part-catalog gate:** every assigned LCSC number must resolve to a catalog
  entry with a non-blank manufacturer, MPN, and description, so each JLCPCB
  BOM match is independently auditable.

## How connectivity is drawn

The sheet is organised left-to-right into labelled blocks (power entry, Teensy +
level shifter → LED strip, sync RX divider + daisy, ID straps / debug, power flags):

- **Power distribution is drawn as visible rails** — a horizontal `+5V_LOGIC` rail
  and `GND` rail with vertical component drops and junctions (the classic ladder); the
  `J1 → F1 → Q_REV → FB` protection/filter chain feeds the rail's left (hub) end.
- **Series/divider paths are wired** — `U1` outputs through the 33 Ω/100 Ω source
  terminators, and the pin-3 RC divider (`R1`/`R2`/`C_SYNC`).
- **Cross-block signals use net labels as ports** — `DATA(_IN)`, `CLK(_IN)`,
  `FRAME_SYNC`, `MASTER_EN`, `SYNC_BUS`, `ID0/1`, `SHIELD` — the conventional way to
  avoid dragging wires around the large Teensy symbol.

It's still a generated **functional layout**; rearrange/beautify freely in Eeschema —
the netlist is what's verified.

## BOM → symbol → footprint

| Ref(s) | Symbol | Footprint | Notes |
|---|---|---|---|
| `U_MCU` | `phantasm:Teensy4.0` | `phantasm:Teensy4.0` (2×14 0.1″ THT) | pad map = top view (component side up), USB end at −X: top row VIN,GND,3V3,23…13 / bottom row GND,0…12; **cut the VIN/VUSB pad before install** — R-ASM-7, see the hand-assembly section below |
| `U1` (A–E) | `74xx:74AHCT125` | `Package_SO:SOIC-14_3.9x8.7mm_P1.27mm` | 4 buffers + power unit |
| `Q_REV` | `Transistor_FET:Q_PMOS_GSD` (AO3401A) | `Package_TO_SOT_SMD:SOT-23` | reverse-polarity protection; pin 3 drain=input, pin 2 source=output, pin 1 gate=GND |
| `F1` | `Device:Fuse` (0.5 A hold) | `Fuse:Fuse_1206_3216Metric` | TLC-NSMD050 resettable fuse; 1 A trip, 0.75 Ω post-trip maximum |
| `C_IN` | `Device:C_Polarized` (100µF) | `Capacitor_THT:CP_Radial_D8.0mm_P3.50mm` | only on-card electrolytic, on +5V_LOGIC; RTV-bond |
| `FB` | `Device:FerriteBead` | `Inductor_SMD:L_1206_3216Metric` | ~600 Ω @100 MHz |
| `R_LF` / `C_LF` | `Device:R` / `Device:C` | 0805 / `C_1206` | bead-LC damper, 22 µF |
| `C_DEC1/2` | `Device:C` (0.1µF) | `Capacitor_SMD:C_0603_1608Metric` | |
| `R_D1/R_D2` | `Device:R` (33Ω) | `Resistor_SMD:R_0805_2012Metric` | DATA/CLK source term |
| `R_S` | `Device:R` (100Ω) | `Resistor_SMD:R_0805_2012Metric`, pads widened to 1.40 mm | SYNC source term; widened land (see the lands note below) |
| `R1/R2` | `Device:R` (10k/15k) | `Resistor_SMD:R_0603_1608Metric`, pads widened to 1.20 mm | sync divider; widened land (see the lands note below) |
| `C_SYNC` | `Device:C` (220pF) | `C_0603` | populated (noise filter) |
| `R_PD` | `Device:R` (10k) | `Resistor_SMD:R_0603_1608Metric`, pads widened to 1.20 mm | master-only bus idle pull-down, widened land; ground-side switched automatically by U1 channel D |
| `R_MEN` | `Device:R` (10k) | `R_0603` | MASTER_EN boot pull-up → 3V3 |
| `D_BUS` | `Device:D_Zener` (Bourns CDSOD323-T05L) | `Diode_SMD:D_SOD-323` with Bourns pad geometry | populated unidirectional 5 V, 1 pF sync-bus TVS; pin 1/cathode on SYNC_BUS, pin 2/anode on GND; exact Bourns land pattern; silkscreen bar marks the cathode end; JLCPCB C1975255 |
| `J1` | `Connector_Generic:Conn_01x02` | `Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical` | +5 V/GND light logic feed, ~1 A; **unkeyed** 0.1″ header — R-PWR-7's keying is unmet on the shipped board (see the deviations note below) |
| `J2` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | strip **signal only**: DI / CI / SIG_GND (no power) |
| `J3A/J3B` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | Belden 8451 daisy |
| `J4` | `Connector_Generic:Conn_01x04` | `PinHeader_1x04_P2.54mm` | debug/serial |
| `H1`–`H4` | — | `MountingHole:MountingHole_2.7mm_M2.5` | four NPTH rotor mounting holes |
| `JP_SHLD/JP_ID0/JP_ID1/JP_ID2` | `Jumper:SolderJumper_2_Open` | `SolderJumper-2_P1.3mm_Open_...` | shield (master only) / ID straps (JP_ID2 read at N=8) |

## Hand assembly (not done by the PCBA house)

The assembly house reflows top-side SMD only; `gen/fab.py` excludes the Teensy,
the connectors, the electrolytic, and the solder jumpers from its BOM/CPL, so
every step below is performed by whoever builds the card.

- **R-ASM-7 — cut the VIN/VUSB pad on every Teensy 4.0. Mandatory, before the
  Teensy is soldered down.** The board feeds Teensy `VIN` from the rotor rail and
  `J4` exposes serial for USB debug. With the pad intact, `VUSB` is tied to `VIN`:
  plugging USB into a powered board back-feeds the live 5 V rotor rail into the
  host's USB `VBUS` (and a host's `VBUS` into the rotor rail when the rail is
  down). Nothing on the card blocks it — `Q_REV` protects the `J1` feed, not the
  USB port — and **no silkscreen or artifact carries this step**, so a build that
  skips it looks correct. Cut the trace between the `VIN` and `VUSB` pads on the
  Teensy's underside; `VIN` is then fed only from the board rail, and the Teensy
  no longer self-powers from USB.
- **`JP_SHLD` is stuffed on the master board only** (R-ASM-4); `JP_ID0`/`JP_ID1`
  strap the segment ID, `JP_ID2` only at N = 8.
- **`C_IN` is RTV-bonded** after soldering (R-PWR-6 / R-MECH-3).

## Notes / deviations from the spec

- **Reverse protection uses one AO3401A P-channel MOSFET** (`Q_REV`, SOT-23), with
  its gate tied directly to GND. It replaces the series Schottky without adding a
  gate resistor or increasing the component count. At a 4.75 V J1 input and 0.15 A,
  a conservative hot calculation uses the fuse's 0.75 Ω post-trip maximum, the
  bead's 0.20 Ω maximum DCR, and twice the MOSFET's 60 mΩ maximum at −4.5 V:
  `VLOGIC = 4.75 − 0.15 × (0.75 + 0.20 + 0.12) = 4.5895 V`. This leaves about
  90 mV above the AHCT125's 4.5 V minimum; verify that J1 itself remains at or above
  4.75 V on the hot, operating rotor because external harness drop is not included.
- **Widened lands on the bench-tuned sync resistors — spec §11.1 is not met on this board.**
  §11.1 mandates the toe-extended KiCad `_HandSolder` land, which keeps the IPC-nominal
  inter-pad gap and adds no copper under the ceramic body. The shipped copper does neither,
  and both parts are reflow-placed. `R1`, `R2` (divider ratio, spec §4.2), `R_PD` (bus idle
  pull-down) and `R_S` (source termination) keep the **stock**
  `Resistor_SMD:R_0603_1608Metric` / `R_0805_2012Metric` footprint id in
  `phantasm.kicad_pcb`, with the pads **widened in place**: the centres stay at the stock ±0.825 mm (0603) / ±0.9125 mm (0805) while the
  width grows from 0.80 → **1.20 mm** (0603) and 1.025 → **1.40 mm** (0805). Toe and heel
  both grow, so the inter-pad gap closes from the IPC-nominal 0.85 mm / 0.80 mm to
  **0.450 mm** (0603) / **0.425 mm** (0805) — the pad inner edge sits 0.225 mm /
  0.2125 mm off the centre-line and copper does reach under the bare ceramic body.
  Three chip lands are live at once: these four, the stock land on every other chip
  passive (`R_MEN`, `R_D1/R_D2`, `R_LF`, the `C_0603`s, `C_LF`), and the Bourns land
  pattern on `D_BUS`.
- **The `_HandSolder` land is generator-only.** `gen/board.py` names
  `R_0603_1608Metric_Pad0.98x0.95mm_HandSolder` /
  `R_0805_2012Metric_Pad1.20x1.40mm_HandSolder` — the §11.1 lands, gap 0.85 / 0.80 mm —
  for those four references and **no committed artifact carries one** — `grep -c HandSolder` is 0 in `phantasm.kicad_sch`,
  `phantasm.kicad_pcb` and `unplaced/phantasm_unplaced.kicad_pcb`. Regenerating the
  schematic would put those ids into it and fail the schematic-parity gate against the
  routed copper. The widening is invisible to that gate in the other direction, because
  the routed board keeps the stock footprint id: KiCad's **Update Footprints from
  Library** would restore the 0.80/1.025 mm pads with no parity difference reported.
- **J1 ships unkeyed — R-PWR-7 is not met on this board.** Both committed artifacts
  carry `Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical`, a plain 0.1″
  header with no key, no shroud and no locking ramp, so nothing mechanically stops the
  +5 V/GND feed going on backwards — the condition R-PWR-7 calls out (“a reversed feed
  destroys the Teensy + '125”). `Q_REV` is series protection on the logic rail, not the
  polarising feature that requirement asks for. `gen/board.py` names the keyed
  `Connector_JST:JST_XA_B02B-XASK-1-A_1x02_P2.50mm_Vertical` (mates XAP-02V-1 +
  SXA-001T-P0.6); it reaches copper only after a re-place and a re-route, because the
  JST body and its 2.50 mm pitch do not fit the header's routed pads. Until then the
  harness carries the polarity marking. The JLC assembly gate cannot catch the
  substitution: `fab.EXCLUDE_FP_SUBSTR` excludes both `PinHeader` and `JST_` as
  hand-soldered, so `gen/tests/test_pcb_lands.py` pins J1's shipped footprint instead.
- **ID straps** use the Teensy's internal pull-ups; the former optional `R_ID0`
  footprint is not required. `D_BUS` is populated on every board. `JP_SHLD` is
  populated **on the master board only**;
  `JP_ID2` is unread at N = 4 and carries the high segment-ID bit at N = 8.
- **LED power is off-board (§2.3).** There is **no `C_BULK` and no `+5V_MAIN` heavy
  rail on the card** — the 1000 µF bulk lives at the strip injection point off-board
  (R-PWR-11), and `J2` carries **signal only** (DI/CI/SIG_GND, no +5 V). `C_IN`
  (≥100 µF) is the card's only electrolytic and sits on the post-bead `+5V_LOGIC` rail
  (R-PWR-3/6, §10).
- **Net naming.** The power chain is `+5V_IN` (J1↔F1), `+5V_RAW` (F1↔Q_REV),
  `+5V_PROT` (Q_REV↔FB), then `+5V_LOGIC` after the bead. `+5V_LOGIC` carries
  C_IN / R_LF / C_DEC / Teensy VIN / U1 Vcc per §10. The strip-return /
  logic-GND star (§R-SI-2) is a single `GND` net in the schematic — the
  load-end star tie is a **layout/harness** concern (SIG_GND meets the heavy
  LED return at the strip GND pin, off-board), not a separate schematic net.
- **Teensy symbol** shows only the **pins this board uses** (VIN, 3V3, GND, D1, D3,
  D5, D11, D13, D21, D22, **D23**); the other ~16 pads are unconnected on this design
  and omitted for readability. Pin **number = the Teensy pad label** (e.g. `11`, `VIN`),
  which matches the generated `phantasm:Teensy4.0` footprint pad names. The footprint
  pad map is the **top view (component side up) with the USB end at −X** — the Teensy
  mounts component-side-up.
- **Mechanical and service access.** The Quilter board retains the existing
  **58.28 × 32 mm** outline and adds four 2.7 mm NPTH M2.5 clearance holes at
  `(3.5, 3.5)`, `(3.5, 28.5)`, `(54.78, 3.5)`, and `(54.78, 28.5)` mm. Each hole
  has a 2.7 mm-radius all-copper routing/zone keepout. The Teensy footprint includes
  a board-envelope 3D model and an 11.5 × 10 mm mating-USB placement keepout; J1 and
  J4 sit below that approach corridor.
- **Identification.** Bottom silkscreen contains the N=4 ID0/ID1 truth table,
  `MASTER = ALL ID OPEN`, the shield instruction, and a large writable board-ID field.

## PCB (`phantasm.kicad_pcb`)

The PCB uses the corrected component-side-up Teensy footprint verified against
PJRC's top-view pinout. The earlier mirrored-footprint warning was obsolete and
has been removed. The completed control-net routing is included in this file.

The committed routed PCB is the source of truth for these facts. Refresh this block with
`python gen/board_metadata.py --write-readme` after an intentional board change.

<!-- BEGIN ROUTED PCB FACTS -->
<!-- Generated by `python gen/board_metadata.py --write-readme`; do not edit. -->
| Routed-board fact | Extracted value |
|---|---|
| Board dimensions | 58.28 × 32 mm |
| Board thickness | 1.6 mm |
| Footprints by side | 32 (F.Cu: 32, B.Cu: 0) |
| Track segments | 339 |
| Vias | 100 |
| Copper pours | 2 (In1.Cu: 1, In2.Cu: 1) |
| Keepout rule areas | 4 (F.Cu: 4, In1.Cu: 4, In2.Cu: 4, B.Cu: 4) |
| Copper layers | 4 (F.Cu, In1.Cu, In2.Cu, B.Cu) |
| Copper thicknesses | F.Cu: 0.035001 mm; In1.Cu: 0.015189 mm; In2.Cu: 0.015189 mm; B.Cu: 0.035001 mm |
| Copper finish | Lead-Free |
<!-- END ROUTED PCB FACTS -->

Connectors are at the **ends** (power/debug `J1`/`J4` at the hub end,
strip/sync `J2`/`J3A`/`J3B` at the far end, R-CON-4). `MASTER_EN` and
`SYNC_PULLDOWN` are fully routed.

### How it was placed & routed

- **Autoplaced + autorouted with Quilter.** Candidate 1 was subsequently corrected
  for mounting, USB access, identifiers, power protection, BOM metadata, and
  signal-integrity placement. That corrected board—not the old unplaced generator
  output—is now the layout source of truth.
- **Ground:** both inner layers have solid `GND` zones, providing an adjacent reference
  plane for traces on both outer layers (R-SI-1).
- **Fast nets:** DATA, CLK, and SYNC routing from the validated input was retained;
  the completed Quilter pass added only the low-rate control routing.
- Critical routing uses at least **0.13 mm trace width** with the JLCPCB
  **4 mil clearance** process limit.

> Routing lives only in `phantasm.kicad_pcb`. **`gen/pcb.py --force` overwrites the
> board and discards routing** — without `--force` it refuses while that file exists.
> Treat `quilter_incremental/` as a historical input snapshot, not as the
> fabrication board.

## Status

- [x] Schematic — complete, netlist verified against spec §10. ERC reported zero
  errors when run by hand; no runner re-checks it, so any schematic edit
  invalidates that result (see Validation)
- [x] Teensy footprint pad map verified for component-side-up mounting
- [x] Corrected Candidate 1 placement and validated routing preserved
- [x] Automatic master-only R_PD circuit added to the schematic and PCB netlist
- [x] Quilter control-net routing imported and verified with a clean DRC

### Layout constraint (R-MECH-6)
**Board width ≤ 35 mm** — mounts along the rotor arm. `PCB_W` is set to **32 mm**
(within the cap, trimmed to the part extent); the packer minimises the length (free)
dimension within that width. Narrowing the width lengthens the board (less room to
pack beside the Teensy). The committed routed dimensions are reported in the generated
facts block above.

## Regenerating

> **`gen/` does not reproduce any committed KiCad artifact.** KiCad wrote
> `phantasm.kicad_sch` and `unplaced/phantasm_unplaced.kicad_pcb` last, and
> `phantasm.kicad_pcb` is routed on top of that schematic. The committed schematic
> carries KiCad's random **v4** uuids while `kicad_common.uid()` emits deterministic
> **v5** ids, so a regeneration renumbers all 28 symbol uuids and dangles every
> `(path …)` link `phantasm.kicad_pcb` holds into the schematic. Even
> uuid-normalised, a fresh `board.py` run differs from the committed schematic by
> roughly 170 lines — KiCad file version, power-symbol annotation (`#PWR` vs
> `#PWR01`), `J1`'s keyed JST footprint, the four spec §11.1 hand-solder lands, and
> some label placement. `gen/` is the design *description* and the source of the
> gates above; it is not a rebuild of what ships.
>
> `board.py` and `pcb.py` refuse to overwrite a committed artifact without
> `--force`. Run the gates freely. Run a `--force` regeneration only when
> deliberately starting a new board revision, and expect to re-annotate, re-place
> and re-route from its output.

Gates — safe against the committed artifacts at any time:

```sh
cd gen
# point at your KiCad install if not the default Windows path:
# export KICAD_SYMBOL_DIR="/path/to/share/kicad/symbols"
# export KICAD_FOOTPRINT_DIR="/path/to/share/kicad/footprints"
python check.py          # gate: exact named-net partition and (ref, pin) nodes
python shorts.py         # union-find short check on the schematic
```

Regenerating for a new revision — **discards the committed design**:

```sh
cd gen
python board.py --force  # ../phantasm.kicad_{sch,sym} + sym-lib-table (.kicad_pro seeded only if absent)
python pcb.py --force    # ../phantasm.kicad_pcb (placed, unrouted) + phantasm.pretty + fp-lib-table
python pcb.py --unplaced --force # ../unplaced/phantasm_unplaced.kicad_pcb (footprints staged below outline, for Quilter)
```

The separate unplaced board for Quilter is configured as **4-layer SIG / GND / GND /
SIG**, 1.6 mm, with 1 oz outer copper, 0.5 oz inner copper, and ENIG. Those generator inputs do not describe the
committed routed board; use the generated facts block above for its construction. The
unplaced stackup is encoded in its file, so Quilter reads it on upload — no need to
hand-enter dielectric/mil values in its UI. Net class is 0.3 mm track / 0.2 mm clearance /
0.6 mm via (well above the 3.5 mil fab minimum).

**Running a future unplaced board through Quilter** — upload the contents of
`unplaced/` together: `phantasm_unplaced.kicad_pcb`,
`phantasm_unplaced.kicad_pro`, and `fp-lib-table`. `gen/pcb.py --unplaced
--force` regenerates the first and the third; the `.kicad_pro` is a captured
artifact (see the file table). There is no unplaced
schematic, and the `.kicad_prl` that the KiCad GUI writes is local only. Run
`python gen/heal_clearance.py` as the final preparation step. Quilter prep:
- **`min_clearance` must be > 0** in the uploaded `.kicad_pro` — Quilter rejects the KiCad
  default of 0 ("min clearance must be > 0"). KiCad **re-zeroes it whenever the project is
  opened in the GUI**, so it must be restored right before upload. `gen/heal_clearance.py`
  restores at least 0.1016 mm in routed projects and preserves the unplaced
  project's 0.2 mm setting.
- The unplaced project intentionally allows Quilter to place unlocked components;
  only explicit mechanical and signal-integrity placements remain locked.
- **Every footprint carries its schematic `(path)`** so Quilter matches board↔schematic
  (groups related parts during placement).
- Select **Preserve copper on internal layers** and preserve the uploaded
  four-layer stackup. The unplaced project starts with 0.60/0.30 mm vias and
  enforces a 0.45/0.20 mm minimum.
- After downloading candidates, run `python gen/analyze_candidates.py <paths>`.
  Candidates with vias below 0.45/0.20 mm are ineligible even if Quilter's DRC
  accepts them. After promoting the chosen board to `phantasm.kicad_pcb`,
  **re-widen the four sync lands by hand** (below), then run `python gen/fab.py`;
  fabrication output repeats the same via gate.

**Required manual step after promoting any re-route: re-widen the four sync lands.**
`unplaced/phantasm_unplaced.kicad_pcb` stages `R1`, `R2`, `R_PD` and `R_S` on the
**stock** library land, so every board that comes back from a placer carries stock pads
and no generator reapplies the widening. Edit the promoted `phantasm.kicad_pcb` in place —
keep the stock footprint id and the stock pad centres, grow only the pad `size`:

| Ref(s) | Footprint id (unchanged) | Pad centres | Pad size stock → shipped |
|---|---|---|---|
| `R1`, `R2`, `R_PD` | `Resistor_SMD:R_0603_1608Metric` | ±0.825 mm | 0.8 × 0.95 → **1.2 × 0.95 mm** |
| `R_S` | `Resistor_SMD:R_0805_2012Metric` | ±0.9125 mm | 1.025 × 1.4 → **1.4 × 1.4 mm** |

`gen/tests/test_pcb_lands.py` (`SHIPPED_CHIP_LANDS`) holds these numbers and fails until
the promoted board carries them; KiCad's DRC and schematic-parity gates cannot see the
difference, because the widened pads keep the stock footprint id.

> `python pcb.py --force` regenerates the placed board and **discards routing**;
> `--unplaced --force` discards the committed KiCad GUI re-save that is the Quilter
> upload input. `python board.py --force` discards the committed schematic.
> Rewriting `phantasm.pretty/Teensy4.0.kicad_mod` — which the routed board
> resolves its Teensy pads against — takes its own `--force-teensy-library`;
> `--force` does not authorize it.

`gen/sexp.py` is a small S-expression parser/serializer; `gen/builder.py` places
stock symbols and emits the `.kicad_sch` (placement transform calibrated against
`kicad-cli` netlist export); `gen/board.py` is the schematic description; `gen/pcb.py`
embeds footprints, assigns pad nets from the netlist, and lays out the board strip.
