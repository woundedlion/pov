# PHANTASM revision 1.1 - routed project and technical reference

This directory contains the routed rev 1.1 project. For the current rev 1.2
Quilter upload and regeneration commands, see the [project index](../README.md).
Paths below are relative to this directory unless a command explicitly runs
from the repository root.

**Routed review candidate: rev 1.1 sync circuit.** Rev 1.2
separates sync transmit from the filtered receive node; see
[Revision 1.2](#revision-12). The routed candidate includes corrected drill
spacing, thermal reliefs, an 8 V TVS, and J2 pin order **DATA / GND / CLK**.
Its cable pinout differs from the earlier DATA / CLK / GND board; repin the
harness before using this candidate. These edited files are not a record of
the previously ordered package; retain that package and its original digests.

KiCad 10 schematic for the per-segment carrier board specified in
[../../../docs/specs/phantasm_pcb_spec.md](../../../docs/specs/phantasm_pcb_spec.md). One identical
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
| `phantasm.kicad_sch` | Schematic — all parts, values, footprints, full §10 connectivity. **Last written by KiCad, not by `../gen/board.py`**, and not reproducible from it — see [Regenerating](#regenerating) |
| `phantasm.kicad_pcb` | Completed routed PCB with validated placement, control routing, planes, mounting, and service clearances |
| `phantasm.kicad_sym` | Project symbol library: custom `Teensy4.0` + `+5V_RAW/+5V_LOGIC` power symbols |
| `phantasm.pretty/` | Project footprint library: generated `Teensy4.0` footprint (2×14 0.1″ THT) |
| `sym-lib-table` / `fp-lib-table` | Register the `phantasm` symbol / footprint libraries |
| `../gen/` | The Python generators (schematic + PCB) — see [Regenerating](#regenerating) |

Open `phantasm.kicad_pro` in KiCad 10. Stock symbols/footprints come from the
standard KiCad libraries; the custom Teensy + power symbols and the Teensy
footprint come from the project `phantasm.kicad_sym` / `phantasm.pretty`.

## DFM defaults and corrected geometry

`../gen/pcb.py` emits the corrected routed board's solder-mask settings: a
**0.10 mm minimum mask web**, no mask bridges within footprints, and vias
**tented on both sides**. Pad mask expansion remains **0 mm**, matching the
accepted fabrication source. `../gen/fab.py` checks these settings before any
export and rejects individual vias that disable tenting.

New projects from `../gen/board.py` also set **0.15 mm silkscreen clearance** and
**0.10 mm solder-mask-to-copper clearance**, with silk over exposed copper an
error. The July 28 JLCDFM reports show the silk-to-pad minimum improving from
0.03 to 0.17 mm and the mask-to-trace minimum from 0.07 to 0.10 mm. Deliberate
schematic or PCB regeneration merges these floors into its output project,
preserving stricter values and unrelated settings. Unplaced generation also
writes or updates its own project with the wider unplaced constraints. Reading
or exporting the accepted board does not change its project. The shared rule
checks and `heal_clearance.py`
preserve its **0.1016 mm hole clearance** floor; the unplaced project's floor
is **0.25 mm**. Both projects require **0.5 mm hole-to-hole clearance**.
The routed board's six affected ground vias are re-spaced, with their connected
trace and zone fills updated; DRC reports no hole-spacing violations.

The corrected placement and routing remain in `phantasm.kicad_pcb`. Start a new
routing job with the current [rev 1.2 project](../README.md). A fresh
`../gen/pcb.py` run creates a placement draft: it reserves space for reference
labels, keeps the back legend between the Teensy's pad rows, and separates
the connector labels from their outlines. It packs every part, including terminals
at the hub and far ends. With `--unplaced`, terminal blocks are locked at
`TERMINAL_EDGE_PLACEMENTS`, the R1/R2/C_SYNC receive filter is locked at
`SYNC_FILTER_PLACEMENTS`. It also locks U_MCU/C_DEC1 shifted 3.5 mm left,
U1/C_DEC2, R_D1/R_D2, C_IN and the four ID/SHLD solder jumpers. Remaining
parts (including R_S/R_PD hand-solder lands) are staged below the outline. It does
not reconstruct routing or the four widened resistor lands. Its clearances
still need checking after placement and routing. The final saved JLCDFM
report retains pad-spacing dangers and annular-ring, mask-expansion, and
silkscreen warnings, so reproducing the accepted board is not a claim of
zero JLCDFM findings. Run KiCad DRC and JLCDFM on every new routed package.

## Validation

The `pcb-tests` CI job runs the generator suite under pinned KiCad 10.0.4,
including schematic generation, committed-board parity and DRC, and
generated-board DRC. Metadata checks also run in CI. ERC and fabrication
exports run locally through `../gen/fab.py` when the package is regenerated.

- **N=8 firmware:** `pio run -e phantasm8` compiles and links the optional
  eight-board profile; this is firmware validation, not rotor qualification.
- **ERC: no error-severity violations** — `../gen/fab.py` runs the pinned
  `kicad-cli sch erc` before producing fabrication outputs, then validates every
  sheet's JSON violation list in `../gen/out/phantasm-erc.json`. A missing or malformed
  report, nonzero tool exit, or any reported violation stops fabrication.
  Power symbols have empty footprints; only `U_MCU` carries the Teensy land.
- **Netlist matches the electrical specification** — `../gen/fab.py` holds the
  netlist it exports to the named-net table in `../gen/check.py`, which also runs
  standalone against the committed schematic: every net in the table must match
  member-for-member, keyed on
  `(ref, pin)` so a connector or IC pinout permutation fails. A named net outside
  the table is reported as a `NOTE` and does not fail the gate. CI enforces the
  same table without KiCad: `../gen/tests/test_check.py` applies it to the pad nets
  of the committed board, so copper that stops matching the spec fails on every
  push.
  The required connections are realized with the correct members
  (logic feed `J1 → F1 → Q_REV → FB → +5V_LOGIC`;
  series terminations `U1 out → R → J2`/bus; the pin-3 divider node ties Teensy D3,
  `R1`/`R2`/`C_SYNC`, plus `U1` ch-C input on rev 1.1; rev 1.2 drives ch-C from
  D4/`SYNC_TX` with `R_TX`; ID0/ID1/ID2 straps; `MASTER_EN`; shield).
- **Copper-connectivity gate:** `../gen/connectivity.py` unions the routed board's
  tracks, vias, pads and pour fills per net and rejects any net whose pads land
  in more than one island. Every other net gate in this list reads pad net
  *attributes*, which survive deleting the copper that realizes them; this one
  walks the copper. It needs no KiCad and runs in CI via
  `../gen/tests/test_connectivity.py`; run it standalone with
  `python ../gen/connectivity.py`.
- **PCB geometry DRC: clean** (`kicad-cli pcb drc`): zero error-severity
  violations and zero unconnected pads.
- **Standard-cost via gate:** `../gen/fab.py` rejects a routed board containing a
  via smaller than 0.45 mm or a drill smaller than 0.20 mm, and rejects any
  different-net via pair with less than 0.15 mm of copper spacing (pad edge to pad edge).
- **Schematic parity gate:** `../gen/fab.py` runs `kicad-cli pcb drc
  --schematic-parity` and rejects any board/schematic difference outside
  `KNOWN_PARITY_ITEMS` (the four mounting holes, which carry no symbol), so
  gerbers for stale copper cannot ship with a BOM exported from a newer
  schematic. `KNOWN_PARITY_WARNING_COUNTS` allows exactly eleven
  `lib_footprint_mismatch` warnings and rejects any other count. Parity items
  are warning severity in KiCad, so this runs separately from the
  error-severity DRC gate.
- **Plot-origin gate:** the gerbers, Excellon drill, and CPL are all exported
  in absolute board coordinates, and `../gen/fab.py` rejects a board carrying a
  non-zero drill/place origin (`aux_axis_origin`) — that would move only the
  origin-relative exports and place every assembled part off-board.
- **Zone-geometry gate:** `../gen/fab.py` rejects a copper pour whose
  `min_thickness` is below 0.13 mm, `thermal_bridge_width` below 0.4 mm,
  or `thermal_gap` below 0.3 mm. Project rules require two resolved spokes.
  The gap and spoke floors include etch margin. KiCad clearance DRC never flags
  these — thermal reliefs are same-net geometry — so a sub-process gap would
  export clean gerbers the fab fills as a solid pour, tying every
  through-hole GND pad to the full plane and starving the hand-soldered
  joints (R-ASM-6).
- **Shipped-land gate:** `../gen/tests/test_pcb_lands.py` pins the pad geometry the
  routed board ships for every chip passive, so restoring a stock land over the
  widened sync-resistor pads fails in CI. It pins as-built geometry, not spec
  §11.1 geometry — see the lands note below. KiCad's `lib_footprint_mismatch` warning count also detects
  that edit — the widened pads keep the stock footprint id. The same file pins
  J1's footprint on both artifacts, which the assembly gate cannot see because
  `EXCLUDE_FP_SUBSTR` excludes every hand-soldered connector spelling.
- **Count-floor ratchets:** `../gen/fab.py` rejects a board carrying fewer than
  99 vias or fewer than 2 copper pours — floors pinned to what the committed
  routed board holds (both counts are in the facts block below; the pours are
  the In1/In2 reference planes, and the mounting-hole keepout rule areas are
  counted separately because they pour no copper). Neither floor updates
  itself: after promoting a re-route, re-measure with `../gen/board_metadata.py`
  and deliberately re-baseline the constants in `../gen/fab.py`.
- **Assembly-metadata gate:** `../gen/fab.py` rejects the assembly outputs
  unless the assembled references exactly match its LCSC assignment table,
  every rotation correction names an assembled part, and every centroid row
  is present, numeric, and on the top (assembly) side.
- **Board-generation gate:** `../gen/tests/test_pcb_generation.py` runs
  `../gen/pcb.py` into a temporary directory and reads the result back through the
  KiCad-free readers: the net table, both inner reference planes, the four
  mounting-hole rule areas, the revision silkscreen, the Teensy footprint
  library, the unmatched-netlist-pin refusal, and a courtyard-overlap check over
  the placed draft. The generator needs KiCad's stock libraries and a
  `kicad-cli` on the pin for the netlist export, so the class skips without them.
- **Layer-name gate:** `../gen/tests/test_pcb_stack.py` rejects a declared layer
  name carrying whitespace on either board. KiCad builds each Gerber's filename
  from the layer name, so an Altium-style alias such as `Ground Layer 1` ships a
  space in the upload zip.
- **Board-revision gate:** `../gen/tests/test_revision.py` reads the revision off
  the bottom silkscreen (`Phantasm Rev 1.1`) and requires the schematic title
  block and the routed board's title block to agree with it. Generated artifacts
  are checked separately against `../gen/builder.py`'s `REVISION` (1.2).
  `pcb.py` rejects a source schematic from another revision, and the electrical
  and assembly gates select their exact requirements by artifact revision.
  The board title block is what KiCad writes into the Gerber X2
  `ProjectId` attribute; without it the gerbers ship `rev?`.
- **Part-catalog gate:** every assigned LCSC number must resolve to a catalog
  entry with a non-blank manufacturer, MPN, and description, so each JLCPCB
  BOM match is independently auditable.
- **Export-content gate:** the upload zip is assembled by filename, so
  `../gen/fab.py` also reads the exported bytes before packaging them: every
  plotted Gerber must draw with an aperture or a filled region (only
  `phantasm-B_Paste.gbp` may be empty — assembly is top-side only), and
  the Excellon files must drill exactly the holes the board carries, counted
  off the board rather than pinned to a number: every via plus every plated pad
  hole in `phantasm-PTH.drl`, the unplated mounting holes in
  `phantasm-NPTH.drl`.
- **Fab-package digest gate:** the exports are stamp-normalized and zipped with
  fixed member metadata, so an unchanged board repackages byte for byte, and
  every run writes `SHA256SUMS.txt` beside the upload zip covering each zipped
  artifact, both assembly CSVs, and the zip itself. `python ../gen/fab.py --verify` re-hashes an
  already-generated package against that manifest and against
  `fab-SHA256SUMS.txt` — the baseline recording what was ordered, kept beside
  the board because `../gen/out/` is gitignored. That baseline is **not yet
  recorded**. Verify the retained package against its manifest and confirm that
  its Gerber zip, BOM, and CPL are the files submitted for the order before
  copying its manifest to `fab-SHA256SUMS.txt`. A fresh export from the current
  board does not establish the bytes used for an earlier order.

## How connectivity is drawn

The sheet is organised left-to-right into labelled blocks (power entry, Teensy +
level shifter → LED strip, sync RX divider + daisy, ID straps / debug, power flags):

- **Power distribution is drawn as visible rails** — a horizontal `+5V_LOGIC` rail
  and `GND` rail with vertical component drops and junctions (the classic ladder); the
  `J1 → F1 → Q_REV → FB` protection/filter chain feeds the rail's left (hub) end.
- **Series/divider paths are wired** — `U1` outputs through the 33 Ω/100 Ω source
  terminators, and the pin-3 RC divider (`R1`/`R2`/`C_SYNC`).
- **Cross-block signals use net labels as ports** — `DATA(_IN)`, `CLK(_IN)`,
  `FRAME_SYNC`, `SYNC_TX` (rev 1.2), `MASTER_EN`, `SYNC_BUS`, `ID0/1`, `SHIELD` — the conventional way to
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
| `R_TX` | `Device:R` (10k) | `Resistor_SMD:R_0603_1608Metric` | rev 1.2 SYNC_TX pull-down |
| `R_MEN` | `Device:R` (10k) | `R_0603` | MASTER_EN boot pull-up → 3V3 |
| `D_BUS` | `Device:D_Zener` (Bourns CDSOD323-T08L) | `Diode_SMD:D_SOD-323` with Bourns pad geometry | populated unidirectional 8 V, 1 pF sync-bus TVS; pin 1/cathode on SYNC_BUS, pin 2/anode on GND; exact Bourns land pattern; silkscreen bar marks the cathode end; JLCPCB C1973344 |
| `J1` | `Connector_Generic:Conn_01x02` | `Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical` | +5 V/GND light logic feed, ~1 A; **unkeyed** 0.1″ header — R-PWR-7's keying is unmet on the shipped board (see the deviations note below) |
| `J2` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | strip **signal only**: DI / SIG_GND / CI (no power) |
| `J3A/J3B` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | Belden 8451 daisy |
| `J4` | `Connector_Generic:Conn_01x04` | `PinHeader_1x04_P2.54mm` | debug/serial |
| `H1`–`H4` | — | `MountingHole:MountingHole_2.7mm_M2.5` | four NPTH rotor mounting holes |
| `JP_SHLD/JP_ID0/JP_ID1/JP_ID2` | `Jumper:SolderJumper_2_Open` | `SolderJumper-2_P1.3mm_Open_...` | shield (master only) / ID straps (JP_ID2 read at N=8) |

## Assembly polarity review

`D_BUS` is polarized: its cathode band must face pad 1 (`SYNC_BUS`), aligned
with the silkscreen bar. On the routed board this is the end **away from
`R_PD`, toward the board edge**; pad 2 (anode) connects to `GND`.
`../gen/fab.py` exports its committed 90-degree placement without a rotation
correction. Check the band against pad 1 in the assembly house's final preview
and on the assembled board; the exported angle alone does not establish the
supplier model's polarity. A reversed `D_BUS` forward-biases during a HIGH
sync pulse and clamps `SYNC_BUS`.

## Hand assembly (not done by the PCBA house)

The assembly house reflows top-side SMD only; `../gen/fab.py` excludes the Teensy,
the connectors, the electrolytic, and the solder jumpers from its BOM/CPL, so
every step below is performed by whoever builds the card.

- **R-ASM-7 — cut the VIN/VUSB pad on every Teensy 4.0. Mandatory, before the
  Teensy is soldered down.** The board feeds Teensy `VIN` from the rotor rail and
  `J4.4` reserves Teensy pin 1 for future UART debug; current diagnostics use
  USB `Serial`, and no `Serial1` driver is enabled. With the pad intact,
  `VUSB` is tied to `VIN`:
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

## Hand-build BOM (the through-hole half)

`../gen/fab.py` emits `phantasm-BOM.csv` for the **assembled SMD** parts only:
`fab.EXCLUDE_FP_SUBSTR` / `EXCLUDE_VAL_SUBSTR` drop every hand-soldered part, and
`fab.PART_BY_LCSC` carries manufacturer + MPN for the assembled half alone. Nothing
below reaches either. **The part numbers are not sourced** — the last column records
what the netlist, the footprint and spec §9 constrain, not a catalogue match; fill it
in before ordering.

| Ref(s) | Qty | What | Footprint / mating requirement | Orderable part |
|---|---|---|---|---|
| `U_MCU` | 1 | Teensy 4.0 (i.MX RT1062) development board | `phantasm:Teensy4.0` — 2×14 0.1″ THT, mounted component-side up, USB end at −X | PJRC **Teensy 4.0**; no distributor SKU pinned |
| `C_IN` | 1 | ≥100 µF radial aluminium electrolytic on `+5V_LOGIC` (spec §9); the card's only electrolytic, RTV-bonded | `Capacitor_THT:CP_Radial_D8.0mm_P3.50mm` — 8.0 mm body, 3.50 mm lead pitch | **unsourced** |
| `J1` | 1 | 2-pin 0.1″ vertical pin header — `+5V_IN` / `GND`, ~1 A | `Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical`; ships **unkeyed** (see the deviations note) | **unsourced** |
| `J2` | 1 | 3-pin 0.1″ vertical pin header — strip signal `DI` / SIG_GND / `CI`, no power | `PinHeader_1x03_P2.54mm_Vertical` | **unsourced** |
| `J3A`, `J3B` | 2 | 3-pin 0.1″ vertical pin headers — SYNC daisy in / out | `PinHeader_1x03_P2.54mm_Vertical`; one Belden 8451 run per link | **unsourced** |
| `J4` | 1 | 4-pin 0.1″ vertical pin header — debug (`+3V3`, `GND`, `MASTER_EN`, `SERIAL1_TX`) | `PinHeader_1x04_P2.54mm_Vertical` | **unsourced** |
| `JP_ID0/1/2`, `JP_SHLD` | 4 | **No part to order** — open solder-bridge pads, closed with solder per board role | `Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm` | — |
| `H1`–`H4` | 4 | M2.5 rotor mounting hardware (screw + nut or standoff) | 2.7 mm NPTH, centres 3.5 mm in from each corner; 5.4 mm square all-copper keepout | **unsourced** |

Off-board and out of this table (spec §2.3 / §9): the 1000 µF `C_BULK` injection bulk at
the strip, the heavy 5 V/GND LED harness, and the Belden 8451 STP for each inter-board run.

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
- **The `_HandSolder` land belongs to rev 1.2.** `../gen/board.py` names
  `R_0603_1608Metric_Pad0.98x0.95mm_HandSolder` /
  `R_0805_2012Metric_Pad1.20x1.40mm_HandSolder` — the §11.1 lands, gap 0.85 / 0.80 mm —
  for those four references and **no rev 1.1 artifact carries one** — `grep -c HandSolder` is 0 in `phantasm.kicad_sch`,
  `phantasm.kicad_pcb`. Regenerating the
  schematic would put those ids into it and fail the schematic-parity gate against the
  routed copper. The routed board keeps the stock footprint id, but restoring stock pads with
  KiCad's **Update Footprints from Library** changes the `lib_footprint_mismatch`
  warning count and fails the land-edit gate.
- **J1 ships unkeyed — R-PWR-7 is not met on this board.** Both committed artifacts
  carry `Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical`, a plain 0.1″
  header with no key, no shroud and no locking ramp, so nothing mechanically stops the
  +5 V/GND feed going on backwards. `Q_REV` is oriented to block that reversed feed from
  the logic rail, so it protects the Teensy and '125. The residual fault is on the return:
  reversal puts +5 V on the board's GND plane. USB ground, J2 pin 2 (SIG_GND), and
  J3A/J3B pin 2 provide return paths through the USB cable, strip ground lead, and
  22 AWG sync conductor. `F1` is only in J1's +5 V leg, so these fault paths are unfused;
  cutting VIN/VUSB per R-ASM-7 does not disconnect them. Verify J1 polarity before
  connecting it on any harnessed card.
  `../gen/board.py` selects the GCT TBC05-02-1-G-G screw terminal for J1.
  It is also unkeyed: R-PWR-7 still requires polarity verification. Its full body
  reservation and 1.3 mm drills apply to newly generated layouts; they do not
  repair the routed board. The assembly gate excludes headers and terminal blocks as
  hand-soldered, so `../gen/tests/test_pcb_lands.py` pins J1's shipped footprint instead.
- **ID straps** use the Teensy's internal pull-ups.
  `D_BUS` is populated on every board. `JP_SHLD` is
  populated **on the master board only**;
  `JP_ID2` is unread at N = 4 and carries the high segment-ID bit at N = 8.
- **LED power is off-board (§2.3).** There is **no `C_BULK` and no `+5V_MAIN` heavy
  rail on the card** — the 1000 µF bulk lives at the strip injection point off-board
  (R-PWR-11), and `J2` carries **signal only** (DI/SIG_GND/CI, no +5 V). `C_IN`
  (≥100 µF) is the card's only electrolytic and sits on the post-bead `+5V_LOGIC` rail
  (R-PWR-3/6, §10).
- **Net naming.** The power chain is `+5V_IN` (J1↔F1), `+5V_RAW` (F1↔Q_REV),
  `+5V_PROT` (Q_REV↔FB), then `+5V_LOGIC` after the bead. `+5V_LOGIC` carries
  C_IN / R_LF / C_DEC / Teensy VIN / U1 Vcc per §10. The strip-return /
  logic-GND star (§R-SI-2) is a single `GND` net in the schematic — the
  load-end star tie is a **layout/harness** concern (SIG_GND meets the heavy
  LED return at the strip GND pin, off-board), not a separate schematic net.
- **Rev 1.1 Teensy symbol** shows only the pins this board uses; its unused
  pads are omitted. The rev 1.2 symbol carries every pad. Pin **number = the Teensy pad label** (e.g. `11`, `VIN`),
  which matches the generated `phantasm:Teensy4.0` footprint pad names. The footprint
  pad map is the **top view (component side up) with the USB end at −X** — the Teensy
  mounts component-side-up.
- **Mechanical and service access.** The Quilter board retains the existing
  **58.28 × 32 mm** outline and adds four 2.7 mm NPTH M2.5 clearance holes at
  `(3.5, 3.5)`, `(3.5, 28.5)`, `(54.78, 3.5)`, and `(54.78, 28.5)` mm. Each hole
  has a 5.4 mm square all-copper routing/zone keepout. The Teensy footprint includes
  a board-envelope 3D model and an 11.5 × 10 mm mating-USB placement keepout; J1 and
  J4 sit below that approach corridor.
- **Identification.** Bottom silkscreen carries the N=4 ID0/ID1 truth table (master
  is the all-open row), the N=8 extension and shield line
  `N8 ID2 OPEN=0-3 GND=4-7; M=OPEN; SHLD=M`, the revision stamp, and a writable
  board-ID field.

## PCB (`phantasm.kicad_pcb`)

The PCB uses a component-side-up Teensy footprint verified against
PJRC's top-view pinout and includes the completed control-net routing.

The committed routed PCB is the source of truth for these facts. Refresh this block with
`python ../gen/board_metadata.py --write-readme` after an intentional board change.

<!-- BEGIN ROUTED PCB FACTS -->
<!-- Generated by `python ../gen/board_metadata.py --write-readme`; do not edit. -->
| Routed-board fact | Extracted value |
|---|---|
| Board dimensions | 58.28 × 32 mm |
| Board thickness | 1.6 mm |
| Footprints by side | 32 (F.Cu: 32, B.Cu: 0) |
| Track segments | 339 |
| Vias | 99 |
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
  **0.1016 mm (4 mil) clearance** process limit.

> The rev 1.1 routing lives in `phantasm.kicad_pcb`. Current generators write
> the separate rev 1.2 project in `../1.2/`.


## Status

- [x] Schematic — complete, netlist verified against spec §10. ERC reported zero
  errors; `../gen/fab.py` re-checks ERC on every fabrication export (see Validation)
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

## Revision 1.2

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

**Firmware compatibility:** build rev 1.2 with `HS_PHANTASM_BOARD_REV=12`;
`platformio.ini` defaults to 11. The rev 1.2 pin map transmits on D4, configures
D3 as input with HYS, and initializes D4 LOW as an output on every board before
enabling the master. Followers keep D4 LOW; only the master emits pulses.
Regenerating the board does not change the firmware build flags or routed copper.

Generating rev 1.2 requires regenerating the schematic before either PCB draft.
The resulting PCB is unrouted and needs placement/routing validation before
fabrication; changing its revision label alone does not convert rev 1.1 copper.

## Regenerating

The generators do not reproduce this routed rev 1.1 project. KiCad wrote its
schematic last; regenerating symbol UUIDs would break the routed board's
schematic links. Preserve these files when working on new revisions.

Current generators write a matching rev 1.2 project to `../1.2/`.
Follow the [project index](../README.md) for generation and Quilter upload
instructions.

To check the accepted rev 1.1 schematic, run from the repository root:

```sh
python hardware/phantasm/gen/check.py hardware/phantasm/1.1/phantasm.kicad_sch
python hardware/phantasm/gen/shorts.py hardware/phantasm/1.1/phantasm.kicad_sch
```
