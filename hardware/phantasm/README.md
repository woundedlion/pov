# PHANTASM Segment Board — KiCad project

KiCad 10 schematic for the per-segment carrier board specified in
[../../docs/phantasm_pcb_spec.md](../../docs/phantasm_pcb_spec.md). One identical
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
| `phantasm.kicad_sch` | Schematic — all parts, values, footprints, full §10 connectivity |
| `phantasm.kicad_pcb` | Completed routed PCB with validated placement, control routing, planes, mounting, and service clearances |
| `quilter_incremental/` | Historical protected input snapshot used for the completed control-net routing |
| `unplaced/phantasm_unplaced.kicad_pcb` | **4-layer** (SIG/GND/GND/SIG) outline + net-assigned footprints **staged below the board, unrouted** — for an autoplacer (Quilter). Stackup encoded in-file. Regenerate: `python gen/pcb.py --unplaced` then `<kicad-python> gen/stackup.py` (both write into `unplaced/`) |
| `phantasm.kicad_sym` | Project symbol library: custom `Teensy4.0` + `+5V_RAW/+5V_LOGIC` power symbols |
| `phantasm.pretty/` | Project footprint library: generated `Teensy4.0` footprint (2×14 0.1″ THT) |
| `sym-lib-table` / `fp-lib-table` | Register the `phantasm` symbol / footprint libraries |
| `gen/` | The Python generators (schematic + PCB) — see [Regenerating](#regenerating) |

Open `phantasm.kicad_pro` in KiCad 10. Stock symbols/footprints come from the
standard KiCad libraries; the custom Teensy + power symbols and the Teensy
footprint come from the project `phantasm.kicad_sym` / `phantasm.pretty`.

## Validation

Both checks run via `kicad-cli` (KiCad 10.0):

- **N=8 firmware:** `pio run -e phantasm8` compiles and links the optional
  eight-board profile; this is firmware validation, not rotor qualification.
- **ERC: 0 errors.** A warning-inclusive KiCad 10 run reports nine
  `lib_symbol_mismatch` notices for embedded copies of stock/custom symbols;
  the exported connectivity is verified separately below.
- **Netlist matches spec §10** — exported with `kicad-cli sch export netlist` and
  diffed against the net table in the spec (see `gen/check.py`), keyed on `(ref, pin)`
  so a connector or IC pinout permutation fails. Every net in §10 is
  realized with the correct members (logic feed `J1 → F1 → Q_REV → FB → +5V_LOGIC`;
  series terminations `U1 out → R → J2`/bus; the pin-3 divider node ties Teensy D3,
  `U1` ch-C input, `R1`/`R2`/`C_SYNC`; ID0/ID1/ID2 straps; `MASTER_EN`; shield).
- **PCB geometry DRC: clean** (`kicad-cli pcb drc`): zero error-severity
  violations and zero unconnected pads.
- **Standard-cost via gate:** `gen/fab.py` rejects a routed board containing a
  via smaller than 0.45 mm with a drill smaller than 0.20 mm.
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
| `U_MCU` | `phantasm:Teensy4.0` | `phantasm:Teensy4.0` (2×14 0.1″ THT) | pad map = top view (component side up), USB end at −X: top row VIN,GND,3V3,23…13 / bottom row GND,0…12 |
| `U1` (A–E) | `74xx:74AHCT125` | `Package_SO:SOIC-14_3.9x8.7mm_P1.27mm` | 4 buffers + power unit |
| `Q_REV` | `Transistor_FET:Q_PMOS_GSD` (AO3401A) | `Package_TO_SOT_SMD:SOT-23` | reverse-polarity protection; pin 3 drain=input, pin 2 source=output, pin 1 gate=GND |
| `F1` | `Device:Fuse` (0.5 A hold) | `Fuse:Fuse_1206_3216Metric` | TLC-NSMD050 resettable fuse; 1 A trip, 0.75 Ω post-trip maximum |
| `C_IN` | `Device:C_Polarized` (100µF) | `Capacitor_THT:CP_Radial_D8.0mm_P3.50mm` | only on-card electrolytic, on +5V_LOGIC; RTV-bond |
| `FB` | `Device:FerriteBead` | `Inductor_SMD:L_1206_3216Metric` | ~600 Ω @100 MHz |
| `R_LF` / `C_LF` | `Device:R` / `Device:C` | 0805 / `C_1206` | bead-LC damper, 22 µF |
| `C_DEC1/2` | `Device:C` (0.1µF) | `Capacitor_SMD:C_0603_1608Metric` | |
| `R_D1/R_D2` | `Device:R` (33Ω) | `Resistor_SMD:R_0805_2012Metric` | DATA/CLK source term |
| `R_S` | `Device:R` (100Ω) | `R_0805_2012Metric_Pad1.20x1.40mm_HandSolder` | SYNC source term; hand-solder land |
| `R1/R2` | `Device:R` (10k/15k) | `R_0603_1608Metric_Pad0.98x0.95mm_HandSolder` | sync divider; hand-solder land |
| `C_SYNC` | `Device:C` (220pF) | `C_0603` | populated (noise filter) |
| `R_PD` | `Device:R` (10k) | `R_0603_1608Metric_Pad0.98x0.95mm_HandSolder` | master-only bus idle pull-down, hand-solder land; ground-side switched automatically by U1 channel D |
| `R_MEN` | `Device:R` (10k) | `R_0603` | MASTER_EN boot pull-up → 3V3 |
| `D_BUS` | `Device:D_TVS` (Bourns CDSOD323-T05L) | `Diode_SMD:D_SOD-323` with Bourns pad geometry | populated 5 V, 1 pF sync-bus TVS; exact Bourns land pattern; JLCPCB C1975255 |
| `J1` | `Connector_Generic:Conn_01x02` | `JST_XA_B02B-XASK-1-A_1x02_P2.50mm_Vertical` | +5 V/GND light logic feed, ~1 A; keyed shrouded 3 A JST XA with a locking ramp + mounting boss (R-PWR-7); mates XAP-02V-1 + SXA-001T-P0.6 contacts |
| `J2` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | strip **signal only**: DI / CI / SIG_GND (no power) |
| `J3A/J3B` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | Belden 8451 daisy |
| `J4` | `Connector_Generic:Conn_01x04` | `PinHeader_1x04_P2.54mm` | debug/serial |
| `H1`–`H4` | — | `MountingHole:MountingHole_2.7mm_M2.5` | four NPTH rotor mounting holes |
| `JP_SHLD/JP_ID0/JP_ID1/JP_ID2` | `Jumper:SolderJumper_2_Open` | `SolderJumper-2_P1.3mm_Open_...` | shield (master only) / ID straps (JP_ID2 read at N=8) |

## Notes / deviations from the spec

- **Reverse protection uses one AO3401A P-channel MOSFET** (`Q_REV`, SOT-23), with
  its gate tied directly to GND. It replaces the series Schottky without adding a
  gate resistor or increasing the component count. At a 4.75 V J1 input and 0.15 A,
  a conservative hot calculation uses the fuse's 0.75 Ω post-trip maximum, the
  bead's 0.20 Ω maximum DCR, and twice the MOSFET's 60 mΩ maximum at −4.5 V:
  `VLOGIC = 4.75 − 0.15 × (0.75 + 0.20 + 0.12) = 4.5895 V`. This leaves about
  90 mV above the AHCT125's 4.5 V minimum; verify that J1 itself remains at or above
  4.75 V on the hot, operating rotor because external harness drop is not included.
- **Hand-solder lands on the bench-tuned sync resistors.** `R1`, `R2` (divider ratio,
  spec §4.2), `R_PD` (bus idle pull-down) and `R_S` (source termination) carry KiCad's
  `_HandSolder` land: the pad centre moves outward with the extra length, so the toe
  grows while the inter-pad gap stays at the IPC-nominal **0.85 mm** (0603) /
  **0.80 mm** (0805) and no copper reaches under the ceramic body. Every other chip
  passive uses the stock IPC-nominal land; `D_BUS` uses the Bourns land pattern.
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
| Copper zones | 6 (F.Cu: 4, In1.Cu: 5, In2.Cu: 5, B.Cu: 4) |
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

- [x] Schematic — complete, zero-error ERC, netlist verified against spec §10
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

The schematic and PCB are generated from the spec by the scripts in `gen/`:

```sh
cd gen
# point at your KiCad install if not the default Windows path:
# export KICAD_SYMBOL_DIR="/path/to/share/kicad/symbols"
# export KICAD_FOOTPRINT_DIR="/path/to/share/kicad/footprints"
python board.py          # ../phantasm.kicad_{sch,sym} + sym-lib-table (.kicad_pro seeded only if absent)
python pcb.py            # ../phantasm.kicad_pcb (placed, unrouted) + phantasm.pretty + fp-lib-table
python pcb.py --unplaced # ../unplaced/phantasm_unplaced.kicad_pcb (footprints staged below outline, for Quilter)
"$KICAD/bin/python" stackup.py  # upgrade unplaced/ board to 4-layer SIG/GND/GND/SIG + stackup, heal min_clearance
python check.py          # gate: exported netlist == spec §10 (by (ref, pin) node)
python shorts.py         # union-find short check on the schematic
```

The separate unplaced board for Quilter is configured as **4-layer SIG / GND / GND /
SIG**, 1.6 mm, with 1 oz outer copper, 0.5 oz inner copper, and ENIG. Those generator inputs do not describe the
committed routed board; use the generated facts block above for its construction. The
unplaced stackup is encoded in its file, so Quilter reads it on upload — no need to
hand-enter dielectric/mil values in its UI. Net class is 0.3 mm track / 0.2 mm clearance /
0.6 mm via (well above the 3.5 mil fab minimum).

**Running a future unplaced board through Quilter** — regenerate and upload the
three `unplaced/phantasm_unplaced.kicad_*` files together. Run
`python gen/heal_clearance.py` as the final preparation step. Quilter prep:
- **`min_clearance` must be > 0** in the uploaded `.kicad_pro` — Quilter rejects the KiCad
  default of 0 ("min clearance must be > 0"). KiCad **re-zeroes it whenever the project is
  opened in the GUI**, so it must be restored right before upload. `gen/heal_clearance.py`
  restores at least 0.1016 mm in routed projects and preserves the unplaced
  project's 0.2 mm setting; `stackup.py` also heals the unplaced project on
  every regeneration.
- The unplaced project intentionally allows Quilter to place unlocked components;
  only explicit mechanical and signal-integrity placements remain locked.
- **Every footprint carries its schematic `(path)`** so Quilter matches board↔schematic
  (groups related parts during placement).
- Select **Preserve copper on internal layers** and preserve the uploaded
  four-layer stackup. The unplaced project starts with 0.60/0.30 mm vias and
  enforces a 0.45/0.20 mm minimum.
- After downloading candidates, run `python gen/analyze_candidates.py <paths>`.
  Candidates with vias below 0.45/0.20 mm are ineligible even if Quilter's DRC
  accepts them. After promoting the chosen board to `phantasm.kicad_pcb`, run
  `python gen/fab.py`; fabrication output repeats the same via gate.

> `python pcb.py --force` regenerates the placed board and **discards routing**; without
> `--force` it refuses to touch the committed `phantasm.kicad_pcb`. `--unplaced` writes a
> *separate* file and needs no flag. `python board.py` guards `phantasm.kicad_sch` the
> same way.

`gen/sexp.py` is a small S-expression parser/serializer; `gen/builder.py` places
stock symbols and emits the `.kicad_sch` (placement transform calibrated against
`kicad-cli` netlist export); `gen/board.py` is the schematic description; `gen/pcb.py`
embeds footprints, assigns pad nets from the netlist, and lays out the board strip.
