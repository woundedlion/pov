# PHANTASM Segment Board — KiCad project

KiCad 10 schematic for the per-segment carrier board specified in
[../../docs/phantasm_pcb_spec.md](../../docs/phantasm_pcb_spec.md). One identical
PCB is built ×4; a solder strap selects each board's role (segment 0 =
master/conductor, 1–3 = flywheel slaves). The built config is **N = 4** (ID0/ID1 via
`JP_ID0` / `JP_ID1`); `JP_ID2` (pin 23) is broken out so a recompiled **N = 8** build
decodes segments 4–7 (R-ID-1). The card is **logic-only (~0.15 A)** — the strips'
4.3 A LED power is injected **off-board** (spec §2.3), so nothing here carries high
current.

Contains the **schematic** (electrical design), a **fully-routed PCB**, and an **unplaced
PCB** variant for autoplacement (Quilter) — see [Files](#files) and [Status](#status).

## Files

| File | What |
|---|---|
| `phantasm.kicad_pro` | Project file |
| `phantasm.kicad_sch` | Schematic — all parts, values, footprints, full §10 connectivity |
| `phantasm.kicad_pcb` | PCB — 53.8 × 32 mm outline, footprints placed + **fully routed** (4-layer SIG/GND/GND/SIG) |
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

- **ERC: 0 violations** (`kicad-cli sch erc --severity-error --severity-warning`).
- **Netlist matches spec §10** — exported with `kicad-cli sch export netlist` and
  diffed against the net table in the spec (see `gen/check.py`). Every net in §10 is
  realized with the correct members (logic feed `J1 → F1 → Q_REV → FB → +5V_LOGIC`;
  series terminations `U1 out → R → J2`/bus; the pin-3 divider node ties Teensy D3,
  `U1` ch-C input, `R1`/`R2`/`C_SYNC`; ID0/ID1/ID2 straps; `MASTER_EN`; shield).
- **PCB DRC: clean** apart from the expected unrouted ratsnest
  (`kicad-cli pcb drc`) — no shorts, clearance, or courtyard-overlap errors.

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
| `U_MCU` | `phantasm:Teensy4.0` | `phantasm:Teensy4.0` (2×14 0.1″ THT) | generated footprint; verify pad map vs your Teensy |
| `U1` (A–E) | `74xx:74AHCT125` | `Package_SO:SOIC-14_3.9x8.7mm_P1.27mm` | 4 buffers + power unit |
| `Q_REV` | `Device:D_Schottky` (B5819W) | `Diode_SMD:D_SOD-123` | series reverse-protect; pin2=anode (in), pin1=cathode (out) |
| `F1` | `Device:Fuse` (0.75 A) | `Fuse:Fuse_1206_3216Metric` | logic-feed overcurrent (~0.15 A draw) |
| `C_IN` | `Device:C_Polarized` (100µF) | `Capacitor_THT:CP_Radial_D8.0mm_P3.50mm` | only on-card electrolytic, on +5V_LOGIC; RTV-bond |
| `FB` | `Device:FerriteBead` | `Inductor_SMD:L_1206_3216Metric` | ~600 Ω @100 MHz |
| `R_LF` / `C_LF` | `Device:R` / `Device:C` | 0805 / `C_1206` | bead-LC damper, 22 µF |
| `C_DEC1/2` | `Device:C` (0.1µF) | `Capacitor_SMD:C_0603_1608Metric` | |
| `R_D1/R_D2` | `Device:R` (33Ω) | `Resistor_SMD:R_0805_2012Metric` | DATA/CLK source term |
| `R_S` | `Device:R` (100Ω) | `R_0805` | SYNC source term |
| `R1/R2` | `Device:R` (10k/15k) | `R_0603` | sync divider |
| `C_SYNC` | `Device:C` (220pF) | `C_0603` | populated (noise filter) |
| `R_PD` | `Device:R` (10k) | `R_0603` | bus idle pull-down |
| `R_MEN` | `Device:R` (10k) | `R_0603` | MASTER_EN boot pull-up → 3V3 |
| `R_ID0` | `Device:R` (10k) | `R_0603` | **DNP** — ID0 pull-up |
| `D_BUS` | `Device:D_TVS` | `Diode_SMD:D_SOD-323` | **DNP** — bus clamp |
| `J1` | `Connector_Generic:Conn_01x02` | `PinHeader_1x02_P2.54mm` | +5 V/GND light logic feed, ~1 A |
| `J2` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | strip **signal only**: DI / CI / SIG_GND (no power) |
| `J3A/J3B` | `Connector_Generic:Conn_01x03` | `PinHeader_1x03_P2.54mm` | Belden 8451 daisy |
| `J4` | `Connector_Generic:Conn_01x04` | `PinHeader_1x04_P2.54mm` | debug/serial |
| `JP_SHLD/JP_ID0/JP_ID1/JP_ID2` | `Jumper:SolderJumper_2_Open` | `SolderJumper-2_P1.3mm_Open_...` | shield (master only) / ID straps (JP_ID2 = N=8) |

## Notes / deviations from the spec

- **Reverse protection is a series Schottky** (`Q_REV`, SOD-123). §R-PWR-7 permits
  "a small series Schottky **or** P-FET" now that the 4.3 A LED current is off the card;
  at ~0.15 A the ideal-diode P-FET (and its gate resistor / SOT pinout ambiguity) is
  unnecessary, so the simpler Schottky is used. Ref kept as `Q_REV` per spec §9/§10.
  The ~0.3 V forward drop at 0.15 A leaves U1 Vcc ≈ 4.7 V, inside the AHCT 4.5–5.5 V
  window; swap to a small SOT-23 P-FET if you want the extra headroom.
- **DNP parts** (`R_ID0`, `D_BUS`) carry the KiCad DNP flag — footprints present,
  not populated, per §R-ASM-4. `JP_SHLD` is populated **on the master board only**;
  `JP_ID2` is left open at N = 4 (stuff for N = 8).
- **LED power is off-board (§2.3).** There is **no `C_BULK` and no `+5V_MAIN` heavy
  rail on the card** — the 1000 µF bulk lives at the strip injection point off-board
  (R-PWR-11), and `J2` carries **signal only** (DI/CI/SIG_GND, no +5 V). `C_IN`
  (≥100 µF) is the card's only electrolytic and sits on the post-bead `+5V_LOGIC` rail
  (R-PWR-3/6, §10).
- **Net naming.** Power chain is `J1 → F1 → Q_REV → FB → +5V_LOGIC`; `+5V_RAW` names
  the F1↔Q_REV node, `+5V_LOGIC` (post-bead) carries C_IN / R_LF / C_DEC / Teensy VIN /
  U1 Vcc per §10. The strip-return / logic-GND star (§R-SI-2) is a single `GND` net in
  the schematic — the load-end star tie is a **layout/harness** concern (SIG_GND meets
  the heavy LED return at the strip GND pin, off-board), not a separate schematic net.
- **Teensy symbol** shows only the **pins this board uses** (VIN, 3V3, GND, D1, D3,
  D5, D11, D13, D21, D22, **D23**); the other ~16 pads are unconnected on this design
  and omitted for readability. Pin **number = the Teensy pad label** (e.g. `11`, `VIN`),
  which matches the generated `phantasm:Teensy4.0` footprint pad names. Verify against
  your actual Teensy footprint before fabricating.

## PCB (`phantasm.kicad_pcb`)

A **fully-routed 4-layer board** (**SIG / GND / GND / SIG**, 1.6 mm): a **53.8 × 32 mm**
`Edge.Cuts` rectangle (trimmed from the 35 mm cap to the part extent — the bottom
strip held only ground pour) with all 29 footprints placed double-sided, connectors at the
**ends** (power/debug `J1`/`J4` at the hub end, strip/sync `J2`/`J3A`/`J3B` at the far
end, R-CON-4), and **every net routed** (350 track segments, 88 vias). `kicad-cli pcb drc`
reports **0 errors, 0 unconnected** (cosmetic silk-overlap warnings from the tight
auto-placement remain — nudge the reference designators in Pcbnew).

### How it was placed & routed

- **Autoplaced + autorouted with Quilter** from `unplaced/phantasm_unplaced.kicad_pcb` (4-layer
  stackup + inner GND planes encoded in-file; connectors pre-locked at the ends). Quilter
  returns several candidates; the chosen one was picked by a fast-net signal-integrity
  comparison (`gen/analyze_candidates.py` — fewest SPI/SYNC vias, shortest fast nets).
- **Ground:** both inner layers are solid `GND` planes, and `GND` is also poured on
  `F.Cu`+`B.Cu` (~77 % of each outer layer), so the `DATA` (F.Cu) and `CLK` (B.Cu) traces
  run with coplanar ground; all four copper layers tie through the board's GND vias (R-SI-1).
- **Fast nets:** `DATA`, `CLK`, `CLK_IN`, `SYNC_BUS`, `FRAME_SYNC` route **via-free**; the
  buffered strip `DATA` output is ~5.4 mm on `F.Cu`. Only `DATA_IN` (Teensy→buffer input)
  takes 2 vias, to cross `MASTER_EN` on `B.Cu` — negligible on the input side.
- Net class `Default` = **0.3 mm track / 0.2 mm clearance / 0.6 mm via**; min clearance
  0.2 mm. **1 oz copper suffices** (R-MECH-5).

### Remaining polish (optional, in Pcbnew)

- **Silk labels:** the auto-placement leaves ~14 reference-designator silk overlaps
  (`silk_overlap` / `silk_over_copper` / clipped-by-edge) — nudge the refdes for a clean
  legend. Cosmetic only; no electrical effect.
- Optional: drop `DATA_IN` to 1 via by rerouting `MASTER_EN` (interactive push-and-shove).

> Routing lives only in `phantasm.kicad_pcb`. **Re-running `gen/pcb.py` overwrites the
> board and discards routing** — don't regenerate the PCB after this point (or route a copy).
> For a fresh Quilter run, upload `unplaced/phantasm_unplaced.kicad_pcb` instead.

## Status

- [x] Schematic — complete, ERC-clean, netlist verified against spec §10
- [x] PCB — placed, net-assigned, **fully routed** (4-layer SIG/GND/GND/SIG, DRC-clean: 0 errors, 0 unconnected)
- [ ] Optional polish — silk-label nudge (above)

### Layout constraint (R-MECH-6)
**Board width ≤ 35 mm** — mounts along the rotor arm. `PCB_W` is set to **32 mm**
(within the cap, trimmed to the part extent); the packer minimises the length (free)
dimension within that width. Narrowing the width lengthens the board (less room to
pack beside the Teensy): the generated outline is **≈58 × 32 mm** at 32 mm wide
(was ≈54 × 35 at the 35 mm cap).

## Regenerating

The schematic and PCB are generated from the spec by the scripts in `gen/`:

```sh
cd gen
# point at your KiCad install if not the default Windows path:
# export KICAD_SYMBOL_DIR="/path/to/share/kicad/symbols"
# export KICAD_FOOTPRINT_DIR="/path/to/share/kicad/footprints"
python board.py          # ../phantasm.kicad_{sch,sym,pro} + sym-lib-table
python pcb.py            # ../phantasm.kicad_pcb (placed, unrouted) + phantasm.pretty + fp-lib-table
python pcb.py --unplaced # ../unplaced/phantasm_unplaced.kicad_pcb (footprints staged below outline, for Quilter)
"$KICAD/bin/python" stackup.py  # upgrade unplaced/ board to 4-layer SIG/GND/GND/SIG + stackup, heal min_clearance
python check.py          # gate: exported netlist == spec §10 (by ref membership)
python shorts.py         # union-find short check on the schematic
```

The unplaced board for Quilter is **4-layer SIG / GND / GND / SIG**, 1.6 mm: 1 oz copper,
~0.21 mm outer prepreg (tight F.Cu→GND coupling), 1.065 mm core, ENIG. Both inner layers
are `GND` planes (R-SI-1). The stackup is encoded in the file, so Quilter reads it on
upload — no need to hand-enter dielectric/mil values in its UI. Net class is 0.3 mm
track / 0.2 mm clearance / 0.6 mm via (well above the 3.5 mil fab minimum).

**Uploading to Quilter** — upload the whole project (`phantasm.kicad_pcb` *placed/routed*
or `unplaced/phantasm_unplaced.kicad_pcb` *for autoplace*, **plus `phantasm.kicad_sch` and
the matching `.kicad_pro`**). **Run `python gen/heal_clearance.py` as the LAST step before
every upload** (see next point). Quilter prep:
- **`min_clearance` must be > 0** in the uploaded `.kicad_pro` — Quilter rejects the KiCad
  default of 0 ("min clearance must be > 0"). KiCad **re-zeroes it whenever the project is
  opened in the GUI**, so it must be restored right before upload. `gen/heal_clearance.py`
  sets it to 0.2 mm in **both** project files (`phantasm.kicad_pro` and the unplaced one);
  `stackup.py` also heals the unplaced pro on every regen. The placed-board upload uses the
  **top-level** `phantasm.kicad_pro` — heal that one too (the common gotcha).
- **Through-hole parts + the tall cap are pre-placed and `locked` on the TOP side**
  (`U_MCU`, `C_IN`, and connectors `J1`–`J4`) so Quilter can't flip them to the bottom —
  keeps the board's back clear of bodies/the 10 mm electrolytic for the ring mount
  (R-MECH). Quilter places only the low-profile SMD (either side). Relock/move in Pcbnew
  to change fixed positions; for a fully bare back, also enable single-sided placement in
  Quilter (costs board area).
- **Every footprint carries its schematic `(path)`** so Quilter matches board↔schematic
  (groups related parts during placement).
- Quilter regenerates its own copper pours, so the inner `GND` zones here are just intent.

> `python pcb.py` (no flag) regenerates the placed board and **discards routing** — only
> rerun it before routing. `--unplaced` writes a *separate* file and never touches
> `phantasm.kicad_pcb`.

`gen/sexp.py` is a small S-expression parser/serializer; `gen/builder.py` places
stock symbols and emits the `.kicad_sch` (placement transform calibrated against
`kicad-cli` netlist export); `gen/board.py` is the schematic description; `gen/pcb.py`
embeds footprints, assigns pad nets from the netlist, and lays out the board strip.
