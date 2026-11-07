# PHANTASM board projects

Use **[rev 1.2](1.2/README.md)** for the current Quilter placement run.
The experimental **[1.3/](1.3/README.md)** project uses differential sync.
Rev 1.2 remains the generator default; select rev 1.3 explicitly.
Open either project in KiCad 10 with its symbol and footprint libraries.

Rev 1.3 uses a different front-only placement: J1/J2/J3A/J3B form a top-edge
bank with north-facing, tangential entries. Every terminal, including sync
GND and SHIELD, reserves a 2.4 mm-wide, 10 mm-long straight approach outside
the board, clear of the module, other bodies and mounting reserves.
All selected parts, nets, outline and mounting centers are retained. Complete
pin legends for all four rev 1.3 connectors are printed together on the back.
The cables must fan left toward the hub and be strain-relieved to the rotor;
actual ferrules, jacket bend radius, clamp locations and swept-envelope fit
remain physical qualification gates. See the [rev 1.3 placement contract](1.3/README.md#board-layout).

| Directory | Contents |
|---|---|
| [1.2/](1.2/README.md) | Current complete rev 1.2 project for Quilter |
| [1.2/routed/](1.2/routed/README.md) | Unvalidated Quilter output; not ordered |
| [1.3/](1.3/README.md) | Experimental differential-sync project with locked placement. Mechanical review, routing, qualification and firmware remain |
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
For rev 1.2, running the PCB generator without `--unplaced` writes a separate
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

The generator emits 0.25 mm minimum pour thickness, 0.5 mm thermal gaps, and
0.5 mm thermal spokes. Before selecting a routed board, rank candidates
with DRC, zone and front-footprint checks. Release validation also refuses
footprints outside the front copper layer.

```sh
python hardware/phantasm/gen/analyze_candidates.py path/to/candidate
```

To check only the zone settings on an individual board, use:

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
and TBC05-04 in `gen/footprints/` are hand-maintained source files.
The PCB generator reads them and derives the footprints in `1.3/phantasm.pretty/`
with the revision's silkscreen clearance and wire-access reservations.

## Validation and release

The [routed rev 1.2 candidate](1.2/routed/README.md) is unvalidated Quilter
output and has **not been ordered**. A routed filename or successful CI run
is not approval for fabrication. CI checks its explicit status and source
hashes; release creation runs the fabrication gates. Generated inputs and
routing candidates can remain in Git while those gates fail.

Use the explicit project and revision together:

```sh
python hardware/phantasm/gen/release.py validate --revision 1.2 --project hardware/phantasm/1.2/routed
python hardware/phantasm/gen/release.py check-candidates
```

Validation snapshots the selected KiCad project, rules, library tables and
local library/model inputs. It checks the R-MECH-6 board width, manufacturing rules, zone/via geometry, mask/origin,
ENIG and copper stackup, copper connectivity, ERC, DRC, schematic parity,
netlist membership, and assembly assignments/metadata. Rev 1.2 does not inherit
rev 1.1's historical footprint-warning allowances; the generator's D_BUS land
is expected, and J1's edge-trimmed silkscreen may report at most once. The first static pass lists
independent failures together; later gates stop on failure. Neither validation
nor release repairs source files. Rev 1.3 release is not supported.

Repair a **separate copy of the whole project**, retaining its library tables
but not its `candidate.json` (the repair tools refuse hash-manifested directories).
Use the existing zone and clearance repair tools on that copy, review the diff,
restore ENIG and inspect the stackup in KiCad. Update the PCB from its schematic
and resolve placement/library differences as needed. Refill and save zones in
KiCad, then run validation against that saved project. Validation and release
refill zones in a disposable copy before all copper gates, retaining both
original and evaluated source hashes. They never change zone settings. Export
uses that evaluated copper without another refill. Check the harness, ferrules,
connector access and Teensy height physically; these are not automated gates.

Once validation passes, choose a new output directory outside the source tree:

```sh
python hardware/phantasm/gen/release.py release --revision 1.2 --project path/to/reviewed-project --output path/to/releases/rev1.2-candidate1
python hardware/phantasm/gen/release.py verify path/to/releases/rev1.2-candidate1
```

Release creation repeats every gate and exports Gerbers, drills, BOM and CPL.
It refuses an existing destination and publishes only after export-content and
digest checks pass. The package contains the exact source snapshot, reports,
file hashes, generator hashes, KiCad/Python versions and fabrication options in
its release manifest. Archive the entire package and retain the printed
**release SHA256** outside it. Verification checks bytes against the manifest;
the independently retained SHA256 identifies the manifest you reviewed.
Source edits require a new release. The older `fab.py` command remains the
legacy rev 1.1 exporter; it does not create this release record.

Before ordering, verify that SHA256, select the recorded ENIG/stackup options
at the fabricator, and review the exact package's assembly preview, rotations,
part substitutions and physical access. Upload only the Gerber archive and
assembly CSVs under the package's output/jlc directory. No command places an
order or marks a candidate ordered automatically.
The rev 1.2 input has no qualified positive overvoltage cutoff. Qualify
slip-ring overshoot and the resulting Teensy VIN/U1 rail before ordering;
passing the package gates does not establish that protection.

After an actual order, create a JSON details file containing vendor, order_id,
ordered_at (ISO date/time), fabrication_options copied from the manifest and
checked against the order, assembly_preview_reviewed=true, and
mechanical_access_reviewed=true. Record it with the reviewed digest:

```sh
python hardware/phantasm/gen/release.py record-order --package path/to/releases/rev1.2-candidate1 --release-sha256 REVIEWED_SHA256 --details path/to/order-details.json --output path/to/order-record.json
```

The recorder verifies the package and options and refuses to overwrite an
existing record. Commit the resulting order record beside the reviewed project
and archive the exact package with it; retain the vendor receipt separately.
The record documents an order already placed, not an authorization to order.
