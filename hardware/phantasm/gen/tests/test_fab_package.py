import contextlib
import csv
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import fab  # noqa: E402
import sexp  # noqa: E402


TITLED_BOARD = sexp.parse_one('(kicad_pcb (title_block (rev "1.1")))')


class ZipMemberTests(unittest.TestCase):
    """Upload archives retain deterministic members and metadata."""

    MEMBER = "phantasm-F_Cu.gtl"
    PAYLOAD = b"X000123Y004567D02*\n" * 4000
    def build_export_zip(self, mtime):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "phantasm-F_Cu.gtl"
            source.write_text("G04 gerber*\n", encoding="utf-8")
            os.utime(source, (mtime, mtime))
            archive = Path(directory) / fab.ARCHIVE
            fab.write_upload_zip(directory, [source.name], archive)
            return archive.read_bytes()

    def test_archive_is_independent_of_source_mtime(self):
        self.assertEqual(self.build_export_zip(946684800), self.build_export_zip(1700000000))

    def test_member_carries_no_host_metadata(self):
        info = fab.zip_member("phantasm-F_Cu.gtl")
        self.assertEqual(info.date_time, (1980, 1, 1, 0, 0, 0))
        self.assertEqual(info.create_system, 3)
        self.assertEqual(info.compress_type, zipfile.ZIP_STORED)


    def build_upload_zip(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (directory / self.MEMBER).write_bytes(self.PAYLOAD)
        path = directory / fab.ARCHIVE
        fab.write_upload_zip(str(directory), [self.MEMBER], str(path))
        return path.read_bytes()

    def test_members_are_stored(self):
        with zipfile.ZipFile(io.BytesIO(self.build_upload_zip())) as archive:
            info = archive.getinfo(self.MEMBER)
            self.assertEqual(info.compress_type, zipfile.ZIP_STORED)
            self.assertEqual(info.compress_size, len(self.PAYLOAD))

    def test_members_round_trip(self):
        with zipfile.ZipFile(io.BytesIO(self.build_upload_zip())) as archive:
            self.assertEqual(archive.namelist(), [self.MEMBER])
            self.assertEqual(archive.read(self.MEMBER), self.PAYLOAD)


class TimestampNormalizationTests(unittest.TestCase):
    """Every artifact's export wall-clock stamp is normalized."""

    # Gerbers are CRLF, Excellon and the job file LF.
    GERBER = ("%TF.GenerationSoftware,KiCad,Pcbnew,10.0.4*%\r\n"
              "%TF.CreationDate,{stamp}-07:00*%\r\n"
              "G04 Created by KiCad (PCBNEW 10.0.4) date {date}*\r\n"
              "D10*\r\nX1000Y2000D02*\r\nM02*\r\n")
    DRILL = ("M48\n; DRILL file KiCad 10.0.4 date {stamp}\n"
             "; #@! TF.CreationDate,{stamp}-07:00\nFMAT,2\nM30\n")
    JOB = ('{{\n  "GeneralSpecs": {{\n'
           '    "CreationDate": "{stamp}-07:00"\n  }}\n}}\n')

    def export(self, stamp, date):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        for name, template in (("phantasm-Top Layer.gtl", self.GERBER),
                               ("phantasm-PTH.drl", self.DRILL),
                               ("phantasm-job.gbrjob", self.JOB)):
            (directory / name).write_bytes(
                template.format(stamp=stamp, date=date).encode("utf-8"))
        rewritten = fab.normalize_fab_timestamps(directory)
        return directory, rewritten

    def test_two_exports_of_one_board_become_identical(self):
        first, rewritten = self.export("2026-08-12T17:25:19", "2026-08-12 17:25:19")
        second, _ = self.export("2019-03-04T01:02:03", "2019-03-04 01:02:03")
        self.assertEqual(len(rewritten), 3)
        for path in sorted(first.iterdir()):
            self.assertEqual(path.read_bytes(),
                             (second / path.name).read_bytes(), path.name)

    def test_line_endings_and_tool_version_survive(self):
        directory, _ = self.export("2026-08-12T17:25:19", "2026-08-12 17:25:19")
        gerber = (directory / "phantasm-Top Layer.gtl").read_bytes()
        self.assertNotIn(b"\n", gerber.replace(b"\r\n", b""))
        self.assertIn(f"%TF.CreationDate,{fab.FAB_TIMESTAMP}*%\r\n".encode(), gerber)
        self.assertIn(b"KiCad (PCBNEW 10.0.4) date", gerber)
        drill = (directory / "phantasm-PTH.drl").read_bytes()
        self.assertNotIn(b"\r", drill)
        self.assertIn(f"date {fab.FAB_TIMESTAMP}\n".encode(), drill)

    def test_an_undecodable_artifact_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (directory / "phantasm-F_Cu.gtl").write_bytes(
            b"%TF.CreationDate,\xff\xfe*%\r\n")
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "cannot read exported artifact"):
            fab.normalize_fab_timestamps(directory)

    def test_an_artifact_no_pattern_reaches_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (directory / "phantasm-F_Cu.gtl").write_bytes(
            b"%TF.Creation-Date,2026-08-12T17:25:19-07:00*%\r\n")
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "phantasm-F_Cu.gtl"):
            fab.normalize_fab_timestamps(directory)

    def test_a_respelled_banner_cannot_report_a_clean_run(self):
        directory, _ = self.export("2026-08-12T17:25:19", "2026-08-12 17:25:19")
        with mock.patch.object(fab, "TIMESTAMP_SUBSTITUTIONS", ()):
            with self.assertRaises(fab.TimestampNormalizationError) as caught:
                fab.normalize_fab_timestamps(directory)
        for name in sorted(path.name for path in directory.iterdir()):
            self.assertIn(name, str(caught.exception))

    def test_a_gerber_carrying_only_one_of_its_two_stamps_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        gerber = self.GERBER.format(stamp="2026-08-12T17:25:19",
                                    date="2026-08-12 17:25:19")
        (directory / "phantasm-F_Cu.gtl").write_bytes(
            gerber.replace("G04 Created by KiCad", "G04 Made by KiCad").encode())
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "gerber banner"):
            fab.normalize_fab_timestamps(directory)

    def test_an_unknown_artifact_type_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (directory / "phantasm-notes.txt").write_bytes(b"G04 gerber*\r\n")
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "phantasm-notes.txt"):
            fab.normalize_fab_timestamps(directory)

    def test_every_zipped_artifact_type_states_its_stamps(self):
        for name in sorted(fab.ZIP_MEMBERS):
            with self.subTest(name=name):
                self.assertTrue(
                    fab.REQUIRED_STAMPS[os.path.splitext(name)[1]])

    def test_normalizing_an_already_stamped_export_stays_clean(self):
        directory, first = self.export("2026-08-12T17:25:19",
                                       "2026-08-12 17:25:19")
        self.assertEqual(fab.normalize_fab_timestamps(directory), first)


class ZipMembershipTests(unittest.TestCase):
    """A full export: GERBER_LAYERS in Protel naming, both drills, the job."""

    EXPORTED = ("phantasm-F_Cu.gtl", "phantasm-In1_Cu.g1", "phantasm-In2_Cu.g2",
                "phantasm-B_Cu.gbl", "phantasm-F_Silkscreen.gto",
                "phantasm-B_Silkscreen.gbo", "phantasm-F_Mask.gts",
                "phantasm-B_Mask.gbs", "phantasm-F_Paste.gtp",
                "phantasm-B_Paste.gbp", "phantasm-Edge_Cuts.gm1",
                "phantasm-PTH.drl", "phantasm-NPTH.drl", "phantasm-job.gbrjob")

    def test_the_fixture_covers_every_required_member(self):
        self.assertEqual(set(self.EXPORTED), fab.ZIP_MEMBERS)

    def test_every_exported_artifact_is_zipped(self):
        self.assertEqual(fab.zip_members(self.EXPORTED + tuple(fab.ZIP_EXCLUDED)),
                         sorted(self.EXPORTED))

    def test_custom_layer_names_are_required_upload_members(self):
        board = sexp.parse_one('(kicad_pcb (layers (4 "In1.Cu" power "GND") '
                               '(6 "In2.Cu" power "Ground") '
                               '(0 "F.Cu" signal "Top.Copper")))')
        expected = (fab.ZIP_MEMBERS - {"phantasm-In1_Cu.g1", "phantasm-In2_Cu.g2",
                                       "phantasm-F_Cu.gtl"}
                    | {"phantasm-GND.g1", "phantasm-Ground.g2", "phantasm-Top_Copper.gtl"})
        self.assertEqual(fab.board_zip_members(board), expected)
        self.assertEqual(fab.zip_members(expected, fab.board_zip_members(board)), sorted(expected))
        with self.assertRaises(fab.UploadPackageError):
            fab.zip_members(self.EXPORTED, fab.board_zip_members(board))

    def test_an_extension_outside_the_allowlist_is_rejected(self):
        with self.assertRaisesRegex(
                fab.UploadPackageError, "phantasm-User_2.gm2"):
            fab.zip_members(self.EXPORTED + ("phantasm-User_2.gm2",))

    def test_an_unexported_layer_is_rejected(self):
        for missing in self.EXPORTED:
            extension = os.path.splitext(missing)[1]
            partial = tuple(n for n in self.EXPORTED
                            if os.path.splitext(n)[1] != extension)
            with self.subTest(extension=extension):
                with self.assertRaisesRegex(fab.UploadPackageError, extension):
                    fab.zip_members(partial)

    def test_missing_npth_drill_is_rejected(self):
        partial = tuple(n for n in self.EXPORTED
                        if n != "phantasm-NPTH.drl")
        with self.assertRaisesRegex(fab.UploadPackageError,
                                    "phantasm-NPTH.drl"):
            fab.zip_members(partial)

    def test_the_upload_zip_itself_is_not_a_member(self):
        self.assertNotIn("phantasm-jlc-gerbers.zip",
                         fab.zip_members(self.EXPORTED
                                         + ("phantasm-jlc-gerbers.zip",)))


FAB_GATES = {
    "kicad_cli": "fixture-cli", "read_board": TITLED_BOARD,
    "validate_plot_origin": None, "validate_via_geometry": 0,
    "validate_solder_mask": None, "validate_zone_geometry": 0,
    "validate_project_rules": 0, "run_drc": (0, 0),
    "run_erc": 0, "run_parity": 0, "validate_netlist_spec": 0,
    "parse_components": {}, "validate_assembled_refs": None,
    "validate_rotation_refs": None, "validate_assembly_metadata": {},
    "validate_part_catalog": None, "normalize_fab_timestamps": [],
    "validate_fab_content": {"plated": 0, "unplated": 0},
}


def patch_fab_run(stack, out, export, **overrides):
    """Stub every `fab.main` gate, export through `export`, and write into `out`.

    `overrides` replace individual gate return values.
    """
    unknown = set(overrides) - set(FAB_GATES)
    assert not unknown, f"not a stubbed gate: {sorted(unknown)}"
    for name, value in {"OUT": str(out), "JLC": str(out / "jlc")}.items():
        stack.enter_context(mock.patch.object(fab, name, value))
    for name, value in {**FAB_GATES, **overrides}.items():
        stack.enter_context(mock.patch.object(fab, name, return_value=value))
    stack.enter_context(mock.patch.object(fab, "run_export", side_effect=export))
    stack.enter_context(contextlib.redirect_stdout(io.StringIO()))


class ZoneRefillTests(unittest.TestCase):
    def test_drc_and_parity_refill_only_when_requested(self):
        for runner in (fab.run_drc, fab.run_parity):
            for options in ({}, {"refill": False}, {"refill": True}):
                with self.subTest(runner=runner.__name__, options=options), \
                        tempfile.TemporaryDirectory() as directory:
                    def fake_run(args, check=True):
                        Path(args[args.index("-o") + 1]).write_text(json.dumps({
                            "violations": [], "unconnected_items": [],
                            "schematic_parity": [],
                        }), encoding="utf-8")
                        self.assertEqual(args.count("--refill-zones"),
                                         int(options.get("refill", False)))
                        self.assertNotIn("--save-board", args)
                        return subprocess.CompletedProcess(args, 0)

                    with mock.patch.object(fab, "kicad_cli", return_value="kicad-cli"), \
                            mock.patch.object(fab, "run", side_effect=fake_run), \
                            mock.patch.object(fab, "require_schematic_parity", return_value=0):
                        runner(Path(directory) / "report.json", **options)

    def test_fabrication_gates_use_the_export_refill_policy(self):
        for refill in (False, True):
            with self.subTest(refill=refill), tempfile.TemporaryDirectory() as directory, \
                    contextlib.ExitStack() as stack:
                patch_fab_run(stack, Path(directory), mock.Mock(side_effect=RuntimeError("export")))
                with self.assertRaisesRegex(RuntimeError, "export"):
                    fab.fabricate(refill_zones=refill)
                fab.run_drc.assert_called_once_with(
                    str(Path(directory) / "phantasm-drc.json"), refill=refill)
                fab.run_parity.assert_called_once_with(
                    str(Path(directory) / "phantasm-parity.json"), refill=refill)


class FabricationWiringTests(unittest.TestCase):
    def fabricate(self, schematic_revision="1.1", board=TITLED_BOARD, comps=None):
        with tempfile.TemporaryDirectory() as directory, contextlib.ExitStack() as stack:
            def export(stage, args):
                target = Path(args[args.index("-o") + 1])
                if stage == "netlist":
                    target.write_text('(export (design (sheet (name "/") '
                                      f'(title_block (rev "{schematic_revision}")))))', encoding="utf-8")
                elif stage == "centroid":
                    target.write_text("Ref,PosX,PosY,Rot,Side\n", encoding="utf-8")
            patch_fab_run(stack, Path(directory), export, read_board=board,
                          parse_components=comps or {})
            fab.fabricate(validate_only=True)

    def test_fabricate_rejects_foreign_schematic_revision(self):
        self.fabricate()
        with self.assertRaisesRegex(SystemExit, "schematic revision does not match"):
            self.fabricate(schematic_revision="1.2")

    def test_fabricate_rejects_assembled_board_exclusion(self):
        comps = {"R1": {"value": "33R", "footprint": "Resistor_SMD:R_0805", "dnp": False}}
        self.fabricate(comps=comps)
        board = sexp.parse_one('(kicad_pcb (title_block (rev "1.1")) '
                               '(footprint "Resistor_SMD:R_0805" '
                               '(property "Reference" "R1") (attr exclude_from_bom)))')
        with self.assertRaisesRegex(SystemExit, "assembled parts the board excludes.*R1"):
            self.fabricate(board=board, comps=comps)


class PackagePromotionTests(unittest.TestCase):
    def test_success_promotes_bom_cpl_archive_and_manifest(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.ExitStack() as stack:
            out = Path(directory)
            stale = out / (fab.STAGE_PREFIX + "stale")
            stale.mkdir()
            (stale / "aborted-export.gbr").write_bytes(b"stranded export")
            jlc = out / "jlc"
            jlc.mkdir()
            (jlc / "previous.txt").write_text("previous package", encoding="utf-8")
            comps = {
                "R_A": {"value": "33R", "footprint": "Resistor:R_0805", "dnp": False},
                "R_B": {"value": "33R", "footprint": "Resistor:R_0805", "dnp": False},
                "R_C": {"value": "47R", "footprint": "Resistor:R_0805", "dnp": False},
                "R_D": {"value": "33R", "footprint": "Resistor:R_0603", "dnp": False},
                "R_E": {"value": "33R", "footprint": "Resistor:R_0805", "dnp": False},
                "U1": {"value": "74AHCT125", "footprint": "Logic:SOIC-14", "dnp": False},
            }
            metadata = {
                ref: {"lcsc": "C17634", "pos_x": str(index + 10),
                      "pos_y": str(index + 20), "side": "top", "rotation": 17.0}
                for index, ref in enumerate(comps)
            }
            metadata["R_E"].update(lcsc="C17408", side="bottom")
            metadata["U1"].update(lcsc="C155176", rotation=180.0)

            def export(stage, args):
                target = Path(args[args.index("-o") + 1])
                if stage == "netlist":
                    target.write_text('(export (design (sheet (name "/") '
                                      '(title_block (rev "1.1")))))', encoding="utf-8")
                elif stage == "centroid":
                    target.write_text("Ref,PosX,PosY,Rot,Side\n", encoding="utf-8")
                elif stage == "gerber":
                    for name in ZipMembershipTests.EXPORTED:
                        (target / name).write_bytes(f"fixture export: {name}\n".encode())

            patch_fab_run(stack, out, export, parse_components=comps,
                          validate_assembly_metadata=metadata)
            fab.main()

            for name in FAB_GATES:
                with self.subTest(gate=name):
                    getattr(fab, name).assert_called()
            fab.validate_via_geometry.assert_called_once_with(
                fab.PCB, min_vias=fab.MIN_BOARD_VIAS, board=TITLED_BOARD)
            fab.normalize_fab_timestamps.assert_called_once()
            staged = Path(fab.normalize_fab_timestamps.call_args.args[0])
            self.assertEqual(staged.parent, out)
            self.assertTrue(staged.name.startswith(fab.STAGE_PREFIX))
            fab.validate_fab_content.assert_called_once_with(str(staged), TITLED_BOARD)

            expected_files = set(ZipMembershipTests.EXPORTED) | {
                "phantasm-BOM.csv", "phantasm-CPL.csv", fab.ARCHIVE, fab.SUMS_FILE}
            self.assertEqual({path.name for path in jlc.iterdir()}, expected_files)
            self.assertFalse((out / "jlc.previous").exists())
            self.assertFalse(any(path.name.startswith(fab.STAGE_PREFIX) for path in out.iterdir()))
            self.assertFalse((out / fab.LOCK_FILE).exists())
            with (jlc / "phantasm-BOM.csv").open(newline="", encoding="utf-8") as stream:
                bom = list(csv.reader(stream))
            self.assertEqual(bom[0], ["Comment", "Designator", "Footprint", "LCSC Part #",
                                      "Manufacturer", "Manufacturer Part Number", "Description"])
            self.assertEqual([row[:4] for row in bom[1:]], [
                ["33R", "R_D", "R_0603", "C17634"],
                ["33R", "R_E", "R_0805", "C17408"],
                ["33R", "R_A,R_B", "R_0805", "C17634"],
                ["47R", "R_C", "R_0805", "C17634"],
                ["74AHCT125", "U1", "SOIC-14", "C155176"],
            ])
            for row in bom[1:]:
                part = fab.PART_BY_LCSC[row[3]]
                self.assertEqual(row[4:], [part["manufacturer"], part["mpn"], part["description"]])
            with (jlc / "phantasm-CPL.csv").open(newline="", encoding="utf-8") as stream:
                cpl = list(csv.reader(stream))
            self.assertEqual(cpl, [
                ["Designator", "Mid X", "Mid Y", "Layer", "Rotation"],
                ["R_A", "10", "20", "top", "17.000000"],
                ["R_B", "11", "21", "top", "17.000000"],
                ["R_C", "12", "22", "top", "17.000000"],
                ["R_D", "13", "23", "top", "17.000000"],
                ["R_E", "14", "24", "bottom", "17.000000"],
                ["U1", "15", "25", "top", "90.000000"],
            ])
            with zipfile.ZipFile(jlc / fab.ARCHIVE) as archive:
                self.assertEqual(archive.namelist(), sorted(ZipMembershipTests.EXPORTED))
                for name in archive.namelist():
                    self.assertEqual(archive.read(name), (jlc / name).read_bytes())
            self.assertEqual(fab.verify_package(str(jlc), None), len(expected_files) - 1)

    def test_invalid_members_preserve_the_previous_package(self):
        for replacement in (None, "unexpected-In1_Cu.g1", "archive-failure", "stale-backup"):
            with self.subTest(replacement=replacement), \
                    tempfile.TemporaryDirectory() as directory, \
                    contextlib.ExitStack() as stack:
                out = Path(directory)
                jlc = out / "jlc"
                jlc.mkdir()
                if replacement == "stale-backup":
                    (out / "jlc.previous").mkdir()
                previous = {fab.ARCHIVE: b"previous archive",
                            fab.SUMS_FILE: b"previous manifest"}
                for name, content in previous.items():
                    (jlc / name).write_bytes(content)

                def export(stage, args):
                    if replacement == "stale-backup" and stage == "gerber":
                        self.fail("stale backup must be diagnosed before gerber export")
                    target = Path(args[args.index("-o") + 1])
                    if stage == "netlist":
                        target.write_text('(export (design (sheet (name "/") '
                                          '(title_block (rev "1.1")))))')
                    elif stage == "centroid":
                        target.write_text("Ref,PosX,PosY,Rot,Side\n")
                    elif stage == "gerber":
                        names = set(ZipMembershipTests.EXPORTED) - {"phantasm-In1_Cu.g1"}
                        if replacement == "archive-failure":
                            names.add("phantasm-In1_Cu.g1")
                        elif replacement:
                            names.add(replacement)
                        for name in names:
                            (target / name).write_text("fixture export")

                patch_fab_run(stack, out, export)
                if replacement == "archive-failure":
                    stack.enter_context(mock.patch.object(
                        fab, "write_upload_zip", side_effect=OSError("injected archive failure")))
                with self.assertRaisesRegex((SystemExit, OSError),
                                            "phantasm-In1_Cu.g1|unexpected-In1_Cu.g1|injected archive failure|recover previous package"):
                    fab.main()
                self.assertEqual({p.name: p.read_bytes() for p in jlc.iterdir()}, previous)
                self.assertFalse((out / fab.LOCK_FILE).exists())


class OutputLockTests(unittest.TestCase):
    def test_held_lock_refuses_without_touching_a_live_stage(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.ExitStack() as stack:
            out = Path(directory)
            live = out / (fab.STAGE_PREFIX + "live")
            live.mkdir()
            (live / "phantasm-F_Cu.gtl").write_text("in progress", encoding="utf-8")
            (out / fab.LOCK_FILE).write_text("1\n", encoding="utf-8")

            def export(stage, args):
                self.fail(f"a refused run must not export {stage}")

            patch_fab_run(stack, out, export)
            with self.assertRaisesRegex(SystemExit, "another fabrication run holds"):
                fab.main()
            self.assertEqual((live / "phantasm-F_Cu.gtl").read_text(encoding="utf-8"),
                             "in progress")
            self.assertEqual((out / fab.LOCK_FILE).read_text(encoding="utf-8"), "1\n")
            self.assertFalse((out / "jlc").exists())

    def test_a_run_started_during_export_keeps_the_live_stage(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.ExitStack() as stack:
            out = Path(directory)
            refusals = []

            def export(stage, args):
                target = Path(args[args.index("-o") + 1])
                if stage == "netlist":
                    target.write_text('(export (design (sheet (name "/") '
                                      '(title_block (rev "1.1")))))', encoding="utf-8")
                elif stage == "centroid":
                    target.write_text("Ref,PosX,PosY,Rot,Side\n", encoding="utf-8")
                elif stage == "gerber":
                    with self.assertRaises(SystemExit) as refused:
                        fab.main()
                    refusals.append(str(refused.exception))
                    for name in ZipMembershipTests.EXPORTED:
                        (target / name).write_bytes(f"fixture export: {name}\n".encode())

            patch_fab_run(stack, out, export)
            fab.main()
            self.assertEqual(len(refusals), 1)
            self.assertIn("another fabrication run holds", refusals[0])
            self.assertEqual({path.name for path in (out / "jlc").iterdir()},
                             set(ZipMembershipTests.EXPORTED)
                             | {"phantasm-BOM.csv", "phantasm-CPL.csv", fab.ARCHIVE, fab.SUMS_FILE})
            self.assertFalse((out / fab.LOCK_FILE).exists())


class PackageManifestTests(unittest.TestCase):
    """Fabrication manifest membership and digests."""

    ARCHIVE = "phantasm-jlc-gerbers.zip"

    def test_manifest_covers_members_assembly_data_and_the_archive(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        payloads = {"phantasm-F_Cu.gtl": b"G04 gerber*",
                    "phantasm-PTH.drl": b"M48 M30",
                    self.ARCHIVE: b"zip bytes"}
        payloads.update((name, name.encode()) for name in fab.ASSEMBLY_MEMBERS)
        for name, data in payloads.items():
            (directory / name).write_bytes(data)
        members = ["phantasm-F_Cu.gtl", "phantasm-PTH.drl"]
        lines = fab.package_manifest(str(directory), members,
                                     self.ARCHIVE).splitlines()
        self.assertEqual(
            lines,
            [f"{hashlib.sha256(payloads[name]).hexdigest()}  {name}"
             for name in members + list(fab.ASSEMBLY_MEMBERS)
             + [self.ARCHIVE]])

    def test_the_assembly_csvs_travel_beside_the_zip(self):
        self.assertEqual(set(fab.ASSEMBLY_MEMBERS),
                         {"phantasm-BOM.csv", "phantasm-CPL.csv"})
        self.assertTrue(set(fab.ASSEMBLY_MEMBERS) <= fab.ZIP_EXCLUDED)
        self.assertFalse(set(fab.ASSEMBLY_MEMBERS) & fab.ZIP_MEMBERS)

    def test_the_manifest_is_not_an_upload_member(self):
        self.assertIn(fab.SUMS_FILE, fab.ZIP_EXCLUDED)
        self.assertNotIn(fab.SUMS_FILE, fab.ZIP_MEMBERS)


class FabContentTests(unittest.TestCase):
    """Zip membership is a filename test; this reads the exported bytes."""

    # An aperture and a flash that uses it, then the same layer plotting nothing.
    GERBER = ("%FSLAX46Y46*%\r\n%MOMM*%\r\n%ADD10C,0.500000*%\r\n"
              "D10*\r\nX1000000Y2000000D03*\r\nM02*\r\n")
    APERTURELESS = "%FSLAX46Y46*%\r\n%MOMM*%\r\nM02*\r\n"
    UNDRAWN = "%FSLAX46Y46*%\r\n%MOMM*%\r\n%ADD10C,0.500000*%\r\nM02*\r\n"
    REGION = ("%FSLAX46Y46*%\r\n%MOMM*%\r\n%LPD*%\r\nG36*\r\n"
              "X0Y0D02*\r\nG01X1000000Y0D01*\r\n"
              "X1000000Y1000000D01*\r\nX0Y1000000D01*\r\n"
              "X0Y0D01*\r\nG37*\r\nM02*\r\n")
    JOB = '{"GeneralSpecs": {}}\n'
    # One via and one plated pad, one unplated mounting hole, one SMD pad.
    BOARD = ('(kicad_pcb (via (at 1 2) (size 0.45) (drill 0.2)) '
             '(footprint "R" (pad "1" thru_hole circle (at 0 0) (drill 1)) '
             '(pad "2" smd roundrect (at 1 0) (drill (offset 0.3 0))) '
             '(pad "" np_thru_hole circle (at 2 0) (drill 2.7))))')

    def drill(self, holes):
        body = "".join(f"X{index}.0Y1.0\n" for index in range(1, holes + 1))
        return "M48\nFMAT,2\nMETRIC\nT1C0.200\n%\nG90\nT1\n" + body + "M30\n"

    def export(self, overrides=None):
        """A staged export of the fixture board, with `overrides` swapped in."""
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        pcb = Path(self.enterContext(tempfile.TemporaryDirectory())) / "content.kicad_pcb"
        pcb.write_text(self.BOARD, encoding="utf-8")
        board = fab.read_board(pcb)
        holes = fab.board_hole_counts(board)
        payloads = {}
        for name in ZipMembershipTests.EXPORTED:
            extension = os.path.splitext(name)[1]
            if extension == ".gbrjob":
                payloads[name] = self.JOB
            elif extension == ".drl":
                payloads[name] = self.drill(holes[fab.DRILL_MEMBERS[name]])
            elif name in fab.EMPTY_GERBER_MEMBERS:
                payloads[name] = self.APERTURELESS
            else:
                payloads[name] = self.GERBER
        payloads.update(overrides or {})
        for name, text in payloads.items():
            (directory / name).write_bytes(text.encode())
        return directory, board

    def validate(self, overrides=None):
        directory, board = self.export(overrides)
        return fab.validate_fab_content(directory, board)

    def test_accepts_an_export_that_carries_the_board(self):
        self.assertEqual(self.validate(), {"plated": 2, "unplated": 1})

    def test_accepts_aperture_free_mask_regions(self):
        self.assertEqual(self.validate({"phantasm-F_Mask.gts": self.REGION,
                                        "phantasm-B_Mask.gbs": self.REGION}),
                         {"plated": 2, "unplated": 1})

    def test_accepts_arc_and_multiple_regions(self):
        arc = ("G36*X0Y0D02*G02X0Y0I1000000J0D01*G37*")
        for body in (arc, arc + arc):
            with self.subTest(body=body):
                self.assertEqual(self.validate({"phantasm-F_Mask.gts": body}),
                                 {"plated": 2, "unplated": 1})

    def test_rejects_missing_or_undrawn_regions_without_apertures(self):
        for body in ("G36*G37*", "G36*X0Y0D02*X1Y1D02*G37*",
                     "G36*X0Y0D02*X1Y1D01*", "X1Y1D01*G37*",
                     "G36*X1Y1D03*G37*", "G36*D01*G37*",
                     "G04 G36*X1Y1D01*G04 G37*", "X1Y1D01*",
                     "G36*G04 X1Y1D01*G37*",
                     "G36*G36*X1Y1D01*G37*G37*",
                     "%AMunused*G36*X1Y1D01*G37*%"):
            with self.subTest(body=body):
                with self.assertRaisesRegex(fab.FabContentError,
                                            "defines no apertures"):
                    self.validate({"phantasm-F_Mask.gts": body})

    def test_rejects_moves_and_comments_as_drawing_operations(self):
        for command in ("X1000000Y2000000D02", "G04 X1000000Y2000000D03",
                        "G04 D01", "G04 D03"):
            with self.subTest(command=command):
                layer = self.GERBER.replace("X1000000Y2000000D03", command)
                with self.assertRaisesRegex(fab.FabContentError, "draws with none"):
                    self.validate({"phantasm-F_Cu.gtl": layer})

    def test_accepts_draw_and_flash_commands(self):
        for command in ("X1000000Y2000000D01", "G01X1000000Y2000000D01",
                        "G02X1000000Y2000000I-100J+200D01", "D03"):
            with self.subTest(command=command):
                layer = self.GERBER.replace("X1000000Y2000000D03", command)
                self.assertEqual(self.validate({"phantasm-F_Cu.gtl": layer}),
                                 {"plated": 2, "unplated": 1})

    def test_accepts_flashes_preceded_by_object_attributes(self):
        layer = self.GERBER.replace(
            "X1000000Y2000000D03*\r\n",
            "%TO.P,J1,1*%\r\n%TO.N,GND*%\r\nX1000000Y2000000D03*\r\n%TD*%\r\n")
        self.assertEqual(self.validate({"phantasm-F_Cu.gtl": layer}),
                         {"plated": 2, "unplated": 1})

    def test_rejects_a_layer_that_plots_nothing(self):
        with self.assertRaisesRegex(fab.FabContentError,
                                    "phantasm-F_Cu.gtl: defines no apertures"):
            self.validate({"phantasm-F_Cu.gtl": self.APERTURELESS})

    def test_rejects_a_layer_that_draws_with_none_of_its_apertures(self):
        with self.assertRaisesRegex(fab.FabContentError, "draws with none"):
            self.validate({"phantasm-In1_Cu.g1": self.UNDRAWN})

    def test_allows_only_the_bottom_stencil_to_be_empty(self):
        self.assertEqual(fab.EMPTY_GERBER_MEMBERS, {"phantasm-B_Paste.gbp"})
        self.assertLess(fab.EMPTY_GERBER_MEMBERS, fab.ZIP_MEMBERS)

    def test_rejects_a_drill_file_short_a_hole(self):
        with self.assertRaisesRegex(
                fab.FabContentError,
                r"phantasm-PTH.drl: drills 1 holes, the board carries 2 plated"):
            self.validate({"phantasm-PTH.drl": self.drill(1)})

    def test_rejects_a_drill_file_the_board_does_not_account_for(self):
        with self.assertRaisesRegex(fab.FabContentError,
                                    "no board hole count covers"):
            self.validate({"phantasm-Extra.drl": self.drill(1)})

    def test_counts_a_slot_as_one_hole(self):
        self.assertEqual(
            len(fab.EXCELLON_HOLE.findall("X1.0Y1.0G85X2.0Y1.0\nX3.0Y1.0\n")),
            2)

    def test_counts_a_route_mode_slot_as_one_hole(self):
        routed = "X3.0Y1.0\nG00X1.0Y1.0\nM15\nG01X2.0Y1.0\nM16\nG05\n"
        self.assertEqual(len(fab.EXCELLON_HOLE.findall(routed)), 2)

    def test_committed_board_holes(self):
        # 99 vias and 45 plated pads; the four mounting holes are unplated.
        # Re-measure after promoting a re-route, as with the via floor.
        self.assertEqual(fab.board_hole_counts(fab.read_board(fab.PCB)),
                         {"plated": 144, "unplated": 4})

    def test_smd_shape_offsets_are_not_drill_holes(self):
        board = sexp.parse_one(self.BOARD)
        self.assertEqual(fab.board_hole_counts(board), {"plated": 2, "unplated": 1})


class PackageVerificationTests(unittest.TestCase):
    """A generated package is re-hashed against its recorded digests."""

    def package(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        for name in (ZipMembershipTests.EXPORTED + (fab.ARCHIVE,)
                     + fab.ASSEMBLY_MEMBERS):
            (directory / name).write_bytes(name.encode())
        members = fab.zip_members(os.listdir(directory))
        (directory / fab.SUMS_FILE).write_text(
            fab.package_manifest(str(directory), members, fab.ARCHIVE),
            encoding="utf-8", newline="\n")
        ordered = Path(self.enterContext(tempfile.TemporaryDirectory()))
        baseline = ordered / "ordered-SHA256SUMS.txt"
        baseline.write_bytes((directory / fab.SUMS_FILE).read_bytes())
        return directory, baseline

    def test_accepts_a_package_matching_manifest_and_baseline(self):
        directory, baseline = self.package()
        self.assertEqual(fab.verify_package(str(directory), str(baseline)),
                         len(fab.ZIP_MEMBERS)
                         + len(fab.ASSEMBLY_MEMBERS) + 1)

    def test_accepts_custom_layer_names_recorded_in_the_manifest(self):
        directory, _ = self.package()
        for old, new in (("phantasm-In1_Cu.g1", "phantasm-GND.g1"),
                         ("phantasm-In2_Cu.g2", "phantasm-Ground.g2")):
            (directory / old).rename(directory / new)
        members = sorted(set(os.listdir(directory)) - fab.ZIP_EXCLUDED)
        (directory / fab.SUMS_FILE).write_text(
            fab.package_manifest(str(directory), members, fab.ARCHIVE), encoding="utf-8")
        self.assertEqual(fab.verify_package(str(directory), None),
                         len(members) + len(fab.ASSEMBLY_MEMBERS) + 1)

    def test_rejects_an_edited_artifact(self):
        directory, baseline = self.package()
        (directory / "phantasm-F_Cu.gtl").write_bytes(b"edited")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "phantasm-F_Cu.gtl"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_an_edited_assembly_csv(self):
        for name in fab.ASSEMBLY_MEMBERS:
            with self.subTest(artifact=name):
                directory, baseline = self.package()
                (directory / name).write_bytes(b"C123456,swapped")
                with self.assertRaisesRegex(fab.PackageVerificationError, name):
                    fab.verify_package(str(directory), str(baseline))

    def test_reports_an_artifact_the_package_lost(self):
        directory, baseline = self.package()
        (directory / "phantasm-CPL.csv").unlink()
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "phantasm-CPL.csv: recorded but not"):
            fab.verify_package(str(directory), str(baseline))

    def test_reports_a_missing_zipped_artifact(self):
        directory, baseline = self.package()
        (directory / "phantasm-F_Cu.gtl").unlink()
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "phantasm-F_Cu.gtl: recorded but not in the package"):
            fab.verify_package(str(directory), str(baseline))

    def test_reports_an_unrecognized_package_file(self):
        directory, baseline = self.package()
        (directory / ".DS_Store").write_bytes(b"metadata")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "contains unrecognized package files: .DS_Store"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_manifest_that_skips_an_artifact(self):
        directory, baseline = self.package()
        manifest = directory / fab.SUMS_FILE
        kept = [line for line in manifest.read_text(encoding="utf-8").splitlines()
                if not line.endswith("phantasm-PTH.drl")]
        manifest.write_text("\n".join(kept) + "\n", encoding="utf-8",
                            newline="\n")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "records no digest for: phantasm-PTH.drl"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_package_that_is_not_what_was_ordered(self):
        directory, baseline = self.package()
        lines = baseline.read_text(encoding="utf-8").splitlines()
        flipped = ("1" if lines[0][0] == "0" else "0") + lines[0][1:]
        baseline.write_text("\n".join([flipped] + lines[1:]) + "\n",
                            encoding="utf-8", newline="\n")
        with self.assertRaisesRegex(fab.PackageVerificationError, "ordered"):
            fab.verify_package(str(directory), str(baseline))

    def test_reports_an_absent_baseline(self):
        directory, baseline = self.package()
        baseline.unlink()
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "no committed digest baseline"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_baseline_missing_an_artifact(self):
        directory, baseline = self.package()
        fab.verify_package(str(directory), str(baseline))
        lines = baseline.read_text(encoding="utf-8").splitlines()
        baseline.write_text("\n".join(line for line in lines
                                      if not line.endswith("phantasm-BOM.csv")) + "\n",
                            encoding="utf-8")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    r"phantasm-BOM.csv: [0-9a-f]{64} is not the ordered absent"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_baseline_with_an_extra_artifact(self):
        directory, baseline = self.package()
        fab.verify_package(str(directory), str(baseline))
        digest = "0" * 64
        baseline.write_text(baseline.read_text(encoding="utf-8")
                            + f"{digest}  obsolete.txt\n", encoding="utf-8")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    f"obsolete.txt: absent is not the ordered {digest}"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_manifest_with_an_uncovered_artifact(self):
        directory, _ = self.package()
        fab.verify_package(str(directory), None)
        manifest = directory / fab.SUMS_FILE
        manifest.write_text(manifest.read_text(encoding="utf-8")
                            + f"{'0' * 64}  extra.txt\n", encoding="utf-8")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "records artifacts the package does not hold: extra.txt"):
            fab.verify_package(str(directory), None)

    def test_checks_the_package_alone_without_a_baseline(self):
        directory, _ = self.package()
        self.assertEqual(fab.verify_package(str(directory), None),
                         len(fab.ZIP_MEMBERS)
                         + len(fab.ASSEMBLY_MEMBERS) + 1)

    def test_rejects_duplicate_manifest_entries(self):
        directory, baseline = self.package()
        manifest = directory / fab.SUMS_FILE
        content = manifest.read_text(encoding="utf-8")
        manifest.write_text(content + content.splitlines()[0] + "\n", encoding="utf-8")
        with self.assertRaisesRegex(fab.PackageVerificationError, "duplicate manifest entry"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_malformed_manifest_line(self):
        directory, baseline = self.package()
        (directory / fab.SUMS_FILE).write_text("not a digest line\n",
                                               encoding="utf-8", newline="\n")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "malformed manifest line"):
            fab.verify_package(str(directory), str(baseline))

    def test_reports_an_ungenerated_package(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "no fab package to verify"):
            fab.verify_package(str(directory / "jlc"), None)


class DirectPackagePromotionTests(unittest.TestCase):
    def make_package(self, root, name, content):
        package = Path(root) / name
        package.mkdir()
        (package / "manifest.txt").write_text(content, encoding="utf-8")
        return str(package)

    def test_promotes_without_a_previous_package(self):
        with tempfile.TemporaryDirectory() as root:
            staged = self.make_package(root, "staged", "new")
            destination = str(Path(root) / "jlc")
            fab.promote_package(staged, destination)
            self.assertEqual((Path(destination) / "manifest.txt").read_text(), "new")
            self.assertFalse(Path(staged).exists())
            self.assertFalse(Path(destination + ".previous").exists())

    def test_replaces_previous_package_and_removes_backup(self):
        with tempfile.TemporaryDirectory() as root:
            staged = self.make_package(root, "staged", "new")
            destination = self.make_package(root, "jlc", "old")
            fab.promote_package(staged, destination)
            self.assertEqual((Path(destination) / "manifest.txt").read_text(), "new")
            self.assertFalse(Path(destination + ".previous").exists())

    def test_failed_swap_restores_previous_package(self):
        with tempfile.TemporaryDirectory() as root:
            staged = self.make_package(root, "staged", "new")
            destination = self.make_package(root, "jlc", "old")
            replace = os.replace
            def fail_swap(source, target):
                if source == staged:
                    raise OSError("swap refused")
                return replace(source, target)
            with mock.patch.object(fab.os, "replace", side_effect=fail_swap):
                with self.assertRaisesRegex(OSError, "swap refused"):
                    fab.promote_package(staged, destination)
            self.assertEqual((Path(destination) / "manifest.txt").read_text(), "old")
            self.assertEqual((Path(staged) / "manifest.txt").read_text(), "new")
            self.assertFalse(Path(destination + ".previous").exists())

    def test_stale_backup_refuses_promotion_without_mutation(self):
        with tempfile.TemporaryDirectory() as root:
            staged = self.make_package(root, "staged", "new")
            destination = self.make_package(root, "jlc", "old")
            backup = self.make_package(root, "jlc.previous", "backup")
            with self.assertRaises(fab.UploadPackageError):
                fab.promote_package(staged, destination)
            for package, content in ((staged, "new"), (destination, "old"), (backup, "backup")):
                self.assertEqual((Path(package) / "manifest.txt").read_text(), content)


if __name__ == "__main__":
    unittest.main()
