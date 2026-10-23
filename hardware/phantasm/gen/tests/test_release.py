"""Candidate provenance and fabrication release contracts."""

import contextlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fab  # noqa: E402
import release  # noqa: E402
from kicad_env import GENERATES, GENERATES_REASON  # noqa: E402


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.project = self.root / "1.2" / "routed"
        self.project.mkdir(parents=True)
        for suffix in (".kicad_pcb", ".kicad_pro", ".kicad_sch"):
            (self.project / ("phantasm" + suffix)).write_text("source", encoding="utf-8")
        self.destination = self.root / "package"

    def status(self):
        release.write_json(self.project / "candidate.json", {
            "revision": "1.2", "status": "unvalidated", "ordered": False,
            "producer": "Quilter", "provenance": "fixture import",
            "source_sha256": release.sources(self.project),
        })

    def test_candidate_requires_unchanged_sources_and_unordered_status(self):
        self.status()
        self.assertFalse(release.check_candidate(self.project)["ordered"])
        (self.project / "phantasm.kicad_pro").write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "source hashes changed"):
            release.check_candidate(self.project)
        self.status()
        status_path = self.project / "candidate.json"
        status = json.loads(status_path.read_text())
        status["ordered"] = True
        release.write_json(status_path, status)
        with self.assertRaisesRegex(ValueError, "ordered=false"):
            release.check_candidate(self.project)

    def test_snapshot_includes_libraries_and_rules_excludes_history(self):
        library = self.project / "phantasm.pretty"
        library.mkdir()
        (library / "part.kicad_mod").write_text("land", encoding="utf-8")
        (self.project / "phantasm.kicad_dru").write_text("rule", encoding="utf-8")
        history = self.project / ".history"
        history.mkdir()
        (history / "old.kicad_pcb").write_text("old", encoding="utf-8")
        hashes = release.snapshot(self.project, self.root / "snapshot")
        self.assertIn("phantasm.pretty/part.kicad_mod", hashes)
        self.assertIn("phantasm.kicad_dru", hashes)
        self.assertNotIn(".history/old.kicad_pcb", hashes)
        self.assertEqual(release.sources(self.root / "snapshot"), hashes)

    def test_real_candidate_reports_finish_and_zone_failures_without_mutation(self):
        before = release.sources(release.CANDIDATE)
        release.check_candidate(release.CANDIDATE)
        with self.assertRaises(ValueError) as caught:
            release.preflight(release.CANDIDATE, "1.2")
        self.assertIn("requires ENIG", str(caught.exception))
        self.assertIn("min_thickness 0.0254", str(caught.exception))
        self.assertEqual(release.sources(release.CANDIDATE), before)

    def test_revision_context_does_not_leak_historical_allowances(self):
        original = fab.KNOWN_PARITY_WARNING_COUNTS
        original_optional = fab.OPTIONAL_PARITY_WARNING_MAX
        with release.project_context(self.project, self.root, "1.2"):
            self.assertEqual(fab.PRO, str(self.project / "phantasm.kicad_pro"))
            self.assertEqual(fab.KNOWN_PARITY_WARNING_COUNTS,
                             {("lib_footprint_mismatch", "D_BUS"): 1})
            self.assertEqual(fab.OPTIONAL_PARITY_WARNING_MAX,
                             {("lib_footprint_mismatch", "J1"): 1})
        self.assertIs(fab.KNOWN_PARITY_WARNING_COUNTS, original)
        self.assertIs(fab.OPTIONAL_PARITY_WARNING_MAX, original_optional)

    def write_package(self):
        jlc = Path(fab.JLC)
        jlc.mkdir()
        for name in fab.ZIP_MEMBERS | set(fab.ASSEMBLY_MEMBERS):
            (jlc / name).write_text(name, encoding="utf-8")
        fab.write_upload_zip(str(jlc), sorted(fab.ZIP_MEMBERS), str(jlc / fab.ARCHIVE))
        (jlc / fab.SUMS_FILE).write_text(
            fab.package_manifest(str(jlc), sorted(fab.ZIP_MEMBERS), fab.ARCHIVE), encoding="utf-8")

    def pipeline(self, stack, fail=None, revision="1.2"):
        events = []
        def preflight(project, revision, check_copper=True):
            events.append("copper" if check_copper else "settings")
            if check_copper:
                self.assertEqual((project / "phantasm.kicad_pcb").read_text(), "refilled")
                if fail == "copper":
                    raise ValueError("injected copper failure")
        def refill(args, check):
            events.append("refill")
            self.assertIn("--refill-zones", args)
            self.assertIn("--save-board", args)
            Path(args[-1]).write_text("refilled", encoding="utf-8")
        expected_revision = revision
        def fabricate(revision, validate_only, refill_zones):
            events.append("fabricate")
            self.assertEqual(revision, expected_revision)
            self.assertFalse(refill_zones)
            self.assertEqual(Path(fab.PCB).read_text(), "refilled")
            if fail == "gate":
                (Path(fab.OUT) / "phantasm-drc.json").write_text("failed DRC")
                raise SystemExit("injected gate failure")
            if fail == "source_change":
                (self.project / "phantasm.kicad_pro").write_text("edit", encoding="utf-8")
            if not validate_only:
                self.write_package()
        for name, replacement in (("preflight", preflight), ("fabrication_options", lambda *args: {"finish": "ENIG"})):
            stack.enter_context(mock.patch.object(release, name, side_effect=replacement))
        stack.enter_context(mock.patch.object(release.subprocess, "run", side_effect=refill))
        stack.enter_context(mock.patch.object(release.subprocess, "check_output", return_value="10.0.4"))
        stack.enter_context(mock.patch.object(fab, "kicad_cli", return_value="fixture-cli"))
        stack.enter_context(mock.patch.object(fab, "fabricate", side_effect=fabricate))
        stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        return events

    def test_release_binds_input_evaluated_copper_and_artifacts(self):
        before = release.sources(self.project)
        with contextlib.ExitStack() as stack:
            events = self.pipeline(stack)
            release.build(self.project, "1.2", self.destination)
        self.assertEqual(events, ["settings", "refill", "copper", "fabricate"])
        record = json.loads((self.destination / "release.json").read_text())
        self.assertEqual(record["source_sha256"], before)
        self.assertNotEqual(record["source_sha256"], record["evaluated_source_sha256"])
        self.assertEqual(record["tools"]["kicad_cli"], "10.0.4")
        self.assertEqual(release.sources(self.project), before)
        self.assertEqual(len(release.verify(self.destination)), 64)
        (self.destination / "output" / "jlc" / "phantasm-BOM.csv").write_text("tampered")
        with self.assertRaisesRegex(ValueError, "contents do not match"):
            release.verify(self.destination)

    def test_manifest_records_finish_stackup_only_when_checked(self):
        for revision in ("1.1", "1.2"):
            with self.subTest(revision=revision), contextlib.ExitStack() as stack:
                self.pipeline(stack, revision=revision)
                destination = self.root / ("package-" + revision)
                release.build(self.project, revision, destination)
                record = json.loads((destination / "release.json").read_text())
                self.assertEqual("finish_stackup" in record["validation"], revision == "1.2")

    def test_failures_never_publish(self):
        for failure in ("copper", "gate", "source_change"):
            with self.subTest(failure=failure), contextlib.ExitStack() as stack:
                self.pipeline(stack, fail=failure)
                with self.assertRaisesRegex((ValueError, SystemExit), "injected|sources changed"):
                    release.build(self.project, "1.2", self.destination)
                self.assertFalse(self.destination.exists())

    def test_failed_validation_retains_reports_and_evaluated_board(self):
        for publish in (False, True):
            with self.subTest(publish=publish), contextlib.ExitStack() as stack:
                self.pipeline(stack, fail="gate")
                stage = self.root / ("publish-stage" if publish else "validate-stage")
                stage.mkdir()
                stack.enter_context(mock.patch.object(release.tempfile, "mkdtemp", return_value=str(stage)))
                stderr = stack.enter_context(contextlib.redirect_stderr(io.StringIO()))
                with self.assertRaisesRegex(SystemExit, "injected gate failure"):
                    release.build(self.project, "1.2", self.destination if publish else None)
                self.assertEqual((stage / "output" / "phantasm-drc.json").read_text(), "failed DRC")
                self.assertEqual((stage / "evaluated" / "phantasm.kicad_pcb").read_text(), "refilled")
                self.assertIn(f"Reports kept in {stage}", stderr.getvalue())
                self.assertFalse(self.destination.exists())

    def test_validation_creates_no_release(self):
        with contextlib.ExitStack() as stack:
            self.pipeline(stack)
            stage = self.root / "validate-stage"
            stage.mkdir()
            stack.enter_context(mock.patch.object(release.tempfile, "mkdtemp", return_value=str(stage)))
            release.build(self.project, "1.2")
        self.assertFalse(stage.exists())
        self.assertFalse(self.destination.exists())

    def test_release_enforces_cross_arm_width(self):
        shipped = Path(__file__).resolve().parents[2] / "1.1"
        for width, accepted in (("32", True), ("35", True), ("35.001", False), ("40", False)):
            with self.subTest(width=width), contextlib.ExitStack() as stack:
                project = self.root / ("width-" + width)
                release.snapshot(shipped, project)
                board = project / "phantasm.kicad_pcb"
                outline = "(end 58.28 32)"
                self.assertEqual(board.read_text(encoding="utf-8").count(outline), 1)
                board.write_text(board.read_text(encoding="utf-8").replace(
                    outline, f"(end 58.28 {width})"), encoding="utf-8")
                stack.enter_context(mock.patch.object(release.subprocess, "run"))
                stack.enter_context(mock.patch.object(release.subprocess, "check_output", return_value="10.0.4"))
                stack.enter_context(mock.patch.object(fab, "kicad_cli", return_value="fixture-cli"))
                stack.enter_context(mock.patch.object(
                    fab, "fabricate",
                    side_effect=lambda revision, validate_only, refill_zones: self.write_package()))
                stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                stderr = stack.enter_context(contextlib.redirect_stderr(io.StringIO()))
                destination = self.root / ("package-" + width)
                arguments = ["release", "--revision", "1.1", "--project", str(project),
                             "--output", str(destination)]
                if accepted:
                    release.main(arguments)
                    record = json.loads((destination / "release.json").read_text())
                    self.assertIn("mechanical_width", record["validation"])
                else:
                    with self.assertRaises(SystemExit) as caught:
                        release.main(arguments)
                    self.assertEqual(caught.exception.code, 1)
                    self.assertIn(f"mechanical_width: R-MECH-6: cross-arm width {width} mm exceeds 35 mm",
                                  stderr.getvalue())
                    self.assertFalse(destination.exists())

    def test_release_refuses_existing_or_nested_destination(self):
        with self.assertRaisesRegex(ValueError, "outside"):
            release.build(self.project, "1.2", self.project / "package")
        self.destination.mkdir()
        with self.assertRaisesRegex(ValueError, "already exists"):
            release.build(self.project, "1.2", self.destination)

    def test_nested_manifest_is_hashed(self):
        nested = self.root / "nested"
        nested.mkdir()
        (nested / "release.json").write_text("nested")
        self.assertIn("nested/release.json", release.inventory(self.root))

    def test_order_record_requires_exact_reviewed_package_and_options(self):
        with contextlib.ExitStack() as stack:
            self.pipeline(stack)
            release.build(self.project, "1.2", self.destination)
        digest = release.verify(self.destination)
        details = {"vendor": "fixture", "order_id": "123", "ordered_at": "2026-10-13",
                   "assembly_preview_reviewed": True, "mechanical_access_reviewed": True,
                   "fabrication_options": {"finish": "ENIG"}}
        with self.assertRaisesRegex(ValueError, "outside the release"):
            release.record_order(self.destination, digest, details, self.destination / "order.json")
        output = self.root / "order.json"
        with self.assertRaisesRegex(ValueError, "reviewed package"):
            release.record_order(self.destination, "bad", details, output)
        self.assertFalse(output.exists())
        with self.assertRaisesRegex(ValueError, "options differ"):
            release.record_order(self.destination, digest, {**details, "fabrication_options": {}}, output)
        release.record_order(self.destination, digest, details, output)
        self.assertEqual(json.loads(output.read_text())["release_sha256"], digest)
        with self.assertRaises(FileExistsError):
            release.record_order(self.destination, digest, details, output)

    @unittest.skipUnless(GENERATES, GENERATES_REASON)
    def test_kicad_refill_preserves_rev_1_2_fabrication_options(self):
        project = self.root / "refilled"
        release.snapshot(release.CANDIDATE.parent, project)
        before = release.fabrication_options(project, "1.2")
        subprocess.run([fab.kicad_cli(), "pcb", "drc", "--refill-zones", "--save-board",
                        "--format", "json", "-o", str(self.root / "drc.json"),
                        str(project / "phantasm.kicad_pcb")],
                       check=True, capture_output=True, text=True)
        self.assertEqual(release.fabrication_options(project, "1.2"), before)

    def test_valid_finish_and_stackup_fixture(self):
        source = release.CANDIDATE.parent
        for path in source.glob("phantasm.kicad_pcb"):
            shutil.copyfile(path, self.project / path.name)
        self.assertEqual(release.fabrication_options(self.project, "1.2")["copper_finish"], "ENIG")
        board = self.project / "phantasm.kicad_pcb"
        original = board.read_text()
        board.write_text(original.replace('(net_name "GND")', '(net_name "VCC")'))
        with self.assertRaisesRegex(ValueError, "GND pours"):
            release.fabrication_options(self.project, "1.2")
        board.write_text(original.replace('(thickness 1.6)', '(thickness 0.8)'))
        with self.assertRaisesRegex(ValueError, "stackup"):
            release.fabrication_options(self.project, "1.2")


if __name__ == "__main__":
    unittest.main()
