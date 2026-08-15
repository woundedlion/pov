import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import heal_clearance
from constraints import UNPLACED_DEFAULT_CLASS, UNPLACED_RULES


class MainTests(unittest.TestCase):
    def test_missing_projects_fails(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            with mock.patch.object(heal_clearance, "OUT", temp_dir):
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr):
                    self.assertEqual(heal_clearance.main(), 1)
        self.assertIn("no uploadable project files", stderr.getvalue())

    def test_malformed_project_fails_cleanly(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "phantasm.kicad_pro"
            project.write_text("{", encoding="utf-8")
            with mock.patch.object(heal_clearance, "OUT", temp_dir):
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr):
                    self.assertEqual(heal_clearance.main(), 1)
        self.assertIn("cannot process phantasm.kicad_pro", stderr.getvalue())

    def test_manifested_snapshot_is_left_untouched(self):
        zeroed = json.dumps(
            {"board": {"design_settings": {"rules": {"min_clearance": 0}}}})
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "phantasm.kicad_pro").write_text(zeroed, encoding="utf-8")
            snapshot = root / "quilter_incremental"
            snapshot.mkdir()
            protected = snapshot / "phantasm.kicad_pro"
            protected.write_text(zeroed, encoding="utf-8")
            (snapshot / "SHA256SUMS.txt").write_text("", encoding="utf-8")

            with mock.patch.object(heal_clearance, "OUT", temp_dir):
                self.assertEqual(
                    heal_clearance.project_files(), [str(root / "phantasm.kicad_pro")])
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(heal_clearance.main(), 0)

            self.assertEqual(protected.read_text(encoding="utf-8"), zeroed)
            healed = json.loads((root / "phantasm.kicad_pro").read_text(encoding="utf-8"))
            self.assertGreater(
                healed["board"]["design_settings"]["rules"]["min_clearance"], 0)


    def healed_bytes(self, source):
        """Heal a one-project tree from raw bytes; return the rewritten bytes."""
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "phantasm.kicad_pro"
            project.write_bytes(source)
            with mock.patch.object(heal_clearance, "OUT", temp_dir):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(heal_clearance.main(), 0)
            return project.read_bytes()

    ZEROED = json.dumps(
        {"board": {"design_settings": {"rules": {"min_clearance": 0}}}},
        indent=2)

    def test_crlf_project_stays_crlf(self):
        healed = self.healed_bytes(
            self.ZEROED.replace("\n", "\r\n").encode("utf-8"))

        self.assertTrue(healed.endswith(b"}\r\n"))
        self.assertNotIn(b"\n", healed.replace(b"\r\n", b""))

    def test_lf_project_stays_lf(self):
        healed = self.healed_bytes(self.ZEROED.encode("utf-8"))

        self.assertTrue(healed.endswith(b"}\n"))
        self.assertNotIn(b"\r", healed)

    def test_unplaced_project_gets_the_unplaced_constraints(self):
        zeroed = json.dumps({
            "board": {"design_settings": {"rules": {"min_clearance": 0}}},
            "net_settings": {"classes": [{"name": "Default"}]}})
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "unplaced").mkdir()
            project = root / "unplaced" / "phantasm_unplaced.kicad_pro"
            project.write_text(zeroed, encoding="utf-8")

            with mock.patch.object(heal_clearance, "OUT", temp_dir):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(heal_clearance.main(), 0)

            healed = json.loads(project.read_text(encoding="utf-8"))
            self.assertEqual(healed["board"]["design_settings"]["rules"],
                             dict(UNPLACED_RULES))
            default = healed["net_settings"]["classes"][0]
            for field, expected in UNPLACED_DEFAULT_CLASS.items():
                self.assertEqual(default[field], expected)


if __name__ == "__main__":
    unittest.main()
