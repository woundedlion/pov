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


if __name__ == "__main__":
    unittest.main()
