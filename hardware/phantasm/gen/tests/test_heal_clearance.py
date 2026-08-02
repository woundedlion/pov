import contextlib
import io
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


if __name__ == "__main__":
    unittest.main()
