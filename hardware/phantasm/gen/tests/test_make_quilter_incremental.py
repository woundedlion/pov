import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import make_quilter_incremental


class SnapshotTests(unittest.TestCase):
    def test_committed_snapshot_matches_manifest(self):
        make_quilter_incremental.verify_snapshot()

    def test_unmanifested_footprint_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot = Path(temp_dir) / "snapshot"
            shutil.copytree(make_quilter_incremental.OUTPUT, snapshot)
            (snapshot / "phantasm.pretty" / "stale.kicad_mod").write_text(
                "stale", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "snapshot manifest mismatch"):
                make_quilter_incremental.verify_snapshot(snapshot)


if __name__ == "__main__":
    unittest.main()
