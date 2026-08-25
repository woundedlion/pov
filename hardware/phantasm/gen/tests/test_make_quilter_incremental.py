import hashlib
import shutil
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import make_quilter_incremental


class SnapshotTests(unittest.TestCase):
    def test_digest_reads_a_crlf_file_as_the_lf_bytes_it_snapshots(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "snapshot.txt"
            path.write_bytes(b"first\r\nsecond\r\n")
            expected = hashlib.sha256(b"first\nsecond\n").hexdigest()
            self.assertEqual(make_quilter_incremental.snapshot_digest(path), expected)

    def test_snapshot_copy_and_manifest_use_lf_bytes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "source.txt"
            snapshot = root / "snapshot"
            snapshot.mkdir()
            source.write_bytes(b"first\r\nsecond\r\n")
            target = snapshot / "payload.txt"
            make_quilter_incremental.copy_snapshot_text(source, target)
            make_quilter_incremental.write_manifest(snapshot)

            payload = b"first\nsecond\n"
            digest = hashlib.sha256(payload).hexdigest().encode("ascii")
            self.assertEqual(target.read_bytes(), payload)
            self.assertEqual(
                (snapshot / "SHA256SUMS.txt").read_bytes(),
                digest + b"  payload.txt\n",
            )
            make_quilter_incremental.verify_snapshot(snapshot)
            sha256sum = shutil.which("sha256sum")
            if sha256sum:
                subprocess.run(
                    [sha256sum, "-c", "SHA256SUMS.txt"], cwd=snapshot,
                    check=True, capture_output=True, text=True,
                )

    def test_committed_snapshot_matches_manifest(self):
        make_quilter_incremental.verify_snapshot()

    def copy_snapshot(self, temp_dir):
        snapshot = Path(temp_dir) / "snapshot"
        shutil.copytree(make_quilter_incremental.OUTPUT, snapshot)
        return snapshot

    def assert_mismatch(self, snapshot, missing, extra):
        with self.assertRaises(RuntimeError) as caught:
            make_quilter_incremental.verify_snapshot(snapshot)
        self.assertEqual(
            str(caught.exception),
            f"snapshot manifest mismatch: missing={missing}, extra={extra}",
        )

    def test_unmanifested_footprint_is_reported_as_extra(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot = self.copy_snapshot(temp_dir)
            (snapshot / "phantasm.pretty" / "stale.kicad_mod").write_text(
                "stale", encoding="utf-8")
            self.assert_mismatch(
                snapshot, [], [str(Path("phantasm.pretty/stale.kicad_mod"))])

    def test_unmanifested_top_level_file_is_reported_as_extra(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot = self.copy_snapshot(temp_dir)
            (snapshot / "phantasm.kicad_dru").write_text("stray", encoding="utf-8")
            self.assert_mismatch(snapshot, [], ["phantasm.kicad_dru"])

    def test_deleted_project_file_is_reported_as_missing(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot = self.copy_snapshot(temp_dir)
            (snapshot / "sym-lib-table").unlink()
            self.assert_mismatch(snapshot, ["sym-lib-table"], [])

    def test_unmanaged_files_are_not_verified(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot = self.copy_snapshot(temp_dir)
            (snapshot / "phantasm.kicad_prl").write_text("local", encoding="utf-8")
            (snapshot / "phantasm.kicad_pcb-bak").write_text("backup", encoding="utf-8")
            (snapshot / "phantasm.pretty" / "fp-info-cache").write_text(
                "cache", encoding="utf-8")
            make_quilter_incremental.verify_snapshot(snapshot)

    def test_committed_snapshot_holds_no_strays(self):
        make_quilter_incremental.require_no_strays(make_quilter_incremental.OUTPUT)

    def test_stray_is_rejected_instead_of_adopted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot = self.copy_snapshot(temp_dir)
            (snapshot / "phantasm.kicad_dru").write_text("stray", encoding="utf-8")
            with self.assertRaises(RuntimeError) as caught:
                make_quilter_incremental.require_no_strays(snapshot)
            self.assertEqual(
                str(caught.exception),
                "snapshot holds files this run does not write: phantasm.kicad_dru",
            )

    def test_unmanaged_stray_is_accepted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot = self.copy_snapshot(temp_dir)
            (snapshot / "phantasm.kicad_prl").write_text("local", encoding="utf-8")
            make_quilter_incremental.require_no_strays(snapshot)

    def test_manifest_regenerates_byte_identically(self):
        manifest = make_quilter_incremental.OUTPUT / "SHA256SUMS.txt"
        rebuilt = "\n".join(
            f"{make_quilter_incremental.snapshot_digest(make_quilter_incremental.OUTPUT / path)}"
            f"  {path.as_posix()}"
            for path in make_quilter_incremental.snapshot_paths(
                make_quilter_incremental.OUTPUT)
        )
        self.assertEqual(manifest.read_text(encoding="utf-8"), rebuilt + "\n")


class MainTests(unittest.TestCase):
    def test_a_snapshot_fault_exits_with_its_message(self):
        with unittest.mock.patch.object(
                make_quilter_incremental, "verify_snapshot",
                side_effect=RuntimeError("snapshot hash mismatch: phantasm.kicad_pcb")), \
                unittest.mock.patch.object(sys, "argv",
                                           ["make_quilter_incremental.py", "--verify"]):
            with self.assertRaises(SystemExit) as caught:
                make_quilter_incremental.main()
        self.assertEqual(str(caught.exception),
                         "snapshot hash mismatch: phantasm.kicad_pcb")


if __name__ == "__main__":
    unittest.main()
