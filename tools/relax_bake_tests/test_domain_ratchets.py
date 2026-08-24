import tempfile
import unittest
from pathlib import Path

from tools import check_domain_ratchets


class ParseFloors(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, name: str, text: str) -> Path:
        path = Path(self.tmp.name) / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_rejects_harness_without_floor(self):
        path = self.write("harness.cpp", "")
        with self.assertRaisesRegex(SystemExit, "no domain coverage floor parsed"):
            check_domain_ratchets.floors(
                path, check_domain_ratchets.HARNESS_FLOOR_RE
            )

    def test_rejects_death_harness_without_floor(self):
        path = self.write("death.h", "")
        with self.assertRaisesRegex(SystemExit, "no domain coverage floor parsed"):
            check_domain_ratchets.death_pins(path)


if __name__ == "__main__":
    unittest.main()
