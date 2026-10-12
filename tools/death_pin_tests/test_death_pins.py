"""Controls for literal death breadcrumb uniqueness."""

import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
import check_death_pins as checker  # noqa: E402


class DeathPinUniqueness(unittest.TestCase):
    def check_fixture(self, source, pin):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "tests").mkdir()
            (root / "core").mkdir()
            (root / "tests/guard_directories.txt").write_text("core\n", encoding="utf-8")
            (root / "tests/test_death.h").write_text(
                '{"ordinary_case", case_fixture, "core/guard.h", "' + pin + '"}',
                encoding="utf-8")
            (root / "core/guard.h").write_text(source, encoding="utf-8")
            with mock.patch.object(sys, "argv", ["check", str(root)]), contextlib.redirect_stdout(io.StringIO()) as out:
                checker.main()
            return out.getvalue()

    def test_non_opleg_duplicate_guard_is_rejected(self):
        source = 'void first() { HS_CHECK(value > 0, "probe: positive"); }'
        pin = "(value > 0) probe: positive"
        self.assertIn("1 literal cases", self.check_fixture(source, pin))
        duplicate = source + 'void second() { HS_CHECK(value > 0, "probe: positive"); }'
        with self.assertRaisesRegex(SystemExit, "ordinary_case: pin identifies 2 guard sites"):
            self.check_fixture(duplicate, pin)

    def test_missing_guard_is_rejected(self):
        pin = "(value > 0) probe: positive"
        self.check_fixture('HS_CHECK(value > 0, "probe: positive");', pin)
        with self.assertRaisesRegex(SystemExit, "pin identifies 0 guard sites"):
            self.check_fixture("", pin)

    def test_formatted_guard_is_unique_for_concrete_and_prefix_pins(self):
        source = 'HS_CHECK(index < count, "index %lu outside [0, %lu)", index, count);'
        for pin in ("(index < count) index 5 outside [0, 2)", "(index < count) index "):
            with self.subTest(pin=pin):
                self.check_fixture(source, pin)
                with self.assertRaisesRegex(SystemExit, "pin identifies 2 guard sites"):
                    self.check_fixture(source + source, pin)

    def test_modulo_condition_is_not_a_printf_conversion(self):
        source = 'HS_CHECK(size % ALIGN == 0, "aligned %lu bytes", size);'
        self.check_fixture(source, "(size % ALIGN == 0) aligned ")
        self.assertFalse(checker.guard_matches("(size % ALIGN == 0) aligned %lu bytes", "(size % OTHER == 0) aligned 4 bytes"))

    def test_direct_cold_failure_is_counted(self):
        source = 'hs::check_fail(HS_CHECK_SITE("false"), "allocation failed");'
        pin = "(false) allocation failed"
        self.check_fixture(source, pin)
        with self.assertRaisesRegex(SystemExit, "pin identifies 2 guard sites"):
            self.check_fixture(source + source.replace("hs::", ""), pin)


if __name__ == "__main__":
    unittest.main()
