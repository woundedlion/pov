"""Controls for literal death breadcrumb uniqueness."""

import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
import check_death_pins as checker


class DeathPinUniqueness(unittest.TestCase):
    def check_fixture(self, source, pin, *, extra_sources=None, table=None):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "tests").mkdir()
            (root / "core").mkdir()
            (root / "tests/guard_directories.txt").write_text("core\n", encoding="utf-8")
            (root / "tests/test_death.h").write_text(
                table if table is not None else
                '{"ordinary_case", case_ordinary_case, "core/guard.h", "' + pin + '"}',
                encoding="utf-8")
            (root / "core/guard.h").write_text(source, encoding="utf-8")
            for name, text in (extra_sources or {}).items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8")
            with mock.patch.object(sys, "argv", ["check", str(root)]), contextlib.redirect_stdout(io.StringIO()) as out:
                checker.main()
            return out.getvalue()

    def test_duplicate_guard_is_rejected(self):
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

    def test_same_named_guard_in_another_directory_is_rejected(self):
        source = 'HS_CHECK(value > 0, "probe: positive");'
        pin = "(value > 0) probe: positive"
        table = '{"ordinary_case", case_ordinary_case, "core/a/guard.h", "' + pin + '"}'
        self.check_fixture("", pin, table=table, extra_sources={"core/a/guard.h": source})
        with self.assertRaisesRegex(SystemExit, "also matches core/b/guard.h"):
            self.check_fixture("", pin, table=table, extra_sources={
                "core/a/guard.h": source, "core/b/guard.h": source})

    def test_malformed_row_is_rejected(self):
        source = 'HS_CHECK(value > 0, "probe: positive");'
        pin = "(value > 0) probe: positive"
        self.check_fixture(source, pin)
        with self.assertRaisesRegex(ValueError, "parsed 0 of 1 death-case rows"):
            self.check_fixture(source, pin, table='{"ordinary_case", case_ordinary_case, "core/guard.h", nullptr}')

    def test_raw_assert_is_rejected(self):
        pin = "(value > 0) probe: positive"
        guard = 'HS_CHECK(value > 0, "probe: positive");\n'
        self.check_fixture(guard + '// assert(x)\nstatic_assert(true);\nauto s = "assert(x)";\n', pin)
        for raw in ("assert(value);\n", "#include <cassert>\n", "#include <assert.h>\n"):
            with self.subTest(raw=raw), self.assertRaisesRegex(SystemExit, r"core/guard.h:2"):
                self.check_fixture(guard + raw, pin)
        with self.assertRaisesRegex(SystemExit, r"core/guard.h:3"):
            self.check_fixture(guard + "auto quote = '\"';\nassert(value);\n", pin)

    def test_empty_table_is_rejected(self):
        source = 'HS_CHECK(value > 0, "probe: positive");'
        pin = "(value > 0) probe: positive"
        self.check_fixture(source, pin)
        with self.assertRaisesRegex(SystemExit, "no literal death pins found"):
            self.check_fixture(source, pin, table="")


class NewGuardsNeedACase(unittest.TestCase):
    GUARD = 'HS_CHECK(value > 0, "probe: positive");\n'
    PIN = '{"probe_case", case_probe_case, "core/guard.h", "(value > 0) probe: positive"}'
    OTHER = '{"other_case", case_other_case, "core/other.h", "(true) other"}'

    def setUp(self):
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)
        self.git("init", "--quiet")
        self.commit({"tests/guard_directories.txt": "core\n",
                     "tests/test_death.h": self.OTHER,
                     "core/other.h": 'HS_CHECK(true, "other");\n'})
        self.base = self.git("rev-parse", "HEAD")

    def git(self, *args):
        return subprocess.run(["git", "-C", str(self.repo), "-c", "user.name=t", "-c", "user.email=t@t",
                               *args], check=True, capture_output=True, text=True).stdout.strip()

    def commit(self, files):
        for name, text in files.items():
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8", newline="\n")
        self.git("add", "-A")
        self.git("commit", "--quiet", "-m", "change")
        return self.git("rev-parse", "HEAD")

    def missing(self, files):
        return checker.new_guards(self.repo, self.base, self.commit(files))

    def test_added_guard_without_a_case_is_reported(self):
        self.assertEqual(self.missing({"core/guard.h": self.GUARD}),
                         ["core/guard.h:1: (value > 0) probe: positive"])

    def test_added_guard_with_a_case_passes(self):
        self.assertEqual(self.missing({"core/guard.h": self.GUARD,
                                       "tests/test_death.h": self.OTHER + "\n" + self.PIN}), [])

    def test_guard_text_already_in_the_base_passes_when_moved(self):
        self.base = self.commit({"core/old.h": self.GUARD})
        (self.repo / "core/old.h").unlink()
        self.assertEqual(self.missing({"core/guard.h": self.GUARD}), [])

    def test_waiver_needs_a_reason(self):
        self.assertEqual(self.missing({"core/guard.h": "// no death case: device-only\n" + self.GUARD}), [])
        self.assertEqual(len(self.missing({"core/guard.h": "// no death case:\n" + self.GUARD})), 1)

    def test_audit_check_needs_no_case(self):
        self.assertEqual(self.missing({"core/guard.h": self.GUARD.replace("HS_CHECK", "HS_AUDIT_CHECK")}), [])

    def test_reworded_guard_counts_as_new(self):
        self.base = self.commit({"core/guard.h": self.GUARD})
        self.assertEqual(len(self.missing({"core/guard.h": self.GUARD.replace("positive", "above zero")})), 1)

    def test_explicit_rename_of_existing_guard_passes(self):
        self.base = self.commit({"core/guard.h": self.GUARD})
        aliases = [("core/guard.h", "(value > 0) probe: positive",
                    "(renamed_value > 0) probe: positive")]
        with mock.patch.object(checker, "GUARD_RENAME_PINS", aliases):
            self.assertEqual(self.missing({"core/guard.h": self.GUARD.replace("value", "renamed_value")}), [])

    def test_command_line_fails_on_a_missing_case(self):
        tip = self.commit({"core/guard.h": self.GUARD})
        with mock.patch.object(sys, "argv", ["check", str(self.repo), "--new-guards", self.base, tip]), \
                contextlib.redirect_stdout(io.StringIO()) as out, \
                self.assertRaisesRegex(SystemExit, "1 new guard"):
            checker.main()
        self.assertIn("core/guard.h:1:", out.getvalue())


class GuardRenamePins(unittest.TestCase):
    def match(self, old, new, before=None, after=None):
        file = 'core/guard.h'
        with mock.patch.object(checker, 'GUARD_RENAME_PINS', [(file, old, new)]):
            return checker.renamed_guards(
                before if before is not None else 'HS_CHECK' + old.replace(') positive', ', "positive");'),
                after if after is not None else 'HS_CHECK' + new.replace(') positive', ', "positive");'), file)

    def test_identifier_rename(self):
        self.assertEqual(self.match('(VALUE > 0) positive', '(value > 0) positive'), {'(value>0)positive'})

    def test_missing_base_guard(self):
        self.assertEqual(self.match('(VALUE > 0) positive', '(value > 0) positive', before=''), set())

    def test_unlisted_guard(self):
        with mock.patch.object(checker, 'GUARD_RENAME_PINS', []):
            self.assertEqual(checker.renamed_guards('HS_CHECK(VALUE > 0, "positive");',
                                                   'HS_CHECK(value > 0, "positive");', 'core/guard.h'), set())

    def test_added_guard(self):
        self.assertEqual(self.match('(VALUE > 0) positive', '(value > 0) positive',
                                   after='HS_CHECK(value > 0, "positive"); HS_CHECK(true, "extra");'), set())

    def test_operator_number_or_message_change(self):
        for new in ('(value >= 0) positive', '(value > 1) positive', '(value > 0) changed'):
            with self.subTest(new=new):
                self.assertEqual(self.match('(VALUE > 0) positive', new,
                                           after='HS_CHECK' + new.replace(') positive', ', "positive");').replace(') changed', ', "changed");')), set())

    def test_inconsistent_identifier_mapping(self):
        self.assertEqual(self.match('(VALUE > 0 && VALUE < 4) positive',
                                    '(lower > 0 && upper < 4) positive'), set())



if __name__ == "__main__":
    unittest.main()
