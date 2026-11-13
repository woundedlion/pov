import contextlib
import io
import subprocess
from unittest import mock
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import doc_coverage as cov  # noqa: E402

ROOT = Path(tempfile.gettempdir()).resolve() / "repo"
PREFIX = ROOT.as_posix()


def warning(path, line, body):
    return f"{PREFIX}/{path}:{line}: warning: {body} is not documented."


def item(name, kind, scope):
    return cov.Undocumented("a.h", 1, name, kind, scope)


class TestParse(unittest.TestCase):
    def test_member_warning_yields_name_kind_and_scope(self):
        (found,) = cov.parse(warning(
            "core/a.h", 12,
            "Member edge_offset(const Point &p, int &axis) const (function) "
            "of struct SDF::WireLattice"), ROOT)
        self.assertEqual(found, cov.Undocumented(
            "core/a.h", 12, "edge_offset(const Point &p, int &axis) const",
            "function", "SDF::WireLattice"))


class TestExempt(unittest.TestCase):
    def test_ordinary_operator_is_not_exempt(self):
        self.assertFalse(cov.exempt(item(
            "operator+(const Vec &b) const", "function", "Vec")))
        self.assertFalse(cov.exempt(item(
            "operator==(const Vec &b) const", "function", "Vec")))


class TestScaffoldMacros(unittest.TestCase):
    def test_macro_undefined_by_an_included_rows_header_is_exempt(self):
        row = cov.Undocumented("t/a.h", 3, "ROW(name)", "macro definition", "a.h")
        other = cov.Undocumented("t/a.h", 5, "KEEP(name)", "macro definition", "a.h")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "t").mkdir()
            (root / "t/rows.h").write_text("ROW(x)\n#undef ROW\n", encoding="utf-8")
            (root / "t/a.h").write_text(
                '#define ROW(name) #name,\n#include "t/rows.h"\n'
                "#define KEEP(name) name\n", encoding="utf-8")
            self.assertEqual(cov.locations([row, other], root), [other])
            (root / "t/a.h").write_text(
                '#define ROW(name) #name,\n#include "rows.h"\n', encoding="utf-8")
            self.assertEqual(cov.locations([row], root), [])
            (root / "t/rows.h").write_text("ROW(x)\n", encoding="utf-8")
            self.assertEqual(cov.locations([row], root), [row])


class TestMain(unittest.TestCase):
    def test_tool_errors_are_distinct_from_coverage_failures(self):
        for error in (OSError("missing executable"), subprocess.CalledProcessError(1, "git")):
            with self.subTest(error=error), \
                 mock.patch.object(cov, "source_files", side_effect=error), \
                 contextlib.redirect_stderr(io.StringIO()) as errors:
                self.assertEqual(cov.main([]), 2)
                self.assertIn("tooling error", errors.getvalue())
        with mock.patch.object(cov, "source_files", return_value=["a.h"]), \
             mock.patch.object(cov, "missing_file_blocks", return_value=[]), \
             mock.patch.object(cov, "run_doxygen", side_effect=FileNotFoundError("doxygen")), \
             contextlib.redirect_stderr(io.StringIO()) as errors:
            self.assertEqual(cov.main([]), 2)
        self.assertIn("doxygen", errors.getvalue())

    def test_complete_input_passes_and_missing_file_block_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "a.h").write_text("/** @file */", encoding="utf-8")
            with mock.patch.object(cov, "source_files", return_value=["a.h"]), \
                 mock.patch.object(cov, "run_doxygen", return_value=""), \
                 contextlib.redirect_stdout(io.StringIO()), \
                 contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(cov.main(["--root", str(root)]), 0)
                (root / "a.h").write_text("", encoding="utf-8")
                self.assertEqual(cov.main(["--root", str(root)]), 1)


class TestSilentDeclarations(unittest.TestCase):
    def test_undocumented_empty_struct_and_concept_are_reported(self):
        self.assertEqual([(i.line, i.name, i.kind) for i in cov.silent_declarations(
            "a.h",
            "struct Tag {};\n"
            "template <typename T>\nconcept C = true;\n"
            "struct Base : std::false_type {};\n")], [
            (1, "Tag", "empty struct"), (3, "C", "concept"),
            (4, "Base", "empty struct")])


if __name__ == "__main__":
    unittest.main()
