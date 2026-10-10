import subprocess
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

    def test_macro_kind_may_contain_spaces(self):
        (found,) = cov.parse(warning(
            "core/r.h", 64,
            "Member HS_RESOLUTIONS(X) (macro definition) of file r.h"), ROOT)
        self.assertEqual((found.kind, found.scope), ("macro definition", "r.h"))

    def test_compound_warning_has_no_scope(self):
        (found,) = cov.parse(warning("core/a.h", 3, "Compound SDF::Plane"),
                             ROOT)
        self.assertEqual(found, cov.Undocumented(
            "core/a.h", 3, "SDF::Plane", "compound", ""))

    def test_windows_separators_are_normalized_to_relative_paths(self):
        (found,) = cov.parse(
            f"{PREFIX}/core/a.h:7: warning: Compound X is not documented."
            .replace("/", "\\"), ROOT)
        self.assertEqual((found.path, found.line), ("core/a.h", 7))

    def test_other_source_warning_keeps_its_text(self):
        (found,) = cov.parse(
            f"{PREFIX}/a.h:1: warning: argument 'x' of command @param is not "
            "found", ROOT)
        self.assertEqual(found, cov.Undocumented(
            "a.h", 1, "argument 'x' of command @param is not found",
            "warning", ""))

    def test_indented_continuation_joins_the_warning(self):
        (found,) = cov.parse(
            f"{PREFIX}/a.h:9: warning: The following parameters of f(int a) "
            "are not documented:\n  parameter 'a'", ROOT)
        self.assertEqual(found.name, "The following parameters of f(int a) "
                         "are not documented: parameter 'a'")

    def test_warnings_without_a_source_location_are_ignored(self):
        self.assertEqual(cov.parse(
            "warning: No output formats selected!", ROOT), [])

    def test_source_warnings_are_never_exempt(self):
        self.assertFalse(cov.exempt(cov.Undocumented(
            "a.h", 1, "operator=(int) =default", "warning", "")))


class TestExempt(unittest.TestCase):
    def test_defaulted_and_deleted_members(self):
        self.assertTrue(cov.exempt(item(
            "operator==(const Axis &) const =default", "function", "Axis")))
        self.assertTrue(cov.exempt(item(
            "Span(const Span &)=delete", "function", "Span")))
        self.assertTrue(cov.exempt(item(
            "operator=(Storage &&other) noexcept", "function", "Storage")))

    def test_ordinary_operator_is_not_exempt(self):
        self.assertFalse(cov.exempt(item(
            "operator+(const Vec &b) const", "function", "Vec")))
        self.assertFalse(cov.exempt(item(
            "operator==(const Vec &b) const", "function", "Vec")))

    def test_container_alias(self):
        self.assertTrue(cov.exempt(item("value_type", "typedef", "Ring")))
        self.assertFalse(cov.exempt(item("value_type", "variable", "Ring")))

    def test_operator_contract_member_only_on_direct_op_models(self):
        self.assertTrue(cov.exempt(item(
            "prepare(const FrameContext &, const Params &p)", "function",
            "Pullback::Interp::Op::WarpVectorNoise")))
        self.assertFalse(cov.exempt(item(
            "prepare(const FrameContext &)", "function",
            "Pullback::Interp::Op::WarpVectorNoise::Inner")))
        self.assertFalse(cov.exempt(item(
            "strength", "variable", "Pullback::Interp::Op::WarpVectorNoise")))

    def test_spec_override_but_not_the_base_spec(self):
        self.assertTrue(cov.exempt(item("HUE", "variable", "AshCloudSpec")))
        self.assertFalse(cov.exempt(item("HUE", "variable", "Pullback::Spec")))
        self.assertFalse(cov.exempt(item("SPEED", "variable", "AshCloudSpec")))

    def test_trait_specialization_type_only(self):
        self.assertTrue(cov.exempt(item(
            "Type", "typedef", "SourcePolicyFor< GridSourceParams, B >")))
        self.assertFalse(cov.exempt(item("Type", "typedef", "SourcePolicyFor")))


class TestLocations(unittest.TestCase):
    def test_instantiations_collapse_to_one_location(self):
        items = [
            cov.Undocumented("a.h", 4, "Params", "typedef", "Effect< A >"),
            cov.Undocumented("a.h", 4, "Params", "typedef", "Effect< B >"),
            cov.Undocumented("a.h", 2, "x", "variable", "S"),
        ]
        self.assertEqual([(i.path, i.line) for i in cov.locations(items)],
                         [("a.h", 2), ("a.h", 4)])

    def test_macro_undefined_in_its_own_file_is_exempt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "a.h").write_text(
                "#define ROW(x) x,\n#define KEEP 1\n#undef ROW\n",
                encoding="utf-8")
            items = [
                cov.Undocumented("a.h", 1, "ROW(x)", "macro definition", "a.h"),
                cov.Undocumented("a.h", 2, "KEEP", "macro definition", "a.h"),
            ]
            self.assertEqual([i.name for i in cov.locations(items, root)],
                             ["KEEP"])

    def test_exempt_items_are_dropped(self):
        self.assertEqual(cov.locations(
            [cov.Undocumented("a.h", 1, "value_type", "typedef", "R")]), [])



class TestDoxyfileSetting(unittest.TestCase):
    def test_continuation_lines_and_append(self):
        text = ("INPUT = core \\\n        effects\nOTHER = x\n"
                "INPUT += targets\n")
        self.assertEqual(cov.doxyfile_setting(text, "INPUT"),
                         ["core", "effects", "targets"])

    def test_reassignment_replaces(self):
        self.assertEqual(cov.doxyfile_setting(
            "EXCLUDE = a\nEXCLUDE = b\n", "EXCLUDE"), ["b"])


class TestFileBlocks(unittest.TestCase):
    def test_files_without_a_file_command_are_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "a.h").write_text("/** @file a.h */\n", encoding="utf-8")
            (root / "b.h").write_text("/** \\file b.h */\n", encoding="utf-8")
            (root / "c.h").write_text("// profile.h\n", encoding="utf-8")
            self.assertEqual(
                cov.missing_file_blocks(root, ["a.h", "b.h", "c.h"]), ["c.h"])

    def test_source_files_apply_patterns_and_excludes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / "Doxyfile").write_text(
                "INPUT = src\nFILE_PATTERNS = *.h\nEXCLUDE = src/gen\n",
                encoding="utf-8")
            (root / "src" / "gen").mkdir(parents=True)
            for name in ("src/a.h", "src/b.md", "src/gen/c.h"):
                (root / name).write_text("x\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "."], check=True)
            self.assertEqual(cov.source_files(root, []), ["src/a.h"])


if __name__ == "__main__":
    unittest.main()
