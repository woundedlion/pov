import os
import sys
import unittest
from pathlib import Path
from unittest import mock

GEN = Path(__file__).resolve().parent.parent
ROUTED = GEN.parent / "phantasm.kicad_pcb"
sys.path.insert(0, str(GEN))

import sexp  # noqa: E402


def atom_kinds(node):
    """Mirror of node's shape with each leaf replaced by "sym" or "str"."""
    if isinstance(node, list):
        return [atom_kinds(child) for child in node]
    return "sym" if isinstance(node, sexp.Sym) else "str"


class StringEscapeTests(unittest.TestCase):
    def test_round_trips_standard_escapes(self):
        source = '(a "line 1\\nline 2\\tindent\\rreturn\\\\slash\\\"quote")'

        parsed = sexp.parse(source)
        dumped = sexp.dumps(parsed[0])

        self.assertEqual(dumped, source)
        self.assertEqual(sexp.parse(dumped), parsed)
        self.assertIn(r"\n", dumped)
        self.assertIn(r"\t", dumped)
        self.assertIn(r"\r", dumped)

    def test_rejects_trailing_backslash(self):
        with self.assertRaisesRegex(ValueError, "trailing backslash"):
            sexp.parse('(a "value\\')

    def test_rejects_unterminated_string(self):
        with self.assertRaisesRegex(ValueError, "unterminated string"):
            sexp.parse('(a "value')

    def test_accepts_kicad_c_hex_octal_and_unknown_escapes(self):
        source = r'(a "\a\b\f\v\x41\101\q")'

        parsed = sexp.parse(source)

        self.assertEqual(parsed[0][1], "\a\b\f\vAA" + r"\q")
        self.assertEqual(sexp.dumps(parsed[0]), source)

    def test_constructed_control_bytes_are_escaped(self):
        node = [sexp.Sym("a"), "\x01\x7f"]

        self.assertEqual(sexp.dumps(node), r'(a "\x01\x7f")')
        self.assertEqual(sexp.parse(sexp.dumps(node))[0], node)

    def test_hex_and_octal_escapes_form_utf8(self):
        source = r'(a "\xC3\xA9 \303\251")'

        parsed = sexp.parse(source)

        self.assertEqual(parsed[0][1], "é é")
        self.assertEqual(sexp.dumps(parsed[0]), source)


class QuotingTests(unittest.TestCase):
    """Sym subclasses str, so parse-tree equality alone cannot see quoting."""

    def test_parse_trees_compare_equal_across_quoting(self):
        self.assertEqual(sexp.parse("(layer F.Cu)"), sexp.parse('(layer "F.Cu")'))

    def test_bare_atoms_parse_as_sym(self):
        (node,) = sexp.parse("(layer F.Cu)")

        self.assertIsInstance(node[0], sexp.Sym)
        self.assertIsInstance(node[1], sexp.Sym)

    def test_quoted_strings_parse_as_plain_str(self):
        (node,) = sexp.parse('(property "Reference" "U1")')

        self.assertIsInstance(node[0], sexp.Sym)
        self.assertNotIsInstance(node[1], sexp.Sym)
        self.assertNotIsInstance(node[2], sexp.Sym)

    def test_dumps_leaves_bare_atoms_unquoted(self):
        self.assertEqual(sexp.dumps(sexp.parse("(layer F.Cu)")[0]), "(layer F.Cu)")

    def test_dumps_quotes_plain_strings(self):
        source = '(property "Reference" "U1")'

        self.assertEqual(sexp.dumps(sexp.parse(source)[0]), source)

    def test_dumps_text_distinguishes_quoting(self):
        self.assertNotEqual(sexp.dumps(sexp.parse("(layer F.Cu)")[0]),
                            sexp.dumps(sexp.parse('(layer "F.Cu")')[0]))

    def test_dumps_quotes_constructed_str_but_not_sym(self):
        node = [sexp.Sym("at"), sexp.Sym("0"), "0"]

        self.assertEqual(sexp.dumps(node), '(at 0 "0")')

    def test_nested_dumps_preserves_quoting(self):
        source = ('(footprint "lib:part"\n'
                  "\t(layer F.Cu)\n"
                  "\t(attr smd)\n"
                  '\t(property "Ref" "U1"\n'
                  "\t\t(at 0 0 90) hide\n"
                  "\t)\n"
                  ")")

        self.assertEqual(sexp.dumps(sexp.parse(source)[0]), source)

    def test_deep_copy_round_trip_preserves_atom_kinds(self):
        source = '(footprint "lib:part" (layer F.Cu) (property "Ref" "U1") hide)'
        original = sexp.parse(source)[0]

        copied = sexp.parse(sexp.dumps(original))[0]

        self.assertEqual(atom_kinds(copied), atom_kinds(original))
        self.assertEqual(atom_kinds(copied),
                         ["sym", "str", ["sym", "sym"], ["sym", "str", "str"], "sym"])
        self.assertEqual(sexp.dumps(copied), sexp.dumps(original))


class ParseTests(unittest.TestCase):
    def test_truncated_list_reports_token_offset(self):
        with self.assertRaisesRegex(ValueError, r"unexpected end of input at token 5"):
            sexp.parse("(root (child)")


class ParseOneTests(unittest.TestCase):
    """The single-document form callers use, so no input raises IndexError."""

    def test_returns_the_only_document(self):
        self.assertEqual(sexp.parse_one("(root a)"), ["root", "a"])

    def test_whitespace_only_input_raises_value_error(self):
        with self.assertRaisesRegex(ValueError, "expected one document, found 0"):
            sexp.parse_one("   ")

    def test_a_second_document_raises_value_error(self):
        with self.assertRaisesRegex(ValueError, "expected one document, found 2"):
            sexp.parse_one("(a) (b)")


class FormattingTests(unittest.TestCase):
    def test_shipped_board_round_trip_preserves_every_byte(self):
        source = ROUTED.read_text(encoding="utf-8")

        dumped = sexp.dumps(sexp.parse(source)[0]) + "\n"

        self.assertEqual(dumped, source)

    def test_edit_preserves_the_surrounding_board_layout(self):
        source = ROUTED.read_text(encoding="utf-8")
        root = sexp.parse(source)[0]
        version = next(child for child in root
                       if isinstance(child, list) and child[0] == "version")
        original_version = str(version[1])
        version[1] = sexp.Sym("99999999")

        dumped = sexp.dumps(root) + "\n"
        changed = [(before, after) for before, after in
                   zip(source.splitlines(), dumped.splitlines())
                   if before != after]

        self.assertEqual(changed, [(f"\t(version {original_version})",
                                    "\t(version 99999999)")])


class FindKicadDataDirTests(unittest.TestCase):
    """Stock data directories must come from the same major as the pinned CLI."""

    def windows_install(self, version):
        return rf"C:\Program Files\KiCad\{version}\share\kicad\symbols"

    def resolve(self, hits):
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(sexp.glob, "glob",
                                  side_effect=lambda p: hits if "Program" in p
                                  and "x86" not in p else []):
            return sexp.find_kicad_data_dir("symbols", "KICAD_SYMBOL_DIR")

    def test_skips_an_install_of_another_major(self):
        newer = self.windows_install(f"{sexp.KICAD_MAJOR + 1}.0")

        self.assertEqual(self.resolve([newer]), "symbols")

    def test_prefers_the_pinned_major_over_a_newer_install(self):
        pinned = self.windows_install(f"{sexp.KICAD_MAJOR}.0")
        newer = self.windows_install(f"{sexp.KICAD_MAJOR + 1}.0")

        self.assertEqual(self.resolve([pinned, newer]), pinned)

    def test_takes_the_newest_minor_of_the_pinned_major(self):
        installs = [self.windows_install(f"{sexp.KICAD_MAJOR}.{minor}")
                    for minor in (0, 3)]

        self.assertEqual(self.resolve(installs), installs[1])

    def test_an_unversioned_prefix_names_no_major_and_is_taken(self):
        self.assertEqual(self.resolve(["/usr/share/kicad/symbols"]),
                         "/usr/share/kicad/symbols")

    def test_env_override_is_not_version_checked(self):
        with mock.patch.dict(os.environ, {"KICAD_SYMBOL_DIR": str(GEN)}):
            self.assertEqual(
                sexp.find_kicad_data_dir("symbols", "KICAD_SYMBOL_DIR"), str(GEN))


if __name__ == "__main__":
    unittest.main()
