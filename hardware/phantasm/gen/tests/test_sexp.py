import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
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

    def test_rejects_unsupported_escape(self):
        with self.assertRaisesRegex(ValueError, r"unsupported escape: \\q"):
            sexp.parse('(a "value\\q")')


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


if __name__ == "__main__":
    unittest.main()
