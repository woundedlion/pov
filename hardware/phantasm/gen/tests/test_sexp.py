import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import sexp  # noqa: E402


class StringEscapeTests(unittest.TestCase):
    def test_round_trips_standard_escapes(self):
        source = '(a "line 1\\nline 2\\tindent\\rreturn\\\\slash\\\"quote")'

        parsed = sexp.parse(source)
        dumped = sexp.dumps(parsed[0])

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


class ParseTests(unittest.TestCase):
    def test_truncated_list_reports_token_offset(self):
        with self.assertRaisesRegex(ValueError, r"unexpected end of input at token 5"):
            sexp.parse("(root (child)")


if __name__ == "__main__":
    unittest.main()
