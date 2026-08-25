import contextlib
import io
import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

from tools import check_coverage


class LineCoverage(unittest.TestCase):
    def test_reads_aggregate_percentage(self):
        document = {"data": [{"totals": {"lines": {"percent": 78.25}}}]}
        self.assertEqual(check_coverage.line_coverage(document), 78.25)

    def test_rejects_missing_summary(self):
        with self.assertRaisesRegex(ValueError, "one data entry"):
            check_coverage.line_coverage({"data": []})

    def test_rejects_impossible_percentage(self):
        document = {"data": [{"totals": {"lines": {"percent": 101}}}]}
        with self.assertRaisesRegex(ValueError, "outside"):
            check_coverage.line_coverage(document)


class Main(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.report = Path(self.tmp.name) / "summary.json"

    def stage(self, document):
        self.report.write_text(json.dumps(document), encoding="utf-8")

    def stage_percentage(self, percent):
        self.stage({"data": [{"totals": {"lines": {"percent": percent}}}]})

    def main_with(self, minimum):
        """Invoke main() against the staged report, capturing its streams."""
        argv = ["check_coverage.py", str(self.report), "--min-lines", minimum]
        self.stdout, self.stderr = io.StringIO(), io.StringIO()
        with unittest.mock.patch.object(sys, "argv", argv), \
                contextlib.redirect_stdout(self.stdout), \
                contextlib.redirect_stderr(self.stderr):
            return check_coverage.main()

    def test_passes_above_the_floor(self):
        self.stage_percentage(78.25)
        self.assertEqual(self.main_with("70"), 0)
        self.assertEqual(self.stdout.getvalue(),
                         "line coverage: 78.25% (minimum 70.00%)\n")

    def test_passes_exactly_on_the_floor(self):
        self.stage_percentage(70.0)
        self.assertEqual(self.main_with("70"), 0)

    def test_fails_below_the_floor(self):
        self.stage_percentage(69.99)
        self.assertEqual(self.main_with("70"), 1)
        self.assertIn("line coverage: 69.99%", self.stdout.getvalue())

    def test_rejects_a_malformed_report(self):
        self.stage({"data": []})
        with self.assertRaises(SystemExit) as raised:
            self.main_with("70")
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("one data entry", self.stderr.getvalue())

    def test_rejects_unparsable_json(self):
        self.report.write_text("{not json", encoding="utf-8")
        with self.assertRaises(SystemExit) as raised:
            self.main_with("70")
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("error: Expecting property name", self.stderr.getvalue())

    def test_rejects_a_missing_report(self):
        with self.assertRaises(SystemExit) as raised:
            self.main_with("70")
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("summary.json", self.stderr.getvalue())

    def test_rejects_a_floor_outside_the_percentage_range(self):
        self.stage_percentage(78.25)
        with self.assertRaises(SystemExit) as raised:
            self.main_with("101")
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("between 0 and 100", self.stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
