import contextlib
import io
import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import check_coverage  # noqa: E402


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


def document_with(files, percent=90.0):
    """An llvm-cov export whose files are (name, count, covered) triples."""
    return {"data": [{
        "totals": {"lines": {"percent": percent}},
        "files": [{"filename": name,
                   "summary": {"lines": {"count": count, "covered": covered}}}
                  for name, count, covered in files],
    }]}


class DirectoryCoverage(unittest.TestCase):
    def test_aggregates_only_files_under_the_directory(self):
        document = document_with([
            ("/w/pov/core/render/scan.h", 100, 90),
            ("/w/pov/core/render/deep/raster.h", 100, 50),
            ("/w/pov/core/mesh/mesh.h", 100, 0),
        ])
        self.assertEqual(
            check_coverage.directory_coverage(document, "core/render"), 70.0)

    def test_rejects_a_directory_the_report_never_names(self):
        document = document_with([("/w/pov/core/mesh/mesh.h", 10, 10)])
        with self.assertRaisesRegex(ValueError, "core/render"):
            check_coverage.directory_coverage(document, "core/render")

    def test_rejects_a_directory_with_no_instrumented_line(self):
        # 0 of 0 lines would read as 100% and satisfy any floor.
        document = document_with([("/w/pov/core/render/scan.h", 0, 0)])
        with self.assertRaisesRegex(ValueError, "no instrumented line"):
            check_coverage.directory_coverage(document, "core/render")

    def test_rejects_a_report_with_no_per_file_summaries(self):
        with self.assertRaisesRegex(ValueError, "per-file summaries"):
            check_coverage.directory_coverage(
                {"data": [{"totals": {}}]}, "core/render")


class Main(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.report = Path(self.tmp.name) / "summary.json"

    def stage(self, document):
        self.report.write_text(json.dumps(document), encoding="utf-8")

    def stage_percentage(self, percent):
        self.stage({"data": [{"totals": {"lines": {"percent": percent}}}]})

    def main_with(self, minimum, *directories):
        """Invoke main() against the staged report, capturing its streams."""
        argv = ["check_coverage.py", str(self.report), "--min-lines", minimum]
        for floor in directories:
            argv += ["--min-directory", floor]
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

    def test_a_directory_below_its_floor_fails_a_passing_aggregate(self):
        self.stage(document_with([
            ("/w/pov/core/render/scan.h", 100, 99),
            ("/w/pov/core/mesh/mesh.h", 100, 10),
        ], percent=90.0))
        self.assertEqual(
            self.main_with("70", "core/render=90", "core/mesh=80"), 1)
        self.assertIn("core/mesh: 10.00% (minimum 80.00%)",
                      self.stdout.getvalue())

    def test_every_directory_on_its_floor_passes(self):
        self.stage(document_with([("/w/pov/core/mesh/mesh.h", 100, 80)]))
        self.assertEqual(self.main_with("70", "core/mesh=80"), 0)

    def test_an_uninstrumented_directory_is_a_tooling_error(self):
        self.stage(document_with([("/w/pov/core/mesh/mesh.h", 0, 0)]))
        with self.assertRaises(SystemExit) as raised:
            self.main_with("70", "core/mesh=80")
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("no instrumented line", self.stderr.getvalue())

    def test_rejects_a_malformed_directory_floor(self):
        self.stage_percentage(78.25)
        with self.assertRaises(SystemExit) as raised:
            self.main_with("70", "core/mesh")
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("<directory>=<percent>", self.stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
