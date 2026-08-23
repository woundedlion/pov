import unittest

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


if __name__ == "__main__":
    unittest.main()
