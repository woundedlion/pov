import subprocess
import sys
import tempfile
import unittest
from pathlib import Path, PurePosixPath

TOOLS = Path(__file__).resolve().parent.parent
FIXTURES = Path(__file__).resolve().parent / "fixtures"
sys.path.insert(0, str(TOOLS))

import docs_check as dc  # noqa: E402


class TestDocumentationChecker(unittest.TestCase):
    def test_valid_links_images_references_and_fences(self):
        text = (FIXTURES / "valid.txt").read_text(encoding="utf-8")
        entries = {
            PurePosixPath("docs/target file.md"),
            PurePosixPath("docs/images"),
            PurePosixPath("docs/images/diagram.svg"),
        }
        self.assertEqual(dc.check_text(PurePosixPath("docs/readme.md"),
                                       text, entries), [])

    def test_unclosed_fence_is_reported_at_its_opening(self):
        text = (FIXTURES / "unclosed_fence.txt").read_text(encoding="utf-8")
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, set())
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].line, 3)
        self.assertIn("unclosed fenced code block", issues[0].message)

    def test_missing_encoded_link_ignores_query_and_fragment(self):
        text = (FIXTURES / "missing_link.txt").read_text(encoding="utf-8")
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, set())
        self.assertEqual(len(issues), 1)
        self.assertIn("docs/missing file.md", issues[0].message)

    def test_undefined_reference_link_is_reported(self):
        text = (FIXTURES / "undefined_reference.txt").read_text(encoding="utf-8")
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, set())
        self.assertEqual(len(issues), 1)
        self.assertIn("undefined reference link [missing]", issues[0].message)

    def test_missing_image_and_reference_targets_are_reported(self):
        text = (FIXTURES / "missing_targets.txt").read_text(encoding="utf-8")
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, set())
        self.assertEqual(len(issues), 2)
        self.assertIn("docs/missing.svg", issues[0].message)
        self.assertIn("docs/missing.md", issues[1].message)

    def test_only_git_tracked_markdown_is_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / "README.md").write_text("[ok](target.txt)\n", encoding="utf-8")
            (root / "target.txt").write_text("ok\n", encoding="utf-8")
            (root / "untracked.md").write_text("[bad](missing.txt)\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "README.md", "target.txt"],
                           check=True)
            markdown, issues = dc.check_repository(root)
        self.assertEqual(markdown, [PurePosixPath("README.md")])
        self.assertEqual(issues, [])


if __name__ == "__main__":
    unittest.main()
