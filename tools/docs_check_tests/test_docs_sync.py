"""Explicit documentation synchronization preserves authored documentation and validates its output."""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path, PurePosixPath
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import docs_check as dc
import docs_sync as ds
from git_test_env import git_run


def entries(*names):
    result = {PurePosixPath(name) for name in names}
    result.update(parent for path in tuple(result) for parent in path.parents
                  if parent != PurePosixPath("."))
    return result


class TreeSync(unittest.TestCase):
    def sync(self, text, paths, checkouts=None):
        with contextlib.redirect_stdout(io.StringIO()):
            return ds.sync_trees(text, paths, checkouts or {})

    def test_new_paths_require_authored_descriptions(self):
        text = "<!-- docs-check: tree exhaustive -->\n```\n```\n"
        with self.assertRaisesRegex(ValueError, "role description for new.h"):
            self.sync(text, entries("new.h"))

    def test_preserves_prose_descriptions_and_spacing_and_is_idempotent(self):
        text = ("# Overview\n\nAuthored prose.\n\n<!-- docs-check: tree exhaustive -->\n```\n"
                "├── core/          Engine description\n"
                "│   ├── gone.h     Removed implementation\n"
                "│   └── live.h     Authored description\n"
                "│                    continuation retained\n│\n"
                "└── docs/          A compact overview\n```\n\nMore prose.\n")
        paths = entries("core/live.h", "docs/guide.md")
        after = self.sync(text, paths)
        self.assertIn("Authored description\n│                    continuation retained\n", after)
        self.assertIn("core/          Engine description", after)
        self.assertIn("docs/          A compact overview", after)
        self.assertNotIn("gone.h", after)
        self.assertNotIn("guide.md", after)
        self.assertTrue(after.startswith("# Overview\n\nAuthored prose.\n\n"))
        self.assertTrue(after.endswith("\n\nMore prose.\n"))
        self.assertEqual(after, self.sync(after, paths))
        self.assertEqual(dc.check_text(PurePosixPath("README.md"), after, paths), [])

    def test_underindented_continuation_fails_without_truncating_text(self):
        text = ("<!-- docs-check: tree exhaustive -->\n```\n"
                "└── core/\n    └── file.h   Keep this description\n"
                "│  under-indented words that matter\n```\n")
        with self.assertRaisesRegex(ValueError, "tree continuation is indented less than its row"):
            self.sync(text, entries("core/file.h"))
        valid = text.replace("│  under-indented", "        under-indented")
        self.assertEqual(self.sync(valid, entries("core/file.h")), valid)


class RepositorySync(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.git("init", "-q")

    def git(self, *args):
        return git_run(self.root, *args).stdout.strip()

    def test_roster_counts_sync_without_rewriting_surrounding_prose(self):
        header = self.root / dc.EFFECT_ROSTER_SOURCE
        header.parent.mkdir(parents=True)
        header.write_text("#define HS_EFFECT_LIST(X) \\\n    X(One) \\\n    X(Two)\n", encoding="utf-8")
        playlist = self.root / dc.PHANTASM_PLAYLIST_SOURCE
        playlist.parent.mkdir(parents=True)
        playlist.write_text("#define HS_PHANTASM_EFFECT_LIST(X) \\\n    X(One, 10)\n", encoding="utf-8")
        readme = self.root / "README.md"
        readme.write_text("The playlist contains 99 effects today.\n", encoding="utf-8")
        self.git("add", "README.md", str(dc.EFFECT_ROSTER_SOURCE), str(dc.PHANTASM_PLAYLIST_SOURCE))
        with contextlib.redirect_stdout(io.StringIO()):
            ds.sync_repository(self.root, {}, {})
        self.assertEqual(readme.read_text(encoding="utf-8"), "The playlist contains 1 effects today.\n")


class AutoCheckout(unittest.TestCase):
    def test_auto_checkout_failure_reports_the_missing_pin(self):
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(ds, "discover_daydream", return_value=None):
                error = io.StringIO()
                with contextlib.redirect_stderr(error):
                    self.assertEqual(dc.main(["--root", directory, "--auto-checkout"]), 2)
                self.assertIn("pinned daydream checkout unavailable", error.getvalue())


if __name__ == "__main__":
    unittest.main()
