"""Build-time synchronization preserves authored documentation and validates its output."""

import contextlib
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path, PurePosixPath
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import docs_check as dc  # noqa: E402
import docs_sync as ds  # noqa: E402


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

    def test_grouped_names_globs_and_implicit_parents_are_preserved(self):
        text = ("<!-- docs-check: tree exhaustive -->\n```\n"
                "├── core/\n│   ├── one.h / two.h   Shared description\n"
                "│   └── *.cpp           Sources\n"
                "└── .github/workflows/  CI\n```\n")
        paths = entries("core/two.h", "core/a.cpp", ".github/workflows/ci.yml")
        after = self.sync(text, paths)
        self.assertIn("two.h           Shared description", after)
        self.assertNotIn("one.h", after)
        self.assertIn("*.cpp", after)
        self.assertNotIn("a.cpp", after)
        self.assertNotIn("── .github/\n", after)
        self.assertEqual(after, self.sync(after, paths))
        self.assertEqual(dc.check_text(PurePosixPath("README.md"), after, paths), [])

    def test_missing_sibling_is_left_untouched(self):
        text = "<!-- docs-check: tree daydream exhaustive -->\n```\n└── old.js   Keep this description\n```\n"
        self.assertEqual(self.sync(text, set()), text)

    def test_multiline_child_lists_keep_layout_without_duplicate_rows(self):
        text = ("<!-- docs-check: tree exhaustive -->\n```\n"
                "└── core/  Helpers (one,\n"
                "               two)\n```\n")
        paths = entries("core/one.h", "core/two.h")
        self.assertEqual(self.sync(text, paths), text)

    def test_broken_fences_and_links_still_fail(self):
        text = "[bad](missing.md)\n<!-- docs-check: tree exhaustive -->\n```\n└── old.h\n"
        after = self.sync(text, entries("new.h"))
        self.assertEqual(after, text)
        issues = dc.check_text(PurePosixPath("README.md"), after, entries("new.h"))
        self.assertTrue(any("unclosed" in issue.message for issue in issues))
        self.assertTrue(any("missing" in issue.message for issue in issues))

    def test_invalid_tree_indentation_fails_without_rewriting(self):
        text = "<!-- docs-check: tree exhaustive -->\n```\n│   └── file.h\n```\n"
        with self.assertRaises(ValueError):
            self.sync(text, entries("file.h"))


class RepositorySync(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.git("init", "-q")

    def git(self, *args):
        return subprocess.check_output(["git", "-C", str(self.root), *args], text=True).strip()

    def test_pinned_tree_ignores_live_and_staged_checkout_changes(self):
        (self.root / "pinned.js").write_text("original", encoding="utf-8")
        self.git("add", "pinned.js")
        tree = self.git("write-tree")
        (self.root / "pinned.js").unlink()
        (self.root / "live.js").write_text("new", encoding="utf-8")
        self.git("add", "pinned.js", "live.js")
        self.assertEqual(dc._tracked_entries(self.root, tree)[1], entries("pinned.js"))
        self.assertEqual(dc._tracked_entries(self.root)[1], entries("live.js"))

    def test_sync_writes_only_when_content_changes(self):
        readme = self.root / "README.md"
        readme.write_text("<!-- docs-check: tree exhaustive -->\n```\n├── old.h  Stale entry\n└── new.h  New source\n```\n", encoding="utf-8")
        (self.root / "new.h").write_text("", encoding="utf-8")
        self.git("add", "README.md", "new.h")
        with contextlib.redirect_stdout(io.StringIO()):
            ds.sync_repository(self.root, {}, {})
        with mock.patch.object(Path, "write_text", side_effect=AssertionError("rewrote unchanged document")):
            ds.sync_repository(self.root, {}, {})

    def test_roster_counts_sync_without_rewriting_surrounding_prose(self):
        header = self.root / dc._EFFECT_ROSTER_SOURCE
        header.parent.mkdir(parents=True)
        header.write_text("#define HS_EFFECT_LIST(X) \\\n    X(One) \\\n    X(Two)\n", encoding="utf-8")
        playlist = self.root / dc._PHANTASM_PLAYLIST_SOURCE
        playlist.parent.mkdir(parents=True)
        playlist.write_text("#define HS_PHANTASM_EFFECT_LIST(X) \\\n    X(One, 10)\n", encoding="utf-8")
        readme = self.root / "README.md"
        readme.write_text("The playlist contains 99 effects today.\n", encoding="utf-8")
        self.git("add", "README.md", str(dc._EFFECT_ROSTER_SOURCE), str(dc._PHANTASM_PLAYLIST_SOURCE))
        with contextlib.redirect_stdout(io.StringIO()):
            ds.sync_repository(self.root, {}, {})
        self.assertEqual(readme.read_text(encoding="utf-8"), "The playlist contains 1 effects today.\n")

    def test_a_correct_spelled_out_count_is_left_alone(self):
        header = self.root / dc._EFFECT_ROSTER_SOURCE
        header.parent.mkdir(parents=True)
        header.write_text("#define HS_EFFECT_LIST(X) \\\n    X(One) \\\n    X(Two)\n", encoding="utf-8")
        playlist = self.root / dc._PHANTASM_PLAYLIST_SOURCE
        playlist.parent.mkdir(parents=True)
        playlist.write_text("#define HS_PHANTASM_EFFECT_LIST(X) \\\n    X(One, 10)\n", encoding="utf-8")
        readme = self.root / "README.md"
        readme.write_text("The playlist contains one effects today.\n", encoding="utf-8")
        self.git("add", "README.md", str(dc._EFFECT_ROSTER_SOURCE), str(dc._PHANTASM_PLAYLIST_SOURCE))
        with contextlib.redirect_stdout(io.StringIO()):
            ds.sync_repository(self.root, {}, {})
        self.assertEqual(readme.read_text(encoding="utf-8"), "The playlist contains one effects today.\n")


if __name__ == "__main__":
    unittest.main()
