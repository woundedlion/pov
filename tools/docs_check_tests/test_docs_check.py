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

    def test_missing_anchors_are_reported(self):
        text = (FIXTURES / "anchors.txt").read_text(encoding="utf-8")
        entries = {PurePosixPath("docs"), PurePosixPath("docs/other.md")}
        anchors = {PurePosixPath("docs/other.md"): {"real-heading"}}
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, entries,
                               anchors)
        self.assertEqual(len(issues), 3)
        self.assertIn("#missing-anchor", issues[0].message)
        self.assertIn("this document", issues[0].message)
        self.assertIn("#ghost", issues[1].message)
        self.assertIn("docs/other.md", issues[1].message)
        self.assertIn("#snakecase-heading", issues[2].message)

    def test_backticked_repo_paths_are_linted(self):
        text = (FIXTURES / "backticked_paths.txt").read_text(encoding="utf-8")
        entries = {
            PurePosixPath("core"),
            PurePosixPath("core/engine"),
            PurePosixPath("core/engine/platform.h"),
            PurePosixPath("tools"),
            PurePosixPath("tools/docs_check.py"),
        }
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, entries)
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].line, 2)
        self.assertIn("core/platform.h", issues[0].message)

    def test_backticked_firmware_linker_and_board_paths_are_linted(self):
        text = (FIXTURES / "native_format_paths.txt").read_text(encoding="utf-8")
        entries = {
            PurePosixPath("targets"),
            PurePosixPath("targets/Phantasm"),
            PurePosixPath("targets/Phantasm/Phantasm.ino"),
            PurePosixPath("tools"),
            PurePosixPath("tools/phantasm.ld"),
            PurePosixPath("hardware"),
            PurePosixPath("hardware/phantasm"),
            PurePosixPath("hardware/phantasm/phantasm.kicad_pcb"),
            PurePosixPath("docs"),
        }
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, entries)
        self.assertEqual([issue.line for issue in issues],
                         [2, 4, 6, 7, 8, 9, 10, 11, 12, 13])
        self.assertIn("targets/Phantasm/Ghost.ino", issues[0].message)
        self.assertIn("tools/ghost.ld", issues[1].message)

    def test_directed_repository_tree_rows_are_linted(self):
        text = (FIXTURES / "tree_fence.txt").read_text(encoding="utf-8")
        entries = {
            PurePosixPath("core"),
            PurePosixPath("core/engine"),
            PurePosixPath("core/engine/platform.h"),
            PurePosixPath("core/engine/memory.h"),
            PurePosixPath("core/engine/memory.cpp"),
            PurePosixPath("tools"),
            PurePosixPath("tools/docs_check_tests"),
            PurePosixPath("justfile"),
        }
        issues = dc.check_text(PurePosixPath("README.md"), text, entries)
        self.assertEqual([issue.line for issue in issues], [7, 8, 27])
        self.assertIn("core/engine/ghost.h", issues[0].message)
        self.assertIn("core/ghosts", issues[1].message)
        self.assertIn("unknown docs-check directive 'trees'", issues[2].message)

    def test_exhaustive_tree_reports_tracked_paths_without_a_row(self):
        text = ("<!-- docs-check: tree exhaustive -->\n"
                "```\n"
                "├── core/                       Engine\n"
                "│   └── engine/                 Machinery\n"
                "│       └── platform.h          Abstraction layer\n"
                "├── docs/                       Design specs\n"
                "└── justfile                    Task runner\n"
                "```\n")
        entries = {
            PurePosixPath("core"),
            PurePosixPath("core/engine"),
            PurePosixPath("core/engine/platform.h"),
            PurePosixPath("core/engine/util.h"),
            PurePosixPath("core/math"),
            PurePosixPath("docs"),
            PurePosixPath("docs/spec.md"),
            PurePosixPath("justfile"),
            PurePosixPath("ruff.toml"),
        }
        issues = dc.check_text(PurePosixPath("README.md"), text, entries)
        # docs/ names no child, so its subtree is elided rather than enumerated.
        self.assertEqual([issue.message for issue in issues], [
            "tree omits tracked path 'core/engine/util.h'",
            "tree omits tracked path 'core/math'",
            "tree omits tracked path 'ruff.toml'",
        ])
        self.assertEqual({issue.line for issue in issues}, {2})

    def test_unmapped_paths_are_exempt_from_the_exhaustive_diff(self):
        text = ("<!-- docs-check: tree exhaustive -->\n"
                "```\n"
                "└── justfile                    Task runner\n"
                "```\n")
        entries = {
            PurePosixPath(".gitignore"),
            PurePosixPath("README.md"),
            PurePosixPath("justfile"),
            PurePosixPath("tests"),
            PurePosixPath("tests/run_tests.cpp"),
        }
        issues = dc.check_text(PurePosixPath("README.md"), text, entries)
        self.assertEqual([issue.message for issue in issues],
                         ["tree omits tracked path 'tests'"])

    def test_checkout_tree_is_validated_against_the_supplied_root(self):
        text = ("<!-- docs-check: tree daydream -->\n"
                "```\n"
                "├── daydream.js                 App entry\n"
                "├── ghost.js                    Stale row\n"
                "└── three.js/                   Optional vendored checkout\n"
                "```\n")
        checkouts = {"daydream": {PurePosixPath("daydream.js")}}
        issues = dc.check_text(PurePosixPath("README.md"), text, set(),
                               checkouts=checkouts)
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].line, 4)
        self.assertIn("ghost.js", issues[0].message)

    def test_exhaustive_checkout_tree_diffs_against_the_checkout(self):
        text = ("<!-- docs-check: tree daydream exhaustive -->\n"
                "```\n"
                "├── daydream.js                 App entry\n"
                "├── tools/                      Standalone geometry tools\n"
                "│   └── color.js                sRGB math\n"
                "└── three.js/                   Optional vendored checkout\n"
                "```\n")
        checkouts = {"daydream": {
            PurePosixPath("daydream.js"),
            PurePosixPath("gui.js"),
            PurePosixPath("tools"),
            PurePosixPath("tools/color.js"),
            PurePosixPath("tools/slider.js"),
        }}
        issues = dc.check_text(PurePosixPath("README.md"), text, set(),
                               checkouts=checkouts)
        # three.js/ is gitignored in the checkout, so no row backs it there.
        self.assertEqual([issue.message for issue in issues], [
            "tree omits tracked path 'gui.js'",
            "tree omits tracked path 'tools/slider.js'",
        ])
        self.assertEqual({issue.line for issue in issues}, {2})

    def test_checkout_tree_without_a_root_is_recorded_as_skipped(self):
        text = ("<!-- docs-check: tree daydream -->\n"
                "```\n"
                "└── ghost.js                    Stale row\n"
                "```\n")
        skipped: set[str] = set()
        issues = dc.check_text(PurePosixPath("README.md"), text, set(),
                               skipped=skipped)
        self.assertEqual(issues, [])
        self.assertEqual(skipped, {"daydream"})

    def test_multi_word_tree_tag_is_parsed_or_reported(self):
        self.assertEqual(dc._tree_directive("tree"),
                         dc.TreeDirective("", False))
        self.assertEqual(dc._tree_directive("tree exhaustive"),
                         dc.TreeDirective("", True))
        self.assertEqual(dc._tree_directive("tree daydream"),
                         dc.TreeDirective("daydream", False))
        self.assertEqual(dc._tree_directive("tree daydream exhaustive"),
                         dc.TreeDirective("daydream", True))
        self.assertIsNone(dc._tree_directive("tree daydream extra"))
        self.assertIsNone(dc._tree_directive("trees"))

    def test_tree_row_indented_past_its_parent_is_reported(self):
        text = ("<!-- docs-check: tree -->\n"
                "```\n"
                "├── core/\n"
                "│   │   ├── platform.h\n"
                "```\n")
        issues = dc.check_text(PurePosixPath("README.md"), text,
                               {PurePosixPath("core")})
        self.assertEqual(len(issues), 1)
        self.assertIn("indented past its parent", issues[0].message)

    def test_absolute_self_repo_blob_links_are_validated(self):
        text = ("[ok](https://github.com/woundedlion/pov/blob/master/docs/real.md)\n"
                "[bad](https://github.com/woundedlion/pov/blob/master/docs/ghost.md)\n"
                "[other](https://github.com/woundedlion/daydream/blob/master/docs/ghost.md)\n"
                "[home](https://github.com/woundedlion/pov)\n")
        entries = {PurePosixPath("docs"), PurePosixPath("docs/real.md")}
        issues = dc.check_text(PurePosixPath("README.md"), text, entries)
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].line, 2)
        self.assertIn("docs/ghost.md", issues[0].message)

    def test_links_above_the_repository_root_are_reported(self):
        text = ("[sibling](../../daydream/daydream.js)\n"
                "[parent](../README.md)\n")
        entries = {PurePosixPath("README.md")}
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, entries)
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].line, 1)
        self.assertIn("escapes the repository root", issues[0].message)

    def test_only_git_tracked_markdown_is_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / "README.md").write_text("[ok](target.txt)\n", encoding="utf-8")
            (root / "target.txt").write_text("ok\n", encoding="utf-8")
            (root / "untracked.md").write_text("[bad](missing.txt)\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "README.md", "target.txt"],
                           check=True)
            markdown, issues, _ = dc.check_repository(root)
        self.assertEqual(markdown, [PurePosixPath("README.md")])
        self.assertEqual(issues, [])

    def test_internal_review_ledger_is_excluded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / "docs").mkdir()
            (root / "tools").mkdir()
            (root / "README.md").write_text("`tools/readme_missing.py`\n",
                                            encoding="utf-8")
            (root / "docs" / "CODE_REVIEW.md").write_text(
                "`tools/ledger_missing.py`\n", encoding="utf-8")
            (root / "tools" / "real.py").write_text("\n", encoding="utf-8")
            subprocess.run([
                "git", "-C", str(root), "add", "README.md",
                "docs/CODE_REVIEW.md", "tools/real.py",
            ], check=True)
            markdown, issues, _ = dc.check_repository(root)
        self.assertEqual(markdown, [PurePosixPath("README.md")])
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].path, "README.md")
        self.assertIn("tools/readme_missing.py", issues[0].message)

    def test_stale_untracked_allowances_are_named(self):
        entries = {PurePosixPath("tools"), PurePosixPath("tools/solids.html")}
        stale = dict(item.split(" ", 1) for item in
                     dc._stale_allowances(entries, {"scripts/run-tests.mjs"}))
        self.assertEqual(stale["tools/solids.html"], "(now tracked)")
        self.assertEqual(stale["tools/solid_codegen.js"], "(uncited)")
        self.assertNotIn("scripts/run-tests.mjs", stale)

    def test_cited_untracked_allowance_is_recorded(self):
        used: set[str] = set()
        issue = dc._path_span_issue(PurePosixPath("README.md"), 1,
                                    "scripts/run-tests.mjs",
                                    {PurePosixPath("scripts")}, used)
        self.assertIsNone(issue)
        self.assertEqual(used, {"scripts/run-tests.mjs"})

    def test_repository_without_markdown_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / "target.txt").write_text("ok\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "target.txt"],
                           check=True)
            self.assertEqual(dc.main(["--root", str(root)]), 2)


if __name__ == "__main__":
    unittest.main()
