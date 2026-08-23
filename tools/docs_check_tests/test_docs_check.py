import contextlib
import io
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
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
        self.assertEqual(len(issues), 4)
        self.assertIn("#missing-anchor", issues[0].message)
        self.assertIn("this document", issues[0].message)
        self.assertIn("#ghost", issues[1].message)
        self.assertIn("docs/other.md", issues[1].message)
        self.assertIn("pins a line number", issues[2].message)
        self.assertIn("#snakecase-heading", issues[3].message)

    def test_backticked_repo_and_untracked_paths_are_linted(self):
        text = (FIXTURES / "backticked_paths.txt").read_text(encoding="utf-8")
        entries = {
            PurePosixPath("core"),
            PurePosixPath("core/engine"),
            PurePosixPath("core/engine/platform.h"),
            PurePosixPath("scripts"),
            PurePosixPath("tools"),
            PurePosixPath("tools/docs_check.py"),
        }
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, entries)
        self.assertEqual(len(issues), 2)
        self.assertEqual(issues[0].line, 2)
        self.assertIn("core/platform.h", issues[0].message)
        self.assertEqual(issues[1].line, 5)
        self.assertIn("scripts/run-tests.mjs", issues[1].message)

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

    def test_prose_child_list_is_gated_in_both_directions(self):
        text = ("<!-- docs-check: tree exhaustive -->\n"
                "```\n"
                "└── core/                       Engine\n"
                "    └── render/                 Per-stage headers (surface,\n"
                "        │                        warp), plus the base (ghost)\n"
                "```\n")
        entries = {
            PurePosixPath("core"),
            PurePosixPath("core/render"),
            PurePosixPath("core/render/surface.h"),
            PurePosixPath("core/render/warp.h"),
            PurePosixPath("core/render/color.h"),
        }
        issues = dc.check_text(PurePosixPath("README.md"), text, entries)
        self.assertEqual([issue.message for issue in issues], [
            "tree omits tracked path 'core/render/color.h'",
            "tree path 'core/render/ghost' does not exist",
        ])

    def test_prose_group_that_is_not_a_stem_list_names_nothing(self):
        text = ("<!-- docs-check: tree exhaustive -->\n"
                "```\n"
                "└── core/                       Engine\n"
                "    └── render/                 Per-stage headers (one per\n"
                "        │                        stage), a work in progress\n"
                "```\n")
        entries = {
            PurePosixPath("core"),
            PurePosixPath("core/render"),
            PurePosixPath("core/render/surface.h"),
        }
        self.assertEqual(
            dc.check_text(PurePosixPath("README.md"), text, entries), [])

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

    def test_line_pinned_repo_links_are_reported(self):
        text = ("[single](../hardware/pov_sync.h#L111)\n"
                "[range](../hardware/pov_sync.h#L405-L451)\n"
                "[blob](https://github.com/woundedlion/pov/blob/master/"
                "hardware/pov_sync.h#L20)\n"
                "[plain](../hardware/pov_sync.h)\n")
        entries = {PurePosixPath("hardware"),
                   PurePosixPath("hardware/pov_sync.h")}
        issues = dc.check_text(PurePosixPath("docs/spec.md"), text, entries)
        self.assertEqual([issue.line for issue in issues], [1, 2, 3])
        for issue in issues:
            self.assertIn("pins a line number", issue.message)

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
        allowances = ("tracked.txt", "uncited.txt", "used.txt")
        with mock.patch.object(dc, "_UNTRACKED_ALLOWED", allowances):
            stale = dict(item.split(" ", 1) for item in dc._stale_allowances(
                {PurePosixPath("tracked.txt")}, {"used.txt"}))
        self.assertEqual(stale["tracked.txt"], "(now tracked)")
        self.assertEqual(stale["uncited.txt"], "(uncited)")
        self.assertNotIn("used.txt", stale)

    def test_cited_untracked_allowance_is_recorded(self):
        used: set[str] = set()
        allowance = "scripts/run-tests.mjs"
        with mock.patch.object(dc, "_UNTRACKED_ALLOWED", (allowance,)):
            issue = dc._path_span_issue(PurePosixPath("README.md"), 1,
                                        allowance,
                                        {PurePosixPath("scripts")}, used)
        self.assertIsNone(issue)
        self.assertEqual(used, {allowance})

    @staticmethod
    def _checkout_fence_repository(root: Path) -> None:
        subprocess.run(["git", "init", "-q", str(root)], check=True)
        (root / "README.md").write_text(
            "<!-- docs-check: tree daydream -->\n"
            "```\n"
            "└── ghost.js                    Stale row\n"
            "```\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(root), "add", "README.md"],
                       check=True)

    def test_unacknowledged_checkout_tree_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._checkout_fence_repository(root)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                status = dc.main(["--root", str(root)])
        self.assertEqual(status, 1)
        self.assertIn("::error::", output.getvalue())
        self.assertIn("daydream", output.getvalue())

    def test_acknowledged_checkout_tree_passes_and_says_so(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._checkout_fence_repository(root)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                status = dc.main(["--root", str(root),
                                  "--skip-checkout", "daydream"])
        self.assertEqual(status, 0)
        self.assertIn("::warning::", output.getvalue())
        self.assertIn("NOT validated: daydream", output.getvalue())

    def test_missing_checkout_root_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._checkout_fence_repository(root)
            absent = root / "no-such-checkout"
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    dc.main(["--root", str(root),
                             "--checkout", f"daydream={absent}"])
        self.assertEqual(raised.exception.code, 2)

    def test_effect_roster_reads_the_macro_continuation(self):
        source = ("#define HS_EFFECT_COUNT_ADD(name) +1\n"
                  "#define HS_EFFECT_LIST(X)         \\\n"
                  "  X(Comets)                       \\\n"
                  "  X(Voronoi)\n"
                  "#define HS_PHANTASM_EFFECT_LIST(X) \\\n"
                  "  X(Comets)\n")
        self.assertEqual(dc.effect_roster(source), {"Comets", "Voronoi"})

    def test_effect_roster_tolerates_whitespace_inside_the_parens(self):
        source = ("#define HS_EFFECT_LIST(X) \\\n"
                  "  X( Comets )             \\\n"
                  "  X(Voronoi)\n")
        self.assertEqual(dc.effect_roster(source), {"Comets", "Voronoi"})

    def test_effect_roster_skips_commented_out_rows(self):
        source = ("#define HS_EFFECT_LIST(X) \\\n"
                  "  // X(Retired)           \\\n"
                  "  /* X(AlsoRetired) */    \\\n"
                  "  X(Voronoi)\n")
        self.assertEqual(dc.effect_roster(source), {"Voronoi"})

    @staticmethod
    def _diagram(count):
        return f"│  │  effects/  ({count} visual algorithms)     │  │\n"

    def test_matching_effects_row_is_clean(self):
        entries = {PurePosixPath("effects"),
                   PurePosixPath("effects/Comets.h"),
                   PurePosixPath("effects/Voronoi.h"),
                   PurePosixPath("effects/shared_palettes.h")}
        row = ("├── effects/  3 headers: one per effect (2) plus the shared\n"
               + self._diagram(2))
        self.assertEqual(
            dc.effects_row_issues(row, entries, {"Comets", "Voronoi"}), [])

    def test_matching_compact_effects_row_is_clean(self):
        entries = {PurePosixPath("effects"),
                   PurePosixPath("effects/Comets.h"),
                   PurePosixPath("effects/PromotedLooks.h")}
        row = ("├── effects/  2 headers covering 3 effects plus shared bases\n"
               + self._diagram(3))
        self.assertEqual(
            dc.effects_row_issues(
                row, entries, {"Comets", "SignalWeave", "KaleidoWave"}), [])

    def test_effects_row_counts_are_checked_against_tree_and_roster(self):
        entries = {PurePosixPath("effects"),
                   PurePosixPath("effects/Comets.h"),
                   PurePosixPath("effects/Voronoi.h"),
                   PurePosixPath("effects/shared_palettes.h")}
        row = ("├── effects/  9 headers: one per effect (7) plus the shared\n"
               + self._diagram(2))
        issues = dc.effects_row_issues(row, entries, {"Comets", "Voronoi"})
        self.assertEqual([issue.line for issue in issues], [1, 1])
        self.assertIn("claims 9 headers, the tracked tree has 3",
                      issues[0].message)
        self.assertIn("claims 7 effects, HS_EFFECT_LIST names 2",
                      issues[1].message)

    def test_diagram_effect_count_is_checked_against_the_roster(self):
        entries = {PurePosixPath("effects"),
                   PurePosixPath("effects/Comets.h")}
        row = ("├── effects/  1 headers covering 1 effects plus shared bases\n"
               + self._diagram(23))
        issues = dc.effects_row_issues(row, entries, {"Comets"})
        self.assertEqual([issue.line for issue in issues], [2])
        self.assertIn("architecture diagram claims 23 effects, "
                      "HS_EFFECT_LIST names 1", issues[0].message)

    def test_unreadable_roster_fails_the_effects_row(self):
        issues = dc.effects_row_issues(
            "├── effects/  1 headers: one per effect (1) plus the shared\n"
            + self._diagram(1),
            {PurePosixPath("effects"), PurePosixPath("effects/Comets.h")},
            None)
        self.assertEqual(len(issues), 2)
        self.assertIn("defines no HS_EFFECT_LIST", issues[0].message)
        self.assertIn("architecture diagram claims 1 effects",
                      issues[1].message)

    def test_deleted_effects_row_is_reported(self):
        issues = dc.effects_row_issues(
            "nothing to see\n" + self._diagram(1), set(), {"Comets"})
        self.assertEqual(len(issues), 1)
        self.assertIn("no effects/ summary row", issues[0].message)

    def test_deleted_diagram_row_is_reported(self):
        issues = dc.effects_row_issues(
            "├── effects/  1 headers covering 1 effects plus shared bases\n",
            {PurePosixPath("effects"), PurePosixPath("effects/Comets.h")},
            {"Comets"})
        self.assertEqual(len(issues), 1)
        self.assertIn("no effects/ architecture-diagram row",
                      issues[0].message)

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
