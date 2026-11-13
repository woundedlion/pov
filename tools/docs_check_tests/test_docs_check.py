import contextlib
import io
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path, PurePosixPath

TOOLS = Path(__file__).resolve().parent.parent
FIXTURES = Path(__file__).resolve().parent / "fixtures"
sys.path.insert(0, str(TOOLS))

import docs_check as dc  # noqa: E402
from git_test_env import git_run  # noqa: E402


class TestDocumentationChecker(unittest.TestCase):
    def test_retired_behavior_wording_is_rejected_until_removed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            git_run(root, 'init', '-q')
            (root / "README.md").write_text(
                "# Contract\n<!-- docs-check: tree exhaustive -->\n```\n"
                "└── engine.h  Engine\n```\n"
                "<!-- docs-check: tree daydream exhaustive -->\n```\n```\n",
                encoding="utf-8")
            source = root / "engine.h"
            source.write_text("// Trips the wire-busy trap.\n", encoding="utf-8")
            (root / "ignored.md").write_text("wire-busy trap\n", encoding="utf-8")
            git_run(root, 'add', 'README.md', 'engine.h')
            arguments = ["--root", str(root), "--skip-checkout", "daydream",
                         "--retired-term", "WIRE-BUSY TRAP"]
            with contextlib.redirect_stdout(io.StringIO()) as output:
                self.assertEqual(dc.main(arguments), 1)
            self.assertIn("engine.h:1: retired wording", output.getvalue())
            source.write_text("// Returns false while busy.\n", encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(dc.main(arguments), 0)

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

    def test_missing_image_and_reference_targets_are_reported(self):
        text = (FIXTURES / "missing_targets.txt").read_text(encoding="utf-8")
        issues = dc.check_text(PurePosixPath("docs/readme.md"), text, set())
        self.assertEqual(len(issues), 3)
        self.assertTrue(any("docs/missing.svg" in issue.message for issue in issues))
        self.assertTrue(any("docs/missing.md" in issue.message for issue in issues))
        self.assertTrue(any("docs/nested-missing.svg" in issue.message for issue in issues))

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
        self.assertEqual(len(issues), 3)
        self.assertEqual(issues[0].line, 1)
        self.assertIn("pins a line number", issues[0].message)
        self.assertEqual(issues[1].line, 2)
        self.assertIn("core/platform.h", issues[1].message)
        self.assertEqual(issues[2].line, 5)
        self.assertIn("scripts/run-tests.mjs", issues[2].message)

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

    def test_stale_untracked_allowances_are_named(self):
        allowances = ("tracked.txt", "uncited.txt", "used.txt")
        used = {dc._cited(dc._UNTRACKED_LIST, "used.txt")}
        with mock.patch.object(dc, "UNTRACKED_ALLOWED", allowances), \
                mock.patch.object(dc, "TREE_UNMAPPED", ()):
            stale = dict(item.split(" ", 1) for item in dc._stale_allowances(
                {PurePosixPath("tracked.txt")}, used))
        self.assertEqual(stale["tracked.txt"], "(now tracked)")
        self.assertEqual(stale["uncited.txt"], "(uncited)")
        self.assertNotIn("used.txt", stale)

    def test_detached_required_tree_directive_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _git_repository(root, {"README.md": "<!-- docs-check: tree exhaustive -->\n"
                "Tree follows.\n"
                "```\n"
                "└── ghost/\n"
                "```\n"})
            required = frozenset({(
                PurePosixPath("README.md"), dc.TreeDirective("", True),
            )})
            with mock.patch.object(dc, "_REQUIRED_TREES", required):
                _, issues, _ = dc.check_repository(root)
        self.assertEqual(len(issues), 1)
        self.assertIn("required docs-check directive 'tree exhaustive'",
                      issues[0].message)

    @staticmethod
    def _checkout_fence_repository(root: Path) -> None:
        _git_repository(root, {"README.md": "<!-- docs-check: tree exhaustive -->\n"
            "```\n"
            "```\n"
            "<!-- docs-check: tree daydream exhaustive -->\n"
            "```\n"
            "└── ghost.js                    Stale row\n"
            "```\n"})

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

    @staticmethod
    def _diagram(count):
        return f"│  │  effects/  ({count} visual algorithms)     │  │\n"

    def test_effects_row_counts_are_checked_against_tree_and_roster(self):
        entries = {PurePosixPath("effects"),
                   PurePosixPath("effects/Comets.h"),
                   PurePosixPath("effects/SphericalHarmonics.h"),
                   PurePosixPath("effects/shared_palettes.h")}
        row = ("├── effects/  9 headers covering 7 effects plus the shared\n"
               + self._diagram(2))
        issues = dc.effects_row_issues(row, entries, {"Comets", "SphericalHarmonics"})
        self.assertEqual([issue.line for issue in issues], [1, 1])
        self.assertIn("claims 9 headers, the tracked tree has 3",
                      issues[0].message)
        self.assertIn("claims 7 effects, HS_EFFECT_LIST names 2",
                      issues[1].message)

    def test_unreferenced_predefined_finds_dead_macro_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "live.h").write_text("#if HS_LIVE\n#endif\n// HS_DEAD\n",
                                          encoding="utf-8")
            dead = dc.unreferenced_predefined(
                root, [PurePosixPath("live.h")],
                {"HS_LIVE", "HS_DEAD", "DOXYGEN"})
        self.assertEqual(dead, {"HS_DEAD"})

    _PLAYLIST_HEADER = ("#define HS_EFFECT_LIST(X) \\\n"
                        "  X(Comets)               \\\n"
                        "  X(SphericalHarmonics)\n"
                        "#define HS_SHADER_PRODUCT_GROUP(X) \\\n"
                        "  X(Comets)                   \\\n"
                        "  X(SphericalHarmonics)\n")
    _PLAYLIST_DEF = ("// show(Commented, 60)\n"
                     "show(Comets, 120)\n"
                     "show(SphericalHarmonics)\n"
                     "exclude(Dynamo) // 96x20 only\n")

    @staticmethod
    def _roster_prose(roster, playlist, sync):
        return (f"The compile-time roster and tests carry {roster} "
                f"firmware-capable effects. The playlist contains {playlist} "
                f"effects across the {sync}-entry roster.\n")

    @staticmethod
    def _effects_prose(products, descriptors):
        return (f"These {products} effects form the product-only group. "
                f"Firmware exposes the {descriptors} promoted fixed "
                "descriptors.")

    _SPEC = PurePosixPath("docs/specs/phantasm_frame_sync_spec.md")
    _EFFECTS = PurePosixPath("docs/effects.md")

    def _roster_issues(self, prose):
        header = self._PLAYLIST_HEADER
        return dc.roster_claim_issues(
            {PurePosixPath("README.md"): prose, self._SPEC: prose,
             self._EFFECTS: self._effects_prose(2, 2)},
            dc.effect_roster(header), dc.phantasm_roster(self._PLAYLIST_DEF),
            dc.shader_product_group(header))

    def test_roster_prose_counts_are_checked_against_the_macros(self):
        issues = self._roster_issues(self._roster_prose(3, 4, 5))
        messages = [issue.message for issue in issues]
        self.assertEqual(len(messages), 4)
        self.assertIn("stated as 3, HS_EFFECT_LIST names 2", messages[0])
        self.assertIn("stated as 4, phantasm_playlist.def names 2",
                      messages[1])
        self.assertIn("stated as 5, phantasm_playlist.def names 2",
                      messages[2])
        self.assertIn("stated as 5, phantasm_playlist.def names 2",
                      messages[3])
        self.assertEqual(issues[3].path, self._SPEC.as_posix())

    _COMPOSED_ROSTER = (
        "### Composed-effect roster\n\n"
        "| Effect ID | Concrete effect | Presets | Legacy source |\n"
        "|---|---|---:|---|\n")

    @staticmethod
    def _composed_row(effect_id, effect, presets):
        return f"| `{effect_id}` | `{effect}` | {presets} | — |\n"

    def _composed_issues(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "effects").mkdir()
            (root / "effects" / "Comets.h").write_text(
                "static constexpr std::array<std::string_view, 2> "
                'PRESET_IDS{\n    "comets", "comets-2"};\n', encoding="utf-8")
            (root / "effects" / "SphericalHarmonics.h").write_text(
                "// static constexpr std::array<std::string_view, 9> "
                "PRESET_IDS{};\n"
                "static constexpr std::array<std::string_view, 1> "
                'PRESET_IDS{"sphericalharmonics"};\n', encoding="utf-8")
            entries = {PurePosixPath("effects", name)
                       for name in ("Comets.h", "SphericalHarmonics.h")}
            return [(issue.line, issue.message)
                    for issue in dc.composed_roster_issues(
                        root, self._COMPOSED_ROSTER + rows, entries,
                        {"Comets", "SphericalHarmonics"})]

    def test_composed_roster_preset_counts_are_checked_against_preset_ids(
            self):
        rows = (self._composed_row("comets", "Comets", 3)
                + self._composed_row("sphericalharmonics", "SphericalHarmonics", 1))
        self.assertEqual(self._composed_issues(rows), [
            (5, "Comets presets stated as 3, PRESET_IDS holds 2")])


def _git_repository(root, files):
    """A git repository at root tracking files, a {path: text} map."""
    root.mkdir(parents=True, exist_ok=True)
    git_run(root, 'init', '-q')
    for name, text in files.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    git_run(root, 'add', *files)


class TestCommentSymbols(unittest.TestCase):
    def test_renamed_symbol_fails_the_repository_check(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _git_repository(root, {
                "README.md": "# Engine\n",
                "core/ring.h": "// Draws with `draw_ring`.\nvoid draw_ring();\n",
                "core/vendor/lib.h": "// `vendor_comment`\nint lib_name;\n",
                "core/use.cpp": "// Calls `lib_name`; not `uncited_name`.\n"
                                "// `comment_only_name`\n",
                "core/other.cpp": "// `comment_only_name` is no code\n",
            })
            arguments = ["--root", str(root), "--skip-checkout", "daydream"]
            with mock.patch.object(dc, "_REQUIRED_TREES", frozenset()), \
                    contextlib.redirect_stdout(io.StringIO()) as output, \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(dc.main(arguments), 1)
            self.assertEqual(
                [line for line in output.getvalue().splitlines()
                 if "comment names" in line],
                ["core/other.cpp:1: comment names `comment_only_name`, which "
                 "no code spells",
                 "core/use.cpp:1: comment names `uncited_name`, which no code "
                 "spells",
                 "core/use.cpp:2: comment names `comment_only_name`, which no "
                 "code spells"])
            for name in ("core/use.cpp", "core/other.cpp"):
                (root / name).write_text("// Calls `lib_name`.\n",
                                         encoding="utf-8")
            with mock.patch.object(dc, "_REQUIRED_TREES", frozenset()):
                self.assertEqual(dc.check_repository(root)[1], [])
            (root / "core/ring.h").write_text(
                "// Draws with `draw_ring`.\nvoid draw_rings();\n",
                encoding="utf-8")
            with mock.patch.object(dc, "_REQUIRED_TREES", frozenset()):
                issues = dc.check_repository(root)[1]
            self.assertEqual([str(issue) for issue in issues], [
                "core/ring.h:1: comment names `draw_ring`, which no code "
                "spells"])

    def test_excluding_every_source_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _git_repository(root, {"core/vendor/lib.h": "// `gone_name`\n"})
            _, entries = dc.tracked_entries(root)
            issues = dc.repository_symbol_issues(root, entries)
        self.assertEqual(len(issues), 1)
        self.assertIn("excluded from the comment symbol scan", issues[0].message)


if __name__ == "__main__":
    unittest.main()
