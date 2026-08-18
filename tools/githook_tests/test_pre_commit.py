#!/usr/bin/env python3
"""Host tests for .githooks/pre-commit path classification.

The hook decides what to run entirely from the staged path list, through two
regexes. A class that stops matching runs no gate and reports nothing, on every
commit that touches only that class -- so the classifiers are extracted as
shell functions and driven here against a table of staged paths, the same way
device_lock.sh and profile_one.sh are tested.

Run:  python -m unittest discover -s tools/githook_tests
"""

import re
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
HOOK = REPO / ".githooks" / "pre-commit"


def shell_function(name: str) -> str:
    source = HOOK.read_text(encoding="utf-8")
    match = re.search(rf"(?ms)^{name}\(\) \{{\n.*?^\}}\n", source)
    if not match:
        raise AssertionError(f"missing shell function: {name}")
    return match.group(0)


def classify(*paths: str) -> set[str]:
    """The gates the hook would run for this staged path list.

    Each path is its own argv entry and the newline-separated list is rebuilt
    inside the shell: a newline embedded in one Windows command-line argument
    does not survive the round trip.
    """
    script = (
        f"{shell_function('hs_paths_native')}\n"
        f"{shell_function('hs_paths_teensy')}\n"
        f"{shell_function('hs_paths_lint')}\n"
        "CHANGED=\"$(printf '%s\\n' \"$@\")\"\n"
        'hs_paths_native "$CHANGED" && echo native\n'
        'hs_paths_teensy "$CHANGED" && echo teensy\n'
        'hs_paths_lint "$CHANGED" && echo lint\n'
        "exit 0\n"
    )
    done = subprocess.run(
        ["bash", "-c", script, "hook-test", *paths],
        capture_output=True, text=True, check=True)
    return set(done.stdout.split())


#: (staged path, gates it must select). One row per class the regexes name.
CASES = (
    ("core/engine/memory.h", {"native", "teensy"}),
    ("core/engine/memory.cpp", {"native", "teensy"}),
    ("core/color/palettes.hpp", {"native", "teensy"}),
    ("core/mesh/mesh.inl", {"native", "teensy"}),
    ("core/vendor/FastNoiseLite.cc", {"native", "teensy"}),
    ("tests/test_mesh.h", {"native", "teensy"}),
    ("targets/Holosphere/Holosphere.ino", {"teensy"}),
    ("platformio.ini", {"teensy"}),
    ("tools/teensy_budgets.json", {"teensy"}),
    ("tools/phantasm.ld", {"teensy"}),
    ("tools/teensy_gate.py", {"teensy", "lint"}),
    # Generators the tree pins committed output against; every one is Python.
    ("tools/teensy_size_trail.py", {"native", "teensy", "lint"}),
    ("scripts/generate_luts.py", {"native", "lint"}),
    ("scripts/generate_reaction_graph.py", {"native", "lint"}),
    ("tools/gen_gamut_lut.py", {"native", "lint"}),
    ("tools/relax_bakes.py", {"native", "lint"}),
    ("CMakeLists.txt", {"native"}),
    ("tests/CMakeLists.txt", {"native"}),
    ("CMakePresets.json", {"native"}),
    ("cmake/toolchain-native-clang.cmake", {"native"}),
    ("tests/check_includes.cmake", {"native"}),
    ("tests/ubsan-ignorelist.txt", {"native"}),
    ("README.md", {"lint"}),
    ("docs/specs/pullback_pipeline_spec.md", {"lint"}),
    ("scripts/wasm_smoke.mjs", {"lint"}),
    ("tools/docs_check.py", {"lint"}),
    ("ruff.toml", {"lint"}),
    ("eslint.config.mjs", {"lint"}),
    ("package.json", {"lint"}),
    ("package-lock.json", {"lint"}),
)

#: Paths that must select no gate at all, so such a commit stays instant.
NEITHER = (
    "docs/screenshots/Voronoi.png",
    "LICENSE",
    "justfile",
    ".github/workflows/ci.yml",
    ".githooks/pre-commit",
)


class Classification(unittest.TestCase):
    def test_each_path_class_selects_its_gates(self):
        for path, expected in CASES:
            with self.subTest(path=path):
                self.assertEqual(classify(path), expected)

    def test_assets_and_ci_only_paths_select_nothing(self):
        for path in NEITHER:
            with self.subTest(path=path):
                self.assertEqual(classify(path), set())

    def test_a_mixed_commit_selects_the_union(self):
        self.assertEqual(
            classify("README.md", "platformio.ini", "CMakeLists.txt"),
            {"native", "teensy", "lint"})

    def test_no_staged_paths_select_nothing(self):
        self.assertEqual(classify(""), set())

    def test_the_extension_anchors_hold(self):
        # Every extension alternative is anchored at end of line, so a path that
        # merely contains one is not a source file.
        for path, expected in (("docs/notes.cmake.md", {"lint"}),
                               ("docs/inl.txt", set()),
                               ("tools/h.json", set()),
                               ("docs/CMakeLists.txt.md", {"lint"}),
                               ("docs/platformio.ini.md", {"lint"})):
            with self.subTest(path=path):
                self.assertEqual(classify(path), expected)


if __name__ == "__main__":
    unittest.main()
