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
        "CHANGED=\"$(printf '%s\\n' \"$@\")\"\n"
        'hs_paths_native "$CHANGED" && echo native\n'
        'hs_paths_teensy "$CHANGED" && echo teensy\n'
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
    ("tools/teensy_gate.py", {"teensy"}),
    # Pins the committed LUTs / trail that CTests and the hook itself read.
    ("tools/teensy_size_trail.py", {"native", "teensy"}),
    ("scripts/generate_luts.py", {"native"}),
    ("tools/gen_gamut_lut.py", {"native"}),
    ("CMakeLists.txt", {"native"}),
    ("tests/CMakeLists.txt", {"native"}),
    ("CMakePresets.json", {"native"}),
    ("cmake/toolchain-native-clang.cmake", {"native"}),
    ("tests/check_includes.cmake", {"native"}),
    ("tests/ubsan-ignorelist.txt", {"native"}),
)

#: Paths that must select neither gate, so a doc-only commit stays fast.
NEITHER = (
    "README.md",
    "docs/specs/shaderball_spec.md",
    "docs/screenshots/Voronoi.png",
    "package.json",
    "LICENSE",
)


class Classification(unittest.TestCase):
    def test_each_path_class_selects_its_gates(self):
        for path, expected in CASES:
            with self.subTest(path=path):
                self.assertEqual(classify(path), expected)

    def test_documentation_and_assets_select_nothing(self):
        for path in NEITHER:
            with self.subTest(path=path):
                self.assertEqual(classify(path), set())

    def test_a_mixed_commit_selects_the_union(self):
        self.assertEqual(
            classify("README.md", "platformio.ini", "CMakeLists.txt"),
            {"native", "teensy"})

    def test_no_staged_paths_select_nothing(self):
        self.assertEqual(classify(""), set())

    def test_the_extension_anchors_hold(self):
        # Every extension alternative is anchored at end of line, so a path that
        # merely contains one is not a source file.
        for path in ("docs/notes.cmake.md", "docs/inl.txt", "tools/h.json",
                     "docs/CMakeLists.txt.md", "docs/platformio.ini.md"):
            with self.subTest(path=path):
                self.assertEqual(classify(path), set())

    def test_a_generator_the_suite_does_not_pin_is_not_native(self):
        # Only the generators a CTest compares committed output against are
        # listed; adding one here without adding it to the hook would pass
        # vacuously, so assert the negative for a generator that is not.
        self.assertEqual(classify("tools/docs_check.py"), set())


if __name__ == "__main__":
    unittest.main()
