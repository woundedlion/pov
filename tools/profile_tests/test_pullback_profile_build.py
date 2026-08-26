"""Tests for pullback profile build provenance."""

import runpy
import subprocess
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "pullback_profile_build.py"
SHA = "0123456789ab"


class FakeEnvironment(dict):
    def Append(self, **values):
        self.appended = values


class PullbackProfileBuild(unittest.TestCase):
    def _run_script(self, status):
        environment = FakeEnvironment(
            PIOENV="profile", PROJECT_DIR="project")
        results = [
            subprocess.CompletedProcess([], 0, f"{SHA}\n", ""),
            subprocess.CompletedProcess([], 0, status, ""),
        ]
        with mock.patch("subprocess.run", side_effect=results) as run:
            runpy.run_path(
                str(SCRIPT),
                init_globals={"env": environment, "Import": lambda _: None},
            )
        return environment.appended["CPPDEFINES"], run.call_args_list

    def test_clean_tree_uses_bare_sha(self):
        defines, calls = self._run_script("")
        self.assertEqual(
            defines, [("HS_PULLBACK_SHORT_SHA", f'\\"{SHA}\\"')])
        self.assertEqual(calls[1].args[0], [
            "git", "-C", "project", "status", "--porcelain=v1",
            "--untracked-files=no",
        ])

    def test_dirty_tree_marks_sha_and_preserves_define_quoting(self):
        defines, _ = self._run_script(" M core/engine/memory.h\n")
        self.assertEqual(
            defines, [("HS_PULLBACK_SHORT_SHA", f'\\"{SHA}-dirty\\"')])


if __name__ == "__main__":
    unittest.main()
