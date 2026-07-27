#!/usr/bin/env python3
"""Host tests for profile_one.sh configuration verification."""

import re
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
PROFILE_ONE = REPO / "tools" / "profile_one.sh"


def shell_function(name):
    source = PROFILE_ONE.read_text(encoding="utf-8")
    match = re.search(rf"(?ms)^{name}\(\) \{{\n.*?^\}}\n", source)
    if not match:
        raise AssertionError(f"missing shell function: {name}")
    return match.group(0)


def verify_log(text, expected="o3"):
    with tempfile.TemporaryDirectory() as directory:
        log = Path(directory) / "capture.log"
        log.write_text(text, encoding="utf-8")
        script = (
            "set -e\n"
            f"{shell_function('verify')}\n"
            'OUT=$1; EFFECT=Fx; TAG=$2; MARKER=""\n'
            "verify\n"
        )
        return subprocess.run(
            ["bash", "-c", script, "profile-test", str(log), expected],
            capture_output=True,
            text=True,
        )


def capture_log(config):
    return (
        f"profile harness: effect=Fx config={config} segments=4 rpm=480 "
        "f_cpu=600000000\n"
        "f 1 w=100 r=90\n"
        "=== profile Fx [288x144] frames 1-1 window=100 us ===\n"
    )


class MainTreeResolution(unittest.TestCase):
    """The build tree must be the main checkout from anywhere, worktrees too."""

    def _main_tree_from(self, script_dir):
        script = f"{shell_function('main_tree')}\nmain_tree\n"
        result = subprocess.run(
            ["bash", "-c", script, str(Path(script_dir) / "profile_one.sh")],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return Path(result.stdout.strip())

    def test_resolves_this_checkout(self):
        self.assertEqual(self._main_tree_from(REPO / "tools").resolve(),
                         Path(subprocess.run(
                             ["git", "-C", str(REPO), "rev-parse",
                              "--path-format=absolute", "--git-common-dir"],
                             capture_output=True, text=True, check=True
                         ).stdout.strip()).parent.resolve())

    def test_worktree_resolves_to_its_main_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "main"
            tree = Path(directory) / "wt"
            subprocess.run(["git", "init", "-q", "-b", "main", str(root)],
                           check=True)
            (root / "seed.txt").write_text("seed\n", encoding="utf-8")
            for args in (["add", "seed.txt"],
                         ["-c", "user.name=t", "-c", "user.email=t@t",
                          "commit", "-qm", "seed"],
                         ["worktree", "add", "-q", str(tree), "-b", "branch"]):
                subprocess.run(["git", "-C", str(root)] + args, check=True)
            self.assertEqual(self._main_tree_from(tree).resolve(),
                             root.resolve())


class ProfileConfigVerification(unittest.TestCase):
    def test_matching_config_passes(self):
        self.assertEqual(verify_log(capture_log("o3")).returncode, 0)

    def test_wrong_config_fails(self):
        result = verify_log(capture_log("ship"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CONFIG TAG MISMATCH", result.stdout)

    def test_missing_config_fails(self):
        result = verify_log(
            "f 1 w=100 r=90\n"
            "=== profile Fx [288x144] frames 1-1 window=100 us ===\n"
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("NO CONFIG TAG", result.stdout)

    def test_platformio_failure_survives_tail_pipeline(self):
        source = PROFILE_ONE.read_text(encoding="utf-8")
        set_line = re.search(r"(?m)^set .+$", source).group(0)
        script = (
            f"{set_line}\n"
            "pio() { return 23; }\n"
            f"{shell_function('flash')}\n"
            'HS_TEENSY_PORT=""; ENV=profile\n'
            "flash\n"
        )
        result = subprocess.run(
            ["bash", "-c", script], capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 23)


if __name__ == "__main__":
    unittest.main()
