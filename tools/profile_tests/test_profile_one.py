#!/usr/bin/env python3
"""Host tests for profile_one.sh configuration verification."""

import hashlib
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


def verify_log(text, expected="o3", provenance=True):
    with tempfile.TemporaryDirectory() as directory:
        artifact = Path(directory) / "firmware.elf"
        artifact.write_bytes(b"profile firmware")
        digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        if provenance:
            text += (
                "profile provenance: version=1\n"
                f"profile provenance: profile_elf_sha256={digest}\n"
                f"profile provenance: artifact_profile_elf={artifact.as_posix()}\n"
            )
        log = Path(directory) / "capture.log"
        log.write_text(text, encoding="utf-8")
        script = (
            "set -e\n"
            f"{shell_function('verify')}\n"
            f"{shell_function('file_sha256')}\n"
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


def attest_toolchains(profile_compiler, phantasm_compiler):
    script = (
        f"{shell_function('assert_matching_toolchains')}\n"
        "elf_compiler() {\n"
        '  if [ "$1" = profile.elf ]; then\n'
        f"    echo '{profile_compiler}'\n"
        "  else\n"
        f"    echo '{phantasm_compiler}'\n"
        "  fi\n"
        "}\n"
        "elf_abi() { echo 'Tag_CPU_name: 7E-M'; }\n"
        "package_fingerprint() { echo 'toolchain 15.2.1'; }\n"
        "PROFILE_ELF=profile.elf\n"
        "PHANTASM_ELF=phantasm.elf\n"
        "PROFILE_BUILD_LOG=profile.log\n"
        "PHANTASM_BUILD_LOG=phantasm.log\n"
        "assert_matching_toolchains\n"
    )
    return subprocess.run(
        ["bash", "-c", script],
        capture_output=True,
        text=True,
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

    def test_missing_provenance_fails(self):
        result = verify_log(capture_log("o3"), provenance=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("NO/INVALID PROFILE PROVENANCE", result.stdout)

    def test_artifact_hash_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            artifact = Path(directory) / "firmware.elf"
            artifact.write_bytes(b"profile firmware")
            text = capture_log("o3") + (
                "profile provenance: version=1\n"
                f"profile provenance: profile_elf_sha256={'0' * 64}\n"
                f"profile provenance: artifact_profile_elf={artifact.as_posix()}\n"
            )
            log = Path(directory) / "capture.log"
            log.write_text(text, encoding="utf-8")
            script = (
                "set -e\n"
                f"{shell_function('verify')}\n"
                f"{shell_function('file_sha256')}\n"
                'OUT=$1; EFFECT=Fx; TAG=o3; MARKER=""\n'
                "verify\n"
            )
            result = subprocess.run(
                ["bash", "-c", script, "profile-test", str(log)],
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PROFILE ARTIFACT HASH MISMATCH", result.stdout)

    def test_platformio_failure_propagates_from_build(self):
        script = (
            "pio() { return 23; }\n"
            f"{shell_function('build_image')}\n"
            "build_image profile $1\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "build.log"
            result = subprocess.run(
                ["bash", "-c", script, "profile-test", str(log)],
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)


class ToolchainAttestation(unittest.TestCase):
    def test_matching_compiler_passes(self):
        result = attest_toolchains("GCC 15.2.1", "GCC 15.2.1")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_compiler_mismatch_fails_before_flash(self):
        result = attest_toolchains("GCC 11.3.1", "GCC 15.2.1")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("compiler mismatch", result.stdout)


if __name__ == "__main__":
    unittest.main()
