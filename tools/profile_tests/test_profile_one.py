#!/usr/bin/env python3
"""Host tests for profile_one.sh configuration verification."""

import hashlib
import os
import re
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
PROFILE_ONE = REPO / "tools" / "profile_one.sh"
EFFECTS = REPO / "effects"


class ProfileSweepTests(unittest.TestCase):
    def test_roster_tool_failures_cannot_report_success(self):
        sweep = (REPO / "tools" / "profile_sweep.sh").as_posix()
        for tool in ("comm", "sort", "sed", "tr"):
            with self.subTest(tool=tool):
                result = subprocess.run(
                    ["bash", "-c", f'{tool}() {{ return 3; }}; export -f {tool}; '
                     'bash "$1" check', "sweep-test", sweep],
                    capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("profile_sweep: cannot", result.stderr)
                self.assertNotIn("effect(s) match", result.stdout)


def shell_function(name):
    source = PROFILE_ONE.read_text(encoding="utf-8")
    match = re.search(rf"(?ms)^{name}\(\) \{{\n.*?^\}}\n", source)
    if not match:
        raise AssertionError(f"missing shell function: {name}")
    return match.group(0)


def verify_log(text, expected="o3", provenance=True, marker=""):
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
            'OUT=$1; EFFECT=Fx; TAG=$2; MARKER=$3\n'
            "verify\n"
        )
        return subprocess.run(
            ["bash", "-c", script, "profile-test", str(log), expected, marker],
            capture_output=True,
            text=True,
            encoding="utf-8",
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
        encoding="utf-8",
    )


class ProfileTreeLock(unittest.TestCase):
    def test_abandoned_empty_lock_is_recovered(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock = root / '.profile-lock'
            lock.mkdir()
            os.utime(lock, (1, 1))
            script = ('. "$1"\n'
                      + 'TREE=$2; SECONDS_ARG=1; acquire_tree_lock\n'
                      + 'rc=$?; [ "$rc" != 0 ] || _hs_break_lock "$TREE_LOCK" "$TREE_TOKEN"; exit "$rc"\n')
            result = subprocess.run(['bash', '-c', script, 'tree-lock-test',
                                     (REPO / 'tools/device_lock.sh').as_posix(), root.as_posix()],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(lock.exists())

    def test_cleanup_continues_after_a_replaced_tree_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / 'cache'
            cache.mkdir()
            (cache / 'large-build').write_text('artifact')
            script = ('set -e\n' + shell_function('cleanup')
                      + 'hs_device_release() { :; }; _hs_break_lock() { return 1; }\n'
                      + 'TREE_LOCK=unused; TREE_TOKEN=old; RETRY_CACHE=$1; cleanup\n')
            result = subprocess.run(['bash', '-c', script, 'cleanup-test', cache.as_posix()],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(cache.exists())

    def test_stale_owner_is_reaped_and_live_owner_is_retained(self):
        for stale in (True, False):
            with self.subTest(stale=stale), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                lock = root / '.profile-lock'
                lock.mkdir()
                deadline = 1 if stale else 9999999999
                (lock / 'info').write_text(f'token=peer\nstarted=1\ndeadline={deadline}\n')
                script = ('. "$1"\n'
                          + 'TREE=$2; SECONDS_ARG=1; acquire_tree_lock\n'
                          + 'rc=$?; [ "$rc" != 0 ] || _hs_break_lock "$TREE_LOCK" "$TREE_TOKEN"; exit "$rc"\n')
                result = subprocess.run(['bash', '-c', script, 'tree-lock-test',
                                         (REPO / 'tools/device_lock.sh').as_posix(), root.as_posix()],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0 if stale else 1, result.stderr)
                self.assertEqual(lock.exists(), not stale)


class ProfileTreeResolution(unittest.TestCase):
    """The build tree follows the checkout containing the invoked script."""

    def _profile_tree_from(self, script_dir):
        script = f"{shell_function('profile_tree')}\nprofile_tree\n"
        result = subprocess.run(
            ["bash", "-c", script, str(Path(script_dir) / "profile_one.sh")],
            capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(result.returncode, 0, result.stderr)
        return Path(result.stdout.strip())

    def test_resolves_this_checkout(self):
        self.assertEqual(self._profile_tree_from(REPO / "tools").resolve(),
                         REPO.resolve())

    def test_worktree_resolves_to_itself(self):
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
            self.assertEqual(self._profile_tree_from(tree).resolve(),
                             tree.resolve())


ARTIFACT_VARS = ("OUT", "PROVENANCE_OUT", "PROFILE_BUILD_LOG",
                 "PHANTASM_BUILD_LOG", "PROFILE_ENVDUMP", "PHANTASM_ENVDUMP",
                 "ATTEST_DIR", "ARTIFACT_BASE")


def derived_paths(profile_out=None):
    """Evaluate profile_one.sh's artifact-path block for one HS_PROFILE_OUT."""
    source = PROFILE_ONE.read_text(encoding="utf-8")
    block = re.search(r"(?ms)^OUT=.*?^ATTEST_DIR=.*?$", source)
    if not block:
        raise AssertionError("missing artifact-path block")
    script = (
        "LOWER=islamicstars; TAG=profile; ENV=profile\n"
        "DEEP_SUFFIX=; MODE_SUFFIX=; MSP_SUFFIX=\n"
        + block.group(0) + "\n"
        + "".join(f'echo "{name}=${name}"\n' for name in ARTIFACT_VARS)
    )
    env = {"PATH": os.environ["PATH"]}
    if profile_out is not None:
        env["HS_PROFILE_OUT"] = profile_out
    result = subprocess.run(["bash", "-c", script], capture_output=True,
                            text=True, encoding="utf-8", env=env)
    if result.returncode != 0:
        raise AssertionError(result.stderr)
    return dict(line.split("=", 1) for line in result.stdout.splitlines())


class ReplayFlags(unittest.TestCase):
    def test_spaced_compact_and_assigned_defines(self):
        source = PROFILE_ONE.read_text(encoding="utf-8")
        block = source[source.index("REPLAY_SUFFIX="):source.index('MSP_FLAGS=""')]
        for suffix in ("", "_AB"):
            for spacing in ("", " "):
                for value in ("", "=1"):
                    flag = f"-D{spacing}HS_MINDSPLATTER_REPLAY{suffix}{value}"
                    with self.subTest(flag=flag):
                        result = subprocess.run(
                            ["bash", "-c", 'EXTRA=$1\n' + block
                             + 'printf "%s" "$REPLAY_SUFFIX"', "test", flag],
                            capture_output=True, text=True, encoding="utf-8")
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertEqual(result.stdout, "_replay" + suffix.lower())


class ArtifactPaths(unittest.TestCase):
    """HS_PROFILE_OUT must move the whole artifact set, not just the log.

    profile_islamic_big.sh renames only the log; artifacts derived from the
    default naming would land on the standard IslamicStars run's.
    """

    def test_default_naming_is_unchanged(self):
        paths = derived_paths()
        self.assertEqual(paths["OUT"], "build/prof/islamicstars_profile.log")
        self.assertEqual(paths["ATTEST_DIR"],
                         "build/prof/attest/islamicstars_profile")

    def test_override_moves_every_artifact(self):
        override = "build/prof/islamicstars_big_ship.log"
        paths = derived_paths(override)
        default = derived_paths()
        self.assertEqual(paths["OUT"], override)
        for name in ARTIFACT_VARS:
            with self.subTest(var=name):
                self.assertIn("islamicstars_big_ship", paths[name])
                self.assertNotEqual(paths[name], default[name])


class ProfileConfigVerification(unittest.TestCase):
    def test_retry_uses_an_isolated_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shared = root / ".pio" / "build_cache"
            env_build = root / ".pio" / "build" / "profile"
            shared.mkdir(parents=True)
            env_build.mkdir(parents=True)
            (shared / "sentinel").write_text("shared\n", encoding="utf-8")
            (env_build / "stale").write_text("stale\n", encoding="utf-8")
            script = (
                f"{shell_function('prepare_retry')}\n"
                "ENV=profile\n"
                f"TMPDIR={root.as_posix()}\n"
                f"cd {root.as_posix()}\n"
                "prepare_retry\n"
                'test ! -e .pio/build/profile/stale\n'
                'test -e .pio/build_cache/sentinel\n'
                'test "$PLATFORMIO_BUILD_CACHE_DIR" != .pio/build_cache\n'
                'test -d "$PLATFORMIO_BUILD_CACHE_DIR"\n'
                'rm -rf "$PLATFORMIO_BUILD_CACHE_DIR"\n'
            )
            result = subprocess.run(["bash", "-c", script], capture_output=True,
                                    text=True, encoding="utf-8")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

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
                encoding="utf-8",
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PROFILE ARTIFACT HASH MISMATCH", result.stdout)

    def test_a_foreign_window_name_is_reported(self):
        # A peer flashing mid-capture splices its board's serial into ours;
        # the log then carries windows the effect under test never rendered.
        text = capture_log("o3") + (
            "=== profile Other [288x144] frames 2-2 window=100 us ===\n")
        result = verify_log(text)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("WINDOW NAME MISMATCH (1/2)", result.stdout)
        self.assertIn("— contention?", result.stdout)

    def test_a_missing_cycler_marker_is_reported(self):
        # A stale image runs the previous build; a cycler that emits no
        # advance marker is the tell that the upload did not take.
        result = verify_log(capture_log("o3"), marker="Preset:")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("NO 'Preset:' MARKER", result.stdout)

    def test_a_present_cycler_marker_passes(self):
        result = verify_log("Preset: 0/3\n" + capture_log("o3"),
                            marker="Preset:")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_a_first_frame_past_the_connect_window_is_reported(self):
        # A freshly flashed board starts at frame 1 and the capture attaches
        # within the connect window, so a high first frame means no reboot.
        text = capture_log("o3").replace("f 1 w=100 r=90", "f 900 w=100 r=90")
        result = verify_log(text)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("STALE IMAGE: first frame 900", result.stdout)

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
                encoding="utf-8",
            )
        self.assertNotEqual(result.returncode, 0)


class ToolchainAttestation(unittest.TestCase):
    def test_envdump_precedes_each_attested_build(self):
        body = shell_function("build_and_attest")
        self.assertLess(
            body.index('pio run -e phantasm -t envdump'),
            body.index('build_image phantasm'),
        )
        self.assertLess(
            body.index('pio run -e "$ENV" -t envdump'),
            body.index('build_image "$ENV"'),
        )

    def test_matching_compiler_passes(self):
        result = attest_toolchains("GCC 15.2.1", "GCC 15.2.1")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_compiler_mismatch_fails_before_flash(self):
        result = attest_toolchains("GCC 11.3.1", "GCC 15.2.1")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("compiler mismatch", result.stdout)


def multi_preset_effects():
    """Effects whose PRESET_IDS holds more than one id."""
    found = set()
    for header in EFFECTS.glob("*.h"):
        text = header.read_text(encoding="utf-8")
        for count in re.findall(
                r"std::array<std::string_view,\s*(\d+)>\s+PRESET_IDS\b", text):
            if int(count) > 1:
                found.add(header.stem)
    return found


def cyclers():
    source = PROFILE_ONE.read_text(encoding="utf-8")
    match = re.search(r'(?m)^CYCLERS="([^"]*)"$', source)
    if not match:
        raise AssertionError("missing CYCLERS list")
    return set(match.group(1).split())


class CyclerRoster(unittest.TestCase):
    """Every multi-preset effect emits an advance marker; the marker guard
    only runs for effects named in CYCLERS."""

    def test_every_multi_preset_effect_is_a_cycler(self):
        presets = multi_preset_effects()
        self.assertTrue(presets, "no multi-preset effect parsed from effects/")
        self.assertEqual(presets - cyclers(), set())


if __name__ == "__main__":
    unittest.main()
