#!/usr/bin/env python3
"""Host tests for profile_one.sh verification, sweep rosters, and tree locks."""

import hashlib
import os
import shlex
import re
import subprocess
import sys
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


def verify_log(text, expected="o3", provenance=True, marker="", *,
               preset="", expected_shape="", msp_marker="",
               artifact_present=True, digest_override=None):
    with tempfile.TemporaryDirectory() as directory:
        artifact = Path(directory) / "firmware.elf"
        artifact.write_bytes(b"profile firmware")
        digest = (hashlib.sha256(artifact.read_bytes()).hexdigest()
                  if digest_override is None else digest_override)
        if not artifact_present:
            artifact.unlink()
        if provenance:
            text += (
                "profile provenance: version=1\n"
                f"profile provenance: profile_elf_sha256={digest}\n"
                f"profile provenance: artifact_profile_elf={artifact.as_posix()}\n"
            )
        log = Path(directory) / "capture.log"
        log.write_text(text, encoding="utf-8")
        script = (
            "set -euo pipefail\n"
            f"{shell_function('verify')}\n"
            f"{shell_function('file_sha256')}\n"
            'OUT=$1; EFFECT=Fx; TAG=$2; MARKER=$3\n'
            'PROFILE_PRESET=$4; MSP_MARKER=$5; HS_PROFILE_EXPECT_SHAPE=$6\n'
            "verify\n"
        )
        return subprocess.run(
            ["bash", "-c", script, "profile-test", str(log), expected, marker,
             preset, msp_marker, expected_shape],
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


def attest_toolchains(profile_compiler, phantasm_compiler,
                      profile_abi="Tag_CPU_name: 7E-M", phantasm_abi="Tag_CPU_name: 7E-M",
                      profile_packages="toolchain 15.2.1", phantasm_packages="toolchain 15.2.1"):
    script = (
        f"{shell_function('assert_matching_toolchains')}\n"
        "elf_compiler() {\n"
        '  if [ "$1" = profile.elf ]; then\n'
        f"    echo '{profile_compiler}'\n"
        "  else\n"
        f"    echo '{phantasm_compiler}'\n"
        "  fi\n"
        "}\n"
        f"elf_abi() {{ if [ \"$1\" = profile.elf ]; then printf '%s\\n' {shlex.quote(profile_abi)}; "
        f"else printf '%s\\n' {shlex.quote(phantasm_abi)}; fi; }}\n"
        f"package_fingerprint() {{ if [ \"$1\" = profile.log ]; then printf '%s\\n' {shlex.quote(profile_packages)}; "
        f"else printf '%s\\n' {shlex.quote(phantasm_packages)}; fi; }}\n"
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
                      + 'TREE=$2; acquire_tree_lock\n'
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
                dead_pid = subprocess.check_output(
                    ['bash', '-c', 'echo $$'], text=True).strip()
                owner = dead_pid if stale else '$$'
                script = ('. "$1"\n'
                          + f"printf 'token=peer\\npid=%s\\nstarted=1\\ndeadline=9999999999\\n' \"{owner}\" > \"$2/.profile-lock/info\"\n"
                          + 'TREE=$2; acquire_tree_lock\n'
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
    """HS_PROFILE_OUT selects paths for the complete capture artifact set."""

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
                "set -e\n"
                f"{shell_function('prepare_retry')}\n"
                "ENV=profile\n"
                f"TMPDIR={root.as_posix()}\n"
                f"cd {root.as_posix()}\n"
                "prepare_retry\n"
                'test ! -e .pio/build/profile/stale\n'
                'test -e .pio/build_cache/sentinel\n'
                'test "$PLATFORMIO_BUILD_CACHE_DIR" != .pio/build_cache\n'
                'test -d "$PLATFORMIO_BUILD_CACHE_DIR"\n'
                f"{shlex.quote(sys.executable)} -c "
                "'import os,sys; sys.exit(os.path.samefile(sys.argv[1], sys.argv[2]))' "
                '"$PLATFORMIO_BUILD_CACHE_DIR" .pio/build_cache\n'
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
        result = verify_log(capture_log("o3"), digest_override="0" * 64)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PROFILE ARTIFACT HASH MISMATCH", result.stdout)

    def test_fixed_preset_requires_the_matching_capture_marker(self):
        result = verify_log("Profile preset: 3/6\n" + capture_log("o3"), preset="3")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for line in ("", "Profile preset: 4/6\n"):
            with self.subTest(line=line):
                result = verify_log(line + capture_log("o3"), preset="3")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("PROFILE PRESET MISMATCH (expected 3)", result.stdout)

    def test_expected_shape_requires_the_matching_spawn_marker(self):
        result = verify_log("Spawning Shape: Cube\n" + capture_log("o3"), expected_shape="Cube")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for line in ("", "Spawning Shape: Sphere\n"):
            with self.subTest(line=line):
                result = verify_log(line + capture_log("o3"), expected_shape="Cube")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("PROFILE SHAPE MISMATCH (expected Cube)", result.stdout)

    def test_mindsplatter_requires_the_selected_instrumentation(self):
        for marker in ("plot render counts particles:", "plot stall: stage=history_vertex"):
            with self.subTest(marker=marker):
                result = verify_log(marker + " 1\n" + capture_log("o3"), msp_marker=marker)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                result = verify_log(capture_log("o3"), msp_marker=marker)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(f"NO '{marker}' INSTRUMENTATION", result.stdout)

    def test_missing_artifact_is_rejected_with_valid_provenance(self):
        result = verify_log(capture_log("o3"))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        result = verify_log(capture_log("o3"), artifact_present=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("MISSING PROFILE ARTIFACT", result.stdout)

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

    def test_platformio_status_propagates_from_build(self):
        for status, expected in ((0, 0), (23, 1)):
            script = (
                f"pio() {{ echo pio-called; return {status}; }}\n"
                f"{shell_function('build_image')}\n"
                "build_image profile $1\n"
            )
            with self.subTest(status=status), tempfile.TemporaryDirectory() as directory:
                log = Path(directory) / "build.log"
                result = subprocess.run(
                    ["bash", "-c", script, "profile-test", str(log)],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                )
                self.assertEqual(result.returncode, expected, result.stderr)
                self.assertIn("pio-called", log.read_text(encoding="utf-8"))


def build_and_attest_calls():
    """Run build_and_attest with stubbed tools; returns (env, step) per call."""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        calls = root / "calls"
        script = (
            "set -e\n"
            f"{shell_function('build_and_attest')}\n"
            f'pio() {{ echo "pio $*" >> {shlex.quote(calls.as_posix())}; }}\n'
            f'build_image() {{ echo "build $1" >> {shlex.quote(calls.as_posix())}; }}\n'
            "assert_matching_toolchains() { :; }\n"
            "_hs_resolve_python() { _HS_LOCK_PYTHON=envdump_filter; }\n"
            "envdump_filter() { cat >/dev/null; }\n"
            "cp() { :; }\n"
            "git() { echo 0000000000000000000000000000000000000000; }\n"
            "file_sha256() { echo abc; }\n"
            "ENV=profile\n"
            f"cd {shlex.quote(root.as_posix())}\n"
            "ATTEST_DIR=attest\n"
            "for name in PHANTASM_ENVDUMP PROFILE_ENVDUMP PHANTASM_BUILD_LOG "
            "PROFILE_BUILD_LOG PHANTASM_ELF PHANTASM_MAP PROFILE_ELF PROFILE_MAP "
            "PROVENANCE_OUT; do printf -v \"$name\" 'attest/%s' \"$name\"; done\n"
            "ARTIFACT_BASE=artifacts/run\n"
            "build_and_attest\n"
        )
        result = subprocess.run(["bash", "-c", script], capture_output=True,
                                text=True, encoding="utf-8")
        if result.returncode != 0:
            raise AssertionError(result.stdout + result.stderr)
        lines = calls.read_text(encoding="utf-8").splitlines()
    steps = []
    for line in lines:
        words = line.split()
        if words[0] == "build":
            steps.append((words[1], "build"))
            continue
        target = words[words.index("-t") + 1] if "-t" in words else "build"
        steps.append((words[words.index("-e") + 1], target))
    return steps


class ToolchainAttestation(unittest.TestCase):
    def test_envdump_precedes_each_attested_build(self):
        steps = build_and_attest_calls()
        for env in ("phantasm", "profile"):
            with self.subTest(env=env):
                self.assertIn((env, "envdump"), steps)
                self.assertIn((env, "build"), steps)
                self.assertLess(steps.index((env, "envdump")),
                                steps.index((env, "build")))

    def test_matching_compiler_passes(self):
        result = attest_toolchains("GCC 15.2.1", "GCC 15.2.1")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_compiler_mismatch_fails_before_flash(self):
        result = attest_toolchains("GCC 11.3.1", "GCC 15.2.1")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("compiler mismatch", result.stdout)


    def test_abi_and_package_failures_refuse_attestation(self):
        for values, message in [
            ({"phantasm_abi": "Tag_CPU_name: other"}, "ARM ABI mismatch"),
            ({"profile_abi": ""}, "no ARM ABI"),
            ({"phantasm_packages": "different-package"}, "PlatformIO package mismatch"),
            ({"profile_packages": ""}, "no package fingerprint"),
        ]:
            with self.subTest(values=values):
                result = attest_toolchains("GCC 15.2.1", "GCC 15.2.1", **values)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stdout)


def multi_preset_effects():
    """Effects with more than one declared preset."""
    found = set()
    for header in EFFECTS.glob("*.h"):
        text = header.read_text(encoding="utf-8")
        for count in re.findall(
                r"std::array<PresetEntry<[^>]*>,\s*(\w+)>\s+PRESETS\b", text):
            if not count.isdecimal():
                value = re.search(r"\b" + re.escape(count) + r"\s*=\s*(\d+)\s*;", text)
                if value is None:
                    raise AssertionError(f"{header.name}: unresolved preset count {count}")
                count = value.group(1)
            if int(count) > 1:
                found.add(header.stem)
        for count in re.findall(
                r"std::array<std::string_view,\s*(\d+)>\s+PRESET_IDS\b", text):
            if int(count) > 1:
                found.add(header.stem)
        for entries in re.findall(
                r"PRESET_IDS\s*=\s*std::to_array<std::string_view>\s*\(\s*\{(.*?)\}\s*\)",
                text, re.DOTALL):
            if len(re.findall(r'"[^"\n]*"', entries)) > 1:
                found.add(header.stem)
        for body in re.findall(
                r"\bPRESETS\s*=\s*\[\]\s*\{(.*?)\}\s*\(\s*\)", text, re.DOTALL):
            registry = re.search(
                r"std::array<PresetEntry<.*?>,\s*std::size\(Solids::(\w+)\)\s*>",
                body, re.DOTALL)
            if registry is None:
                raise AssertionError(f"{header.name}: unresolved lambda preset count")
            declarations = (REPO / "core/mesh/solids.h").read_text(encoding="utf-8")
            entries = re.search(
                r"\b" + re.escape(registry.group(1)) + r"\s*\[\s*\]\s*=\s*\{(.*?)\};",
                declarations, re.DOTALL)
            if entries is None:
                raise AssertionError(f"{header.name}: unresolved preset registry {registry.group(1)}")
            if len(re.findall(r'\{\s*"', entries.group(1))) > 1:
                found.add(header.stem)
    return found


def cyclers():
    source = PROFILE_ONE.read_text(encoding="utf-8")
    match = re.search(r'(?m)^CYCLERS="([^"]*)"$', source)
    if not match:
        raise AssertionError("missing CYCLERS list")
    return set(match.group(1).split())


class CyclerRoster(unittest.TestCase):
    """Every multi-preset effect is in CYCLERS, so its marker is verified."""

    def test_every_multi_preset_effect_is_a_cycler(self):
        presets = multi_preset_effects()
        self.assertTrue(presets, "no multi-preset effect parsed from effects/")
        self.assertIn("HyperLattice", presets)
        self.assertTrue({"Comets", "DreamBalls", "MeshFeedback", "MindSplatter",
                         "ShapeShifter", "IslamicStars"}.issubset(presets))
        self.assertNotIn("Fishbowl", presets)
        self.assertEqual(presets - cyclers(), set())


if __name__ == "__main__":
    unittest.main()
