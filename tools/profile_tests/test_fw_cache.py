#!/usr/bin/env python3
"""Host tests for the firmware image cache and profile_one.sh's use of it."""

import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import git_run, isolated_env  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FW_CACHE = REPO / "tools" / "fw_cache.sh"
PROFILE_ONE = REPO / "tools" / "profile_one.sh"

FAKE_PIO = r"""#!/bin/bash
if [ "$1" = --version ]; then echo "PlatformIO Core, version ${FAKE_PIO_VERSION:-6.1.19}"; exit 0; fi
echo "$*" >> "$FAKE_PIO_CALLS"
env=$3
case " $* " in *" -t envdump "*)
  printf "{\n  'ENV': {'PATH': 'host'},\n  'PIOENV': '%s',\n}\n" "$env"
  exit 0;;
esac
[ -z "${FAKE_PIO_FAIL:-}" ] || exit 1
out=${PLATFORMIO_BUILD_DIR:-.pio/build}/$env
mkdir -p "$out"
printf 'elf %s %s %s\n' "$env" "$PLATFORMIO_BUILD_FLAGS" "$(cat source.txt)" > "$out/firmware.elf"
printf 'hex %s\n' "$env" > "$out/firmware.hex"
printf 'map %s\n' "$env" > "$out/firmware.map"
echo "PLATFORM: Teensy (5.2.0)"
echo " - framework-arduinoteensy @ 1.162.0"
echo " - toolchain-gccarmnoneeabi-teensy @ 1.150201.0"
echo "========== [SUCCESS] $env"
"""

FAKE_READELF = r"""#!/bin/bash
case "$1" in
  -p) echo "  [     0]  GCC: (GNU) 15.2.1";;
  -A) echo '  Tag_CPU_name: "7E-M"';;
esac
"""


class FirmwareCacheFixture(unittest.TestCase):
    """A committed checkout, a stub PlatformIO and a private cache root."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        self.repo = root / "repo"
        self.cache = root / "cache"
        self.core = root / "core"
        self.calls = root / "pio_calls"
        bin_dir = root / "bin"
        bin_dir.mkdir()
        for name, body in (("pio", FAKE_PIO), ("readelf", FAKE_READELF)):
            path = bin_dir / name
            path.write_text(body, encoding="utf-8", newline="\n")
            path.chmod(0o755)
        self.readelf = bin_dir / "readelf"
        self.manifest(self.core / "platforms" / "teensy" / "platform.json", "5.2.0")
        self.manifest(self.core / "packages" / "toolchain-gccarmnoneeabi-teensy"
                      / "package.json", "1.150201.0")
        self.env = isolated_env()
        self.env = {key: value for key, value in self.env.items()
                    if not key.startswith(("PLATFORMIO_", "HS_FW_", "HS_DEVICE_",
                                           "HS_TEENSY_", "HS_PROFILE_"))}
        self.env.update(PATH=str(bin_dir) + os.pathsep + self.env["PATH"],
                        PLATFORMIO_CORE_DIR=self.core.as_posix(),
                        HS_FW_CACHE_DIR=self.cache.as_posix(),
                        HS_ARM_READELF=self.readelf.as_posix(),
                        FAKE_PIO_CALLS=self.calls.as_posix())

        self.repo.mkdir()
        self.git("init", "--quiet", "-b", "main")
        (self.repo / ".gitignore").write_text(
            "/.pio/\n/build/\n/.profile-lock/\n/.profile-lock.guard\n", encoding="utf-8")
        identity = self.repo / "core" / "platform" / "effect_identity.h"
        identity.parent.mkdir(parents=True)
        identity.write_text('  X(Fx, "fx", "Fx", "fx", "fx")\n', encoding="utf-8")
        (self.repo / "tools").mkdir()
        shutil.copyfile(REPO / "tools" / "profile_envdump.py",
                        self.repo / "tools" / "profile_envdump.py")
        self.commit_source("one")

    @staticmethod
    def manifest(path, version):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('{\n  "name": "x",\n  "version": "%s",\n  "deps": {\n'
                        '    "version": "~0.0.1"\n  }\n}\n' % version, encoding="utf-8")

    def git(self, *args):
        return git_run(self.repo, *args, env=self.env)

    def commit_source(self, text):
        (self.repo / "source.txt").write_text(text + "\n", encoding="utf-8")
        self.git("add", "-A")
        self.git("commit", "--quiet", "--allow-empty", "-m", text)

    def pio_builds(self):
        if not self.calls.exists():
            return []
        return [line for line in self.calls.read_text().splitlines()
                if "-t envdump" not in line]

    def entries(self):
        return sorted(p.name for p in self.cache.iterdir()) if self.cache.exists() else []

    def fw(self, script, **extra):
        """Run bash with fw_cache.sh sourced in the fixture checkout."""
        return subprocess.run(
            ["bash", "-c", 'set -u; . "$1"; ' + script, "fw-test", FW_CACHE.as_posix()],
            cwd=self.repo, env=dict(self.env, **extra),
            capture_output=True, text=True, encoding="utf-8")

    def key(self, env="profile", **extra):
        result = self.fw(f'hs_fw_cache_key {env} && echo "$HS_FW_KEY"', **extra)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.strip()


class CacheKey(FirmwareCacheFixture):
    def test_the_key_covers_tree_env_flags_and_toolchain(self):
        base = self.key(PLATFORMIO_BUILD_FLAGS="-D A -D B")
        self.assertRegex(base, r"^[0-9a-f]{32}$")
        self.assertEqual(self.key(PLATFORMIO_BUILD_FLAGS="  -D A\t-D  B "), base)
        self.assertNotEqual(self.key(PLATFORMIO_BUILD_FLAGS="-D A"), base)
        self.assertNotEqual(self.key("phantasm", PLATFORMIO_BUILD_FLAGS="-D A -D B"), base)
        self.assertNotEqual(self.key(PLATFORMIO_BUILD_FLAGS="-D A -D B",
                                     PLATFORMIO_BUILD_UNFLAGS="-Os"), base)
        self.assertNotEqual(self.key(PLATFORMIO_BUILD_FLAGS="-D A -D B",
                                     FAKE_PIO_VERSION="6.1.20"), base)
        self.manifest(self.core / "packages" / "toolchain-gccarmnoneeabi-teensy"
                      / "package.json", "1.150202.0")
        self.assertNotEqual(self.key(PLATFORMIO_BUILD_FLAGS="-D A -D B"), base)

    def test_a_new_tree_changes_the_key_and_an_empty_commit_does_not(self):
        base = self.key()
        self.git("commit", "--quiet", "--allow-empty", "-m", "empty")
        self.assertEqual(self.key(), base)
        self.commit_source("two")
        self.assertNotEqual(self.key(), base)

    def test_uncommitted_or_untracked_changes_refuse_a_key(self):
        for change in ("modified", "untracked"):
            with self.subTest(change=change):
                path = self.repo / ("source.txt" if change == "modified" else "new.txt")
                original = path.read_bytes() if path.exists() else None
                path.write_text("dirty\n", encoding="utf-8")
                result = self.fw('hs_fw_cache_key profile')
                self.assertEqual(result.returncode, 3, result.stderr)
                self.assertIn("uncommitted changes", result.stderr)
                if original is None:
                    path.unlink()
                else:
                    path.write_bytes(original)


class CacheStore(FirmwareCacheFixture):
    def build_and_store(self, **extra):
        return self.fw(
            'hs_fw_cache_lookup profile; rc=$?; [ "$rc" = 1 ] || exit 9\n'
            'mkdir -p build; pio run -e profile >build/out.log || exit 8\n'
            'hs_fw_cache_store profile .pio/build/profile build/out.log && echo "$HS_FW_ENTRY"',
            **extra)

    def test_a_miss_is_stored_and_then_hit(self):
        stored = self.build_and_store()
        self.assertEqual(stored.returncode, 0, stored.stdout + stored.stderr)
        entry = Path(stored.stdout.strip())
        meta = (entry / "meta").read_text(encoding="utf-8")
        head = self.git("rev-parse", "HEAD").stdout.strip()
        elf = entry / "firmware.elf"
        self.assertIn(f"key {entry.name}\n", meta)
        self.assertIn(f"commit {head}\n", meta)
        self.assertIn(f"elf_sha256 {hashlib.sha256(elf.read_bytes()).hexdigest()}\n", meta)
        self.assertTrue((entry / "firmware.hex").is_file())
        hit = self.fw('hs_fw_cache_lookup profile && echo "$HS_FW_ENTRY"')
        self.assertEqual(hit.returncode, 0, hit.stderr)
        self.assertEqual(Path(hit.stdout.strip()), entry)

    def test_a_corrupted_entry_is_a_miss(self):
        stored = self.build_and_store()
        self.assertEqual(stored.returncode, 0, stored.stdout + stored.stderr)
        self.assertEqual(self.fw("hs_fw_cache_lookup profile").returncode, 0)
        (Path(stored.stdout.strip()) / "firmware.elf").write_text("tampered\n", encoding="utf-8")
        self.assertEqual(self.fw("hs_fw_cache_lookup profile").returncode, 1)

    def test_a_tree_change_during_the_build_is_not_stored(self):
        result = self.fw(
            'hs_fw_cache_lookup profile\n'
            'mkdir -p build; pio run -e profile >build/out.log\n'
            'echo three > source.txt; git commit -qam three\n'
            'hs_fw_cache_store profile .pio/build/profile build/out.log')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("changed during the profile build", result.stderr)
        self.assertEqual(self.entries(), [])

    def test_a_disabled_cache_neither_reuses_nor_stores(self):
        self.assertEqual(self.build_and_store().returncode, 0)
        self.assertEqual(self.fw("hs_fw_cache_lookup profile", HS_FW_CACHE="0").returncode, 1)
        self.commit_source("two")
        result = self.build_and_store(HS_FW_CACHE="0")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(self.entries()), 1)

    def test_prune_keeps_the_newest_entries(self):
        for text in ("two", "three", "four"):
            self.commit_source(text)
            self.assertEqual(self.build_and_store(HS_FW_CACHE_KEEP="2").returncode, 0)
        self.assertEqual(len(self.entries()), 2)
        self.assertEqual(self.fw("hs_fw_cache_lookup profile").returncode, 0)

    def test_build_reuses_a_stored_image_without_platformio(self):
        first = self.fw('hs_fw_cache_build phantasm && echo "IMAGE=$HS_FW_IMAGE"')
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertIn("fw_cache: phantasm stored", first.stdout)
        self.assertEqual(len(self.pio_builds()), 1)
        second = self.fw('hs_fw_cache_build phantasm && echo "IMAGE=$HS_FW_IMAGE"')
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("fw_cache: phantasm reused", second.stdout)
        self.assertIn("[SUCCESS] phantasm", second.stdout)
        self.assertEqual(len(self.pio_builds()), 1)
        self.assertEqual(first.stdout.split("IMAGE=")[1], second.stdout.split("IMAGE=")[1])

    def test_build_of_an_uncommitted_checkout_flashes_the_build_dir(self):
        (self.repo / "source.txt").write_text("dirty\n", encoding="utf-8")
        result = self.fw('hs_fw_cache_build phantasm && echo "IMAGE=$HS_FW_IMAGE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("IMAGE=.pio/build/phantasm", result.stdout)
        self.assertEqual(self.entries(), [])

    def test_each_flag_set_builds_in_its_own_directory(self):
        script = 'hs_fw_cache_build profile >/dev/null && echo "$HS_FW_IMAGE"'
        images = {}
        for flags in ("", "-D A", "-D B", " -D  A "):
            result = self.fw(script, PLATFORMIO_BUILD_FLAGS=flags, HS_FW_CACHE="0",
                             HS_PIO_JOBS="4")
            self.assertEqual(result.returncode, 0, result.stderr)
            images[flags] = result.stdout.strip()
        self.assertEqual(images[""], ".pio/build/profile")
        self.assertRegex(images["-D A"], r"^\.pio/flags/[0-9a-f]{12}/profile$")
        self.assertNotEqual(images["-D A"], images["-D B"])
        self.assertEqual(images[" -D  A "], images["-D A"])
        for image in images.values():
            self.assertTrue((self.repo / image / "firmware.elf").is_file(), image)
        self.assertTrue(all(" -j 4" in call for call in self.pio_builds()))
        self.assertEqual(self.fw(script, PLATFORMIO_BUILD_FLAGS="-D C", HS_FW_CACHE="0",
                                 HS_FW_BUILD_DIRS_KEEP="1").returncode, 0)
        self.assertEqual(len(list((self.repo / ".pio" / "flags").iterdir())), 1)
        explicit = self.fw(script, PLATFORMIO_BUILD_DIR="out", HS_FW_CACHE="0")
        self.assertEqual(explicit.stdout.strip(), "out/profile")

    def test_a_failed_build_is_not_stored(self):
        result = self.fw("hs_fw_cache_build phantasm", FAKE_PIO_FAIL="1")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.entries(), [])

    def test_put_stores_only_a_build_of_the_current_tree(self):
        shutil.copyfile(FW_CACHE, self.repo / "tools" / "fw_cache.sh")
        self.commit_source("script")
        tree = self.git("rev-parse", "HEAD^{tree}").stdout.strip()
        self.fw("mkdir -p build; pio run -e phantasm > build/out.log")

        def put():
            return subprocess.run(
                ["bash", "tools/fw_cache.sh", "put", "phantasm", tree, "build/out.log"],
                cwd=self.repo, env=self.env, capture_output=True, text=True)

        stored = put()
        self.assertEqual(stored.returncode, 0, stored.stderr)
        self.assertEqual(Path(stored.stdout.strip()).parent, self.cache)
        self.assertEqual(self.fw("hs_fw_cache_lookup phantasm").returncode, 0)
        self.commit_source("two")
        refused = put()
        self.assertEqual(refused.returncode, 1, refused.stderr)
        self.assertIn(f"no longer at tree {tree}", refused.stderr)
        self.assertEqual(len(self.entries()), 1)


class ProfileOneCache(FirmwareCacheFixture):
    def profile(self, *options, window="16"):
        return subprocess.run(
            ["bash", PROFILE_ONE.as_posix(), *options, "fx", "profile", "20", window],
            cwd=self.repo, env=dict(self.env, HS_PROFILE_TREE=self.repo.as_posix()),
            capture_output=True, text=True, encoding="utf-8")

    def provenance(self):
        text = (self.repo / "build" / "prof" / "fx_ship.provenance").read_text(encoding="utf-8")
        return dict(line.removeprefix("profile provenance: ").split("=", 1)
                    for line in text.splitlines())

    def test_a_second_build_reuses_both_images_with_the_same_provenance(self):
        first = self.profile("--build-only")
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        built = self.provenance()
        self.assertEqual(len(self.pio_builds()), 2)
        self.assertRegex(built["profile_image"], r"^built:[0-9a-f]{32}$")
        self.assertRegex(built["phantasm_image"], r"^built:[0-9a-f]{32}$")

        second = self.profile("--build-only")
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        reused = self.provenance()
        self.assertEqual(len(self.pio_builds()), 2)
        self.assertEqual(reused["profile_image"], "cache:" + built["profile_image"][6:])
        self.assertEqual(reused["phantasm_image"], "cache:" + built["phantasm_image"][6:])
        for field in ("version", "source_sha", "compiler", "profile_elf_sha256",
                      "phantasm_elf_sha256", "profile_envdump_sha256",
                      "phantasm_envdump_sha256", "artifact_dir"):
            with self.subTest(field=field):
                self.assertEqual(reused[field], built[field])
        self.assertEqual(reused["source_sha"], self.git("rev-parse", "HEAD").stdout.strip())
        artifact = Path(self.repo / reused["artifact_profile_elf"])
        self.assertEqual(hashlib.sha256(artifact.read_bytes()).hexdigest(),
                         reused["profile_elf_sha256"])
        self.assertIn("-D HS_PROFILE_TARGET=Fx -D HS_PROFILE_WINDOW=16",
                      artifact.read_text(encoding="utf-8"))

    def test_dry_run_reports_misses_then_hits_and_builds_nothing(self):
        before = self.profile("--dry-run")
        self.assertEqual(before.returncode, 0, before.stderr)
        self.assertIn("dry-run: phantasm miss", before.stdout)
        self.assertIn("dry-run: profile miss", before.stdout)
        self.assertEqual(self.pio_builds(), [])
        self.assertFalse((self.repo / ".profile-lock").exists())
        self.assertEqual(self.profile("--build-only").returncode, 0)
        after = self.profile("--dry-run")
        self.assertIn("dry-run: phantasm hit", after.stdout)
        self.assertIn("dry-run: profile hit", after.stdout)
        self.assertIn("dry-run: profile miss", self.profile("--dry-run", window="32").stdout)

    def test_an_uncommitted_checkout_builds_without_caching(self):
        (self.repo / "source.txt").write_text("dirty\n", encoding="utf-8")
        self.assertIn("dry-run: profile uncached (uncommitted changes)",
                      self.profile("--dry-run").stdout)
        result = self.profile("--build-only")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.provenance()["profile_image"], "built")
        self.assertEqual(self.entries(), [])
        self.assertEqual(self.profile("--build-only").returncode, 0)
        self.assertEqual(len(self.pio_builds()), 4)

    def test_elf_accepts_only_its_own_checkout_env_and_flags(self):
        self.assertEqual(self.profile("--build-only").returncode, 0)
        elf = self.cache / self.provenance()["profile_image"].split(":", 1)[1] / "firmware.elf"
        accepted = self.profile("--dry-run", "--elf", elf.as_posix())
        self.assertEqual(accepted.returncode, 0, accepted.stdout + accepted.stderr)
        self.assertIn("--elf accepted", accepted.stdout)

        used = self.profile("--build-only", "--elf", elf.as_posix())
        self.assertEqual(used.returncode, 0, used.stdout + used.stderr)
        self.assertEqual(self.provenance()["profile_image"],
                         "prebuilt:" + elf.parent.name)
        self.assertEqual(len(self.pio_builds()), 2)

        refused = self.profile("--dry-run", "--elf", elf.as_posix(), window="32")
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("--elf refused", refused.stdout)
        refused = self.profile("--build-only", "--elf", elf.as_posix(), window="32")
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("PROFILE PROVENANCE FAIL", refused.stdout)

        build_elf = next((self.repo / ".pio" / "flags").glob("*/profile/firmware.elf"))
        self.assertNotEqual(self.profile("--dry-run", "--elf", build_elf.as_posix()).returncode, 0)

        self.commit_source("two")
        self.assertNotEqual(self.profile("--dry-run", "--elf", elf.as_posix()).returncode, 0)


if __name__ == "__main__":
    unittest.main()
