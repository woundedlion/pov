"""Render comparison distinguishes build failures from pixel differences."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import git_run, isolated_env

REPO = Path(__file__).resolve().parents[2]


class RenderDiffTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.env = isolated_env()
        git_run(self.root, "init", "--quiet", env=self.env)
        git_run(self.root, "config", "core.autocrlf", "false", env=self.env)
        (self.root / "tests").mkdir()
        (self.root / "tests" / "test_render_hashes.h").write_text("fixture\n")
        git_run(self.root, "add", "--", "tests/test_render_hashes.h", env=self.env)
        git_run(self.root, "commit", "--quiet", "-m", "fixture", env=self.env)
        (self.root / "stubs").mkdir()
        stub = self.root / "stubs" / "cmake"
        stub.write_text("""#!/usr/bin/env bash
set -eu
if [ "$1" = --build ]; then
  stage=build
else
  stage=configure
fi
if [ "${FAIL_STAGE:-}" = "$stage" ]; then
  echo "fixture $stage diagnostic"
  exit 1
fi
if [ "$stage" = configure ]; then
  while [ "$1" != -B ]; do shift; done
  mkdir -p "$2"
  touch "$2/CMakeCache.txt"
  cksum < tests/test_render_hashes.h > "$2/source-id"
else
  mkdir -p "$2/tests"
  cat > "$2/tests/$4" <<'RUNNER'
#!/usr/bin/env bash
if [ "$1" = --list ]; then
  printf '  plot_parity\\n  render_hashes\\n'
else
  echo "render-hash fixture p0=$(cut -d' ' -f1 "$(dirname "$0")/../source-id")"
fi
RUNNER
  chmod +x "$2/tests/$4"
fi
""", encoding="utf-8", newline="\n")
        stub.chmod(0o755)
        self.env["HS_RENDER_DIFF_DIR"] = str(self.root / "base")

    def compare(self, *args):
        return subprocess.run(
            ["bash", "-c", 'export PATH="$PWD/stubs:$PATH"; exec bash "$@"',
             "fixture", str(REPO / "tools" / "render_diff.sh"), *args, "HEAD"],
            cwd=self.root, env=self.env, capture_output=True, text=True, timeout=30)

    def test_matching_renders_succeed(self):
        for args in ((), ("--report",)):
            with self.subTest(args=args):
                result = self.compare(*args)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("unchanged" if not args else "no pixel changes", result.stdout)

    def test_unowned_checkouts_keep_local_work(self):
        sibling = self.root / "other"
        git_run(self.root, "worktree", "add", "--quiet", "--detach", str(sibling),
                env=self.env)
        for checkout in (self.root, sibling):
            tracked = checkout / "tests/test_render_hashes.h"
            untracked = checkout / "untracked"
            tracked.write_text("local edit")
            untracked.write_text("local file")
            self.env["HS_RENDER_DIFF_DIR"] = str(checkout)
            for args in ((), ("--report",)):
                with self.subTest(checkout=checkout, args=args):
                    result = self.compare(*args)
                    self.assertEqual(result.returncode, 0 if args else 3,
                                     result.stdout + result.stderr)
                    self.assertIn("refusing to reuse", result.stderr)
                    self.assertEqual(tracked.read_text(), "local edit")
                    self.assertEqual(untracked.read_text(), "local file")

    def test_changed_renders_fail_and_report_the_changed_preset(self):
        (self.root / "tests" / "test_render_hashes.h").write_text("changed fixture\n")
        result = self.compare()
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("rendered pixels changed", result.stderr)
        report = self.compare("--report")
        self.assertEqual(report.returncode, 0, report.stdout + report.stderr)
        self.assertIn("1 of 1 renders differ", report.stdout)
        self.assertIn("  fixture p0", report.stdout)

    def test_configure_and_build_failures_print_diagnostics(self):
        for stage in ("configure", "build"):
            for args in ((), ("--report",)):
                with self.subTest(stage=stage, args=args):
                    self.env["FAIL_STAGE"] = stage
                    result = self.compare(*args)
                    self.assertEqual(result.returncode, 0 if args else 3,
                                     result.stdout + result.stderr)
                    self.assertIn(f"fixture {stage} diagnostic", result.stderr)


if __name__ == "__main__":
    unittest.main()
