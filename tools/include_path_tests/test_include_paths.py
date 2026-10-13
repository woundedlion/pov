"""Contract checks for canonical first-party include paths."""

import contextlib
import io
import subprocess
from unittest import mock
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import include_paths  # noqa: E402


class IncludePathsTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        self.files = ["core/math/vector.h", "effects/sample.h"]
        (self.root / "core/math").mkdir(parents=True)
        (self.root / "effects").mkdir()
        (self.root / "core/math/vector.h").write_text("#pragma once\n", encoding="utf-8")

    def check(self, source):
        (self.root / "effects/sample.h").write_text(source, encoding="utf-8")
        return include_paths.check_includes(self.root, self.files)

    def test_canonical_external_and_generated_includes_are_accepted(self):
        self.assertEqual(self.check('#include "core/math/vector.h"\n'
                                    '#include "Arduino.h"\n'
                                    '#include "generated_fixture.h"\n'), [])

    def test_forbidden_dependency_directions(self):
        rules = {
            "core/": ("effects/", "targets/", "workbench/", "hardware/", "tests/", "tools/"),
            "hardware/": ("effects/", "targets/", "workbench/", "tests/"),
            "effects/": ("targets/", "workbench/", "hardware/", "tests/"),
            "workbench/": ("effects/", "targets/"),
        }
        for prefix, forbidden in rules.items():
            for target_prefix in forbidden:
                with self.subTest(source=prefix, target=target_prefix):
                    source, target = prefix + "fixture.h", target_prefix + "fixture.h"
                    for name in (source, target):
                        (self.root / name).parent.mkdir(parents=True, exist_ok=True)
                        (self.root / name).write_text("#pragma once\n", encoding="utf-8")
                    self.assertEqual(include_paths.check_includes(self.root, [source, target]), [])
                    (self.root / source).write_text(f'#include "{target}"\n', encoding="utf-8")
                    self.assertEqual(include_paths.check_includes(self.root, [source, target]), [
                        f'{source}:1: forbidden include "{target}" from {prefix}'])

    def test_legacy_effect_and_vendor_dependency_exemptions(self):
        target = "effects/sample.h"
        self.check("#pragma once\n")
        for source in ("core/engine/effects_legacy.h", "core/vendor/upstream.h"):
            (self.root / source).parent.mkdir(parents=True, exist_ok=True)
            (self.root / source).write_text(f'#include "{target}"\n', encoding="utf-8")
            self.assertEqual(include_paths.check_includes(self.root, [source, target]), [])

    def test_project_config_obeys_core_dependency_direction(self):
        source = include_paths.CONFIG_HEADER
        target = "effects/sample.h"
        self.check("#pragma once\n")
        (self.root / source).parent.mkdir(parents=True, exist_ok=True)
        (self.root / source).write_text('#include "core/math/vector.h"\n', encoding="utf-8")
        files = [*self.files, source]
        self.assertEqual(include_paths.check_includes(self.root, files), [])
        (self.root / source).write_text(f'#include "{target}"\n', encoding="utf-8")
        self.assertEqual(include_paths.check_includes(self.root, files), [
            f'{source}:1: forbidden include "{target}" from core/'])

    def test_missing_tracked_source_is_diagnosed(self):
        diagnostics = include_paths.check_includes(self.root, self.files)
        self.assertEqual(len(diagnostics), 1)
        self.assertTrue(diagnostics[0].startswith("effects/sample.h:1: cannot read source"))

    def test_core_relative_alias_is_rejected(self):
        self.assertEqual(self.check('#include "math/vector.h"\n'), [
            'effects/sample.h:1: include "core/math/vector.h" instead of "math/vector.h"'])

    def test_parent_relative_alias_is_rejected(self):
        self.assertEqual(self.check('#include "../core/math/vector.h"\n'), [
            'effects/sample.h:1: include "core/math/vector.h" instead of "../core/math/vector.h"'])

    def test_comments_do_not_count_as_includes(self):
        self.assertEqual(self.check('/*\n#include "math/vector.h"\n*/\n'
                                    '// #include "math/vector.h"\n'), [])

    def test_vendor_body_is_exempt_but_project_config_is_checked(self):
        vendor = self.root / "core/vendor"
        vendor.mkdir()
        files = [*self.files, "core/vendor/upstream.h", include_paths.CONFIG_HEADER]
        (vendor / "upstream.h").write_text('#include "math/vector.h"\n', encoding="utf-8")
        (self.root / include_paths.CONFIG_HEADER).write_text(
            '#include "math/vector.h"\n', encoding="utf-8")
        self.check('#include "core/math/vector.h"\n')
        self.assertEqual(include_paths.check_includes(self.root, files), [
            f'{include_paths.CONFIG_HEADER}:1: include "core/math/vector.h" instead of "math/vector.h"'])


class TestMain(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        self.files = ["core/math/vector.h", "effects/sample.h"]
        for name in self.files:
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("#pragma once\n", encoding="utf-8")

    def run_main(self, files):
        with mock.patch.object(include_paths, "git_nul_list", return_value=files), \
                contextlib.redirect_stdout(io.StringIO()) as out, \
                contextlib.redirect_stderr(io.StringIO()) as err:
            result = include_paths.main(["--root", str(self.root)])
        return result, out.getvalue(), err.getvalue()

    def test_clean_and_bad_include_exit_codes(self):
        source = self.root / "effects/sample.h"
        source.write_text('#include "core/math/vector.h"\n', encoding="utf-8")
        self.assertEqual(self.run_main(self.files)[0], 0)
        source.write_text('#include "math/vector.h"\n', encoding="utf-8")
        result, out, _ = self.run_main(self.files)
        self.assertEqual(result, 1)
        self.assertIn('instead of "math/vector.h"', out)

    def test_empty_or_non_source_listing_is_a_tooling_error(self):
        for files in ([], ["README.md"]):
            with self.subTest(files=files):
                result, out, err = self.run_main(files)
                self.assertEqual(result, 2)
                self.assertNotIn("PASS", out)
                self.assertIn("no tracked C/C++ sources", err)

    def test_failed_or_undecodable_listing_is_a_tooling_error(self):
        for error in (subprocess.CalledProcessError(128, "git"),
                      UnicodeDecodeError("utf-8", b"\xff", 0, 1, "invalid"),
                      FileNotFoundError("git")):
            with self.subTest(error=type(error).__name__), \
                    mock.patch.object(include_paths, "git_nul_list", side_effect=error), \
                    contextlib.redirect_stderr(io.StringIO()) as err:
                self.assertEqual(include_paths.main(["--root", str(self.root)]), 2)
                self.assertIn("tooling error", err.getvalue())


if __name__ == "__main__":
    unittest.main()
