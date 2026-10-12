"""Contract checks for canonical first-party include paths."""

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
        for source in ("core/engine/effects_legacy.h", include_paths.CONFIG_HEADER):
            (self.root / source).parent.mkdir(parents=True, exist_ok=True)
            (self.root / source).write_text(f'#include "{target}"\n', encoding="utf-8")
            self.assertEqual(include_paths.check_includes(self.root, [source, target]), [])

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


if __name__ == "__main__":
    unittest.main()
