import contextlib
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import license_check as lc  # noqa: E402

POLYFORM_HEADER = ("/*\n"
                   f" * {lc.NOTICE}\n"
                   f" * {lc.POLYFORM}\n"
                   " */\n")
RESERVED_HEADER = ("/*\n"
                   f" * {lc.NOTICE}\n"
                   f" * {lc.RESERVED} No redistribution or use without explicit\n"
                   " * permission.\n"
                   " */\n")
LICENSE_HEADING = ("Exceptions outside `effects/`, each carrying its own terms "
                   "in its own header:\n\n")


class TestLicenseExceptions(unittest.TestCase):
    def test_repository_license_matches_executable_exceptions(self):
        license_text = (TOOLS.parent / "LICENSE").read_text(encoding="utf-8")
        self.assertEqual(lc.license_exception_paths(license_text),
                         set(lc.EXCEPTIONS))

    def test_paths_must_match_in_both_directions(self):
        license_text = LICENSE_HEADING + "- `licensed.h` — its own terms\n"
        with mock.patch.object(lc, "EXCEPTIONS", {"configured.h": lc.RESERVED}):
            issues = lc.license_exception_issues(license_text)
        self.assertEqual(len(issues), 2)
        self.assertIn("licensed.h", issues[0])
        self.assertIn("configured.h", issues[1])


class TestExpectedMarker(unittest.TestCase):
    def test_effects_are_all_rights_reserved(self):
        self.assertEqual(lc.expected_marker("effects/Voronoi.h"), lc.RESERVED)

    def test_everything_else_is_polyform(self):
        for path in ("core/math/3dmath.h", "hardware/pov_sync.h",
                     "tests/run_tests.cpp", "targets/wasm/wasm.cpp"):
            self.assertEqual(lc.expected_marker(path), lc.POLYFORM, path)

    def test_named_exceptions_carry_their_own_terms(self):
        self.assertEqual(lc.expected_marker("core/engine/effects_legacy.h"),
                         lc.RESERVED)
        self.assertEqual(lc.expected_marker("core/math/projections.h"),
                         lc.MIT_GRANT)
        self.assertEqual(lc.expected_marker("core/vendor/FastNoiseLite.h"),
                         lc.MIT_TITLE)
        self.assertEqual(lc.expected_marker("workbench/shader/shader_host.h"),
                         lc.RESERVED)

    def test_the_first_party_file_under_core_vendor_is_polyform(self):
        self.assertEqual(
            lc.expected_marker("core/vendor/FastNoiseLite_config.h"),
            lc.POLYFORM)


class TestHeaderIssue(unittest.TestCase):
    def test_matching_headers_pass(self):
        self.assertIsNone(lc.header_issue("core/math/3dmath.h", POLYFORM_HEADER))
        self.assertIsNone(lc.header_issue("effects/Voronoi.h", RESERVED_HEADER))

    def test_reserving_rights_LICENSE_grants_is_reported(self):
        issue = lc.header_issue("core/math/noise_field.h", RESERVED_HEADER)
        self.assertIn(lc.POLYFORM, issue)

    def test_granting_rights_LICENSE_reserves_is_reported(self):
        issue = lc.header_issue("effects/Voronoi.h", POLYFORM_HEADER)
        self.assertIn(lc.RESERVED, issue)

    def test_a_header_carrying_both_markers_is_reported(self):
        both = POLYFORM_HEADER.replace(" */\n", f" * {lc.RESERVED}\n */\n")
        for path in ("core/math/noise_field.h", "effects/Voronoi.h"):
            self.assertIsNotNone(lc.header_issue(path, both), path)

    def test_a_third_party_header_carrying_first_party_terms_is_reported(self):
        head = (f"// {lc.MIT_TITLE}\n// {lc.MIT_GRANT}, free of charge\n"
                f"// {lc.POLYFORM}\n// {lc.RESERVED}\n")
        issue = lc.header_issue("core/vendor/FastNoiseLite.h", head)
        self.assertIn(lc.POLYFORM, issue)

    def test_a_first_party_header_carrying_an_MIT_grant_is_reported(self):
        head = POLYFORM_HEADER.replace(" */\n", f" * {lc.MIT_GRANT}\n */\n")
        issue = lc.header_issue("core/math/noise_field.h", head)
        self.assertIn(lc.MIT_GRANT, issue)

    def test_every_marker_a_path_can_expect_has_contradictions(self):
        markers = set(lc.EXCEPTIONS.values()) | {lc.POLYFORM, lc.RESERVED}
        self.assertEqual(markers - set(lc.CONTRADICTIONS), set())

    def test_a_file_with_no_notice_is_reported(self):
        issue = lc.header_issue("core/math/3dmath.h", "#pragma once\n")
        self.assertIn("no copyright notice", issue)

    def test_a_third_party_header_stripped_of_its_banner_is_reported(self):
        issue = lc.header_issue("core/vendor/FastNoiseLite.h",
                                "#pragma once\n")
        self.assertIn(lc.MIT_TITLE, issue)

    def test_a_third_party_header_needs_no_first_party_notice(self):
        head = f"// {lc.MIT_TITLE}\n// {lc.MIT_GRANT}, free of charge\n"
        self.assertIsNone(lc.header_issue("core/vendor/FastNoiseLite.h", head))

    def test_the_alternate_spelling_of_the_license_name_passes(self):
        head = POLYFORM_HEADER.replace("PolyForm", "Polyform")
        self.assertIsNone(lc.header_issue("tests/x.h", head))

    def test_a_header_past_the_scanned_window_is_reported(self):
        head = "// " + "x" * lc.HEAD_BYTES + "\n" + POLYFORM_HEADER
        self.assertIsNotNone(lc.header_issue("core/math/3dmath.h",
                                             head[:lc.HEAD_BYTES]))


class TestMain(unittest.TestCase):
    def test_main_checks_tracked_sources_and_returns_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            source = root / "sample.h"
            source.write_text(POLYFORM_HEADER, encoding="utf-8")
            (root / "LICENSE").write_text(LICENSE_HEADING, encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "sample.h"],
                           check=True)
            with mock.patch.object(lc, "EXCEPTIONS", {}):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(lc.main(["--root", str(root)]), 0)

                source.write_text("#pragma once\n", encoding="utf-8")
                with contextlib.redirect_stdout(io.StringIO()), \
                        contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(lc.main(["--root", str(root)]), 1)

    def test_stale_exception_is_a_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            source = root / "sample.h"
            source.write_text(POLYFORM_HEADER, encoding="utf-8")
            (root / "LICENSE").write_text(
                LICENSE_HEADING + "- `retired.h` — retired terms\n",
                encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "sample.h"],
                           check=True)
            exceptions = {"retired.h": lc.RESERVED}
            with mock.patch.object(lc, "EXCEPTIONS", exceptions), \
                    contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(lc.main(["--root", str(root)]), 1)


if __name__ == "__main__":
    unittest.main()
