import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import kicad_common  # noqa: E402
import sexp  # noqa: E402


class FindKicadCliTests(unittest.TestCase):
    """Every tool that shells out to KiCad resolves the binary through here."""

    def test_env_override_wins_when_it_exists(self):
        with tempfile.TemporaryDirectory() as directory:
            cli = Path(directory) / "kicad-cli"
            cli.touch()
            with mock.patch.dict(os.environ, {"KICAD_CLI": str(cli)}), \
                    mock.patch.object(kicad_common, "kicad_cli_major",
                                      return_value=sexp.KICAD_MAJOR):
                self.assertEqual(kicad_common.find_kicad_cli(), str(cli))

    def test_env_override_rejects_another_major(self):
        with tempfile.TemporaryDirectory() as directory:
            cli = Path(directory) / "kicad-cli"
            cli.touch()
            with mock.patch.dict(os.environ, {"KICAD_CLI": str(cli)}), \
                    mock.patch.object(kicad_common, "kicad_cli_major",
                                      return_value=sexp.KICAD_MAJOR + 1), \
                    self.assertRaisesRegex(SystemExit, "fab gates require"):
                kicad_common.find_kicad_cli()

    def test_falls_back_to_the_path_name(self):
        with mock.patch.dict(os.environ, {"KICAD_CLI": "missing-kicad-cli"}), \
                mock.patch.object(kicad_common.glob, "glob", return_value=[]):
            self.assertEqual(kicad_common.find_kicad_cli(), "kicad-cli")

    def windows_install(self, version):
        return rf"C:\Program Files\KiCad\{version}\bin\kicad-cli.exe"

    def resolve(self, installs):
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(kicad_common.glob, "glob",
                                  side_effect=lambda p: installs if "Program" in p
                                  and "x86" not in p else []):
            return kicad_common.find_kicad_cli()

    def test_prefers_the_pinned_major_over_a_newer_install(self):
        pinned = self.windows_install(f"{sexp.KICAD_MAJOR}.0")
        newer = self.windows_install(f"{sexp.KICAD_MAJOR + 1}.0")

        self.assertEqual(self.resolve([pinned, newer]), pinned)

    def test_takes_the_newest_minor_of_the_pinned_major(self):
        installs = [self.windows_install(f"{sexp.KICAD_MAJOR}.{minor}")
                    for minor in (0, 3)]

        self.assertEqual(self.resolve(installs), installs[1])

    def test_exits_when_no_install_is_the_pinned_major(self):
        newer = self.windows_install(f"{sexp.KICAD_MAJOR + 1}.0")

        with self.assertRaises(SystemExit) as caught:
            self.resolve([newer])
        self.assertIn(f"KiCad {sexp.KICAD_MAJOR}", str(caught.exception))
        self.assertIn(newer, str(caught.exception))
        self.assertIn("KICAD_CLI", str(caught.exception))

    def test_an_unversioned_path_is_asked_for_its_version(self):
        completed = subprocess.CompletedProcess(
            [], 0, stdout=f"{sexp.KICAD_MAJOR}.0.1\n")
        with mock.patch.object(kicad_common.subprocess, "run",
                               return_value=completed) as run:
            self.assertEqual(self.resolve(["/usr/bin/kicad-cli"]),
                             "/usr/bin/kicad-cli")
        self.assertEqual(run.call_args.args[0],
                         ["/usr/bin/kicad-cli", "--version"])

    def test_an_unversioned_path_reporting_another_major_is_refused(self):
        completed = subprocess.CompletedProcess(
            [], 0, stdout=f"{sexp.KICAD_MAJOR + 1}.0.0\n")
        with mock.patch.object(kicad_common.subprocess, "run",
                               return_value=completed), \
                self.assertRaises(SystemExit):
            self.resolve(["/usr/bin/kicad-cli"])


class KicadCliTests(unittest.TestCase):
    """Resolution is deferred to first use: find_kicad_cli() exits when no
    install is the pinned major, and an exit during import takes down every
    test in the importing suite, KiCad-dependent or not."""

    def setUp(self):
        patcher = mock.patch.object(kicad_common, "_KCLI", None)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_resolves_once_and_memoizes(self):
        with mock.patch.object(kicad_common, "find_kicad_cli",
                               return_value="kicad-cli") as find:
            self.assertEqual(kicad_common.kicad_cli(), "kicad-cli")
            self.assertEqual(kicad_common.kicad_cli(), "kicad-cli")
        find.assert_called_once_with()

    def load_fresh(self, name):
        """Execute a gate module into a throwaway namespace."""
        spec = importlib.util.spec_from_file_location(name, GEN / f"{name}.py")
        spec.loader.exec_module(importlib.util.module_from_spec(spec))

    def test_importing_a_gate_module_resolves_nothing(self):
        with mock.patch.object(kicad_common, "find_kicad_cli",
                               side_effect=SystemExit):
            for name in ("analyze_candidates", "check", "fab", "pcb"):
                with self.subTest(module=name):
                    self.load_fresh(name)


class ExportNetlistTests(unittest.TestCase):
    """The netlist gates run from a shell; a failed export must read as one."""

    def export(self, error):
        with mock.patch.object(kicad_common.subprocess, "run",
                               side_effect=error), \
                self.assertRaises(SystemExit) as caught:
            kicad_common.export_netlist("kicad-cli", "phantasm.kicad_sch")
        return str(caught.exception)

    def test_reports_a_missing_kicad_cli(self):
        message = self.export(FileNotFoundError())
        self.assertIn("kicad-cli not found: kicad-cli", message)
        self.assertIn("KICAD_CLI", message)

    def test_reports_a_failed_export(self):
        message = self.export(subprocess.CalledProcessError(
            2, [], stderr="schematic is broken\n"))
        self.assertIn("kicad-cli exited 2", message)
        self.assertIn("phantasm.kicad_sch", message)


class RequireWritableTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.path = Path(directory.name) / "phantasm.kicad_pcb"

    def refusal(self, **kwargs):
        with self.assertRaises(SystemExit) as caught:
            kicad_common.require_writable(self.path, False, **kwargs)
        return str(caught.exception)

    def test_allows_an_absent_target(self):
        self.assertIsNone(kicad_common.require_writable(self.path, False))

    def test_allows_a_forced_overwrite(self):
        self.path.write_text("routed", encoding="utf-8")

        self.assertIsNone(kicad_common.require_writable(self.path, True))

    def test_refuses_an_existing_target(self):
        self.path.write_text("routed", encoding="utf-8")

        message = self.refusal()
        self.assertIn(f"refusing to overwrite {self.path}", message)
        self.assertIn("routing, vias, silk, hand edits", message)
        self.assertIn("--force", message)

    def test_states_the_caller_reason(self):
        self.path.write_text("routed", encoding="utf-8")

        message = self.refusal(reason="It is the fabrication source.")
        self.assertIn("It is the fabrication source.", message)
        self.assertNotIn("routing, vias, silk, hand edits", message)

    def test_names_the_authorizing_flag(self):
        self.path.write_text("routed", encoding="utf-8")

        message = self.refusal(flag="--force-teensy-library")
        self.assertIn("Re-run with --force-teensy-library", message)


if __name__ == "__main__":
    unittest.main()
