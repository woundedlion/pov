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


class FindKicadCliTests(unittest.TestCase):
    def resolve(self, reported, installs=None):
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(kicad_common.glob, "glob",
                                  return_value=installs or []), \
                mock.patch.object(kicad_common.subprocess, "run",
                                  return_value=subprocess.CompletedProcess(
                                      [], 0, stdout=reported)) as run:
            result = kicad_common.find_kicad_cli()
        self.assertEqual(run.call_args.args[0][-1], "--version")
        return result

    def test_exact_release_on_path(self):
        self.assertEqual(self.resolve(kicad_common.KICAD_VERSION), "kicad-cli")

    def test_adjacent_patch_is_rejected(self):
        with self.assertRaises(SystemExit):
            self.resolve("10.0.5")

    def test_windows_directory_does_not_substitute_for_reported_version(self):
        cli = r"C:\Program Files\KiCad\10.0\bin\kicad-cli.exe"
        with self.assertRaises(SystemExit):
            self.resolve("10.0.5", [cli])
        self.assertEqual(self.resolve(kicad_common.KICAD_VERSION, [cli]), cli)

    def test_env_override_requires_exact_release(self):
        with tempfile.TemporaryDirectory() as directory:
            cli = Path(directory) / "kicad-cli"
            cli.touch()
            with mock.patch.dict(os.environ, {"KICAD_CLI": str(cli)}), \
                    mock.patch.object(kicad_common, "kicad_cli_version",
                                      return_value="10.0.5"), \
                    self.assertRaises(SystemExit):
                kicad_common.find_kicad_cli()

    def test_unreadable_release(self):
        with mock.patch.object(kicad_common.subprocess, "run",
                               side_effect=OSError):
            self.assertIsNone(kicad_common.kicad_cli_version("missing"))

    def test_major_only_release_is_rejected(self):
        with self.assertRaises(SystemExit):
            self.resolve("10.0")


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


class AtomicWriteTests(unittest.TestCase):
    def test_failed_replace_keeps_original(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "board.kicad_sch"
            path.write_text("original", encoding="utf-8")
            with mock.patch.object(kicad_common.os, "replace", side_effect=OSError("failed")):
                with self.assertRaises(OSError):
                    kicad_common.atomic_write_text(path, "replacement")
            self.assertEqual(path.read_text(encoding="utf-8"), "original")
            self.assertEqual(list(Path(directory).iterdir()), [path])
            kicad_common.atomic_write_text(path, "complete\n")
            self.assertEqual(path.read_bytes(), b"complete\n")
