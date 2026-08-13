import json
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest import mock

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import analyze_candidates  # noqa: E402
import kicad_common  # noqa: E402


class CandidateBoardTests(unittest.TestCase):
    def test_accepts_any_single_board_filename(self):
        with tempfile.TemporaryDirectory() as directory:
            board = Path(directory) / "phantasm.kicad_pcb"
            board.touch()

            self.assertEqual(
                analyze_candidates.candidate_board(directory), str(board))

    def test_rejects_missing_board(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "0 found"):
                analyze_candidates.candidate_board(directory)

    def test_rejects_ambiguous_boards(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "a.kicad_pcb").touch()
            Path(directory, "b.kicad_pcb").touch()

            with self.assertRaisesRegex(ValueError, "2 found"):
                analyze_candidates.candidate_board(directory)


class ParserTests(unittest.TestCase):
    def test_parentheses_inside_strings_do_not_split_blocks(self):
        source = (
            '(kicad_pcb '
            '(footprint "A" (descr "contains ) and ( characters")) '
            '(footprint "B"))'
        )

        found = analyze_candidates.blocks(source, "footprint")

        self.assertEqual([str(block[1]) for block in found], ["A", "B"])

    def test_uses_shared_kicad_cli_discovery(self):
        self.assertEqual(analyze_candidates.KCLI,
                         kicad_common.find_kicad_cli())


SYNTHETIC_BOARD = """
(kicad_pcb
  (segment (start 0 0) (end 10 0) (width 0.25) (layer "F.Cu") (net "/DATA"))
  (via (at 3 4) (size 0.6) (drill 0.3) (layers "F.Cu" "B.Cu") (net "GND"))
  (footprint "Resistor_SMD:R_0402"
    (at 10 20)
    (property "Reference" "R1" (at 0 0 0))
    (property "Value" "10k" (at 0 0 0)))
)
"""


class AnalyzeTests(unittest.TestCase):
    def analyze_source(self, source):
        with tempfile.TemporaryDirectory() as directory:
            board = Path(directory) / "phantasm.kicad_pcb"
            board.write_text(source, encoding="utf-8")
            return analyze_candidates.analyze(str(board))

    def test_reads_segment_via_and_footprint(self):
        r = self.analyze_source(SYNTHETIC_BOARD)

        self.assertEqual(r["nseg"], 1)
        self.assertEqual(r["nvia"], 1)
        self.assertAlmostEqual(r["total_len"], 10.0)
        # net_of strips the hierarchical path prefix
        self.assertAlmostEqual(r["netlen"]["DATA"], 10.0)
        self.assertEqual(r["netseg"]["DATA"], 1)
        self.assertEqual(r["netlayers"]["DATA"], {"F.Cu"})
        self.assertEqual(r["netvias"]["GND"], 1)
        self.assertEqual(r["gnd_vias"], 1)
        self.assertEqual(r["small_vias"], 0)
        self.assertEqual(r["pos"]["R1"], (10.0, 20.0))

    def test_flags_undersized_via(self):
        source = SYNTHETIC_BOARD.replace("(size 0.6) (drill 0.3)",
                                         "(size 0.4) (drill 0.15)")

        self.assertEqual(self.analyze_source(source)["small_vias"], 1)

    def test_reads_the_id_first_net_form(self):
        source = (SYNTHETIC_BOARD.replace('(net "/DATA")', '(net 9 "/DATA")')
                  .replace('(net "GND")', '(net 1 "GND")'))
        r = self.analyze_source(source)

        self.assertAlmostEqual(r["netlen"]["DATA"], 10.0)
        self.assertAlmostEqual(r["crit_len"], 10.0)
        self.assertEqual(r["netvias"]["GND"], 1)

    def test_rejects_a_board_whose_nets_are_ids_only(self):
        source = SYNTHETIC_BOARD.replace('(net "/DATA")', "(net 9)")

        with self.assertRaisesRegex(ValueError, "no critical net"):
            self.analyze_source(source)


class ResolveKicadCliTests(unittest.TestCase):
    def test_absolute_install_path_is_used_as_is(self):
        with tempfile.TemporaryDirectory() as directory:
            cli = Path(directory) / "kicad-cli"
            cli.touch()

            with mock.patch.object(analyze_candidates, "KCLI", str(cli)):
                self.assertEqual(analyze_candidates.resolve_kicad_cli(), str(cli))

    def test_bare_name_resolves_through_path(self):
        with mock.patch.object(analyze_candidates, "KCLI", "kicad-cli"), \
                mock.patch.object(analyze_candidates.shutil, "which",
                                  return_value="/opt/homebrew/bin/kicad-cli") as which:
            self.assertEqual(analyze_candidates.resolve_kicad_cli(),
                             "/opt/homebrew/bin/kicad-cli")
        which.assert_called_once_with("kicad-cli")

    def test_unresolvable_name_is_none(self):
        with mock.patch.object(analyze_candidates, "KCLI", "kicad-cli"), \
                mock.patch.object(analyze_candidates.shutil, "which",
                                  return_value=None):
            self.assertIsNone(analyze_candidates.resolve_kicad_cli())

    def test_drc_gate_reports_missing_only_when_unresolvable(self):
        with mock.patch.object(analyze_candidates, "resolve_kicad_cli",
                               return_value=None):
            self.assertEqual(analyze_candidates.run_drc("board.kicad_pcb")["status"],
                             analyze_candidates.DRC_MISSING)


class RunDrcReportTests(unittest.TestCase):
    """The gate reads kicad-cli's JSON report, so a shape change can't read clean."""

    def run_drc(self, report_text):
        def write_report(cmd, **kwargs):
            Path(cmd[cmd.index("-o") + 1]).write_text(report_text,
                                                      encoding="utf-8")
            return subprocess.CompletedProcess(cmd, 0)

        with mock.patch.object(analyze_candidates, "resolve_kicad_cli",
                               return_value="kicad-cli"), \
                mock.patch.object(analyze_candidates.subprocess, "run",
                                  side_effect=write_report) as run:
            result = analyze_candidates.run_drc("board.kicad_pcb")
        self.assertIn("json", run.call_args.args[0])
        return result

    def test_counts_violations_by_rule_and_unconnected_items(self):
        report = self.run_drc(json.dumps({
            "violations": [{"type": "clearance"}, {"type": "clearance"},
                           {"type": "shorting_items"}],
            "unconnected_items": [{"type": "unconnected_items"}],
        }))

        self.assertEqual(report["status"], analyze_candidates.DRC_OK)
        self.assertEqual(report["errors"], 3)
        self.assertEqual(report["unconnected"], 1)
        self.assertEqual(report["by_type"],
                         Counter({"clearance": 2, "shorting_items": 1}))
        self.assertEqual(analyze_candidates.real_faults(report["by_type"]), 1)

    def test_clean_report_reports_no_faults(self):
        report = self.run_drc(json.dumps({"violations": [],
                                          "unconnected_items": []}))

        self.assertEqual(report["status"], analyze_candidates.DRC_OK)
        self.assertEqual((report["errors"], report["unconnected"]), (0, 0))

    def test_report_without_sections_fails_rather_than_reading_clean(self):
        self.assertEqual(self.run_drc(json.dumps({"violations": []}))["status"],
                         analyze_candidates.DRC_FAILED)


if __name__ == "__main__":
    unittest.main()
