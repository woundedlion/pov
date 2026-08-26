import json
import math
import subprocess
import sys
import tempfile
import unittest
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

    def test_default_discovery_matches_quilter_folder_names(self):
        with tempfile.TemporaryDirectory() as directory:
            candidates = Path(directory) / "candidates"
            expected = candidates / "Candidate 1"
            expected.mkdir(parents=True)
            (candidates / "Candidate_2").mkdir()

            with mock.patch.object(analyze_candidates, "PROJ", directory):
                self.assertEqual(analyze_candidates.default_candidates(),
                                 [str(expected)])

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
        self.assertEqual(analyze_candidates.kicad_cli, kicad_common.kicad_cli)


SYNTHETIC_BOARD = """
(kicad_pcb
  (segment (start 0 0) (end 10 0) (width 0.25) (layer "F.Cu") (net "/DATA"))
  (segment (start 0 1) (end 1 1) (width 0.25) (layer "F.Cu") (net "/CLK"))
  (segment (start 0 2) (end 1 2) (width 0.25) (layer "F.Cu") (net "/DATA_IN"))
  (segment (start 0 3) (end 1 3) (width 0.25) (layer "F.Cu") (net "/CLK_IN"))
  (segment (start 0 4) (end 1 4) (width 0.25) (layer "F.Cu") (net "/DATA_SRC"))
  (segment (start 0 5) (end 1 5) (width 0.25) (layer "F.Cu") (net "/CLK_SRC"))
  (segment (start 0 6) (end 1 6) (width 0.25) (layer "F.Cu") (net "/SYNC_BUS"))
  (segment (start 0 7) (end 1 7) (width 0.25) (layer "F.Cu") (net "/FRAME_SYNC"))
  (segment (start 0 8) (end 1 8) (width 0.25) (layer "F.Cu") (net "/SYNC_SRC"))
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

        self.assertEqual(r["nseg"], len(analyze_candidates.CRIT))
        self.assertEqual(r["nvia"], 1)
        self.assertAlmostEqual(r["total_len"], 18.0)
        # net_name drops the leading slash
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
        self.assertAlmostEqual(r["crit_len"], 18.0)
        self.assertEqual(r["netvias"]["GND"], 1)

    def test_rejects_a_board_whose_nets_are_ids_only(self):
        source = SYNTHETIC_BOARD
        for net in analyze_candidates.CRIT:
            source = source.replace(f'(net "/{net}")', "(net 9)")

        with self.assertRaisesRegex(ValueError, "no critical net"):
            self.analyze_source(source)

    def test_rejects_a_missing_critical_net(self):
        source = SYNTHETIC_BOARD.replace(
            '  (segment (start 0 0) (end 10 0) (width 0.25) '
            '(layer "F.Cu") (net "/DATA"))\n', '')

        with self.assertRaisesRegex(ValueError, "unrouted critical nets.*DATA"):
            self.analyze_source(source)

    def test_rejects_a_zero_length_critical_net(self):
        source = SYNTHETIC_BOARD.replace("(end 10 0)", "(end 0 0)")

        with self.assertRaisesRegex(ValueError, "unrouted critical nets.*DATA"):
            self.analyze_source(source)


class ScoreTests(unittest.TestCase):
    def test_fully_routed_candidates_keep_the_shorter_copper_ordering(self):
        base = {
            "spi_vias": 0,
            "sync_vias": 0,
            "ergo": {
                "decap_u1": 0.0,
                "divider": 0.0,
                "term_j2": 0.0,
            },
        }

        short = analyze_candidates.score({**base, "crit_len": 90.0})
        long = analyze_candidates.score({**base, "crit_len": 180.0})

        self.assertEqual(short, (9.0, 10.0))
        self.assertEqual(long, (8.0, 10.0))
        self.assertGreater(short[0], long[0])


class MainTests(unittest.TestCase):
    def test_unrouted_candidate_stops_before_skew_or_score_reports(self):
        source = SYNTHETIC_BOARD.replace("(end 10 0)", "(end 0 0)")
        with tempfile.TemporaryDirectory() as directory:
            board = Path(directory) / "Candidate 1.kicad_pcb"
            board.write_text(source, encoding="utf-8")

            with mock.patch("builtins.print") as emit, \
                    mock.patch.object(analyze_candidates, "run_drc") as run_drc:
                result = analyze_candidates.main([str(board)])

        output = "\n".join(" ".join(map(str, call.args))
                           for call in emit.call_args_list)
        self.assertEqual(result, 1)
        self.assertIn("unrouted critical nets", output)
        self.assertNotIn("PAIR SKEW", output)
        self.assertNotIn("COMPOSITE SCORE", output)
        run_drc.assert_not_called()

    def test_two_runs_of_one_label_stop_rather_than_rank_one(self):
        with tempfile.TemporaryDirectory() as directory:
            boards = []
            for run in ("runA", "runB"):
                folder = Path(directory) / run / "Candidate 1"
                folder.mkdir(parents=True)
                board = folder / "phantasm.kicad_pcb"
                board.write_text(SYNTHETIC_BOARD, encoding="utf-8")
                boards.append(str(folder))

            with mock.patch("builtins.print") as emit, \
                    mock.patch.object(analyze_candidates, "run_drc") as run_drc:
                result = analyze_candidates.main(boards)

        output = "\n".join(" ".join(map(str, call.args))
                           for call in emit.call_args_list)
        self.assertEqual(result, 1)
        self.assertIn("names two boards", output)
        for board in boards:
            self.assertIn(board, output)
        self.assertEqual(run_drc.call_count, 1)


class ClosestSpacingTests(unittest.TestCase):
    NAN = float("nan")

    def test_returns_the_smaller_spacing(self):
        self.assertEqual(analyze_candidates.closest(5.0, 3.0), 3.0)

    def test_a_missing_part_forfeits_whichever_one_vanished(self):
        self.assertTrue(math.isnan(analyze_candidates.closest(self.NAN, 3.0)))
        self.assertTrue(math.isnan(analyze_candidates.closest(3.0, self.NAN)))


class ResolveKicadCliTests(unittest.TestCase):
    def test_absolute_install_path_is_used_as_is(self):
        with tempfile.TemporaryDirectory() as directory:
            cli = Path(directory) / "kicad-cli"
            cli.touch()

            with mock.patch.object(analyze_candidates, "kicad_cli",
                                   return_value=str(cli)):
                self.assertEqual(analyze_candidates.resolve_kicad_cli(), str(cli))

    def test_bare_name_resolves_through_path(self):
        with mock.patch.object(analyze_candidates, "kicad_cli",
                               return_value="kicad-cli"), \
                mock.patch.object(analyze_candidates.shutil, "which",
                                  return_value="/opt/homebrew/bin/kicad-cli") as which:
            self.assertEqual(analyze_candidates.resolve_kicad_cli(),
                             "/opt/homebrew/bin/kicad-cli")
        which.assert_called_once_with("kicad-cli")

    def test_unresolvable_name_is_none(self):
        with mock.patch.object(analyze_candidates, "kicad_cli",
                               return_value="kicad-cli"), \
                mock.patch.object(analyze_candidates.shutil, "which",
                                  return_value=None):
            self.assertIsNone(analyze_candidates.resolve_kicad_cli())

    def test_an_unusable_kicad_is_ungated_not_fatal(self):
        with mock.patch.object(analyze_candidates, "kicad_cli",
                               side_effect=SystemExit("KiCad 9, gates need 10")):
            self.assertIsNone(analyze_candidates.resolve_kicad_cli())

    def test_drc_gate_reports_missing_only_when_unresolvable(self):
        with mock.patch.object(analyze_candidates, "resolve_kicad_cli",
                               return_value=None):
            self.assertEqual(analyze_candidates.run_drc("board.kicad_pcb")["status"],
                             analyze_candidates.DRC_MISSING)


def zone_clearance():
    """A missing-antipad clearance error: a via against a pour."""
    return {"type": "clearance",
            "items": [{"description": "Via [/DATA] on F.Cu - B.Cu"},
                      {"description": "Zone [GND] on In1.Cu"}]}


def track_clearance():
    """A clearance error between two different-net tracks -- no refill clears it."""
    return {"type": "clearance",
            "items": [{"description": "Track [/DATA] on F.Cu, length 3.4 mm"},
                      {"description": "Track [/CLK] on F.Cu, length 2.1 mm"}]}


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

    def test_counts_violations_and_unconnected_items(self):
        report = self.run_drc(json.dumps({
            "violations": [zone_clearance(), zone_clearance(),
                           {"type": "shorting_items", "items": [
                               {"description": "Track [/DATA] on F.Cu"},
                               {"description": "Track [/CLK] on F.Cu"}]}],
            "unconnected_items": [{"type": "unconnected_items"}],
        }))

        self.assertEqual(report["status"], analyze_candidates.DRC_OK)
        self.assertEqual(report["errors"], 3)
        self.assertEqual(report["unconnected"], 1)
        self.assertEqual(report["real"], 1)

    def test_track_to_track_clearance_counts_as_a_real_fault(self):
        report = self.run_drc(json.dumps({
            "violations": [zone_clearance(), track_clearance()],
            "unconnected_items": [],
        }))

        self.assertEqual(report["errors"], 2)
        self.assertEqual(report["real"], 1)

    def test_zone_clearance_is_refill_fixable(self):
        self.assertTrue(analyze_candidates.refill_fixable(zone_clearance()))

    def test_shorting_zone_is_not_refill_fixable(self):
        self.assertFalse(analyze_candidates.refill_fixable(
            {"type": "shorting_items",
             "items": [{"description": "Zone [GND] on In1.Cu"},
                       {"description": "Track [/DATA] on F.Cu"}]}))

    def test_clearance_without_items_is_a_real_fault(self):
        self.assertFalse(analyze_candidates.refill_fixable({"type": "clearance"}))
        self.assertFalse(analyze_candidates.refill_fixable("clearance"))

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
