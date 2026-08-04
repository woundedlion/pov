import sys
import tempfile
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import analyze_candidates  # noqa: E402
import fab  # noqa: E402


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
        self.assertEqual(analyze_candidates.KCLI, fab.find_kicad_cli())


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


if __name__ == "__main__":
    unittest.main()
