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


if __name__ == "__main__":
    unittest.main()
