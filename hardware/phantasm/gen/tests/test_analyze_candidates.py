import sys
import tempfile
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import analyze_candidates  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
