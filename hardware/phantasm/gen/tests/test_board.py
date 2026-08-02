import sys
import unittest
import unittest.mock
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))


class BoardEntryPointTests(unittest.TestCase):
    def test_import_ignores_host_arguments(self):
        with unittest.mock.patch.object(sys, "argv", ["host", "--unknown"]):
            import board

        self.assertTrue(callable(board.main))


if __name__ == "__main__":
    unittest.main()
