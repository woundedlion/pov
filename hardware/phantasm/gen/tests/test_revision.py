import re
import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import builder  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F  # noqa: E402

REVISIONS = ("1.1", "1.2")

SILK = re.compile(r"^Phantasm Rev (\S+)$")


def title_block_rev(path):
    root = sexp.parse(path.read_text(encoding="utf-8"))[0]
    for block in F(root, "title_block"):
        rev = sexp.val(block, "rev")
        if rev:
            return str(rev[0])
    return None


def silk_revisions(board):
    root = sexp.parse(board.read_text(encoding="utf-8"))[0]
    found = []
    for node in root:
        if not isinstance(node, list) or not node or node[0] != "gr_text":
            continue
        match = SILK.match(str(node[1]))
        if match:
            found.append(match.group(1))
    return found


class RevisionTests(unittest.TestCase):
    """One board revision, spelled in three places a reader or a fab reads.

    The silkscreen is what a built board carries; the schematic title block
    labels the sheet; the routed board's title block is what KiCad writes into
    the Gerber X2 ProjectId attribute, which reads `rev?` when it is absent.
    """

    def test_revision_labels(self):
        for revision in REVISIONS:
            with self.subTest(revision=revision):
                project = GEN.parent / revision
                board = project / "phantasm.kicad_pcb"
                self.assertEqual(silk_revisions(board), [revision])
                self.assertEqual(title_block_rev(board), revision)
                self.assertEqual(title_block_rev(project / "phantasm.kicad_sch"), revision)

    def test_generator_targets_revision_1_2(self):
        self.assertEqual(builder.REVISION, "1.2")

    def test_the_board_generator_stamps_the_silk_from_builder(self):
        self.assertEqual(pcb.SILK_REVISION,
                         f"Phantasm Rev {builder.REVISION}")


if __name__ == "__main__":
    unittest.main()
