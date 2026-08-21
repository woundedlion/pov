import re
import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import builder  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F  # noqa: E402

PROJ = GEN.parent
ROUTED = PROJ / "phantasm.kicad_pcb"
SCH = PROJ / "phantasm.kicad_sch"

SILK = re.compile(r"^Phantasm Rev (\S+)$")


def title_block_rev(path):
    root = sexp.parse(path.read_text(encoding="utf-8"))[0]
    for block in F(root, "title_block"):
        rev = sexp.val(block, "rev")
        if rev:
            return str(rev[0])
    return None


def silk_revisions():
    root = sexp.parse(ROUTED.read_text(encoding="utf-8"))[0]
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

    def setUp(self):
        revisions = silk_revisions()
        self.assertEqual(len(revisions), 1, revisions)
        self.revision = revisions[0]

    def test_routed_board_title_block_carries_the_silk_revision(self):
        self.assertEqual(title_block_rev(ROUTED), self.revision)

    def test_schematic_title_block_carries_the_silk_revision(self):
        self.assertEqual(title_block_rev(SCH), self.revision)

    def test_the_schematic_generator_emits_the_silk_revision(self):
        self.assertEqual(builder.REVISION, self.revision)


if __name__ == "__main__":
    unittest.main()
