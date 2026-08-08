#!/usr/bin/env python3
"""Host tests for the relax-bake header generator (tools/relax_bakes.py)."""

import sys
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import relax_bakes  # noqa: E402


def make_dump(name, iterations, verts, topology_hash):
    """Build a RELAX_BAKE block whose declared output hash matches its bits."""
    words = []
    for x, y, z in verts:
        words.extend((x, y, z))
    out = relax_bakes.vertex_hash(words)
    lines = [
        f"RELAX_BAKE_BEGIN {name} {iterations} {len(verts)} 2 6 "
        f"{topology_hash:08x} {out:08x}"
    ]
    for x, y, z in verts:
        lines.append(f"RELAX_BAKE_DATA {x:08x} {y:08x} {z:08x}")
    lines.append("RELAX_BAKE_END")
    return "\n".join(lines), out


class ParseDump(unittest.TestCase):
    def test_parses_block_fields_and_bits(self):
        dump, out = make_dump("foo", 100, [(1, 2, 3), (4, 5, 6)], 0xABCD1234)
        bakes = relax_bakes.parse_dump(dump)
        self.assertEqual(len(bakes), 1)
        b = bakes[0]
        self.assertEqual(b["name"], "foo")
        self.assertEqual(b["iterations"], 100)
        self.assertEqual(b["vertices"], 2)
        self.assertEqual(b["topology_hash"], 0xABCD1234)
        self.assertEqual(b["output_hash"], out)
        self.assertEqual(b["bits"], [1, 2, 3, 4, 5, 6])

    def test_dedupes_identical_repeat(self):
        dump, _ = make_dump("foo", 100, [(1, 2, 3)], 0x1)
        bakes = relax_bakes.parse_dump(dump + "\n" + dump)
        self.assertEqual(len(bakes), 1)

    def test_preserves_first_seen_order(self):
        a, _ = make_dump("a", 8, [(1, 1, 1)], 0x1)
        b, _ = make_dump("b", 8, [(2, 2, 2)], 0x2)
        bakes = relax_bakes.parse_dump(b + "\n" + a)
        self.assertEqual([x["name"] for x in bakes], ["b", "a"])

    def test_rejects_output_hash_mismatch(self):
        dump, _ = make_dump("foo", 100, [(1, 2, 3)], 0x1)
        dump = dump.replace(
            "RELAX_BAKE_DATA 00000001 00000002 00000003",
            "RELAX_BAKE_DATA deadbeef 00000002 00000003",
        )
        with self.assertRaises(ValueError):
            relax_bakes.parse_dump(dump)

    def test_rejects_nondeterministic_duplicate(self):
        a, _ = make_dump("foo", 8, [(1, 1, 1)], 0x1)
        b, _ = make_dump("foo", 8, [(9, 9, 9)], 0x1)
        with self.assertRaises(ValueError):
            relax_bakes.parse_dump(a + "\n" + b)

    def test_empty_dump_raises(self):
        with self.assertRaises(ValueError):
            relax_bakes.parse_dump("no blocks here\n")

    def test_rejects_unterminated_block(self):
        # Harness died mid-block: the partial bake must not be dropped.
        dump, _ = make_dump("foo", 8, [(1, 1, 1), (2, 2, 2)], 0x1)
        truncated = dump.split("RELAX_BAKE_DATA 00000002")[0].rstrip()
        with self.assertRaises(ValueError):
            relax_bakes.parse_dump(truncated)

    def test_rejects_unterminated_block_after_a_complete_one(self):
        good, _ = make_dump("a", 8, [(1, 1, 1)], 0x1)
        bad, _ = make_dump("b", 8, [(2, 2, 2)], 0x2)
        head = bad.split("RELAX_BAKE_END")[0].rstrip()
        with self.assertRaises(ValueError):
            relax_bakes.parse_dump(good + "\n" + head)

    def test_rejects_unterminated_block_before_the_next_begin(self):
        a, _ = make_dump("a", 8, [(1, 1, 1)], 0x1)
        b, _ = make_dump("b", 8, [(2, 2, 2)], 0x2)
        head = a.split("RELAX_BAKE_END")[0].rstrip()
        with self.assertRaises(ValueError):
            relax_bakes.parse_dump(head + "\n" + b)

    def test_rejects_dump_starting_mid_block(self):
        dump, _ = make_dump("foo", 8, [(1, 1, 1)], 0x1)
        with self.assertRaises(ValueError):
            relax_bakes.parse_dump(dump.split("\n", 1)[1])


class EmitHeader(unittest.TestCase):
    def test_emits_named_struct_and_bits(self):
        dump, out = make_dump("foo_bar", 100, [(0x11, 0x22, 0x33)], 0xC0FFEE)
        bakes = relax_bakes.parse_dump(dump)
        header = relax_bakes.emit_header(bakes)
        self.assertIn("namespace Solids {\nnamespace RelaxBakes {", header)
        self.assertIn('#include "mesh/conway.h"', header)
        self.assertIn(
            "inline const uint32_t foo_bar_bits[] "
            "HS_PROGMEM_UNIQUE(foo_bar_bits) = {",
            header,
        )
        self.assertIn('.name = "foo_bar", .vertex_bits = foo_bar_bits,', header)
        self.assertIn(
            ".vertex_count = 1, .face_count = 2, .index_count = 6, "
            ".iterations = 100,",
            header,
        )
        self.assertIn(f".topology_hash = 0x{0xC0FFEE:08x}u,", header)
        self.assertIn(f".output_hash = 0x{out:08x}u}};", header)


if __name__ == "__main__":
    unittest.main()
