#!/usr/bin/env python3
"""Host tests for the PC-sampler dump symbolizer (tools/pc_profile.py)."""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import pc_profile as pc  # noqa: E402

HEADER = ("pc profile: version=1 hz=10000 bin_bytes=16 capacity=4096 "
          "defer_counts=48 counts_per_us=24 elapsed_ms=1000")
COUNTERS = ("pc profile: ticks=1000 idle=500 deferred=60 missed=5 handler=3 "
            "recorded=430 dropped=5 bins=6 service_avg=110 service_max=300")
LATENCY = "pc profile latency: 400 300 200 50 20 10 5 5 5 3 1 1"
BINS = {0x00001000: 50, 0x00001010: 30, 0x00001020: 40, 0x60002000: 200,
        0x60002040: 100, 0x20000100: 10}

NM = """\
00000000 A __absolute
00001000 00000020 T fast_fn()
00001028 00000010 t helper(int)
00001028 00000010 W folded_alias()
00001040 W no_size_fn
00001060 00000020 T last_itcm
20001000 00000004 B g_var
60002000 00000040 W Plot::PlanarChords<288, 144>::rasterize_run<Sink, Lambda>(Sink&)
60002041 00000010 T thumb_bit_symbol
6000cbc0 00000100 R some_rodata
"""


def dump_lines(bins=BINS, header=HEADER, counters=COUNTERS, end=True):
    lines = ["profile harness: effect=Fx config=ship", "f 1 w=100 r=90",
             header, counters, LATENCY]
    lines += [f"pc {address:08x} {count}" for address, count in bins.items()]
    if end:
        lines.append("pc profile end")
    lines.append("=== profile Fx [288x144] frames 1-32 window=100 us ===")
    return lines


def rows_by_name(rows):
    return {row.name: row for row in rows}


class ParseDump(unittest.TestCase):
    def test_reads_header_counters_latency_and_bins(self):
        dump = pc.parse_dump(dump_lines())
        self.assertEqual(dump.header["bin_bytes"], 16)
        self.assertEqual(dump.header["elapsed_ms"], 1000)
        self.assertEqual(dump.counters["recorded"], 430)
        self.assertEqual(dump.latency[:3], [400, 300, 200])
        self.assertEqual(dump.bins, BINS)
        self.assertEqual(dump.render, 500)

    def test_rejects_incomplete_or_inconsistent_dumps(self):
        cases = {
            "no pc profile dump": [line for line in dump_lines()
                                   if not line.startswith("pc ")],
            "no end line": dump_lines(end=False),
            "more than one": dump_lines() + dump_lines(),
            "bins sum to 400": dump_lines(bins={**BINS, 0x1000: 20}),
            "5 pc bins listed": dump_lines(
                bins={k: v for k, v in BINS.items() if k != 0x20000100},
                counters=COUNTERS.replace("recorded=430", "recorded=420")),
            "fields": dump_lines(counters=COUNTERS + " extra=1"),
            "unsupported pc profile version": dump_lines(
                header=HEADER.replace("version=1", "version=2")),
        }
        for message, lines in cases.items():
            with self.subTest(message=message):
                with self.assertRaisesRegex(pc.DumpError, message):
                    pc.parse_dump(lines)

    def test_rejects_a_duplicate_bin(self):
        lines = dump_lines()
        lines.insert(lines.index("pc profile end"), "pc 00001000 1")
        with self.assertRaisesRegex(pc.DumpError, "duplicate pc bin 00001000"):
            pc.parse_dump(lines)

    def test_elf_comes_from_provenance(self):
        lines = dump_lines() + [
            "profile provenance: artifact_profile_elf=build/prof/a/profile.elf"]
        self.assertEqual(pc.elf_from_log(lines), "build/prof/a/profile.elf")
        self.assertIsNone(pc.elf_from_log(dump_lines()))


class Symbols(unittest.TestCase):
    def test_parse_nm_keeps_text_symbols_and_clears_thumb_bit(self):
        symbols = pc.parse_nm(NM)
        self.assertEqual([s[3] for s in symbols], [
            "fast_fn()", "helper(int)", "folded_alias()", "no_size_fn",
            "last_itcm",
            "Plot::PlanarChords<288, 144>::rasterize_run<Sink, Lambda>(Sink&)",
            "thumb_bit_symbol"])
        self.assertEqual(symbols[3], (0x1040, None, "W", "no_size_fn"))
        self.assertEqual(symbols[6], (0x60002040, 0x10, "T", "thumb_bit_symbol"))

    def test_build_functions_folds_aliases_and_infers_sizes(self):
        functions = {f.address: f for f in pc.build_functions(pc.parse_nm(NM))}
        self.assertEqual(sorted(functions), [0x1000, 0x1028, 0x1040, 0x1060,
                                             0x60002000, 0x60002040])
        folded = functions[0x1028]
        self.assertEqual((folded.name, folded.aliases, folded.size),
                         ("folded_alias()", 1, 0x10))
        self.assertEqual(functions[0x1040].size, 0x20)
        self.assertEqual(functions[0x1000].aliases, 0)

    def test_regions_split_at_the_itcm_and_flexspi_bounds(self):
        for address, region in ((0x0, "ITCM"), (0x7FFFF, "ITCM"),
                                (0x80000, "other"), (0x20000000, "other"),
                                (0x5FFFFFFF, "other"), (0x60000000, "flash"),
                                (0x7FFFFFFF, "flash"), (0x80000000, "other")):
            with self.subTest(address=hex(address)):
                self.assertEqual(pc.region_of(address), region)


class Attribution(unittest.TestCase):
    def rows(self):
        dump = pc.parse_dump(dump_lines())
        functions = pc.build_functions(pc.parse_nm(NM))
        return dump, pc.attribute(dump, pc.Symbolizer(functions))

    def test_straddling_bins_split_by_bytes(self):
        _, rows = self.rows()
        named = rows_by_name(rows)
        self.assertEqual(named["fast_fn()"].samples, 80)
        self.assertEqual(named["folded_alias()"].samples, 20)
        self.assertEqual(named["<no symbol: ITCM>"].samples, 20)
        self.assertEqual(named["<no symbol: other>"].samples, 10)
        self.assertEqual(named["thumb_bit_symbol"].samples, 100)
        self.assertEqual(sum(row.samples for row in rows), 430)

    def test_rows_rank_by_samples_with_regions(self):
        _, rows = self.rows()
        self.assertEqual([(r.name[:22], r.region, r.samples) for r in rows[:3]], [
            ("Plot::PlanarChords<288", "flash", 200),
            ("thumb_bit_symbol", "flash", 100),
            ("fast_fn()", "ITCM", 80)])
        self.assertEqual(pc.region_totals(rows),
                         {"ITCM": 120, "flash": 300, "other": 10})

    def test_hot_ranges_merge_contiguous_bins(self):
        _, rows = self.rows()
        fast = rows_by_name(rows)["fast_fn()"]
        self.assertEqual(pc.hot_ranges(fast, 16, 4), [(0x1000, 0x1020, 80)])

    def test_report_states_shares_regions_and_ranges(self):
        dump, rows = self.rows()
        text = pc.report(dump, rows, top=8, ranges=1, range_rows=2,
                         name_width=30, lines={0x60002000: "chords.h:491"})
        self.assertIn("render samples 500 (50.0% of all): binned 430, "
                      "deferred 60, missed 5, dropped 5, handler-mode 3", text)
        self.assertIn("regions (% of render): ITCM 24.0%  flash 60.0%  "
                      "other 2.0%  deferred 13.0%  dropped 1.0%", text)
        self.assertIn("regions (% of binned): ITCM 27.9%  flash 69.8%  "
                      "other 2.3%", text)
        table = text[text.index("   pct  samples region   size  function") + 1:]
        self.assertEqual(table[:3], [
            " 40.0%      200 flash      64  Plot::PlanarChords<288, 144...",
            " 20.0%      100 flash      16  thumb_bit_symbol",
            " 16.0%       80 ITCM       32  fast_fn()"])
        self.assertIn(" 13.0%       65 -           0  "
                      "<deferred: higher-priority ...", text)
        self.assertIn("  1.0%        5 -           0  <dropped: histogram full>",
                      text)
        self.assertIn("    +0x0000..+0x000f      200  40.0%  chords.h:491", text)
        self.assertTrue(any(line.startswith("WARNING: 5 samples") for line in text))


class CommandLine(unittest.TestCase):
    def run_main(self, lines, *extra):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "capture.log"
            log.write_text("\n".join(lines) + "\n", encoding="utf-8")
            nm = Path(directory) / "nm.txt"
            nm.write_text(NM, encoding="utf-8")
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = pc.main([str(log), "--nm-file", str(nm), *extra])
            return code, out.getvalue(), err.getvalue()

    def test_symbolizes_a_capture_from_an_nm_listing(self):
        code, out, err = self.run_main(dump_lines(), "--top", "2")
        self.assertEqual(code, 0, err)
        table = out.split("hot ranges")[0]
        self.assertIn("thumb_bit_symbol", table)
        self.assertNotIn("fast_fn()", table)

    def test_a_capture_without_a_dump_fails(self):
        code, out, err = self.run_main(dump_lines()[:2])
        self.assertEqual(code, 1)
        self.assertIn("no pc profile dump", err)


if __name__ == "__main__":
    unittest.main()
