#!/usr/bin/env python3
"""Host tests for the hot-code cache layout planner and gate (tools/hot_layout.py).

Run:  python -m unittest discover -s tools/hot_layout_tests
"""

import os
import struct
import subprocess
import sys
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import hot_layout as hl  # noqa: E402

PHANTASM = hl.Origin(r".pio\build\phantasm\src\targets\Phantasm\Phantasm.ino.cpp.o")
LIBM = hl.Origin("C:/toolchain/hard/libm.a", "libm_a-sf_sin.o")
WAY = hl.CACHE_WAY_BYTES


MAP_TEXT = """\
Discarded input sections

 .text.hot._ZN4Cold4goneEv
                0x00000000       0x40 x.o

Linker script and memory map

.text.code      0x60001400    0x9000
 .startup       0x60001400      0x20c .pio\\build\\phantasm\\libFrameworkArduino.a(startup.c.o)
                0x600015e4                ResetHandler
                0x60008000                . = ALIGN (0x20)
                0x60008000                hs_hot_layout_start = .
 .text.hot._ZN8RingSpinILi288ELi144EE10draw_frameEv
                0x60008000      0x120 .pio\\build\\phantasm\\src\\targets\\Phantasm\\Phantasm.ino.cpp.o
                0x60008000                RingSpin<288, 144>::draw_frame()
 .text.sinf     0x60008120       0x80 C:/toolchain/hard\\libm.a(libm_a-sf_sin.o)
 .text.empty    0x600081a0        0x0 x.o
.text.itcm      0x00000000     0x1000 load address 0x601b4124
 .text._ZN6Canvas5flushEv
                0x00000020       0x40 .pio\\build\\phantasm\\src\\targets\\Phantasm\\Phantasm.ino.cpp.o
.text.progmem   0x600f3580       0x10
 .rodata.table  0x600f3580       0x10 x.o
"""


class ParseMap(unittest.TestCase):
    def setUp(self):
        self.sections, self.symbols = hl.parse_map(MAP_TEXT)

    def test_reads_code_sections_after_the_memory_map_marker(self):
        names = [(s.output, s.name) for s in self.sections]
        self.assertEqual(names, [
            (".text.code", ".startup"),
            (".text.code", ".text.hot._ZN8RingSpinILi288ELi144EE10draw_frameEv"),
            (".text.code", ".text.sinf"),
            (".text.itcm", ".text._ZN6Canvas5flushEv")])

    def test_wrapped_line_carries_address_size_and_object(self):
        hot = self.sections[1]
        self.assertEqual((hot.address, hot.size), (0x60008000, 0x120))
        self.assertEqual(hot.origin, hl.Origin(
            r".pio\build\phantasm\src\targets\Phantasm\Phantasm.ino.cpp.o"))

    def test_archive_member_splits_into_archive_and_member(self):
        sinf = self.sections[2]
        self.assertEqual(sinf.origin.member, "libm_a-sf_sin.o")
        self.assertEqual(sinf.origin.archive, "libm.a")

    def test_symbol_assignments(self):
        self.assertEqual(self.symbols["hs_hot_layout_start"], 0x60008000)

    def test_rejects_a_file_without_the_memory_map(self):
        with self.assertRaises(hl.LayoutError):
            hl.parse_map("Discarded input sections\n")


class ParseAlignments(unittest.TestCase):
    def test_object_and_archive_sections(self):
        obj = (
            "x.o:     file format elf32-littlearm\n\nSections:\n"
            "Idx Name          Size      VMA       LMA       File off  Algn  Flags\n"
            "  4 .text.hot._ZN1A1fEv 00000014  00000000  00000000  0000b7dc  2**2"
            "  CONTENTS, ALLOC, LOAD, READONLY, CODE\n"
            "  5 .text.hot._ZN1A1gEv 00000040  00000000  00000000  0000b7f0  2**5"
            "  CONTENTS, ALLOC, LOAD, READONLY, CODE\n")
        self.assertEqual(hl.parse_section_alignments(obj, "x.o"), {
            (hl.Origin("x.o"), ".text.hot._ZN1A1fEv"): 4,
            (hl.Origin("x.o"), ".text.hot._ZN1A1gEv"): 32})
        archive = (
            "In archive libm.a:\n\nlibm_a-sf_sin.o:     file format elf32-littlearm\n\n"
            "Sections:\nIdx Name Size VMA LMA File off Algn Flags\n"
            "  3 .text.sinf    00000190  00000000  00000000  00000034  2**3"
            "  CONTENTS, ALLOC, LOAD, RELOC, READONLY, CODE\n")
        self.assertEqual(hl.parse_section_alignments(archive, "libm.a"),
                         {(hl.Origin("libm.a", "libm_a-sf_sin.o"), ".text.sinf"): 8})


class ParseFunctions(unittest.TestCase):
    def test_clears_the_thumb_bit_and_keeps_the_largest_alias(self):
        text = (
            "   1: 60008001   288 FUNC    GLOBAL DEFAULT    2 _ZN1A4drawEv\n"
            "   2: 60008001     0 FUNC    GLOBAL DEFAULT    2 alias\n"
            "   3: 20000000     4 OBJECT  GLOBAL DEFAULT    3 data\n")
        self.assertEqual(hl.parse_functions(text),
                         [hl.Function("_ZN1A4drawEv", 0x60008000, 288)])


class ParseReferences(unittest.TestCase):
    def test_literals_require_exact_function_starts_but_branches_allow_interiors(self):
        text = ("60008000:\t.word\t0x00000105\n"
                "60008004:\t.word\t0x00000101\n"
                "60008008:\tbl\t104 <callee+0x4>\n")
        code, _ = hl.parse_references(text, {0x100})
        self.assertEqual(code, [(0x60008004, 0x100), (0x60008008, 0x104)])

    def test_branches_literals_and_move_pairs_split_code_from_data(self):
        text = (
            "60008000 <_ZN1A4drawEv>:\n"
            "60008000:\tbl\t60009000 <_ZN1B3runEv>\n"
            "60008004:\tbne.w\t60009100 <_ZN1C3runEv>\n"
            "60008008:\tcbz\tr3, 6000a000 <_ZN1D3runEv>\n"
            "6000800a:\tmovw\tr3, #45057\t@ 0xb001\n"
            "6000800e:\tmovt\tr3, #24576\t@ 0x6000\n"
            "60008012:\tldr\tr0, [pc, #4]\t@ 60008018 <_ZN1A4drawEv+0x18>\n"
            "60008018:\t.word\t0x6000c001\n"
            "6000801c:\t.word\t0x20000ec0\n")
        code, data = hl.parse_references(text)
        self.assertEqual(code, [
            (0x60008000, 0x60009000), (0x60008004, 0x60009100),
            (0x60008008, 0x6000a000), (0x6000800e, 0x6000b000),
            (0x60008018, 0x6000c000)])
        self.assertEqual(data, [(0x6000801c, 0x20000ec0)])


class DataObjects(unittest.TestCase):
    def test_read_words_takes_aligned_words_from_the_loaded_section(self):
        elf = bytearray(0xB0)
        struct.pack_into("<I", elf, 0x20, 0x60)
        struct.pack_into("<HH", elf, 0x2E, 40, 2)
        struct.pack_into("<4I", elf, 0x40, 1, 2, 3, 4)
        struct.pack_into("<IIIIII", elf, 0x60 + 40, 0, 1, 2, 0x60001000, 0x40, 16)
        self.assertEqual(hl.read_words(bytes(elf), [(0x60001004, 8), (0x60002000, 4)]),
                         {0x60001004: [2, 3]})

    def _callees(self, words):
        program = _program([
            ("_ZN4Beta4drawEv", 0x60010000, 0x40, ".text.hot._ZN4Beta4drawEv", PHANTASM),
            ("_ZN4Beta6targetEv", 0x60020000, 0x40, ".text.hot._ZN4Beta6targetEv", PHANTASM),
        ], [], data=[(0x60010004, 0x600f0000)], objects={0x600f0000: (4 * len(words), words)})
        return program.edges.get(0x60010000, set())

    def test_pointer_table_reaches_its_functions(self):
        self.assertEqual(self._callees([0x60020001, 0]), {0x60020000})

    def test_function_address_in_a_data_table_is_no_edge(self):
        self.assertEqual(self._callees([0x3F800000, 0x40000000, 0x12345678, 0x60020001]), set())


class EffectNames(unittest.TestCase):
    def test_mangled_identifiers_are_length_prefixed(self):
        found = hl.mangled_identifiers(
            "_ZN4Plot4Mesh16plot_steady_edgeILi288ELi144EZN13HopfFibrationILi288E")
        self.assertTrue({"Plot", "Mesh", "plot_steady_edge", "HopfFibration"} <= found)
        self.assertNotIn("HopfFibrationILi288E", found)

    def test_effect_list_macro_body(self):
        text = ("#define OTHER(X) X(Nope, 1)\n"
                "#define HS_PHANTASM_EFFECT_LIST(X)  \\\n"
                "  X(LiquidWindows, 60)              \\\n"
                "  X(DreamBalls, hs_preset_window_seconds<DreamBalls>())\n"
                "  X(After, 1)\n")
        self.assertEqual(hl.parse_effect_list(text, "HS_PHANTASM_EFFECT_LIST"),
                         ["LiquidWindows", "DreamBalls"])

    def test_missing_macro_raises(self):
        with self.assertRaises(hl.LayoutError):
            hl.parse_effect_list("", "HS_PHANTASM_EFFECT_LIST")


def _program(functions, references, data=(), objects=None):
    """Program with one input section per function, kind from its section name."""
    sections = []
    for name, address, size, section, origin in functions:
        output = hl.ITCM_OUTPUT if address < hl.ITCM_END else hl.FLASH_OUTPUT
        sections.append(hl.Section(output, section, origin, address, size))
    return hl.Program([hl.Function(f[0], f[1], f[2]) for f in functions],
                      sections, references, data, objects)


class HotSets(unittest.TestCase):
    def setUp(self):
        self.program = _program([
            ("_ZN5Alpha4drawEv", 0x60010000, 0x100, ".text.hot._ZN5Alpha4drawEv", PHANTASM),
            ("___ZN6Shared4walkEv_veneer", 0x60020000, 8, ".text.stub", PHANTASM),
            ("_ZN6Shared4walkEv", 0x00001000, 0x200, ".text._ZN6Shared4walkEv", PHANTASM),
            ("_ZN6Shared6kernelEv", 0x60030000, 0x400, ".text.hot._ZN6Shared6kernelEv", PHANTASM),
            ("_ZN6Shared4coldEv", 0x60040000, 0x80, ".text.unlikely._ZN6Shared4coldEv", PHANTASM),
            ("_ZN6Shared4deepEv", 0x60050000, 0x40, ".text.hot._ZN6Shared4deepEv", PHANTASM),
            ("_ZN6Shared5colderEv", 0x60058000, 0x40, ".text.unlikely._ZN6Shared5colderEv", PHANTASM),
            ("_ZN6Shared5tableEv", 0x60090000, 0x40, ".text.hot._ZN6Shared5tableEv", PHANTASM),
            ("_ZN5Alpha4initEv", 0x60060000, 0x40, ".text._ZN5Alpha4initEv", PHANTASM),
            ("_ZN5Setup4onlyEv", 0x60070000, 0x40, ".text.hot._ZN5Setup4onlyEv", PHANTASM),
            ("sinf", 0x60080000, 0x80, ".text.sinf", LIBM),
        ], [
            (0x60010004, 0x60020000),   # draw -> veneer
            (0x60020004, 0x00001000),   # veneer literal -> walk in ITCM
            (0x00001004, 0x60030000),   # walk -> kernel
            (0x60030004, 0x60040000),   # kernel -> cold
            (0x60040004, 0x60050000),   # cold -> deep
            (0x60040008, 0x60058000),   # cold -> colder
            (0x60030008, 0x60080000),   # kernel -> sinf
            (0x60060004, 0x60070000),   # init -> setup-only hot code
        ], data=[(0x6003000c, 0x600f0004)],   # kernel -> inside a pointer table
           objects={0x600f0000: (8, [0x20000000, 0x60090001])})

    def test_kinds(self):
        kinds = {self.program.by_address[a].name: k for a, k in self.program.kind.items()}
        self.assertEqual(kinds["_ZN5Alpha4drawEv"], hl.HOT)
        self.assertEqual(kinds["___ZN6Shared4walkEv_veneer"], hl.VENEER)
        self.assertEqual(kinds["_ZN6Shared4walkEv"], hl.ITCM)
        self.assertEqual(kinds["_ZN6Shared4coldEv"], hl.COLD)
        self.assertEqual(kinds["_ZN5Alpha4initEv"], hl.SETUP)
        self.assertEqual(kinds["sinf"], hl.LIBRARY)

    def test_closure_ranks_by_depth_then_call_order_and_skips_setup_code(self):
        hot = self.program.hot_set("Alpha", [])
        ranked = [(level, key[1]) for level, key in hot.ranked]
        self.assertEqual(ranked, [
            (0, ".text.hot._ZN5Alpha4drawEv"),
            (2, ".text.hot._ZN6Shared6kernelEv"),
            (3, ".text.unlikely._ZN6Shared4coldEv"),
            (3, ".text.sinf"),
            (3, ".text.hot._ZN6Shared5tableEv"),
            (4, ".text.hot._ZN6Shared4deepEv")])

    def test_gated_set_skips_members_that_would_overflow_the_budget(self):
        budget = hl.GATED_BYTES
        try:
            hl.GATED_BYTES = 0x100 + 0x400 + 0x7f
            hot = self.program.hot_set("Alpha", [])
        finally:
            hl.GATED_BYTES = budget
        self.assertEqual([key[1] for key in hot.gated], [
            ".text.hot._ZN5Alpha4drawEv", ".text.hot._ZN6Shared6kernelEv",
            ".text.hot._ZN6Shared5tableEv"])
        self.assertEqual([key[1] for key in hot.overflow], [
            ".text.unlikely._ZN6Shared4coldEv", ".text.sinf", ".text.hot._ZN6Shared4deepEv"])

    def test_cold_code_is_placed_only_when_hot_code_calls_it(self):
        names = {self.program.by_key[key].name for _, key in self.program.hot_set("Alpha", []).ranked}
        self.assertIn(".text.unlikely._ZN6Shared4coldEv", names)
        self.assertNotIn(".text.unlikely._ZN6Shared5colderEv", names)

    def test_effect_without_code_has_no_hot_set(self):
        self.assertEqual(self.program.hot_sets(["Alpha", "Missing"]).keys(), {"Alpha"})


class Gate(unittest.TestCase):
    def test_three_lines_in_one_set_fail_and_name_the_functions(self):
        program = _program([
            (f"_ZN4Beta{len(n)}{n}Ev", 0x60010000 + i * WAY, 64, f".text.hot._ZN4Beta{len(n)}{n}Ev",
             PHANTASM) for i, n in enumerate(("aa", "bb", "cc"))], [])
        hot_sets = program.hot_sets(["Beta"])
        report, = hl.evaluate(program, hot_sets)
        self.assertEqual((report.worst, report.over_sets, report.passed), (3, 2, False))
        self.assertEqual(report.offenders, ["_ZN4Beta2aaEv", "_ZN4Beta2bbEv", "_ZN4Beta2ccEv"])
        lines = hl.format_reports("phantasm", [report], [])
        self.assertIn("FAIL", lines[0])
        self.assertTrue(any("_ZN4Beta2ccEv" in line for line in lines))

    def test_two_lines_per_set_pass(self):
        program = _program([
            (f"_ZN4Beta{len(n)}{n}Ev", 0x60010000 + i * WAY, 64, f".text.hot._ZN4Beta{len(n)}{n}Ev",
             PHANTASM) for i, n in enumerate(("aa", "bb"))], [])
        report, = hl.evaluate(program, program.hot_sets(["Beta"]))
        self.assertTrue(report.passed)
        self.assertEqual(report.worst, 2)

    def test_missing_expected_effect_fails(self):
        lines = hl.format_reports("profile", [], ["HopfFibration"])
        self.assertIn("FAIL", lines[0])
        self.assertIn("HopfFibration", lines[1])


def _worst(program, keys):
    """Most distinct lines of `keys` mapping to one cache set."""
    per_set = {}
    for key in keys:
        section = program.by_key[key]
        for line in hl.lines_of(section.address, section.size):
            per_set.setdefault(line % hl.CACHE_SETS, set()).add(line)
    return max(len(lines) for lines in per_set.values())


class Planner(unittest.TestCase):
    def _relocated(self, program, order, start):
        at = start
        for key in order:
            at = (at + 3) // 4 * 4
            program.by_key[key] = replace(program.by_key[key], address=at)
            at += program.by_key[key].size
        return program

    def _shared_program(self):
        functions = []
        address = 0x60010000
        for effect, count in (("Alpha", 9), ("Beta", 9), ("Gamma", 5)):
            for index in range(count):
                name = f"_ZN{len(effect)}{effect}2k{index}Ev"
                functions.append((name, address, 1500, ".text.hot." + name, PHANTASM))
                address += WAY
        shared = "_ZN6Shared4coreEv"
        functions.append((shared, address, 4000, ".text.hot." + shared, PHANTASM))
        references = [(f[1] + 4, address) for f in functions[:-1]]
        return _program(functions, references)

    def test_shared_and_private_code_fit_two_lines_per_set(self):
        program = self._shared_program()
        hot_sets = program.hot_sets(["Alpha", "Beta", "Gamma"])
        self.assertGreater(max(_worst(program, h.gated) for h in hot_sets.values()), 2)
        alignments = {s.key: 4 for s in program.sections}
        order = hl.plan(program, hot_sets, alignments, 0x60008010)
        self.assertLessEqual({k for h in hot_sets.values() for k in h.gated}, set(order))
        relocated = self._relocated(program, order, 0x60008010)
        self.assertEqual({e: _worst(relocated, h.gated) <= 2 for e, h in hot_sets.items()},
                         {"Alpha": True, "Beta": True, "Gamma": True})

    def test_plan_is_deterministic(self):
        script = (
            "from test_hot_layout import Planner, hl\n"
            "program = Planner()._shared_program()\n"
            "hot_sets = program.hot_sets(['Alpha', 'Beta', 'Gamma'])\n"
            "alignments = {s.key: 4 for s in program.sections}\n"
            "print(hl.plan(program, hot_sets, alignments, 0x60008010))\n"
        )
        orders = []
        for seed in ("1", "2"):
            result = subprocess.run(
                [sys.executable, "-c", script], cwd=Path(__file__).parent,
                env=dict(os.environ, PYTHONHASHSEED=seed),
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            orders.append(result.stdout)
        self.assertEqual(orders[0], orders[1])

    def test_fragment_names_archive_members_and_checks_the_end(self):
        program = _program([
            ("_ZN5Alpha4drawEv", 0x60010000, 0x100, ".text.hot._ZN5Alpha4drawEv", PHANTASM),
            ("sinf", 0x60080000, 0x80, ".text.sinf", LIBM)], [])
        order = [s.key for s in program.sections]
        fragment = hl.render_fragment(program, order, 0x60008180)
        self.assertIn("\n*Phantasm.ino.cpp.o(.text.hot._ZN5Alpha4drawEv)\n", fragment)
        self.assertIn("\n*libm.a:libm_a-sf_sin.o(.text.sinf)\n", fragment)
        self.assertIn('ASSERT(ABSOLUTE(.) == 0x60008180, "hot code layout differs from its plan");',
                      fragment)
        self.assertEqual(hl.planned_end(program, order, {k: 4 for k in order}, 0x60008000),
                         0x60008180)


class CheckVerdict(unittest.TestCase):
    def test_missing_and_overbound_effects_fail_with_a_valid_companion(self):
        def program(count):
            return _program([(f"_ZN5Alpha2f{i}Ev", 0x60010000 + i * WAY, 32,
                              f".text.hot._ZN5Alpha2f{i}Ev", PHANTASM)
                             for i in range(count)], [])
        for image, expected, accepted in ((program(1), ["Alpha"], True),
                                           (program(1), ["Missing"], False),
                                           (program(3), ["Alpha"], False)):
            with self.subTest(expected=expected, accepted=accepted), \
                    mock.patch.object(hl, "load_program", return_value=(image, {})):
                lines, passed = hl.check(None, "firmware.elf", "firmware.map", ["Alpha"], expected, "fixture")
                self.assertEqual(passed, accepted)
                self.assertIn("PASS" if accepted else "FAIL", lines[0])


if __name__ == "__main__":
    unittest.main()
