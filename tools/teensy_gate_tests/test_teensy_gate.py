#!/usr/bin/env python3
"""Host tests for the Teensy 4 size/layout gate (tools/teensy_gate.py).

The gate must gate itself: a size/layout check that can never fail is worse than
none (permanent false-green). So every layout invariant and region ceiling is
proven to FAIL on a deliberately-broken fixture and PASS on the good one. These
are plain stdlib `unittest` tests — no ARM toolchain, no PlatformIO — so they run
in milliseconds in the existing CI test lane.

Run:  python -m unittest discover -s tools/teensy_gate_tests
"""

import contextlib
import copy
import dataclasses
import io
import json
import os
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import teensy_gate as tg          # noqa: E402
import teensy_warnings as tw      # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures"
REAL_DIR = FIX / "real"
BUDGETS = tg.load_budgets(TOOLS / "teensy_budgets.json")

# Verbatim toolchain output the synthetic fixtures cannot stand in for. All are
# committed, so a missing one is a deleted fixture, not an optional capture.
REAL_CAPTURES = (
    "cold_env_section.txt",
    "holosphere_readelf_secs.txt",
    "holosphere_readelf_syms.txt",
    "holosphere_size_a.txt",
    "holosphere_teensy_size.txt",
    "phantasm_readelf_secs.txt",
    "phantasm_readelf_syms.txt",
    "phantasm_size_a.txt",
    "phantasm_teensy_size.txt",
    "verbose_build_log.txt",
    "warm_env_section.txt",
)


def setUpModule():
    """Fail the suite when a real capture is gone, rather than skipping it."""
    missing = [name for name in REAL_CAPTURES
               if not (REAL_DIR / name).exists()]
    if missing:
        raise AssertionError(
            f"missing real captures in {REAL_DIR}: {', '.join(missing)}")


def _read(name):
    return (FIX / name).read_text(encoding="utf-8")


def _load_gate_extra():
    """Load tools/teensy_gate_extra.py (SCons glue) outside PlatformIO.

    The module runs `Import("env")` and `env.subst(...)` at import time, so it
    cannot be plain-imported; inject no-op SCons hooks and a stub env, then exec.
    """
    class _Env:
        def __init__(self):
            self.post_actions = []
            self.dependencies = []

        def subst(self, s):
            return str(TOOLS.parent) if s == "$PROJECT_DIR" else s

        def AddPostAction(self, target, action):
            self.post_actions.append((target, action))

        def Depends(self, target, dependency):
            self.dependencies.append((target, dependency))

    src = (TOOLS / "teensy_gate_extra.py").read_text(encoding="utf-8")
    mod = types.ModuleType("teensy_gate_extra")
    mod.__file__ = str(TOOLS / "teensy_gate_extra.py")
    mod.__dict__["Import"] = lambda *a, **k: None
    mod.__dict__["env"] = _Env()
    exec(compile(src, mod.__file__, "exec"), mod.__dict__)
    return mod


def _eval(env, teensy_size_file, syms_file, secs_file="good_readelf_secs.txt"):
    sizes = tg.parse_teensy_size(_read(teensy_size_file))
    symbols = tg.parse_readelf_symbols(_read(syms_file))
    sections = tg.parse_readelf_sections(_read(secs_file))
    return tg.evaluate(env, BUDGETS[env], sizes, symbols, sections)


def _codes(result):
    return sorted(v.code for v in result.violations)


def _size_a(itcm, dtcm, ocram, flash):
    return (
        "elf:\nsection size addr\n"
        f".text.itcm 0x{itcm:x} 0x0\n"
        f".bss 0x{dtcm:x} 0x20000000\n"
        f".bss.dma 0x{ocram:x} 0x20200000\n"
        f".text.progmem 0x{flash:x} 0x60000000\n")


def _invalid_size_a_cases():
    valid = _size_a(0x10000, 0x40000, 0x70000, 0x20000)
    return {
        "empty": "",
        "malformed": valid + ".broken not-a-size 0x60010000\n",
        "metadata-only": (
            "elf:\nsection size addr\n"
            ".ARM.attributes 0x30 0x0\n"
            ".comment 0x67 0x0\n"),
        "missing-flash": _size_a(0x10000, 0x40000, 0x70000, 0),
        "missing-itcm": _size_a(0, 0x40000, 0x70000, 0x20000),
        "missing-dtcm": _size_a(0x10000, 0, 0x70000, 0x20000),
        "missing-ocram": _size_a(0x10000, 0x40000, 0, 0x20000),
    }


#: teensy_size output whose component blobs carry no `name:bytes` pair (an `=`
#: separator variant). Every region is still reported, so a parser that summed
#: the pairs it found would report `used` 0 and clear every ceiling.
UNPARSEABLE_TEENSY_SIZE = (
    "teensy_size: Memory Usage on Teensy 4.0:\n"
    "teensy_size:   FLASH: code=158788, data=13684, headers=8460"
    "   free for files: 1883136\n"
    "teensy_size:   RAM1: variables=351280, code=62240, padding=30496"
    "   free for local variables: 68256\n"
    "teensy_size:   RAM2: variables=497920   free for malloc/new: 26368\n")


class TestAddressClassifier(unittest.TestCase):
    """The load-bearing replacement for nm — the easiest place for an off-by-one."""

    def test_region_buckets(self):
        self.assertEqual(tg.region_for_address(0x00000000), "ITCM")
        self.assertEqual(tg.region_for_address(0x20000040), "DTCM")
        self.assertEqual(tg.region_for_address(0x20200000), "OCRAM")
        self.assertEqual(tg.region_for_address(0x60010000), "FLASH")

    def test_half_open_boundaries(self):
        # lo inclusive, hi exclusive — guard the exact edges.
        self.assertEqual(tg.region_for_address(0x20000000), "DTCM")   # DTCM lo
        self.assertEqual(tg.region_for_address(0x2007FFFF), "DTCM")   # DTCM hi-1
        self.assertEqual(tg.region_for_address(0x20080000), "OTHER")  # past DTCM
        self.assertEqual(tg.region_for_address(0x201FFFFF), "OTHER")  # gap before OCRAM
        self.assertEqual(tg.region_for_address(0x20200000), "OCRAM")  # OCRAM lo
        self.assertEqual(tg.region_for_address(0x5FFFFFFF), "OTHER")  # just below flash
        self.assertEqual(tg.region_for_address(0x60000000), "FLASH")  # flash lo


class TestParsers(unittest.TestCase):
    def test_parse_teensy_size(self):
        sizes = tg.parse_teensy_size(_read("good_teensy_size.txt"))
        self.assertEqual(sizes["flash"]["used"], 158788 + 13684 + 8460)
        self.assertEqual(sizes["flash"]["free"], 1883136)
        self.assertEqual(sizes["ram1"]["used"], 351280 + 62240 + 30496)
        self.assertEqual(sizes["ram1"]["free"], 68256)
        self.assertEqual(sizes["ram2"]["used"], 497920)
        self.assertEqual(sizes["ram2"]["free"], 26368)

    def test_parse_teensy_size_single_space_before_free(self):
        # The "free for ..." separator width is not contractual; a single-space
        # variant must still parse all three regions, not silently yield
        # "region-missing". Internal single spaces in the component blob (", data:")
        # must NOT be mistaken for the separator.
        text = (
            "teensy_size: FLASH: code:62788, data:13684, headers:8460 free for files: 1979136\n"
            "teensy_size: RAM1: variables:305152, code:62240, padding:30496 free for local variables: 88512\n"
            "teensy_size: RAM2: variables:497920 free for malloc/new: 26368\n"
        )
        sizes = tg.parse_teensy_size(text)
        self.assertEqual(set(sizes), {"flash", "ram1", "ram2"})
        self.assertEqual(sizes["flash"]["used"], 62788 + 13684 + 8460)
        self.assertEqual(sizes["flash"]["free"], 1979136)
        self.assertEqual(sizes["ram1"]["used"], 305152 + 62240 + 30496)
        self.assertEqual(sizes["ram1"]["free"], 88512)
        self.assertEqual(sizes["ram2"]["used"], 497920)
        self.assertEqual(sizes["ram2"]["free"], 26368)

    def test_parse_teensy_size_rejects_an_unparseable_component_blob(self):
        # A blob that yields no pair sums to `used` 0, which is under every
        # ceiling: the parse must fail loud rather than hand the gate zeros.
        with self.assertRaises(tg.TeensySizeFormatError) as caught:
            tg.parse_teensy_size(UNPARSEABLE_TEENSY_SIZE)
        self.assertIn("FLASH", str(caught.exception))

    def test_parse_teensy_size_rejects_an_empty_component_blob(self):
        text = ("teensy_size:   FLASH: code:158788   free for files: 1883136\n"
                "teensy_size:   RAM2:    free for malloc/new: 26368\n")
        with self.assertRaises(tg.TeensySizeFormatError) as caught:
            tg.parse_teensy_size(text)
        self.assertIn("RAM2", str(caught.exception))

    def test_parse_teensy_size_rejects_multiple_environment_blocks(self):
        text = (_read("broken_over_cap_teensy_size.txt") + "\n" +
                _read("good_teensy_size.txt"))
        with self.assertRaises(tg.TeensySizeFormatError) as caught:
            tg.parse_teensy_size(text)
        self.assertIn("FLASH", str(caught.exception))
        self.assertIn("more than once", str(caught.exception))

    def test_parse_readelf_symbols_uses_real_mangled_names(self):
        syms = {s.name: s for s in tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))}
        self.assertIn("_ZN6Effect8buffer_aE", syms)              # class static, mangled
        self.assertIn("_ZL18global_arena_block", syms)           # file-scope static -> _ZL mangled LOCAL
        self.assertEqual(syms["_ZL18global_arena_block"].bind, "LOCAL")
        self.assertEqual(syms["_ZN6Effect8buffer_aE"].value, 0x20200000)
        self.assertEqual(syms["_ZN6Effect8buffer_aE"].region, "OCRAM")
        self.assertEqual(syms["_ZL18global_arena_block"].size, 305152)

class TestGoodBuildPasses(unittest.TestCase):
    def test_holosphere_good_build_clean(self):
        result = _eval("holosphere", "good_teensy_size.txt", "good_readelf_syms.txt")
        self.assertTrue(result.passed, msg=_codes(result))

    def test_phantasm_shares_layout_and_passes(self):
        result = _eval("phantasm", "good_teensy_size.txt", "good_readelf_syms.txt")
        self.assertTrue(result.passed, msg=_codes(result))


class TestLayoutInvariantsFail(unittest.TestCase):
    """One deliberately-broken fixture per invariant; each must turn the gate red."""

    def test_framebuffer_dropped_dmamem_lands_in_dtcm(self):
        result = _eval("holosphere", "good_teensy_size.txt",
                       "broken_framebuffer_dtcm_syms.txt")
        self.assertFalse(result.passed)
        self.assertIn("symbol-wrong-region", _codes(result))
        self.assertTrue(any("framebuffer" in v.message for v in result.violations))

    def test_reaction_graph_dropped_const_lands_in_ram(self):
        # phantasm env: it instantiates reaction-diffusion effects so the K-NN
        # table is present + budgeted (Holosphere GC's it out).
        result = _eval("phantasm", "good_teensy_size.txt",
                       "broken_reaction_graph_ram_syms.txt")
        self.assertFalse(result.passed)
        self.assertIn("symbol-wrong-region", _codes(result))
        self.assertTrue(any("reaction_graph" in v.message for v in result.violations))

    def test_dma_tx_buffer_dropped_dmamem_lands_in_dtcm(self):
        # phantasm env: the segment LED controller's eDMA TX buffers must stay in
        # OCRAM; a vague-linkage DMAMEM drop strands them in DTCM.
        result = _eval("phantasm", "good_teensy_size.txt",
                       "broken_dma_tx_dtcm_syms.txt")
        self.assertFalse(result.passed)
        self.assertIn("symbol-wrong-region", _codes(result))
        self.assertTrue(any("dma_tx_buffer" in v.message for v in result.violations))

    def test_arena_8mb_test_build_leak_trips_magnitude(self):
        result = _eval("holosphere", "good_teensy_size.txt",
                       "broken_arena_8mb_syms.txt")
        self.assertFalse(result.passed)
        self.assertIn("symbol-too-large", _codes(result))

    def test_renamed_symbol_fails_loud_not_silent(self):
        result = _eval("holosphere", "good_teensy_size.txt",
                       "broken_missing_symbol_syms.txt")
        self.assertFalse(result.passed)
        self.assertIn("symbol-not-found", _codes(result))

    def test_und_reference_row_does_not_trip_magnitude_or_region(self):
        # A same-named UND reference row (size 0, ndx UND, null address) alongside
        # the real definition must NOT spuriously fire symbol-too-small (0 < floor)
        # or symbol-wrong-region (addr 0 -> ITCM, not DTCM). Only defined rows count.
        sizes = tg.parse_teensy_size(_read("good_teensy_size.txt"))
        sections = tg.parse_readelf_sections(_read("good_readelf_secs.txt"))
        symbols = tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))
        symbols.append(tg.Symbol(
            num=9999, value=0, size=0, type="NOTYPE", bind="GLOBAL",
            ndx="UND", name="_ZL18global_arena_block"))
        result = tg.evaluate("holosphere", BUDGETS["holosphere"], sizes,
                             symbols, sections)
        self.assertTrue(result.passed, msg=_codes(result))
        self.assertNotIn("symbol-too-small", _codes(result))
        self.assertNotIn("symbol-wrong-region", _codes(result))


class TestRegionCeilingsFail(unittest.TestCase):
    def test_negative_free_reports_headroom_violation(self):
        sizes = tg.parse_teensy_size(_read("broken_negative_free_teensy_size.txt"))
        self.assertEqual(sizes["ram1"]["free"], -1024)

        symbols = tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))
        sections = tg.parse_readelf_sections(_read("good_readelf_secs.txt"))
        result = tg.evaluate("holosphere", BUDGETS["holosphere"], sizes,
                             symbols, sections)
        self.assertIn("headroom-below-floor", _codes(result))
        self.assertNotIn("region-missing", _codes(result))
        self.assertTrue(any("RAM1 free-for-local-variables -1,024 B" in v.message
                            for v in result.violations))

    def test_over_cap_trips_every_region(self):
        # Good symbols (layout fine) + over-cap totals: only region checks fire.
        result = _eval("holosphere", "broken_over_cap_teensy_size.txt",
                       "good_readelf_syms.txt")
        self.assertFalse(result.passed)
        codes = _codes(result)
        self.assertEqual(codes.count("region-over-budget"), 2)   # FLASH + RAM2
        self.assertIn("headroom-below-floor", codes)             # DTCM stack room


class TestEmptyBudgetFails(unittest.TestCase):
    """A budget with nothing to check must FAIL, and the shipped budgets must
    never be one: both loops in evaluate() are vacuous without their key."""

    def test_wholly_empty_budget_fails(self):
        result = tg.evaluate("x", {}, {}, [])
        self.assertFalse(result.passed)
        self.assertEqual(_codes(result), ["budget-empty"])

    def test_singular_key_typo_fails(self):
        budget = {"region": {"ram1": {"max_bytes": 1}},
                  "symbol": {"arena": {"name": "_ZL18global_arena_block"}}}
        result = tg.evaluate("x", budget, {}, [])
        self.assertEqual(_codes(result), ["budget-empty"])

    def test_either_key_alone_is_enough(self):
        regions_only = tg.evaluate(
            "x", {"regions": {"ram1": {"max_bytes": 1 << 20}}},
            {"ram1": {"used": 0, "free": 0}}, [])
        self.assertTrue(regions_only.passed, msg=_codes(regions_only))
        symbols_only = tg.evaluate(
            "x", {"symbols": {"arena": {"name": "a"}}}, {},
            [tg.Symbol(1, 0x20000000, 8, "OBJECT", "GLOBAL", "3", "a")])
        self.assertTrue(symbols_only.passed, msg=_codes(symbols_only))

    def test_every_shipped_budget_declares_regions_and_symbols(self):
        for env, budget in BUDGETS.items():
            with self.subTest(env=env):
                self.assertTrue(budget.get("regions"))
                self.assertTrue(budget.get("symbols"))
                self.assertNotIn("budget-empty",
                                 _codes(tg.evaluate(env, budget, {}, [])))


class TestBudgetSchema(unittest.TestCase):
    """Every ceiling is an optional `.get()`, so a misspelled key removes its
    check and the gate reports PASS. load_budgets() must reject the file first."""

    def _load(self, budgets):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "budgets.json"
            path.write_text(json.dumps(budgets), encoding="utf-8")
            return tg.load_budgets(path)

    def test_shipped_budgets_validate_unchanged(self):
        self.assertEqual(tg.validate_budgets(copy.deepcopy(BUDGETS)), BUDGETS)

    def _assert_typo_disables(self, env, path, real, wrong, code, sizes, symbols):
        """The intact budget rejects this build; the typo'd one PASSes; the
        schema now rejects the typo'd budgets file."""
        intact = copy.deepcopy(BUDGETS)
        self.assertEqual(_codes(tg.evaluate(env, intact[env], sizes, symbols)),
                         [code])
        typoed = copy.deepcopy(BUDGETS)
        spec = typoed[env]
        for key in path:
            spec = spec[key]
        spec[wrong] = spec.pop(real)
        self.assertTrue(tg.evaluate(env, typoed[env], sizes, symbols).passed)
        with self.assertRaises(tg.BudgetSchemaError) as ctx:
            self._load(typoed)
        self.assertIn(wrong, str(ctx.exception))

    def test_components_typo_drops_the_itcm_bank_ceiling(self):
        # variables 300,000 + the 12,288 floor -> 10 DTCM banks -> 196,608 B
        # ITCM ceiling; code 200,000 is over it and nothing regional fires.
        text = _read("good_teensy_size.txt").replace(
            "RAM1: variables:351280, code:62240, padding:30496"
            "   free for local variables: 68256",
            "RAM1: variables:300000, code:200000, padding:1000"
            "   free for local variables: 40288")
        self._assert_typo_disables(
            "phantasm", ("regions", "ram1"), "components", "component",
            "component-over-derived-ceiling", tg.parse_teensy_size(text),
            tg.parse_readelf_symbols(_read("good_readelf_syms.txt")))

    def test_max_bytes_typo_drops_the_flash_ceiling(self):
        text = _read("good_teensy_size.txt").replace(
            "FLASH: code:158788", "FLASH: code:1958788")
        self._assert_typo_disables(
            "phantasm", ("regions", "flash"), "max_bytes", "max_byte",
            "region-over-budget", tg.parse_teensy_size(text),
            tg.parse_readelf_symbols(_read("good_readelf_syms.txt")))

    def test_region_typo_drops_the_ocram_placement_invariant(self):
        # Only framebuffer_a loses DMAMEM, so exactly one invariant is in play.
        symbols = [dataclasses.replace(s, value=0x20060000)
                   if s.name == "_ZN6Effect8buffer_aE" else s
                   for s in tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))]
        self._assert_typo_disables(
            "holosphere", ("symbols", "framebuffer_a"), "region", "regoin",
            "symbol-wrong-region",
            tg.parse_teensy_size(_read("good_teensy_size.txt")), symbols)

    def _assert_deletion_disables(self, env, path, key, code, sizes, symbols):
        """The intact budget rejects this build; the one with `key` deleted
        PASSes; the schema now rejects the stripped budgets file."""
        intact = copy.deepcopy(BUDGETS)
        self.assertEqual(_codes(tg.evaluate(env, intact[env], sizes, symbols)),
                         [code])
        stripped = copy.deepcopy(BUDGETS)
        spec = stripped[env]
        for step in path:
            spec = spec[step]
        del spec[key]
        self.assertTrue(tg.evaluate(env, stripped[env], sizes, symbols).passed)
        with self.assertRaises(tg.BudgetSchemaError) as ctx:
            self._load(stripped)
        self.assertIn(key, str(ctx.exception))

    def test_deleting_max_bytes_drops_the_flash_ceiling(self):
        text = _read("good_teensy_size.txt").replace(
            "FLASH: code:158788", "FLASH: code:1958788")
        self._assert_deletion_disables(
            "phantasm", ("regions", "flash"), "max_bytes",
            "region-over-budget", tg.parse_teensy_size(text),
            tg.parse_readelf_symbols(_read("good_readelf_syms.txt")))

    def test_deleting_region_drops_the_ocram_placement_invariant(self):
        # Only framebuffer_a loses DMAMEM, so exactly one invariant is in play.
        symbols = [dataclasses.replace(s, value=0x20060000)
                   if s.name == "_ZN6Effect8buffer_aE" else s
                   for s in tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))]
        self._assert_deletion_disables(
            "holosphere", ("symbols", "framebuffer_a"), "region",
            "symbol-wrong-region",
            tg.parse_teensy_size(_read("good_teensy_size.txt")), symbols)

    def test_gutted_region_and_symbol_objects_are_rejected(self):
        for path in (("phantasm", "regions", "flash"),
                     ("phantasm", "regions", "ram1"),
                     ("phantasm", "regions", "ram2"),
                     ("holosphere", "regions", "flash"),
                     ("holosphere", "regions", "ram1"),
                     ("holosphere", "regions", "ram2"),
                     ("phantasm", "symbols", "arena"),
                     ("phantasm", "symbols", "reaction_graph"),
                     ("phantasm", "symbols", "dma_tx_buffer"),
                     ("holosphere", "symbols", "framebuffer_a"),
                     ("holosphere", "symbols", "framebuffer_b")):
            with self.subTest(path=path):
                budgets = copy.deepcopy(BUDGETS)
                spec = budgets
                for step in path[:-1]:
                    spec = spec[step]
                spec[path[-1]] = {}
                with self.assertRaises(tg.BudgetSchemaError):
                    self._load(budgets)

    def test_gutted_component_object_is_rejected(self):
        # Both component ceilings are optional individually; declaring neither
        # retires the ITCM bank-boundary ratchet with nothing else observing it.
        text = _read("good_teensy_size.txt").replace(
            "RAM1: variables:351280, code:62240, padding:30496"
            "   free for local variables: 68256",
            "RAM1: variables:300000, code:200000, padding:1000"
            "   free for local variables: 40288")
        sizes = tg.parse_teensy_size(text)
        symbols = tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))
        intact = copy.deepcopy(BUDGETS)
        self.assertEqual(
            _codes(tg.evaluate("phantasm", intact["phantasm"], sizes, symbols)),
            ["component-over-derived-ceiling"])
        gutted = copy.deepcopy(BUDGETS)
        gutted["phantasm"]["regions"]["ram1"]["components"]["code"] = {}
        self.assertTrue(
            tg.evaluate("phantasm", gutted["phantasm"], sizes, symbols).passed)
        with self.assertRaises(tg.BudgetSchemaError) as ctx:
            self._load(gutted)
        self.assertIn("max_banks_from_stack_floor", str(ctx.exception))

    def test_unknown_key_rejected_at_every_level(self):
        for path in (("phantasm",),
                     ("phantasm", "regions", "ram1"),
                     ("phantasm", "regions", "ram1", "components", "code"),
                     ("phantasm", "regions", "ram1", "components", "code",
                      "max_banks_from_stack_floor"),
                     ("phantasm", "symbols", "arena")):
            with self.subTest(level=path[-1]):
                budgets = copy.deepcopy(BUDGETS)
                spec = budgets
                for key in path:
                    spec = spec[key]
                spec["bogus"] = 1
                with self.assertRaises(tg.BudgetSchemaError):
                    self._load(budgets)

    def test_required_keys_enforced(self):
        no_name = copy.deepcopy(BUDGETS)
        del no_name["phantasm"]["symbols"]["arena"]["name"]
        with self.assertRaises(tg.BudgetSchemaError):
            self._load(no_name)
        no_bank = copy.deepcopy(BUDGETS)
        del (no_bank["phantasm"]["regions"]["ram1"]["components"]["code"]
             ["max_banks_from_stack_floor"]["bank_bytes"])
        with self.assertRaises(tg.BudgetSchemaError):
            self._load(no_bank)
        no_headroom = copy.deepcopy(BUDGETS)
        del (no_headroom["phantasm"]["regions"]["ram1"]["components"]["code"]
             ["max_banks_from_stack_floor"]["min_headroom_bytes"])
        with self.assertRaises(tg.BudgetSchemaError):
            self._load(no_headroom)

    def test_ram2_free_floor_is_required(self):
        # RAM2's max_bytes is the OCRAM size itself, so the free floor is the
        # region's only reachable constraint.
        for env in BUDGETS:
            with self.subTest(env=env):
                stripped = copy.deepcopy(BUDGETS)
                del stripped[env]["regions"]["ram2"]["free_min_bytes"]
                with self.assertRaises(tg.BudgetSchemaError) as ctx:
                    self._load(stripped)
                self.assertIn("free_min_bytes", str(ctx.exception))

    def test_boundary_headroom_must_be_a_non_negative_integer(self):
        for value in (-1, 1.5, True):
            with self.subTest(value=value):
                budgets = copy.deepcopy(BUDGETS)
                derived = (budgets["phantasm"]["regions"]["ram1"]
                           ["components"]["code"]
                           ["max_banks_from_stack_floor"])
                derived["min_headroom_bytes"] = value
                with self.assertRaises(tg.BudgetSchemaError):
                    self._load(budgets)

    def test_deleting_any_shipped_region_is_rejected(self):
        for env, budget in BUDGETS.items():
            for region in budget["regions"]:
                with self.subTest(env=env, region=region):
                    dropped = copy.deepcopy(BUDGETS)
                    del dropped[env]["regions"][region]
                    # The truncated budget still evaluates clean: nothing but
                    # the schema notices that the ceiling is gone.
                    self.assertTrue(tg.evaluate(
                        env, dropped[env],
                        tg.parse_teensy_size(_read("good_teensy_size.txt")),
                        tg.parse_readelf_symbols(
                            _read("good_readelf_syms.txt"))).passed)
                    with self.assertRaises(tg.BudgetSchemaError) as ctx:
                        self._load(dropped)
                    self.assertIn(region, str(ctx.exception))

    def test_deleting_any_shipped_layout_symbol_is_rejected(self):
        for env, budget in BUDGETS.items():
            for key in budget["symbols"]:
                with self.subTest(env=env, symbol=key):
                    dropped = copy.deepcopy(BUDGETS)
                    del dropped[env]["symbols"][key]
                    with self.assertRaises(tg.BudgetSchemaError) as ctx:
                        self._load(dropped)
                    self.assertIn(key, str(ctx.exception))

    def test_non_object_spec_rejected(self):
        budgets = copy.deepcopy(BUDGETS)
        budgets["phantasm"]["regions"] = [1, 2]
        with self.assertRaises(tg.BudgetSchemaError):
            self._load(budgets)

    def test_cli_reports_schema_error_as_cannot_run(self):
        budgets = copy.deepcopy(BUDGETS)
        budgets["phantasm"]["regions"]["ram1"]["max_byte"] = 1
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "budgets.json"
            path.write_text(json.dumps(budgets), encoding="utf-8")
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                code = tg.main(["--env", "phantasm", "--budgets", str(path),
                                "--teensy-size", str(FIX / "good_teensy_size.txt"),
                                "--readelf-syms", str(FIX / "good_readelf_syms.txt")])
        self.assertEqual(code, 2)
        self.assertIn("max_byte", err.getvalue())


class TestComponentCeilings(unittest.TestCase):
    """Per-component ceilings: the static max_bytes form, plus the shared
    fail-loud rules (missing component, size -A fallback)."""

    _STATIC_BUDGET = {"regions": {"ram1": {
        "components": {"code": {"max_bytes": 100000}}}}}

    def test_good_build_passes_component_ceiling(self):
        result = _eval("phantasm", "good_teensy_size.txt", "good_readelf_syms.txt")
        self.assertTrue(result.passed, msg=_codes(result))

    def test_code_component_over_static_ceiling_fails(self):
        # Static max_bytes form: only the component fires, nothing regional.
        sizes = tg.parse_teensy_size(_read("broken_component_over_teensy_size.txt"))
        result = tg.evaluate("phantasm", self._STATIC_BUDGET, sizes, [], {})
        self.assertFalse(result.passed)
        self.assertEqual(_codes(result), ["component-over-budget"])
        self.assertTrue(any("'code'" in v.message for v in result.violations))

    def test_code_component_under_static_ceiling_passes(self):
        sizes = tg.parse_teensy_size(_read("good_teensy_size.txt"))
        result = tg.evaluate("phantasm", self._STATIC_BUDGET, sizes, [], {})
        self.assertTrue(result.passed, msg=_codes(result))

    def test_missing_code_component_fails_loud_not_silent(self):
        result = _eval("phantasm", "broken_component_missing_teensy_size.txt",
                       "good_readelf_syms.txt")
        self.assertFalse(result.passed)
        self.assertIn("component-missing", _codes(result))

    def test_component_ceiling_is_opt_in_per_target(self):
        # Holosphere configures no components key; the same over-component
        # figures pass its (region-only) budget untouched.
        result = _eval("holosphere", "broken_component_over_teensy_size.txt",
                       "good_readelf_syms.txt")
        self.assertTrue(result.passed, msg=_codes(result))

    def test_declares_components_is_per_target(self):
        self.assertTrue(tg.declares_components(BUDGETS["phantasm"]))
        self.assertFalse(tg.declares_components(BUDGETS["holosphere"]))

    def test_size_a_fallback_reports_component_missing(self):
        # The `size -A` fallback synthesizes region totals without components,
        # so a component-ceiling target must fail loud there, not false-green.
        sizes = tg.fallback_sizes_from_size_a(_size_a(0x10000, 0x40000, 0x70000,
                                                      0x20000))
        symbols = tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))
        result = tg.evaluate("phantasm", BUDGETS["phantasm"], sizes, symbols, {})
        self.assertIn("component-missing", _codes(result))

    def test_only_teensy_size_reports_region_components(self):
        # The one shape difference between the two size producers, which every
        # per-component ceiling and the uncalibrated-pass exit code turn on.
        parsed = tg.parse_teensy_size(_read("good_teensy_size.txt"))
        fallback = tg.fallback_sizes_from_size_a(_read("good_size_a.txt"))
        self.assertTrue(parsed and fallback)
        self.assertTrue(all("components" in r for r in parsed.values()))
        self.assertFalse(any("components" in r for r in fallback.values()))


class TestDerivedComponentCeiling(unittest.TestCase):
    """The stack-floor-derived ITCM code ceiling: DTCM reserves
    ceil((variables + free_min_bytes) / bank) FlexRAM banks and code may fill
    the remaining banks. phantasm's ram1.code uses this form."""

    FLOOR = 12288      # phantasm ram1.free_min_bytes
    BANK = 32768

    @staticmethod
    def _ts(variables, code, free):
        return (
            "teensy_size: Memory Usage on Teensy 4.0:\n"
            "teensy_size:   FLASH: code:158788, data:13684, headers:8460"
            "   free for files: 1883136\n"
            f"teensy_size:   RAM1: variables:{variables}, code:{code}, "
            f"padding:1000   free for local variables: {free}\n"
            "teensy_size:   RAM2: variables:497920   free for malloc/new: 26368\n")

    @staticmethod
    def _budget():
        # Region + component checks only; the layout symbols have their own
        # test class (TestLayoutInvariantsFail).
        budget = copy.deepcopy(BUDGETS["phantasm"])
        budget.pop("symbols", None)
        (budget["regions"]["ram1"]["components"]["code"]
         ["max_banks_from_stack_floor"]["min_headroom_bytes"]) = 0
        return budget

    def _eval_ts(self, text, budget=None):
        sizes = tg.parse_teensy_size(text)
        return tg.evaluate("phantasm", budget or self._budget(), sizes, [], {})

    def test_ceiling_derivation(self):
        # variables 312,704 + floor 12,288 = 324,992 -> 10 DTCM banks ->
        # 6 ITCM banks -> 196,608 B ceiling; code well under -> pass.
        result = self._eval_ts(self._ts(312704, 150200, 14976))
        self.assertTrue(result.passed, msg=_codes(result))
        self.assertTrue(any("196,608" in n for n in result.notes))

    def test_code_over_derived_ceiling_fails_naming_the_floor(self):
        # variables 300,000 + floor -> still 10 DTCM banks (ceiling 196,608);
        # region total stays under its cap so only the derived check fires.
        result = self._eval_ts(self._ts(300000, 200000, 27680))
        self.assertFalse(result.passed)
        self.assertIn("component-over-derived-ceiling", _codes(result))
        self.assertNotIn("region-over-budget", _codes(result))
        msg = next(v.message for v in result.violations
                   if v.code == "component-over-derived-ceiling")
        self.assertIn("196,608", msg)                      # the derived ceiling
        self.assertIn("3,392", msg)                        # the overage
        self.assertIn("stack floor", msg)                  # the binding constraint
        self.assertIn("12,288", msg)

    def test_code_exactly_at_ceiling_passes_one_over_fails(self):
        at = self._eval_ts(self._ts(312704, 196608, 14976))
        self.assertTrue(at.passed, msg=_codes(at))
        over = self._eval_ts(self._ts(312704, 196609, 14976))
        self.assertIn("component-over-derived-ceiling", _codes(over))

    def test_shipping_budget_reserves_boundary_headroom(self):
        budget = copy.deepcopy(BUDGETS["phantasm"])
        budget.pop("symbols", None)
        derived = (budget["regions"]["ram1"]["components"]["code"]
                   ["max_banks_from_stack_floor"])
        self.assertEqual(derived["min_headroom_bytes"], 3072)
        at = self._eval_ts(self._ts(312704, 193536, 14976), budget)
        self.assertTrue(at.passed, msg=_codes(at))
        short = self._eval_ts(self._ts(312704, 193537, 14976), budget)
        self.assertIn("component-over-derived-ceiling", _codes(short))
        message = next(v.message for v in short.violations
                       if v.code == "component-over-derived-ceiling")
        self.assertIn("3,072 B of boundary headroom", message)

    def test_variables_crossing_a_bank_shrinks_the_ceiling(self):
        # variables 315,392 + 12,288 = 327,680 = exactly 10 banks -> ceiling
        # stays 196,608; one more variable byte forces an 11th DTCM bank and the
        # same code now violates the 163,840 B ceiling.
        ok = self._eval_ts(self._ts(315392, 180000, 12288))
        self.assertTrue(ok.passed, msg=_codes(ok))
        squeezed = self._eval_ts(self._ts(315393, 180000, 12287 + self.BANK))
        codes = _codes(squeezed)
        self.assertIn("component-over-derived-ceiling", codes)
        msg = next(v.message for v in squeezed.violations
                   if v.code == "component-over-derived-ceiling")
        self.assertIn("163,840", msg)

    def test_missing_variables_component_fails_loud_not_silent(self):
        text = (
            "teensy_size: Memory Usage on Teensy 4.0:\n"
            "teensy_size:   FLASH: code:158788, data:13684, headers:8460"
            "   free for files: 1883136\n"
            "teensy_size:   RAM1: code:150200, padding:1000"
            "   free for local variables: 14976\n"
            "teensy_size:   RAM2: variables:497920   free for malloc/new: 26368\n")
        result = self._eval_ts(text)
        self.assertIn("component-missing", _codes(result))
        self.assertTrue(any("'variables'" in v.message
                            for v in result.violations))

    def test_missing_free_min_bytes_fails_loud_not_silent(self):
        budget = self._budget()
        del budget["regions"]["ram1"]["free_min_bytes"]
        result = self._eval_ts(self._ts(312704, 150200, 14976), budget)
        self.assertIn("component-floor-missing", _codes(result))

    def test_informational_note_reports_growth_headroom(self):
        # Measured, ceiling, remaining, and next-bank-boundary distance must all
        # appear in the report so intra-bank growth is visible without failing.
        result = self._eval_ts(self._ts(312704, 150200, 14976))
        note = next(n for n in result.notes if "derived ceiling" in n)
        self.assertIn("150,200", note)                     # measured
        self.assertIn("196,608", note)                     # ceiling
        self.assertIn("46,408", note)                      # remaining
        self.assertIn("13,640", note)                      # to next 32 KiB boundary
        self.assertIn(note, tg.render_report(result))

    def test_note_present_even_when_over(self):
        result = self._eval_ts(self._ts(300000, 200000, 27680))
        self.assertTrue(any("derived ceiling" in n for n in result.notes))


class TestWarningRatchet(unittest.TestCase):
    def test_normalize_strips_line_and_col_keeps_identity(self):
        a = tw.normalize("core/effects/Foo.h:120:7: warning: unused variable 'x' [-Wunused-variable]")
        b = tw.normalize("core/effects/Foo.h:998:3: warning: unused variable 'x' [-Wunused-variable]")
        self.assertEqual(a, b)  # line/col stripped -> stable across unrelated edits
        self.assertIn("unused variable 'x'", a)

    def test_absolute_path_relativized_to_first_party(self):
        got = tw.normalize(
            "/home/runner/work/Holosphere/Holosphere/hardware/dma_led.h:42:5: "
            "warning: comparison is always true [-Wtype-limits]")
        self.assertTrue(got.startswith("hardware/dma_led.h:"))

    def test_fileless_diagnostics_are_normalized(self):
        lines = [
            '<command-line>: warning: "HS_PROFILE" redefined',
            "cc1plus: warning: command-line option '-Wmissing-prototypes' is valid for C",
            "ld: warning: firmware.elf has a LOAD segment with RWX permissions",
        ]
        self.assertEqual(tw.extract_warnings("\n".join(lines)), {
            '<command-line>: warning: "HS_PROFILE" redefined',
            "cc1plus: warning: command-line option '-Wmissing-prototypes' is valid for C",
            "ld: warning: firmware.elf has a LOAD segment with RWX permissions",
        })

    def test_library_warning_excluded(self):
        self.assertIsNone(tw.normalize(
            "/root/.platformio/lib/FastLED/FastLED.h:9:1: warning: foo [-Wbar]"))

    def test_library_path_with_first_party_segment_excluded(self):
        # A vendored lib whose nested dir reuses a first-party name (effects/)
        # must NOT collapse to a first-party key and pollute the baseline.
        self.assertIsNone(tw.normalize(
            "/root/.platformio/lib/SomeLib/effects/reverb.h:5:1: warning: w [-Wx]"))
        self.assertIsNone(tw.normalize(
            "/x/.pio/libdeps/teensy40/Foo/lib/core/bar.h:2:1: warning: w [-Wy]"))
        # No `lib/` segment: libdeps/ itself has to carry the exclusion.
        self.assertIsNone(tw.normalize(
            "/x/.pio/libdeps/teensy40/Foo/src/effects/baz.h:2:1: warning: w [-Wz]"))

    def test_nested_paths_do_not_alias_to_one_key(self):
        # A nested targets/.../effects/Foo.h and a top-level effects/Foo.h are
        # distinct files; relativizing to the repo root must keep them apart so a
        # new warning in one cannot be masked by a baseline entry from the other.
        root = "/home/runner/work/Holosphere/Holosphere/"
        nested = tw.normalize(
            root + "targets/Phantasm/effects/Foo.h:7:1: warning: w [-Wx]")
        top = tw.normalize(root + "effects/Foo.h:7:1: warning: w [-Wx]")
        self.assertTrue(nested.startswith("targets/Phantasm/effects/Foo.h:"))
        self.assertTrue(top.startswith("effects/Foo.h:"))
        self.assertNotEqual(nested, top)

    def test_new_warning_fails_but_reorder_passes(self):
        baseline = {"core/a.h: warning: w1 [-Wx]", "core/b.h: warning: w2 [-Wy]"}
        # Same set, different order -> no new warnings (set-based).
        same = {"core/b.h: warning: w2 [-Wy]", "core/a.h: warning: w1 [-Wx]"}
        self.assertEqual(same - baseline, set())
        # A genuinely new warning -> flagged.
        added = baseline | {"core/c.h: warning: w3 [-Wz]"}
        self.assertEqual(added - baseline, {"core/c.h: warning: w3 [-Wz]"})

    def test_extract_warnings_dedups_and_filters(self):
        log = "\n".join([
            "core/effects/Foo.h:1:1: warning: dup [-Wd]",
            "core/effects/Foo.h:9:9: warning: dup [-Wd]",          # same after normalize
            "/x/.platformio/packages/framework-arduinoteensy/cores/teensy4/usb.c:5:1: warning: lib [-Wl]",
            "hardware/pov_single.h:3:2: warning: real [-Wr]",
        ])
        got = tw.extract_warnings(log)
        self.assertEqual(got, {
            "core/effects/Foo.h: warning: dup [-Wd]",
            "hardware/pov_single.h: warning: real [-Wr]",
        })


def _run_ratchet(log_text, *extra, envs=None):
    """Run the ratchet over `log_text` against an empty baseline; return its exit.

    `envs` is the environment set the build was asked to produce, written to a
    throwaway platformio.ini. It defaults to the environments `log_text` itself
    contains, which isolates a test from the repo's real environment list.
    """
    with tempfile.TemporaryDirectory() as d:
        log = Path(d) / "build.log"
        log.write_text(log_text, encoding="utf-8")
        base = Path(d) / "baseline.txt"
        base.write_text("", encoding="utf-8")
        if envs is None:
            envs = [s.name for s in tw.parse_env_sections(log_text)]
        ini = Path(d) / "platformio.ini"
        ini.write_text("".join(f"[env:{e}]\n" for e in envs), encoding="utf-8")
        return tw.main(["--build-log", str(log), "--baseline", str(base),
                        "--platformio-ini", str(ini), *extra])


def _banner(env, *sources):
    """A `Processing <env> (...)` banner declaring exactly these first-party TUs."""
    terms = ", ".join(["-<*>"] + [f"+<{s}>" for s in sources])
    return (f"Processing {env} (board: teensy40; build_src_filter: {terms}; "
            f"platform: teensy@5.2.0; framework: arduino)")


class TestWarningRatchetCaptureEvidence(unittest.TestCase):
    """A broken capture must not read as today's healthy green."""

    PIO_LINE = "Compiling .pio/build/phantasm/targets/Phantasm/Phantasm.ino.cpp.o"
    PIO_BANNER = _banner("phantasm", "targets/Phantasm/Phantasm.ino.cpp")
    VERBOSE_LINE = ("arm-none-eabi-g++ -o .pio/build/phantasm/core/engine/memory.cpp.o "
                    "-c core/engine/memory.cpp")
    # `pio run -v`: `-c` is a bare flag, the source is the LAST argument.
    VERBOSE_SOURCE_LAST = (
        "arm-none-eabi-g++ -o .pio/build/phantasm/src/core/engine/memory.cpp.o "
        "-c -std=gnu++20 -fno-exceptions -O3 -DPLATFORMIO=60119 "
        "-I. -Icore -Ieffects -Ihardware core/engine/memory.cpp")

    def _run(self, log_text):
        return _run_ratchet(log_text)

    def test_pio_step_line_counts_as_first_party(self):
        self.assertEqual(tw.count_first_party_compiles(self.PIO_LINE), 1)

    def test_verbose_invocation_counts_as_first_party(self):
        self.assertEqual(tw.count_first_party_compiles(self.VERBOSE_LINE), 1)

    def test_third_party_compiles_are_not_evidence(self):
        log = "\n".join([
            "Compiling .pio/build/phantasm/FrameworkArduino/usb.c.o",
            "arm-none-eabi-gcc -o x.o -c /root/.platformio/packages/f/cores/t4/usb.c",
        ])
        self.assertEqual(tw.count_first_party_compiles(log), 0)

    def test_empty_log_fails(self):
        self.assertEqual(self._run(""), 1)

    def test_log_with_no_compiles_fails_even_with_zero_warnings(self):
        self.assertEqual(self._run("Environment    Status    Duration\nphantasm  SUCCESS\n"), 1)

    def test_log_with_first_party_compiles_and_no_new_warnings_passes(self):
        self.assertEqual(self._run(self.PIO_BANNER + "\n" + self.PIO_LINE + "\n"), 0)

    def test_capture_evidence_does_not_mask_a_new_warning(self):
        log = (self.PIO_BANNER + "\n" + self.PIO_LINE
               + "\ncore/render/sdf.h:9:1: warning: novel [-Wnovel]\n")
        self.assertEqual(self._run(log), 1)

    def test_source_last_invocation_counts_as_first_party(self):
        self.assertEqual(tw.count_first_party_compiles(self.VERBOSE_SOURCE_LAST), 1)
        banner = _banner("phantasm", "core/engine/memory.cpp")
        self.assertEqual(self._run(banner + "\n" + self.VERBOSE_SOURCE_LAST + "\n"), 0)

    def test_include_flags_alone_are_not_evidence(self):
        # -Icore/-Ieffects/-Ihardware ride on every invocation, third-party ones
        # included; only the positional source may vote.
        line = ("arm-none-eabi-g++ -o .pio/build/phantasm/lib999/FastLED/noise.cpp.o "
                "-c -std=gnu++20 -I. -Icore -Ieffects -Ihardware "
                ".pio/libdeps/phantasm/FastLED/src/noise.cpp")
        self.assertEqual(tw.count_first_party_compiles(line), 0)

    def test_link_line_naming_first_party_objects_is_not_evidence(self):
        # No -c: linking cannot emit a compile warning.
        line = ("arm-none-eabi-g++ -o .pio/build/phantasm/firmware.elf -T imxrt1062.ld "
                ".pio/build/phantasm/src/core/engine/memory.cpp.o")
        self.assertEqual(tw.count_first_party_compiles(line), 0)

    def test_preprocess_only_invocation_is_not_evidence(self):
        # PlatformIO's .ino -> .cpp preprocess pass runs -E, not -c.
        line = ('arm-none-eabi-g++ -o "/w/pov/targets/Phantasm/Phantasm.ino.cpp" '
                '-x c++ -fpreprocessed -dD -E "/tmp/tmpm5mgtz1h"')
        self.assertEqual(tw.count_first_party_compiles(line), 0)


class TestColdCaptureAudit(unittest.TestCase):
    """A partially cached build must FAIL, not pass on a shrunken warning set.

    The old guard only demanded ONE first-party compile, so a build serving 12 of
    18 TUs from PlatformIO's object cache read as green while two thirds of the
    warning surface was never compiled. The expectation is derived from
    `build_src_filter` in PlatformIO's own banner, so it tracks a new TU or a new
    environment with nothing to keep in sync.
    """

    TUS = ("core/engine/memory.cpp", "core/spatial/reaction_graph.cpp",
           "targets/Phantasm/Phantasm.ino.cpp")

    def _compile(self, env, source):
        return (f"arm-none-eabi-g++ -o .pio/build/{env}/src/{source}.o -c "
                f"-std=gnu++20 -O3 -I. -Icore -Ieffects -Ihardware {source}")

    def _cache_hit(self, env, source):
        return f"Retrieved `.pio/build/{env}/src/{source}.o' from cache"

    def _log(self, envs, compiled, cached=()):
        lines = []
        for env in envs:
            lines.append(_banner(env, *self.TUS))
            lines += [self._compile(env, s) for s in compiled]
            lines += [self._cache_hit(env, s) for s in cached]
        return "\n".join(lines) + "\n"

    def test_expectation_is_derived_from_the_banner(self):
        section, = tw.parse_env_sections(_banner("phantasm", *self.TUS))
        self.assertEqual(section.name, "phantasm")
        self.assertEqual(tw.declared_first_party_sources(section), set(self.TUS))

    def test_exact_count_across_every_environment_passes(self):
        log = self._log(("holosphere", "phantasm", "profile"), self.TUS)
        self.assertEqual(_run_ratchet(log), 0)

    def test_short_count_from_the_object_cache_fails(self):
        # The reproduced failure: only the sketch recompiles, the shared core TUs
        # come from build_cache_dir.
        log = self._log(("holosphere", "phantasm", "profile"),
                        self.TUS[2:], cached=self.TUS[:2])
        self.assertEqual(_run_ratchet(log), 1)

    def test_short_count_in_a_single_environment_fails(self):
        # One fully-cold env cannot vouch for another: the audit is per-env.
        log = (self._log(("holosphere",), self.TUS)
               + self._log(("phantasm",), self.TUS[:2]))
        self.assertEqual(_run_ratchet(log), 1)

    def test_short_count_diagnostic_names_the_cache_and_the_missing_tus(self):
        log = self._log(("phantasm",), self.TUS[2:], cached=self.TUS[:2])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(_run_ratchet(log, "--github"), 1)
        out = buf.getvalue()
        self.assertIn("::error::", out)
        self.assertIn("2 of 3 first-party translation unit(s)", out)
        self.assertIn("object cache", out)
        for tu in self.TUS[:2]:
            self.assertIn(tu, out)

    def test_incremental_build_with_no_cache_hits_fails_and_says_so(self):
        log = self._log(("phantasm",), self.TUS[2:])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(_run_ratchet(log), 1)
        self.assertIn("incremental build", buf.getvalue())

    def test_a_new_tu_raises_the_bar_with_no_second_edit(self):
        # Adding a TU to build_src_filter moves the expectation by itself.
        grown = self.TUS + ("core/engine/newtu.cpp",)
        banner = _banner("phantasm", *grown)
        compiles = "\n".join(self._compile("phantasm", s) for s in self.TUS)
        self.assertEqual(_run_ratchet(banner + "\n" + compiles + "\n"), 1)
        self.assertEqual(_run_ratchet(
            banner + "\n" + compiles + "\n"
            + self._compile("phantasm", "core/engine/newtu.cpp") + "\n"), 0)

    def test_pio_step_line_shape_also_satisfies_the_expectation(self):
        # Without -v PlatformIO names the OBJECT; the `.o` suffix must not stop it
        # matching the source build_src_filter declares.
        log = (_banner("phantasm", *self.TUS) + "\n"
               + "\n".join(f"Compiling .pio/build/phantasm/src/{s}.o"
                           for s in self.TUS) + "\n")
        self.assertEqual(_run_ratchet(log), 0)

    def test_missing_banner_fails_rather_than_falling_back(self):
        # If PlatformIO's banner format ever changes the expectation cannot be
        # derived; that must go red, not silently revert to "at least one compile".
        self.assertEqual(_run_ratchet(self._compile("phantasm", self.TUS[0]) + "\n"), 1)

    def test_banner_without_build_src_filter_fails(self):
        log = ("Processing phantasm (board: teensy40; platform: teensy@5.2.0)\n"
               + self._compile("phantasm", self.TUS[0]) + "\n")
        self.assertEqual(_run_ratchet(log), 1)

    def test_first_party_glob_in_the_filter_fails_loudly(self):
        # A glob is not countable from the log; the tool must say so, not guess.
        section, = tw.parse_env_sections(_banner("phantasm", "core/engine/*.cpp"))
        with self.assertRaises(tw.CaptureError):
            tw.declared_first_party_sources(section)

    def test_third_party_filter_terms_are_not_expected_tus(self):
        section, = tw.parse_env_sections(
            _banner("phantasm", "core/engine/memory.cpp", "lib/Foo/foo.cpp"))
        self.assertEqual(tw.declared_first_party_sources(section),
                         {"core/engine/memory.cpp"})

    def test_exclusion_term_removes_a_declared_tu(self):
        header = _banner("phantasm", *self.TUS).replace(
            "; platform:", " -<core/engine/memory.cpp>; platform:")
        section, = tw.parse_env_sections(header)
        self.assertEqual(tw.declared_first_party_sources(section), set(self.TUS[1:]))

    def test_update_baseline_is_also_gated_on_a_cold_capture(self):
        # Regenerating from a partial build would silently drop warnings.
        log = self._log(("phantasm",), self.TUS[2:], cached=self.TUS[:2])
        self.assertEqual(_run_ratchet(log, "--update-baseline"), 1)

    def test_update_baseline_refuses_an_env_narrowed_capture(self):
        # A whole-file rewrite from one environment's warnings would delete the
        # warnings every other environment carries; the file stays untouched.
        log = self._log(("phantasm",), self.TUS)
        with tempfile.TemporaryDirectory() as d:
            build_log = Path(d) / "build.log"
            build_log.write_text(log, encoding="utf-8")
            base = Path(d) / "baseline.txt"
            base.write_text("core/effects/Foo.h: warning: w [-Wx]\n",
                            encoding="utf-8")
            ini = Path(d) / "platformio.ini"
            ini.write_text("[env:phantasm]\n[env:holosphere]\n",
                           encoding="utf-8")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                code = tw.main(["--build-log", str(build_log),
                                "--baseline", str(base),
                                "--platformio-ini", str(ini),
                                "--env", "phantasm", "--update-baseline"])
            self.assertEqual(code, 1)
            self.assertIn("refusing to rewrite", buf.getvalue())
            self.assertIn("holosphere", buf.getvalue())
            self.assertIn("Foo.h", base.read_text(encoding="utf-8"))

    def test_update_baseline_accepts_a_capture_of_every_environment(self):
        log = self._log(("holosphere", "phantasm"), self.TUS)
        self.assertEqual(_run_ratchet(log, "--update-baseline"), 0)


class TestExpectedEnvironmentSet(unittest.TestCase):
    """The audited environments must be the ones the build was asked to produce.

    Sections come from banners, so a `pio run` over six environments that dies in
    the first prints ONE banner: the other five are absent rather than short, and
    the per-environment coldness audit has nothing to complain about.
    """

    TU = "core/engine/memory.cpp"
    CI_ENVS = ("holosphere", "holosphere_dma", "phantasm", "phantasm8",
               "profile", "profile_o3")

    def _cold_env(self, env):
        return (_banner(env, self.TU) + "\n"
                + f"arm-none-eabi-g++ -o .pio/build/{env}/src/{self.TU}.o -c "
                  f"-std=gnu++20 -O3 -I. -Icore {self.TU}\n")

    def test_truncated_run_fails(self):
        self.assertEqual(
            _run_ratchet(self._cold_env(self.CI_ENVS[0]), envs=self.CI_ENVS), 1)

    def test_truncated_run_diagnostic_names_the_absent_environments(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(_run_ratchet(self._cold_env(self.CI_ENVS[0]),
                                          "--github", envs=self.CI_ENVS), 1)
        out = buf.getvalue()
        self.assertIn("::error::", out)
        self.assertIn("5 of 6 expected environment(s)", out)
        for env in self.CI_ENVS[1:]:
            self.assertIn(env, out)

    def test_every_expected_environment_present_passes(self):
        log = "".join(self._cold_env(e) for e in self.CI_ENVS)
        self.assertEqual(_run_ratchet(log, envs=self.CI_ENVS), 0)

    def test_env_flag_narrows_the_expectation(self):
        # A deliberate subset build states its own set instead of platformio.ini's.
        self.assertEqual(
            _run_ratchet(self._cold_env(self.CI_ENVS[0]), "--env", self.CI_ENVS[0],
                         envs=self.CI_ENVS), 0)

    def test_expectation_is_read_from_the_repo_platformio_ini(self):
        envs = tw.declared_environments(TOOLS.parent / "platformio.ini")
        self.assertNotIn("env", envs)          # the shared [env] base section
        self.assertLessEqual({"holosphere", "phantasm", "profile"}, set(envs))

    def test_missing_platformio_ini_fails_rather_than_skipping_the_check(self):
        with self.assertRaises(tw.CaptureError):
            tw.declared_environments(TOOLS / "no_such_platformio.ini")

    def test_ini_without_environments_fails(self):
        with tempfile.TemporaryDirectory() as d:
            ini = Path(d) / "platformio.ini"
            ini.write_text("[platformio]\nsrc_dir = .\n[env]\n", encoding="utf-8")
            with self.assertRaises(tw.CaptureError):
                tw.declared_environments(ini)


class TestRealColdVersusWarmCapture(unittest.TestCase):
    """Verbatim `pio run -v` sections from a real cold and a real warm build.

    fixtures/real/{cold,warm}_env_section.txt are the `holosphere` sections of two
    consecutive runs of the CI command: the first with `.pio/build_cache` deleted,
    the second reusing it. Only a real capture pins PlatformIO's banner text and
    SCons's `Retrieved … from cache` line, which the derived expectation reads.
    """

    COLD = (REAL_DIR / "cold_env_section.txt").read_text(encoding="utf-8")
    WARM = (REAL_DIR / "warm_env_section.txt").read_text(encoding="utf-8")
    TUS = {"core/engine/memory.cpp", "core/spatial/reaction_graph.cpp",
           "targets/Holosphere/Holosphere.ino.cpp"}

    def test_real_banner_declares_three_first_party_tus(self):
        section, = tw.parse_env_sections(self.COLD)
        self.assertEqual(section.name, "holosphere")
        self.assertEqual(tw.declared_first_party_sources(section), self.TUS)

    def test_real_cold_section_compiles_every_declared_tu(self):
        audit = tw.audit_capture(self.COLD)
        self.assertEqual(audit.missing, ())
        self.assertEqual(audit.cache_hits, 0)
        self.assertEqual(_run_ratchet(self.COLD), 0)

    def test_real_warm_section_is_short_and_fails(self):
        audit = tw.audit_capture(self.WARM)
        self.assertEqual(audit.missing_count, 2)
        self.assertEqual(audit.first_party_cache_hits, 2)
        self.assertEqual(_run_ratchet(self.WARM), 1)


class TestRealVerboseCapture(unittest.TestCase):
    """Capture evidence against REAL `pio run -v` lines, both CI and Windows.

    fixtures/real/verbose_build_log.txt holds verbatim invocations in a fixed
    order: three first-party then two third-party from a Windows build
    (backslash paths), then one of each from the Linux CI runner (forward
    slashes). Only a real capture pins the argument order `-v` emits — `-c` is a
    bare flag and the source trails the whole flag list.
    """

    LINES = (REAL_DIR / "verbose_build_log.txt").read_text(
        encoding="utf-8").splitlines()
    FIRST_PARTY = LINES[:3] + LINES[5:6]
    THIRD_PARTY = LINES[3:5] + LINES[6:]

    def test_every_line_yields_exactly_one_source(self):
        for line in self.LINES:
            self.assertEqual(len(tw.compiled_paths(line)), 1, line[-60:])

    def test_real_first_party_invocations_are_evidence(self):
        for line in self.FIRST_PARTY:
            self.assertEqual(tw.count_first_party_compiles(line), 1, line[-60:])

    def test_real_third_party_invocations_are_not_evidence(self):
        for line in self.THIRD_PARTY:
            self.assertEqual(tw.count_first_party_compiles(line), 0, line[-60:])
        self.assertEqual(
            tw.count_first_party_compiles("\n".join(self.THIRD_PARTY)), 0)

    def test_whole_real_log_counts_first_party_only(self):
        self.assertEqual(
            tw.count_first_party_compiles("\n".join(self.LINES)), 4)


class TestRealCapture(unittest.TestCase):
    """Parse REAL toolchain output from the two shipping images.

    These exercise the toolchain quirks the synthetic fixtures can't: teensy_size
    on stderr, readelf printing large sizes in HEX (0x4a800), and the arena's real
    _ZL-mangled internal-linkage name. Captured from `pio run -e holosphere
    -e phantasm` on the toolchain platformio.ini pins, so holosphere is the
    96x20 canvas the env ships.
    """

    def test_real_holosphere_build_passes_the_calibrated_gate(self):
        sizes = tg.parse_teensy_size((REAL_DIR / "holosphere_teensy_size.txt").read_text())
        syms = tg.parse_readelf_symbols((REAL_DIR / "holosphere_readelf_syms.txt").read_text())
        secs = tg.parse_readelf_sections((REAL_DIR / "holosphere_readelf_secs.txt").read_text())
        result = tg.evaluate("holosphere", BUDGETS["holosphere"], sizes, syms, secs)
        self.assertTrue(result.passed, msg=_codes(result))
        arena = next(s for s in syms if s.name == "_ZL18global_arena_block")
        self.assertEqual(arena.region, "DTCM")
        self.assertEqual(arena.size, 305152)  # 0x4a800, parsed from the hex Size column

    def test_real_phantasm_reaction_graph_is_in_flash(self):
        syms = tg.parse_readelf_symbols((REAL_DIR / "phantasm_readelf_syms.txt").read_text())
        rg = next(s for s in syms if s.name == "_ZN13ReactionGraph9neighborsE")
        self.assertEqual(rg.region, "FLASH")
        self.assertEqual(rg.size, 92160)

    def test_real_phantasm_dma_tx_buffer_is_in_ocram(self):
        # Spelled from the budget, so a linkage-name drift between the shipped
        # budget and the captured ELF fails here instead of passing on a name
        # only this test still knows.
        name = BUDGETS["phantasm"]["symbols"]["dma_tx_buffer"]["name"]
        syms = tg.parse_readelf_symbols((REAL_DIR / "phantasm_readelf_syms.txt").read_text())
        led = next(s for s in syms if s.name == name)
        self.assertEqual(led.region, "OCRAM")

    def test_real_phantasm_build_passes_the_calibrated_gate(self):
        sizes = tg.parse_teensy_size((REAL_DIR / "phantasm_teensy_size.txt").read_text())
        syms = tg.parse_readelf_symbols((REAL_DIR / "phantasm_readelf_syms.txt").read_text())
        secs = tg.parse_readelf_sections((REAL_DIR / "phantasm_readelf_secs.txt").read_text())
        result = tg.evaluate("phantasm", BUDGETS["phantasm"], sizes, syms, secs)
        self.assertTrue(result.passed, msg=_codes(result))
        # The bank-derived ITCM ceiling is phantasm-only and unmeasurable from
        # `size -A`; a capture without the code component would pass every other
        # check while leaving that branch unexercised.
        self.assertIn("code", sizes["ram1"]["components"])

    def test_real_phantasm_calibration_prose_matches_the_gate(self):
        sizes = tg.parse_teensy_size(
            (REAL_DIR / "phantasm_teensy_size.txt").read_text())
        ram1 = sizes["ram1"]
        components = ram1["components"]
        self.assertEqual(
            (components["variables"], components["code"], ram1["free"]),
            (313760, 193496, 13920))

        result = tg.evaluate("phantasm", BUDGETS["phantasm"], sizes, [], {})
        note = next(line for line in result.notes
                    if "'code' derived ceiling" in line)
        remaining_text = note.split("remaining ", 1)[1].split(" B;", 1)[0]
        remaining = int(remaining_text.replace(",", ""))
        allocated_code = components["code"] + components["padding"]

        budgets_text = " ".join(
            line.strip().removeprefix("//").removeprefix("#").strip()
            for line in (TOOLS / "teensy_budgets.json")
            .read_text().splitlines())
        self.assertIn(
            f"it uses {ram1['used']:,} B of RAM1 "
            f"(variables {components['variables']:,} + code "
            f"{components['code']:,} rounded up to a whole 32 KiB FlexRAM "
            f"bank, {allocated_code:,}) with {ram1['free']:,} B free for locals",
            budgets_text)

        ci_text = " ".join(
            line.strip().removeprefix("//").removeprefix("#").strip()
            for line in (TOOLS.parent / ".github/workflows/ci.yml")
            .read_text().splitlines())
        self.assertIn(
            f"free-for-locals {ram1['free']:,} B over a measured 12 KiB floor, "
            f"and ITCM code {remaining:,} B under its bank-derived ceiling",
            ci_text)

    def test_real_phantasm_exception_index_is_routed_to_flash(self):
        # tools/phantasm.ld keeps .ARM.exidx out of the FlexRAM banks; in ITCM it
        # would consume the headroom the derived code ceiling hands out.
        secs = tg.parse_readelf_sections((REAL_DIR / "phantasm_readelf_secs.txt").read_text())
        addr = next(a for name, a in secs.values() if name == ".ARM.exidx")
        self.assertEqual(tg.region_for_address(addr), "FLASH")


class TestSizeAFallback(unittest.TestCase):
    """The `--size-a` fallback path (no teensy_size): VMA bucketing + 0x80000-ram1
    free-headroom arithmetic, end-to-end through main() and evaluate()."""

    def test_free_headroom_is_0x80000_minus_used(self):
        sizes = tg.fallback_sizes_from_size_a(_read("good_size_a.txt"))
        # ITCM 0xf320 + 0x2a00 occupies three whole FlexRAM banks; DTCM is the
        # raw 0x1a00 + 0x53c00.
        expected_ram1 = 3 * tg.FLEXRAM_BANK_BYTES + 0x1a00 + 0x53c00
        self.assertEqual(sizes["ram1"]["used"], expected_ram1)
        self.assertEqual(sizes["ram1"]["free"], 0x80000 - expected_ram1)
        self.assertEqual(sizes["ram2"]["free"], 0x80000 - 0x79900)

    def test_itcm_is_charged_by_whole_flexram_banks(self):
        # phantasm.ld rounds .text.itcm up to a bank before DTCM gets the rest,
        # so one byte over a boundary costs a whole bank of RAM1 headroom.
        under = tg.fallback_sizes_from_size_a(
            _size_a(tg.FLEXRAM_BANK_BYTES, 0x40000, 0x70000, 0x20000))
        over = tg.fallback_sizes_from_size_a(
            _size_a(tg.FLEXRAM_BANK_BYTES + 1, 0x40000, 0x70000, 0x20000))
        self.assertEqual(under["ram1"]["used"],
                         tg.FLEXRAM_BANK_BYTES + 0x40000)
        self.assertEqual(over["ram1"]["used"],
                         under["ram1"]["used"] + tg.FLEXRAM_BANK_BYTES)

    def test_main_size_a_fallback_passes_a_fitting_build(self):
        rc, out = self._run_main_size_a(
            _size_a(0x10000, 0x40000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "holosphere")
        self.assertEqual(rc, tg.EXIT_UNCALIBRATED_PASS, msg=out)
        self.assertIn("PASS", out)

    def test_main_size_a_fallback_free_arithmetic_trips_headroom_floor(self):
        # DTCM large enough that 0x80000 - ram1 drops below the 32,768 B floor,
        # proving the fallback `free` figure actually drives the gate decision.
        # ram1 = 0x10000 + 0x70000 = 0x80000 -> free 0; floor 32768 -> violation.
        rc, out = self._run_main_size_a(
            _size_a(0x10000, 0x70000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "holosphere")
        self.assertEqual(rc, 1, msg=out)
        self.assertIn("free-for-local-variables", out)

    def test_fallback_pass_is_marked_advisory_not_calibrated(self):
        # A fallback PASS must carry the advisory note so it is never mistaken
        # for a teensy_size-calibrated verdict.
        rc, out = self._run_main_size_a(
            _size_a(0x10000, 0x40000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "holosphere")
        self.assertEqual(rc, tg.EXIT_UNCALIBRATED_PASS, msg=out)
        self.assertIn("PASS", out)
        self.assertIn("ADVISORY", out)
        self.assertIn("not calibrated", out.lower())

    def test_fallback_pass_exit_code_differs_from_a_calibrated_pass(self):
        # The exit status is what an automated caller reads: an uncalibrated
        # PASS must not be reportable as a real one, and must stay distinct from
        # the cannot-run code so a guess is not confused with a tooling break.
        rc_advisory, out_advisory = self._run_main_size_a(
            _size_a(0x10000, 0x40000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "holosphere")
        rc_real, out_real = self._run_main_teensy_size(
            "good_teensy_size.txt", "good_readelf_syms.txt", "holosphere")
        self.assertEqual(rc_real, 0, msg=out_real)
        self.assertIn("PASS", out_real)
        self.assertNotIn("ADVISORY", out_real)
        self.assertEqual(rc_advisory, tg.EXIT_UNCALIBRATED_PASS, msg=out_advisory)
        self.assertNotIn(tg.EXIT_UNCALIBRATED_PASS, (0, 1, 2))

    def test_main_rejects_a_component_env_under_the_fallback(self):
        # phantasm's ram1 code ceiling is unmeasurable from `size -A`; the missing
        # tool must read as cannot-run, not as a component-missing violation.
        rc, out = self._run_main_size_a(
            _size_a(0x10000, 0x40000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "phantasm")
        self.assertEqual(rc, 2, msg=out)
        self.assertIn("per-component ceilings", out)
        self.assertNotIn("component-missing", out)

    def test_main_rejects_invalid_size_a_output_as_tooling_error(self):
        for name, text in _invalid_size_a_cases().items():
            with self.subTest(name=name):
                rc, out = self._run_main_size_a(
                    text, "good_readelf_syms.txt", "holosphere")
                self.assertEqual(rc, 2, msg=out)
                self.assertIn("invalid `size -A` output", out)
                self.assertIn("tooling/format error", out)

    def test_real_size_a_captures_parse_through_the_fallback(self):
        # Real `arm-none-eabi-size -A` output opens with a `<path> :` header whose
        # path starts with '.'; the synthetic fixtures' header does not, so only a
        # real capture proves the fallback is reachable at all.
        for name in ("holosphere_size_a.txt", "phantasm_size_a.txt"):
            with self.subTest(capture=name):
                sizes = tg.fallback_sizes_from_size_a(
                    (REAL_DIR / name).read_text(encoding="utf-8"))
                self.assertEqual(set(sizes), {"flash", "ram1", "ram2"})
                for region, measured in sizes.items():
                    self.assertGreater(measured["used"], 0, msg=region)

    def test_real_holosphere_capture_reaches_an_advisory_verdict(self):
        rc, out = self._run_main(
            "--size-a",
            (REAL_DIR / "holosphere_size_a.txt").read_text(encoding="utf-8"),
            (REAL_DIR / "holosphere_readelf_syms.txt").read_text(encoding="utf-8"),
            "holosphere")
        self.assertEqual(rc, tg.EXIT_UNCALIBRATED_PASS, msg=out)
        self.assertIn("PASS", out)
        self.assertIn("ADVISORY", out)

    def _run_main_size_a(self, size_a_text, syms_fixture, env):
        return self._run_main("--size-a", size_a_text, _read(syms_fixture), env)

    def _run_main_teensy_size(self, teensy_size_fixture, syms_fixture, env):
        return self._run_main("--teensy-size", _read(teensy_size_fixture),
                              _read(syms_fixture), env)

    def _run_main(self, sizes_flag, sizes_text, syms_text, env):
        with tempfile.TemporaryDirectory() as d:
            sizes = Path(d) / "sizes.txt"
            sizes.write_text(sizes_text, encoding="utf-8")
            syms = Path(d) / "syms.txt"
            syms.write_text(syms_text, encoding="utf-8")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = tg.main([
                    "--env", env,
                    "--budgets", str(TOOLS / "teensy_budgets.json"),
                    sizes_flag, str(sizes),
                    "--readelf-syms", str(syms),
                ])
            return rc, buf.getvalue()


class TestNonUtf8Captures(unittest.TestCase):
    """Both gates answer by exit code, and a decode error has none.

    A Windows `pio run -v 2>&1 | tee` interleaves cp1252 bytes into the stream,
    so a capture that is not valid UTF-8 must still produce a verdict.
    """

    # RIGHT SINGLE QUOTATION MARK in cp1252; not a valid UTF-8 sequence.
    CP1252 = b"don\x92t"

    def test_gate_reads_a_capture_with_a_cp1252_byte(self):
        with tempfile.TemporaryDirectory() as d:
            sizes = Path(d) / "sizes.txt"
            sizes.write_bytes(_read("good_teensy_size.txt").encode("utf-8")
                              + b"note: " + self.CP1252 + b" care\n")
            syms = Path(d) / "syms.txt"
            syms.write_bytes(_read("good_readelf_syms.txt").encode("utf-8"))
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = tg.main(["--env", "holosphere",
                              "--budgets", str(TOOLS / "teensy_budgets.json"),
                              "--teensy-size", str(sizes),
                              "--readelf-syms", str(syms)])
            self.assertEqual(rc, 0, msg=buf.getvalue())

    def test_ratchet_reads_a_build_log_with_a_cp1252_byte(self):
        log_text = (_banner("phantasm", "core/engine/memory.cpp") + "\n"
                    "arm-none-eabi-g++ -o .pio/build/phantasm/core/engine/"
                    "memory.cpp.o -c core/engine/memory.cpp\n")
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "build.log"
            log.write_bytes(log_text.encode("utf-8")
                            + b"note: " + self.CP1252 + b" care\n")
            base = Path(d) / "baseline.txt"
            base.write_text("", encoding="utf-8")
            ini = Path(d) / "platformio.ini"
            ini.write_text("[env:phantasm]\n", encoding="utf-8")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = tw.main(["--build-log", str(log), "--baseline", str(base),
                              "--platformio-ini", str(ini)])
            self.assertEqual(rc, 0, msg=buf.getvalue())


class TestStripJsoncComments(unittest.TestCase):
    """The bespoke JSONC comment stripper guarding the budgets file."""

    def test_strips_line_and_block_comments(self):
        text = '{\n  // line\n  "a": 1, /* block */ "b": 2\n}'
        self.assertEqual(tg._strip_jsonc_comments(text),
                         '{\n  \n  "a": 1,  "b": 2\n}')

    def test_keeps_comment_sequences_inside_strings(self):
        # // and /* inside a JSON string value must survive untouched (a path,
        # URL, or note must not be mangled).
        text = '{"url": "http://x/y", "note": "a /* not a */ comment"}'
        self.assertEqual(tg._strip_jsonc_comments(text), text)

    def test_escaped_quote_does_not_end_string(self):
        text = r'{"s": "a\"// still in string", "n": 1}'
        self.assertEqual(tg._strip_jsonc_comments(text), text)

    def test_block_comment_with_newlines(self):
        text = '{\n/* multi\nline */\n"a": 1}'
        self.assertEqual(tg._strip_jsonc_comments(text), '{\n\n"a": 1}')

    def test_unterminated_block_comment_raises(self):
        with self.assertRaises(ValueError):
            tg._strip_jsonc_comments('{"a":1} /* never closed')

    def test_round_trips_through_json_load(self):
        import json
        cfg = json.loads(tg._strip_jsonc_comments(
            '{\n  // c\n  "x": 1, /* y */ "s": "10 // 2"\n}'))
        self.assertEqual(cfg, {"x": 1, "s": "10 // 2"})


class TestToolingFailureExits(unittest.TestCase):
    """A tooling break exits 2 (cannot-run), never 1 (budget violation).

    An uncaught exception exits 1, which every status-reading caller reads as a
    size regression — the misattribution the gate's exit contract exists to
    prevent.
    """

    def _run(self, **over):
        argv = {"--budgets": str(TOOLS / "teensy_budgets.json"),
                "--teensy-size": str(FIX / "good_teensy_size.txt"),
                "--readelf-syms": str(FIX / "good_readelf_syms.txt")}
        argv.update(over)
        err, out = io.StringIO(), io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(out):
            rc = tg.main(["--env", "holosphere"]
                         + [a for pair in argv.items() for a in pair
                            if pair[1] is not None])
        return rc, err.getvalue() + out.getvalue()

    def test_malformed_budgets_json_is_cannot_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "budgets.json"
            path.write_text('{ "holosphere": { ', encoding="utf-8")
            rc, text = self._run(**{"--budgets": str(path)})
        self.assertEqual(rc, 2, msg=text)

    def test_unterminated_jsonc_comment_is_cannot_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "budgets.json"
            path.write_text('{ /* never closed\n"a": 1}', encoding="utf-8")
            rc, text = self._run(**{"--budgets": str(path)})
        self.assertEqual(rc, 2, msg=text)

    def test_missing_budgets_file_is_cannot_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            rc, text = self._run(**{"--budgets": str(Path(tmp) / "no.json")})
        self.assertEqual(rc, 2, msg=text)

    def test_unlisted_env_is_cannot_run(self):
        # Same verdict tools/teensy_gate_extra.py gives: an env with no budget
        # entry is a budgets-file error, not a size-budget violation.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "budgets.json"
            path.write_text('{ "phantasm": {} }', encoding="utf-8")
            rc, text = self._run(**{"--budgets": str(path)})
        self.assertEqual(rc, 2, msg=text)

    def test_unparseable_teensy_size_blob_is_cannot_run(self):
        # Zeroed region totals would otherwise render PASS on an unmeasured image.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "size.txt"
            path.write_text(UNPARSEABLE_TEENSY_SIZE, encoding="utf-8")
            rc, text = self._run(**{"--teensy-size": str(path)})
        self.assertEqual(rc, 2, msg=text)
        self.assertIn("invalid teensy_size output", text)
        self.assertIn("tooling/format error", text)
        self.assertNotIn("PASS", text)

    def test_missing_capture_file_is_cannot_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            gone = str(Path(tmp) / "gone.txt")
            for flag in ("--teensy-size", "--readelf-syms"):
                with self.subTest(flag=flag):
                    rc, text = self._run(**{flag: gone})
                    self.assertEqual(rc, 2, msg=text)
            rc, text = self._run(**{"--size-a": gone, "--teensy-size": None})
        self.assertEqual(rc, 2, msg=text)


class TestGateExtra(unittest.TestCase):
    """The PlatformIO post-build glue: toolchain discovery + exit(2) guards."""

    def setUp(self):
        self.ge = _load_gate_extra()

    def test_gate_runs_after_the_elf_links(self):
        self.assertEqual([t for t, _ in self.ge.env.post_actions],
                         ["$BUILD_DIR/${PROGNAME}.elf"])
        self.assertIs(self.ge.env.post_actions[0][1], self.ge.run_gate)

    def test_elf_depends_on_gate_inputs(self):
        # Without these edges a warm build dir links nothing when an input
        # changes, and the post-action never runs.
        elf = self.ge.env.post_actions[0][0]
        deps = [d for target, group in self.ge.env.dependencies
                if target == elf for d in group]
        for name in ("teensy_budgets.json", "teensy_gate.py",
                     "teensy_gate_extra.py", "phantasm.ld"):
            self.assertIn(str(TOOLS / name), deps)

    def test_declared_gate_sources_exist(self):
        for path in self.ge.GATE_SOURCES:
            self.assertTrue(Path(path).is_file(), f"missing {path}")

    def test_tool_output_replaces_undecodable_bytes(self):
        seen = {}

        def run(*args, **kwargs):
            seen.update(kwargs)
            return subprocess.CompletedProcess(args, 0, "ok", "")

        with mock_patch(self.ge.subprocess, "run", run):
            self.assertEqual(self.ge._run(["tool"]), "ok")
        self.assertEqual(seen["encoding"], "utf-8")
        self.assertEqual(seen["errors"], "replace")

    def test_tool_derives_sibling_arm_tools(self):
        self.assertEqual(self.ge._tool("/opt/arm/bin/arm-none-eabi-gcc", "size"),
                         "/opt/arm/bin/arm-none-eabi-size")
        self.assertEqual(self.ge._tool("C:/x/arm-none-eabi-gcc.exe", "readelf"),
                         "C:/x/arm-none-eabi-readelf.exe")

    def test_tool_falls_back_to_path_lookup(self):
        self.assertEqual(self.ge._tool("clang", "size"), "arm-none-eabi-size")

    def test_candidates_probe_the_tool_package_before_path(self):
        # teensy_size is not on PATH; it ships in the tool-teensy package, so the
        # package paths must be probed first and the bare names only last.
        cands = self.ge._teensy_size_candidates(self._env("/pkg/tool-teensy"))
        self.assertEqual(cands[0], os.path.join("/pkg/tool-teensy", "teensy_size"))
        self.assertEqual(cands[-2:], ["teensy_size", "teensy_size.exe"])

    def test_candidates_include_the_default_core_dir(self):
        # No PlatformIO platform object (plain SCons env): the gate still probes
        # the default core packages dir rather than giving up on PATH alone.
        cands = self.ge._teensy_size_candidates(self._env())
        expected = os.path.join(os.path.expanduser("~"), ".platformio",
                                "packages", "tool-teensy", "teensy_size.exe")
        self.assertIn(expected, cands)

    def test_candidates_surface_platform_lookup_failures(self):
        env = self._env()
        env.PioPlatform = lambda: types.SimpleNamespace(
            get_package_dir=lambda name: (_ for _ in ()).throw(RuntimeError(name)))
        with self.assertRaisesRegex(RuntimeError, "tool-teensy"):
            self.ge._teensy_size_candidates(env)

    def test_find_teensy_size_resolves_from_the_tool_package(self):
        # Only the packaged binary exists; a PATH-only probe would return None
        # and drop the gate to the uncalibrated fallback.
        packaged = os.path.join("/pkg/tool-teensy", "teensy_size.exe")

        def _probe(args, **kwargs):
            if args[0] != packaged:
                raise OSError("not found")
            return types.SimpleNamespace(stdout="", stderr="teensy_size: usage")

        with mock_patch(subprocess, "run", _probe):
            self.assertEqual(self.ge._find_teensy_size(self._env("/pkg/tool-teensy")),
                             packaged)

    def test_find_teensy_size_accepts_self_identifying_on_path(self):
        def _probe(args, **kwargs):
            if args[0] not in ("teensy_size", "teensy_size.exe"):
                raise OSError("not found")
            return types.SimpleNamespace(stdout="usage: teensy_size <elf>", stderr="")

        with mock_patch(subprocess, "run", _probe):
            self.assertEqual(self.ge._find_teensy_size(self._env()), "teensy_size")

    def test_find_teensy_size_rejects_foreign_binary(self):
        probe = types.SimpleNamespace(stdout="usage: other-tool", stderr="")
        with mock_patch(subprocess, "run", lambda *a, **k: probe):
            self.assertIsNone(self.ge._find_teensy_size(self._env()))

    def test_find_teensy_size_none_when_absent(self):
        def _raise(*a, **k):
            raise OSError("not found")
        with mock_patch(subprocess, "run", _raise):
            self.assertIsNone(self.ge._find_teensy_size(self._env()))

    def test_toolchain_oserror_exits_2(self):
        # A tool step raising OSError is a build/tooling break -> exit(2), never a
        # size-budget "violation".
        self.ge._find_teensy_size = lambda env: None
        def _boom(*a, **k):
            raise OSError("no such tool")
        self.ge._run = _boom
        rc, out = self._run_gate("holosphere")
        self.assertEqual(rc, 2)
        self.assertIn("toolchain step failed", out)

    def test_empty_regions_exits_2(self):
        # Tool output the parser no longer recognizes (no FLASH/RAM1/RAM2) is a
        # format break -> exit(2), not a region-missing "violation". teensy_gate
        # is the shared module, so restore parse_teensy_size after the patch.
        self.ge._find_teensy_size = lambda env: "teensy_size"
        self.ge._run = lambda *a, **k: ""
        with mock_patch(self.ge.teensy_gate, "parse_teensy_size", lambda text: {}):
            rc, out = self._run_gate("holosphere")
        self.assertEqual(rc, 2)
        self.assertIn("parsed no FLASH/RAM1/RAM2 regions", out)

    def test_unparseable_teensy_size_blob_exits_2(self):
        self.ge._find_teensy_size = lambda env: "teensy_size"

        def _run(args, check=True):
            if args[0] == "teensy_size":
                return UNPARSEABLE_TEENSY_SIZE
            if "-sW" in args:
                return _read("good_readelf_syms.txt")
            return _read("good_readelf_secs.txt")

        self.ge._run = _run
        rc, out = self._run_gate("holosphere")
        self.assertEqual(rc, 2, msg=out)
        self.assertIn("invalid teensy_size output", out)
        self.assertNotIn("PASS", out)

    def test_size_a_fallback_pass_exits_advisory_not_zero(self):
        # CI accepts the build on this status, so a bucketed guess must not
        # report as a calibrated PASS.
        self.ge._find_teensy_size = lambda env: None

        def _run(args, check=True):
            if "-A" in args:
                return _size_a(0x10000, 0x40000, 0x70000, 0x20000)
            if "-sW" in args:
                return _read("good_readelf_syms.txt")
            return _read("good_readelf_secs.txt")

        self.ge._run = _run
        rc, out = self._run_gate("holosphere")
        self.assertEqual(rc, tg.EXIT_UNCALIBRATED_PASS, msg=out)
        self.assertIn("using `size -A` fallback", out)
        self.assertIn("ADVISORY", out)
        self.assertIn("UNCALIBRATED", out)

    def test_size_a_fallback_on_a_component_env_exits_2(self):
        # A missing teensy_size is a tooling break for phantasm too: its ram1
        # code ceiling has no measurement under the fallback, so the gate must
        # not red-flag it as a renamed field or a code-size regression.
        self.ge._find_teensy_size = lambda env: None

        def _run(args, check=True):
            if "-A" in args:
                return _size_a(0x10000, 0x40000, 0x70000, 0x20000)
            if "-sW" in args:
                return _read("good_readelf_syms.txt")
            return _read("good_readelf_secs.txt")

        self.ge._run = _run
        rc, out = self._run_gate("phantasm")
        self.assertEqual(rc, 2, msg=out)
        self.assertIn("per-component ceilings", out)
        self.assertNotIn("component-missing", out)

    def test_size_a_fallback_rejects_invalid_output_as_tooling_error(self):
        self.ge._find_teensy_size = lambda env: None
        for name, text in _invalid_size_a_cases().items():
            with self.subTest(name=name):
                self.ge._run = lambda *args, output=text, **kw: output
                rc, out = self._run_gate("holosphere")
                self.assertEqual(rc, 2, msg=out)
                self.assertIn("invalid `size -A` output", out)
                self.assertIn("tooling/format error", out)

    def _env(self, package_dir=None):
        """Stub SCons env; `package_dir` adds the PlatformIO platform object."""
        class _Env(dict):
            def subst(self, s):
                return self.get(s.lstrip("$"), s)
        env = _Env(PIOENV="holosphere", CC="/opt/arm/bin/arm-none-eabi-gcc")
        if package_dir is not None:
            env.PioPlatform = lambda: types.SimpleNamespace(
                get_package_dir=lambda name: package_dir)
        return env

    def _run_gate(self, pioenv):
        env = self._env()
        env["PIOENV"] = pioenv
        target = [str(TOOLS / "build" / "firmware.elf")]
        buf = io.StringIO()
        rc = 0
        try:
            with contextlib.redirect_stdout(buf):
                self.ge.run_gate(None, target, env)
        except SystemExit as exc:
            rc = exc.code
        return rc, buf.getvalue()


@contextlib.contextmanager
def mock_patch(obj, attr, value):
    """Temporarily set obj.attr = value (stdlib-only stand-in for mock.patch)."""
    orig = getattr(obj, attr)
    setattr(obj, attr, value)
    try:
        yield
    finally:
        setattr(obj, attr, orig)


if __name__ == "__main__":
    unittest.main()
