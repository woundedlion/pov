#!/usr/bin/env python3
"""Host tests for the Teensy 4 size/layout gate (docs/teensy_ci_gate_spec.md §9.1).

The gate must gate itself: a size/layout check that can never fail is worse than
none (permanent false-green). So every layout invariant and region ceiling is
proven to FAIL on a deliberately-broken fixture and PASS on the good one. These
are plain stdlib `unittest` tests — no ARM toolchain, no PlatformIO — so they run
in milliseconds in the existing CI test lane.

Run:  python -m unittest discover -s tools/teensy_gate_tests
"""

import contextlib
import io
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


def _read(name):
    return (FIX / name).read_text(encoding="utf-8")


def _load_gate_extra():
    """Load tools/teensy_gate_extra.py (SCons glue) outside PlatformIO.

    The module runs `Import("env")` and `env.subst(...)` at import time, so it
    cannot be plain-imported; inject no-op SCons hooks and a stub env, then exec.
    """
    class _Env:
        def subst(self, s):
            return str(TOOLS.parent) if s == "$PROJECT_DIR" else s

        def AddPostAction(self, *a, **k):
            pass

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
            "teensy_size: RAM1: variables:343040, code:62240, padding:30496 free for local variables: 88512\n"
            "teensy_size: RAM2: variables:497920 free for malloc/new: 26368\n"
        )
        sizes = tg.parse_teensy_size(text)
        self.assertEqual(set(sizes), {"flash", "ram1", "ram2"})
        self.assertEqual(sizes["flash"]["used"], 62788 + 13684 + 8460)
        self.assertEqual(sizes["flash"]["free"], 1979136)
        self.assertEqual(sizes["ram1"]["used"], 343040 + 62240 + 30496)
        self.assertEqual(sizes["ram1"]["free"], 88512)
        self.assertEqual(sizes["ram2"]["used"], 497920)
        self.assertEqual(sizes["ram2"]["free"], 26368)

    def test_parse_readelf_symbols_uses_real_mangled_names(self):
        syms = {s.name: s for s in tg.parse_readelf_symbols(_read("good_readelf_syms.txt"))}
        self.assertIn("_ZN6Effect8buffer_aE", syms)              # class static, mangled
        self.assertIn("_ZL18global_arena_block", syms)           # file-scope static -> _ZL mangled LOCAL
        self.assertEqual(syms["_ZL18global_arena_block"].bind, "LOCAL")
        self.assertEqual(syms["_ZN6Effect8buffer_aE"].value, 0x20200000)
        self.assertEqual(syms["_ZN6Effect8buffer_aE"].region, "OCRAM")
        self.assertEqual(syms["_ZL18global_arena_block"].size, 343040)

    def test_parse_size_a_bucketing(self):
        totals = tg.region_totals_from_size_a(_read("good_size_a.txt"))
        # .ARM.attributes/.comment sit at VMA 0 but are non-allocated: they must be
        # dropped, not bucketed into ITCM (whose base is also 0x0).
        self.assertEqual(totals["ITCM"], 0xf320 + 0x2a00)   # .text.itcm + .text.code
        self.assertEqual(totals["FLASH"], 0x2cc00 + 0x8)    # .text.progmem + .ARM.exidx
        self.assertEqual(totals["DTCM"], 0x1a00 + 0x53c00)
        self.assertEqual(totals["OCRAM"], 0x79900)


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
        # table is present + budgeted (Holosphere GC's it out, §8).
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
    def test_over_cap_trips_every_region(self):
        # Good symbols (layout fine) + over-cap totals: only region checks fire.
        result = _eval("holosphere", "broken_over_cap_teensy_size.txt",
                       "good_readelf_syms.txt")
        self.assertFalse(result.passed)
        codes = _codes(result)
        self.assertEqual(codes.count("region-over-budget"), 2)   # FLASH + RAM2
        self.assertIn("headroom-below-floor", codes)             # DTCM stack room


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
        # Same set, different order -> no new warnings (set-based, §7.2).
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


class TestRealCapture(unittest.TestCase):
    """Parse REAL arm-gcc 11.3.1 output from a fitting Holosphere build (§9.1).

    These exercise the toolchain quirks the synthetic fixtures can't: teensy_size
    on stderr, readelf printing large sizes in HEX (0x53c00), and the arena's real
    _ZL-mangled internal-linkage name. Captured by the Phase-0 spike.
    """

    @unittest.skipUnless((REAL_DIR / "holosphere_readelf_syms.txt").exists(),
                         "real captures not present")
    def test_real_holosphere_build_passes_the_calibrated_gate(self):
        sizes = tg.parse_teensy_size((REAL_DIR / "holosphere_teensy_size.txt").read_text())
        syms = tg.parse_readelf_symbols((REAL_DIR / "holosphere_readelf_syms.txt").read_text())
        secs = tg.parse_readelf_sections((REAL_DIR / "holosphere_readelf_secs.txt").read_text())
        result = tg.evaluate("holosphere", BUDGETS["holosphere"], sizes, syms, secs)
        self.assertTrue(result.passed, msg=_codes(result))
        arena = next(s for s in syms if s.name == "_ZL18global_arena_block")
        self.assertEqual(arena.region, "DTCM")
        self.assertEqual(arena.size, 343040)  # 0x53c00, parsed from the hex Size column

    @unittest.skipUnless((REAL_DIR / "phantasm_readelf_syms.txt").exists(),
                         "real captures not present")
    def test_real_phantasm_reaction_graph_is_in_flash(self):
        syms = tg.parse_readelf_symbols((REAL_DIR / "phantasm_readelf_syms.txt").read_text())
        rg = next(s for s in syms if s.name == "_ZN13ReactionGraph9neighborsE")
        self.assertEqual(rg.region, "FLASH")
        self.assertEqual(rg.size, 92160)

    @unittest.skipUnless((REAL_DIR / "phantasm_readelf_syms.txt").exists(),
                         "real captures not present")
    def test_real_phantasm_dma_tx_buffer_is_in_ocram(self):
        syms = tg.parse_readelf_symbols((REAL_DIR / "phantasm_readelf_syms.txt").read_text())
        led = next(s for s in syms
                   if s.name == "_ZN12POVSegmentedILi288ELi4ELi480EE14ledController_E")
        self.assertEqual(led.region, "OCRAM")


class TestSizeAFallback(unittest.TestCase):
    """The `--size-a` fallback path (no teensy_size): VMA bucketing + 0x80000-ram1
    free-headroom arithmetic, end-to-end through main() and evaluate()."""

    def _size_a(self, itcm, dtcm, ocram, flash):
        # One section per region at its base VMA, sizes in hex.
        return (
            "elf:\nsection size addr\n"
            f".text.itcm 0x{itcm:x} 0x0\n"
            f".bss 0x{dtcm:x} 0x20000000\n"
            f".bss.dma 0x{ocram:x} 0x20200000\n"
            f".text.progmem 0x{flash:x} 0x60000000\n")

    def test_free_headroom_is_0x80000_minus_used(self):
        totals = tg.region_totals_from_size_a(_read("good_size_a.txt"))
        ram1 = totals.get("ITCM", 0) + totals.get("DTCM", 0)
        # Replicate main()'s fallback synthesis and assert the free arithmetic.
        sizes = {
            "ram1": {"used": ram1, "free": 0x80000 - ram1},
            "ram2": {"used": totals.get("OCRAM", 0),
                     "free": 0x80000 - totals.get("OCRAM", 0)},
        }
        self.assertEqual(sizes["ram1"]["free"], 0x80000 - ram1)
        self.assertEqual(sizes["ram2"]["free"], 0x80000 - totals["OCRAM"])

    def test_main_size_a_fallback_passes_a_fitting_build(self):
        rc, out = self._run_main_size_a(
            self._size_a(0x10000, 0x40000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "holosphere")
        self.assertEqual(rc, 0, msg=out)
        self.assertIn("PASS", out)

    def test_main_size_a_fallback_free_arithmetic_trips_headroom_floor(self):
        # DTCM large enough that 0x80000 - ram1 drops below the 32,768 B floor,
        # proving the fallback `free` figure actually drives the gate decision.
        # ram1 = 0x10000 + 0x70000 = 0x80000 -> free 0; floor 32768 -> violation.
        rc, out = self._run_main_size_a(
            self._size_a(0x10000, 0x70000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "holosphere")
        self.assertEqual(rc, 1, msg=out)
        self.assertIn("free-for-local-variables", out)

    def test_fallback_pass_is_marked_advisory_not_calibrated(self):
        # A fallback PASS must carry the advisory note so it is never mistaken
        # for a teensy_size-calibrated verdict.
        rc, out = self._run_main_size_a(
            self._size_a(0x10000, 0x40000, 0x70000, 0x20000),
            "good_readelf_syms.txt", "holosphere")
        self.assertEqual(rc, 0, msg=out)
        self.assertIn("PASS", out)
        self.assertIn("ADVISORY", out)
        self.assertIn("not calibrated", out.lower())

    def _run_main_size_a(self, size_a_text, syms_fixture, env):
        with tempfile.TemporaryDirectory() as d:
            sa = Path(d) / "size_a.txt"
            sa.write_text(size_a_text, encoding="utf-8")
            syms = Path(d) / "syms.txt"
            syms.write_text(_read(syms_fixture), encoding="utf-8")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = tg.main([
                    "--env", env,
                    "--budgets", str(TOOLS / "teensy_budgets.json"),
                    "--size-a", str(sa),
                    "--readelf-syms", str(syms),
                ])
            return rc, buf.getvalue()


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


class TestGateExtra(unittest.TestCase):
    """The PlatformIO post-build glue: toolchain discovery + exit(2) guards."""

    def setUp(self):
        self.ge = _load_gate_extra()

    def test_tool_derives_sibling_arm_tools(self):
        self.assertEqual(self.ge._tool("/opt/arm/bin/arm-none-eabi-gcc", "size"),
                         "/opt/arm/bin/arm-none-eabi-size")
        self.assertEqual(self.ge._tool("C:/x/arm-none-eabi-gcc.exe", "readelf"),
                         "C:/x/arm-none-eabi-readelf.exe")

    def test_tool_falls_back_to_path_lookup(self):
        self.assertEqual(self.ge._tool("clang", "size"), "arm-none-eabi-size")

    def test_find_teensy_size_accepts_self_identifying(self):
        probe = types.SimpleNamespace(stdout="usage: teensy_size <elf>", stderr="")
        with mock_patch(subprocess, "run", lambda *a, **k: probe):
            self.assertEqual(self.ge._find_teensy_size(), "teensy_size")

    def test_find_teensy_size_rejects_foreign_binary(self):
        probe = types.SimpleNamespace(stdout="usage: other-tool", stderr="")
        with mock_patch(subprocess, "run", lambda *a, **k: probe):
            self.assertIsNone(self.ge._find_teensy_size())

    def test_find_teensy_size_none_when_absent(self):
        def _raise(*a, **k):
            raise OSError("not found")
        with mock_patch(subprocess, "run", _raise):
            self.assertIsNone(self.ge._find_teensy_size())

    def test_toolchain_oserror_exits_2(self):
        # A tool step raising OSError is a build/tooling break -> exit(2), never a
        # size-budget "violation".
        self.ge._find_teensy_size = lambda: None
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
        self.ge._find_teensy_size = lambda: "teensy_size"
        self.ge._run = lambda *a, **k: ""
        with mock_patch(self.ge.teensy_gate, "parse_teensy_size", lambda text: {}):
            rc, out = self._run_gate("holosphere")
        self.assertEqual(rc, 2)
        self.assertIn("parsed no FLASH/RAM1/RAM2 regions", out)

    def _run_gate(self, pioenv):
        class _Env(dict):
            def subst(self, s):
                return self.get(s.lstrip("$"), s)
        env = _Env(PIOENV=pioenv, CC="/opt/arm/bin/arm-none-eabi-gcc")
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
