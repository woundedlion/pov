#!/usr/bin/env python3
"""Host tests for the Teensy 4 zero-warning gate.

No ARM toolchain or PlatformIO is required.

Run:  python -m unittest discover -s tools/teensy_gate_tests
"""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import teensy_warnings as tw      # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures"
REAL_DIR = FIX / "real"

# Verbatim toolchain output the synthetic fixtures cannot stand in for. All are
# committed, so a missing one is a deleted fixture, not an optional capture.
REAL_CAPTURES = (
    "cold_env_section.txt",
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


class TestWarningGate(unittest.TestCase):
    def test_toolchain_warning_uses_innermost_first_party_inline_frame(self):
        context = (
            "In file included from ./effects/Voronoi.h:3:\n"
            "    inlined from 'nearest' at ./core/spatial/kd_tree.h:182:14,\n"
            "    inlined from 'classify' at ./effects/Voronoi.h:295:28:\n")
        toolchain = "/x/.platformio/packages/toolchain-gccarmnoneeabi-teensy/include/bits/stl_algo.h"
        diagnostic = ":1818:32: warning: array subscript 16 is outside array bounds [-Warray-bounds=]"
        self.assertEqual(tw.extract_warnings(context + toolchain + diagnostic), {
            "core/spatial/kd_tree.h: warning: array subscript 16 is outside array bounds [-Warray-bounds=]"})
        framework = "/x/.platformio/packages/framework-arduinoteensy/cores/teensy4/usb.c"
        self.assertEqual(tw.extract_warnings(context + framework + diagnostic), set())
        self.assertEqual(tw.extract_warnings(context + framework + diagnostic + "\n"
                                            + toolchain + diagnostic), set())

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

    def test_checkout_parent_names_do_not_hide_first_party_warnings(self):
        for parent in ("/work/packages/repo", "/home/user/lib/Holosphere"):
            warning = parent + "/hardware/dma_led.h:42:5: warning: unused variable"
            self.assertEqual(tw.normalize(warning),
                             "hardware/dma_led.h: warning: unused variable")

    def test_fileless_diagnostics_are_normalized(self):
        lines = [
            '<command-line>: warning: "HS_PROFILE" redefined',
            "cc1plus.exe: warning: command-line option '-Wmissing-prototypes' is valid for C",
            "/opt/toolchain/bin/ld: warning: firmware.elf has a LOAD segment with RWX permissions",
            r"C:\toolchain\bin\ld.exe: warning: firmware.elf has a LOAD segment with RWX permissions",
        ]
        self.assertEqual(tw.extract_warnings("\n".join(lines)), {
            '<command-line>: warning: "HS_PROFILE" redefined',
            "cc1plus: warning: command-line option '-Wmissing-prototypes' is valid for C",
            "ld: warning: firmware.elf has a LOAD segment with RWX permissions",
        })

        self.assertIsNone(tw.normalize("cc1.exe: warning: C++ option ignored"))

    def test_library_warning_excluded(self):
        self.assertIsNone(tw.normalize(
            "/root/.platformio/lib/FastLED/FastLED.h:9:1: warning: foo [-Wbar]"))

    def test_library_path_with_first_party_segment_excluded(self):
        # A vendored lib whose nested dir reuses a first-party name (effects/)
        # remains excluded from the first-party warning set.
        self.assertIsNone(tw.normalize(
            "/root/.platformio/lib/SomeLib/effects/reverb.h:5:1: warning: w [-Wx]"))
        self.assertIsNone(tw.normalize(
            "/x/.pio/libdeps/teensy40/Foo/lib/core/bar.h:2:1: warning: w [-Wy]"))
        self.assertIsNone(tw.normalize(
            "/x/.pio/libdeps/teensy40/Foo/src/effects/baz.h:2:1: warning: w [-Wz]"))
        # No toolchain root: a third-party dir under the first-party anchor
        # carries the exclusion.
        self.assertIsNone(tw.normalize(
            "/home/runner/work/Holosphere/Holosphere/effects/lib/Foo/x.h:2:1: "
            "warning: w [-Wx]"))

    def test_nested_paths_do_not_alias_to_one_key(self):
        # A nested targets/.../effects/Foo.h and a top-level effects/Foo.h are
        # distinct files in the normalized warning set.
        root = "/home/runner/work/Holosphere/Holosphere/"
        nested = tw.normalize(
            root + "targets/Phantasm/effects/Foo.h:7:1: warning: w [-Wx]")
        top = tw.normalize(root + "effects/Foo.h:7:1: warning: w [-Wx]")
        self.assertTrue(nested.startswith("targets/Phantasm/effects/Foo.h:"))
        self.assertTrue(top.startswith("effects/Foo.h:"))
        self.assertNotEqual(nested, top)

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


def _run_warning_gate(log_text, *extra, envs=None):
    """Run the zero-warning gate over `log_text`; return its exit.

    `envs` is the environment set the build was asked to produce, written to a
    throwaway platformio.ini. It defaults to the environments `log_text` itself
    contains, which isolates a test from the repo's real environment list.
    """
    with tempfile.TemporaryDirectory() as d:
        log = Path(d) / "build.log"
        log.write_text(log_text, encoding="utf-8")
        if envs is None:
            envs = [s.name for s in tw.parse_env_sections(log_text)]
        ini = Path(d) / "platformio.ini"
        ini.write_text("".join(f"[env:{e}]\n" for e in envs), encoding="utf-8")
        return tw.main(["--build-log", str(log),
                        "--platformio-ini", str(ini), *extra])


def _banner(env, *sources):
    """A `Processing <env> (...)` banner declaring exactly these first-party TUs."""
    terms = ", ".join(["-<*>"] + [f"+<{s}>" for s in sources])
    return (f"Processing {env} (board: teensy40; build_src_filter: {terms}; "
            f"platform: teensy@5.2.0; framework: arduino)")


class TestWarningGateCaptureEvidence(unittest.TestCase):
    """A broken capture must not read as today's healthy green."""

    PIO_LINE = "Compiling .pio/build/phantasm/targets/Phantasm/Phantasm.ino.cpp.o"
    PIO_BANNER = _banner("phantasm", "targets/Phantasm/Phantasm.ino.cpp")
    VERBOSE_LINE = ("arm-none-eabi-g++ -o .pio/build/phantasm/core/memory.cpp.o "
                    "-c core/memory.cpp")
    # `pio run -v`: `-c` is a bare flag, the source is the LAST argument.
    VERBOSE_SOURCE_LAST = (
        "arm-none-eabi-g++ -o .pio/build/phantasm/src/core/memory.cpp.o "
        "-c -std=gnu++20 -fno-exceptions -O3 -DPLATFORMIO=60119 "
        "-I. -Icore -Ieffects -Ihardware core/memory.cpp")

    def _run(self, log_text):
        return _run_warning_gate(log_text)

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
        banner = _banner("phantasm", "core/memory.cpp")
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
                ".pio/build/phantasm/src/core/memory.cpp.o")
        self.assertEqual(tw.count_first_party_compiles(line), 0)

    def test_preprocess_only_invocation_is_not_evidence(self):
        # PlatformIO's .ino -> .cpp preprocess pass runs -E, not -c.
        line = ('arm-none-eabi-g++ -o "/w/pov/targets/Phantasm/Phantasm.ino.cpp" '
                '-x c++ -fpreprocessed -dD -E "/tmp/tmpm5mgtz1h"')
        self.assertEqual(tw.count_first_party_compiles(line), 0)


class TestColdCaptureAudit(unittest.TestCase):
    """A partially cached build must FAIL, not pass on a shrunken warning set.

    Every expected first-party translation unit must compile. The set is derived from
    `build_src_filter` in PlatformIO's own banner.
    """

    TUS = ("core/memory.cpp", "core/engine/static_storage.cpp",
           "core/spatial/reaction_graph.cpp",
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
        self.assertEqual(_run_warning_gate(log), 0)

    def test_short_count_from_the_object_cache_fails(self):
        # Cached shared core TUs do not satisfy the per-environment cold-build gate.
        log = self._log(("holosphere", "phantasm", "profile"),
                        self.TUS[2:], cached=self.TUS[:2])
        self.assertEqual(_run_warning_gate(log), 1)

    def test_short_count_in_a_single_environment_fails(self):
        # One fully-cold env cannot vouch for another: the audit is per-env.
        log = (self._log(("holosphere",), self.TUS)
               + self._log(("phantasm",), self.TUS[:2]))
        self.assertEqual(_run_warning_gate(log), 1)

    def test_short_count_diagnostic_names_the_cache_and_the_missing_tus(self):
        log = self._log(("phantasm",), self.TUS[2:], cached=self.TUS[:2])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(_run_warning_gate(log, "--github"), 1)
        out = buf.getvalue()
        self.assertIn("::error::", out)
        self.assertIn("2 of 4 first-party translation unit(s)", out)
        self.assertIn("object cache", out)
        for tu in self.TUS[:2]:
            self.assertIn(tu, out)

    def test_incremental_build_with_no_cache_hits_fails_and_says_so(self):
        log = self._log(("phantasm",), self.TUS[2:])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(_run_warning_gate(log), 1)
        self.assertIn("incremental build", buf.getvalue())

    def test_a_new_tu_raises_the_bar_with_no_second_edit(self):
        # Adding a TU to build_src_filter moves the expectation by itself.
        grown = self.TUS + ("core/engine/newtu.cpp",)
        banner = _banner("phantasm", *grown)
        compiles = "\n".join(self._compile("phantasm", s) for s in self.TUS)
        self.assertEqual(_run_warning_gate(banner + "\n" + compiles + "\n"), 1)
        self.assertEqual(_run_warning_gate(
            banner + "\n" + compiles + "\n"
            + self._compile("phantasm", "core/engine/newtu.cpp") + "\n"), 0)

    def test_pio_step_line_shape_also_satisfies_the_expectation(self):
        # Without -v PlatformIO names the OBJECT; the `.o` suffix must not stop it
        # matching the source build_src_filter declares.
        log = (_banner("phantasm", *self.TUS) + "\n"
               + "\n".join(f"Compiling .pio/build/phantasm/src/{s}.o"
                           for s in self.TUS) + "\n")
        self.assertEqual(_run_warning_gate(log), 0)

    def test_missing_banner_fails_rather_than_falling_back(self):
        # If PlatformIO's banner format ever changes the expectation cannot be
        # derived; that must go red, not silently revert to "at least one compile".
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = _run_warning_gate(self._compile("phantasm", self.TUS[0]) + "\n",
                                   envs=("phantasm",))
        self.assertEqual(rc, 1)
        self.assertIn("no `Processing <env> (...)` banner", buf.getvalue())

    def test_banner_without_build_src_filter_fails(self):
        log = ("Processing phantasm (board: teensy40; platform: teensy@5.2.0)\n"
               + self._compile("phantasm", self.TUS[0]) + "\n")
        self.assertEqual(_run_warning_gate(log), 1)

    def test_first_party_glob_in_the_filter_fails_loudly(self):
        # A glob is not countable from the log; the tool must say so, not guess.
        section, = tw.parse_env_sections(_banner("phantasm", "core/engine/*.cpp"))
        with self.assertRaises(tw.CaptureError):
            tw.declared_first_party_sources(section)

    def test_third_party_filter_terms_are_not_expected_tus(self):
        section, = tw.parse_env_sections(
            _banner("phantasm", "core/memory.cpp", "lib/Foo/foo.cpp"))
        self.assertEqual(tw.declared_first_party_sources(section),
                         {"core/memory.cpp"})

    def test_exclusion_term_removes_a_declared_tu(self):
        header = _banner("phantasm", *self.TUS).replace(
            "; platform:", " -<core/memory.cpp>; platform:")
        section, = tw.parse_env_sections(header)
        self.assertEqual(tw.declared_first_party_sources(section), set(self.TUS[1:]))


class TestExpectedEnvironmentSet(unittest.TestCase):
    """The audited environments must be the ones the build was asked to produce.

    Sections come from banners, so a `pio run` over six environments that dies in
    the first prints ONE banner: the other five are absent rather than short, and
    the per-environment coldness audit has nothing to complain about.
    """

    TU = "core/memory.cpp"
    ENVS = ("holosphere", "holosphere_dma", "phantasm", "phantasm8",
            "profile", "profile_o3")

    def _cold_env(self, env):
        return (_banner(env, self.TU) + "\n"
                + f"arm-none-eabi-g++ -o .pio/build/{env}/src/{self.TU}.o -c "
                  f"-std=gnu++20 -O3 -I. -Icore {self.TU}\n")

    def test_truncated_run_fails(self):
        self.assertEqual(
            _run_warning_gate(self._cold_env(self.ENVS[0]), envs=self.ENVS), 1)

    def test_truncated_run_diagnostic_names_the_absent_environments(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(_run_warning_gate(self._cold_env(self.ENVS[0]),
                                          "--github", envs=self.ENVS), 1)
        out = buf.getvalue()
        self.assertIn("::error::", out)
        self.assertIn("5 of 6 expected environment(s)", out)
        for env in self.ENVS[1:]:
            self.assertIn(env, out)

    def test_every_expected_environment_present_passes(self):
        log = "".join(self._cold_env(e) for e in self.ENVS)
        self.assertEqual(_run_warning_gate(log, envs=self.ENVS), 0)

    def test_env_flag_narrows_the_expectation(self):
        # A deliberate subset build states its own set instead of platformio.ini's.
        self.assertEqual(
            _run_warning_gate(self._cold_env(self.ENVS[0]), "--env", self.ENVS[0],
                         envs=self.ENVS), 0)

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
    """Historical `pio run -v` sections from before static_storage.cpp split out.

    fixtures/real/{cold,warm}_env_section.txt are the `holosphere` sections of two
    consecutive runs of the CI command: the first with `.pio/build_cache` deleted,
    the second reusing it. Only a real capture pins PlatformIO's banner text and
    SCons's `Retrieved … from cache` line, which the derived expectation reads.
    """

    COLD = (REAL_DIR / "cold_env_section.txt").read_text(encoding="utf-8")
    WARM = (REAL_DIR / "warm_env_section.txt").read_text(encoding="utf-8")
    TUS = {"core/engine/memory.cpp", "core/spatial/reaction_graph.cpp",
           "targets/Holosphere/Holosphere.ino.cpp"}

    def test_historical_banner_declares_three_first_party_tus(self):
        section, = tw.parse_env_sections(self.COLD)
        self.assertEqual(section.name, "holosphere")
        self.assertEqual(tw.declared_first_party_sources(section), self.TUS)

    def test_real_cold_section_compiles_every_declared_tu(self):
        audit = tw.audit_capture(self.COLD)
        self.assertEqual(audit.missing, ())
        self.assertEqual(audit.cache_hits, 0)
        self.assertEqual(_run_warning_gate(self.COLD), 0)

    def test_real_warm_section_is_short_and_fails(self):
        audit = tw.audit_capture(self.WARM)
        self.assertEqual(audit.missing_count, 2)
        self.assertEqual(audit.first_party_cache_hits, 2)
        self.assertEqual(_run_warning_gate(self.WARM), 1)


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


class TestNonUtf8Captures(unittest.TestCase):
    """The warning gate answers by exit code, and a decode error has none.

    A Windows `pio run -v 2>&1 | tee` interleaves cp1252 bytes into the stream,
    so a capture that is not valid UTF-8 must still produce a verdict.
    """

    # RIGHT SINGLE QUOTATION MARK in cp1252; not a valid UTF-8 sequence.
    CP1252 = b"don\x92t"

    def test_warning_gate_reads_a_build_log_with_a_cp1252_byte(self):
        log_text = (_banner("phantasm", "core/memory.cpp") + "\n"
                    "arm-none-eabi-g++ -o .pio/build/phantasm/core/engine/"
                    "memory.cpp.o -c core/memory.cpp\n")
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "build.log"
            log.write_bytes(log_text.encode("utf-8")
                            + b"note: " + self.CP1252 + b" care\n")
            ini = Path(d) / "platformio.ini"
            ini.write_text("[env:phantasm]\n", encoding="utf-8")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = tw.main(["--build-log", str(log),
                              "--platformio-ini", str(ini)])
            self.assertEqual(rc, 0, msg=buf.getvalue())
