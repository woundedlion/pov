#!/usr/bin/env python3
"""Host tests for the PlatformIO build hooks (platformio.ini extra_scripts).

Each hook is executed against a fake SCons construction environment and its
build effects asserted directly. No ARM toolchain or PlatformIO is required.

Run:  python -m unittest discover -s tools/teensy_hook_tests
"""

import configparser
import os
import re
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TOOLS = REPO / "tools"
PIO_INI = REPO / "platformio.ini"

sys.path.insert(0, str(TOOLS))

import teensy_gate as tg   # noqa: E402
from scons_hook_test_env import FakeEnv, load_hook  # noqa: E402

GATE_SCRIPT = "post:tools/teensy_gate_extra.py"

# Hooks every environment needs.
REQUIRED_SCRIPTS = (
    "pre:tools/teensy_pre.py",
    "pre:tools/teensy_map.py",
    "post:tools/teensy_isystem.py",
)


def _abs(*parts):
    """An absolute-looking path with this platform's separator."""
    return os.path.join(os.sep + parts[0], *parts[1:])


LIBDEPS = _abs("work", "Holosphere", ".pio", "libdeps", "phantasm", "FastLED", "src")
FRAMEWORK = _abs("home", "runner", ".platformio", "packages",
                 "framework-arduinoteensy", "cores", "teensy4")
RELOCATED = _abs("opt", "pio-core", "packages", "framework-arduinoteensy",
                 "libraries", "SPI")
FIRST_PARTY = [".", "core", "effects", "hardware"]
COMPILE_COMMANDS = {"CCCOM": "$CC -c", "CXXCOM": "$CXX -c", "ASPPCOM": "$CC -x"}


def _pio_config():
    cfg = configparser.RawConfigParser()
    cfg.read(PIO_INI, encoding="utf-8")
    return cfg


def _pio_envs():
    return [s.split(":", 1)[1] for s in _pio_config().sections()
            if s.startswith("env:")]


_INTERP = re.compile(r"\$\{([^.}]+)\.([^}]+)\}")


def _option_lines(cfg, section, option, depth=0):
    """Resolve inherited options and ${section.option} references without merging lists."""
    if depth > 8:
        raise RecursionError(f"[{section}] {option} interpolates cyclically")
    sections = [("env", depth), (section, depth)] if section.startswith("env:") else [(section, depth)]
    while sections:
        section, inherited_depth = sections.pop()
        if inherited_depth > 8:
            raise RecursionError(f"[{section}] extends cyclically")
        if cfg.has_option(section, option):
            break
        parents = cfg.get(section, "extends", fallback="")
        sections.extend((parent.strip(), inherited_depth + 1)
                        for parent in re.split(r"[,\n]", parents) if parent.strip())
    else:
        return []
    out = []
    for line in cfg.get(section, option, fallback="").splitlines():
        line = line.strip()
        if not line:
            continue
        expanded = _INTERP.sub(
            lambda ref: " ".join(_option_lines(cfg, ref.group(1), ref.group(2), depth + 1)),
            line)
        out.extend(expanded.split())
    return out


class TestOptionInterpolation(unittest.TestCase):
    def test_inline_source_references_include_shared_and_local_filters(self):
        cfg = configparser.RawConfigParser()
        cfg.read_dict({
            "env": {"build_src_filter": "-<*>\n+<core/memory.cpp>"},
            "shared": {"sources": "${env.build_src_filter} +<core/engine/static_storage.cpp>"},
            "env:fixture": {"build_src_filter": "${shared.sources} +<targets/Profile/Profile.ino.cpp>"},
        })
        self.assertEqual(_option_lines(cfg, "env:fixture", "build_src_filter"), [
            "-<*>", "+<core/memory.cpp>", "+<core/engine/static_storage.cpp>",
            "+<targets/Profile/Profile.ino.cpp>"])

    def test_inline_script_references_include_shared_and_local_hooks(self):
        cfg = configparser.RawConfigParser()
        cfg.read_dict({
            "env": {"extra_scripts": "pre:tools/teensy_pre.py\npre:tools/teensy_map.py"},
            "env:fixture": {"extra_scripts": "${env.extra_scripts} pre:tools/teensy_nano.py"},
        })
        self.assertEqual(_option_lines(cfg, "env:fixture", "extra_scripts"), [
            "pre:tools/teensy_pre.py", "pre:tools/teensy_map.py", "pre:tools/teensy_nano.py"])


class TestOptionInheritance(unittest.TestCase):
    def test_inherited_option_overrides_env_and_last_parent_wins(self):
        cfg = configparser.RawConfigParser()
        cfg.read_dict({
            "env": {"extra_scripts": "post:gate.py"},
            "first": {"extra_scripts": "pre:first.py"},
            "second": {"extra_scripts": "pre:second.py"},
            "shared": {"extends": "first, second"},
            "env:fixture": {"extends": "shared"},
        })
        self.assertEqual(_option_lines(cfg, "env:fixture", "extra_scripts"),
                         ["pre:second.py"])
        cfg.set("env:fixture", "extra_scripts", "pre:local.py")
        self.assertEqual(_option_lines(cfg, "env:fixture", "extra_scripts"),
                         ["pre:local.py"])


class TestExtraScriptsExist(unittest.TestCase):
    """Every hook platformio.ini names must be a real file."""

    def test_every_extra_script_resolves(self):
        cfg = _pio_config()
        seen = 0
        for section in cfg.sections():
            for script in _option_lines(cfg, section, "extra_scripts"):
                rel = script.split(":", 1)[-1]
                self.assertTrue((REPO / rel).is_file(),
                                f"[{section}] extra_scripts names missing {rel}")
                seen += 1
        self.assertGreater(seen, 0, "parsed no extra_scripts from platformio.ini")


class TestBudgetedEnvsWireTheGate(unittest.TestCase):
    """The size/layout gate hook is wired to exactly the budgeted envs."""

    def test_every_budgeted_env_resolves_the_gate_hook(self):
        budgets = tg.load_budgets(TOOLS / "teensy_budgets.json")
        self.assertTrue(budgets, "teensy_budgets.json names no environments")
        cfg = _pio_config()
        for name in budgets:
            with self.subTest(env=name):
                self.assertIn(f"env:{name}", cfg.sections(),
                              f"budgeted env '{name}' has no platformio.ini section")
                self.assertIn(GATE_SCRIPT,
                              _option_lines(cfg, f"env:{name}", "extra_scripts"),
                              f"budgeted env '{name}' does not wire {GATE_SCRIPT}")

    def test_gate_hook_never_runs_without_a_budget(self):
        # teensy_gate_extra.py looks its env up in teensy_budgets.json by name.
        budgets = tg.load_budgets(TOOLS / "teensy_budgets.json")
        cfg = _pio_config()
        for name in _pio_envs():
            if GATE_SCRIPT in _option_lines(cfg, f"env:{name}", "extra_scripts"):
                self.assertIn(name, budgets,
                              f"env '{name}' wires the gate but has no budget entry")


class TestRequiredHooksReachEveryEnv(unittest.TestCase):
    """Every env carries the required hooks, including envs that re-type
    extra_scripts (list options do not merge)."""

    def test_every_base_hook_but_the_gate_reaches_every_env(self):
        # Derived from [env]; the gate hook is the one line envs may drop.
        cfg = _pio_config()
        base = [s for s in _option_lines(cfg, "env", "extra_scripts")
                if s != GATE_SCRIPT]
        self.assertTrue(base, "[env] declares no extra_scripts besides the gate")
        for name in _pio_envs():
            resolved = _option_lines(cfg, f"env:{name}", "extra_scripts")
            for script in base:
                with self.subTest(env=name, script=script):
                    self.assertIn(script, resolved,
                                  f"env '{name}' re-types extra_scripts and "
                                  f"drops {script}, which [env] declares")

    def test_shipping_family_resolves_nano_specs(self):
        cfg = _pio_config()
        matched = set()
        for name in _pio_envs():
            sources = " ".join(_option_lines(cfg, f"env:{name}", "build_src_filter"))
            if any(f"targets/{target}/" in sources for target in ("Phantasm", "Profile", "Bench")):
                matched.add(name)
                with self.subTest(env=name):
                    self.assertIn("pre:tools/teensy_nano.py",
                                  _option_lines(cfg, f"env:{name}", "extra_scripts"))

        self.assertTrue({"phantasm", "phantasm8", "profile", "profile_o3", "bench"} <= matched,
                        f"shipping environments missing from source selection: {matched}")

    def test_base_env_declares_every_required_hook(self):
        declared = _option_lines(_pio_config(), "env", "extra_scripts")
        for script in REQUIRED_SCRIPTS:
            self.assertIn(script, declared, f"[env] does not declare {script}")

    def test_no_env_wires_a_hook_twice(self):
        # ${env.extra_scripts} plus a hand-listed copy double-appends its flags.
        cfg = _pio_config()
        for name in _pio_envs():
            resolved = _option_lines(cfg, f"env:{name}", "extra_scripts")
            self.assertEqual(len(resolved), len(set(resolved)),
                             f"env '{name}' lists a hook more than once: {resolved}")


class TestHooksRegisterNoBuildActions(unittest.TestCase):
    """The non-gate hooks only mutate build flags, compile commands and sketch
    handling."""

    def test_firmware_hooks_register_no_build_actions(self):
        cfg = _pio_config()
        for name in _pio_envs():
            with self.subTest(env=name):
                projenv = FakeEnv(CPPPATH=list(FIRST_PARTY))
                env = FakeEnv(
                    PIOENV=name, PROJECT_DIR=str(REPO),
                    BUILD_DIR=os.path.join(".pio", "build", name),
                    CPPPATH=FIRST_PARTY + [FRAMEWORK],
                    PROJECT_PACKAGES_DIR=_abs("home", "runner", ".platformio", "packages"),
                    PROJECT_LIBDEPS_DIR=_abs("work", "Holosphere", ".pio", "libdeps"),
                    ENV={"PATH": ""}, **COMPILE_COMMANDS)
                hooks = [s.split(":", 1)[1]
                         for s in _option_lines(cfg, f"env:{name}", "extra_scripts")
                         if s != GATE_SCRIPT]
                self.assertTrue(hooks)
                for rel in hooks:
                    self.assertEqual(Path(rel).parent, Path("tools"), rel)
                    load_hook(Path(rel).name, env=env, projenv=projenv,
                              COMMAND_LINE_TARGETS=[])
                for build_env in (env, projenv):
                    self.assertEqual(build_env.pre_actions, [])
                    self.assertEqual(build_env.post_actions, [])
                self.assertEqual(set(env.methods), {"FindInoNodes", "ConvertInoToCpp"})


class TestSketchSelection(unittest.TestCase):
    """teensy_pre.py selects each env's sketch."""

    def _run(self, pioenv):
        env = FakeEnv(PIOENV=pioenv, PROJECT_DIR=str(REPO))
        return load_hook("teensy_pre.py", env=env), env

    def test_cpp_only_flag_is_removed_before_shared_flag_processing(self):
        for flags in (["-O3 -fno-use-cxa-atexit", "-D FLAG=-fno-use-cxa-atexit"],
                      "-O3 -fno-use-cxa-atexit -fno-exceptions"):
            with self.subTest(flags=flags):
                env = FakeEnv(PIOENV="holosphere", PROJECT_DIR=str(REPO),
                              BUILD_FLAGS=flags, CXXFLAGS=["-std=gnu++20"])
                load_hook("teensy_pre.py", env=env)
                shared = " ".join(env["BUILD_FLAGS"]).split()
                self.assertNotIn("-fno-use-cxa-atexit", shared)
                self.assertIn("-O3", shared)
                self.assertEqual(env["CXXFLAGS"], ["-std=gnu++20", "-fno-use-cxa-atexit"])
                if isinstance(flags, list):
                    self.assertIn("FLAG=-fno-use-cxa-atexit", shared)

    def test_mapped_sketches_exist(self):
        mod, _ = self._run("phantasm")
        for name, rel in mod.SKETCH.items():
            self.assertTrue((REPO / rel).is_file(), f"{name} -> missing {rel}")

    def test_find_ino_nodes_returns_only_this_envs_sketch(self):
        for pioenv, leaf in (("holosphere", "Holosphere.ino"),
                             ("phantasm8", "Phantasm.ino"),
                             ("profile_o3", "Profile.ino")):
            with self.subTest(env=pioenv):
                _, env = self._run(pioenv)
                nodes = env.methods["FindInoNodes"](env)
                self.assertEqual(len(nodes), 1)
                self.assertEqual(os.path.basename(nodes[0]), leaf)
                self.assertTrue(os.path.isabs(nodes[0]))

    def test_each_mapping_names_the_sketch_that_env_compiles(self):
        # A mapping to the other target's sketch links the wrong setup()/loop()
        # and still builds green.
        mod, _ = self._run("phantasm")
        cfg = _pio_config()
        for name in _pio_envs():
            with self.subTest(env=name):
                compiled = re.findall(
                    r"\+<([^>]+\.ino)\.cpp>",
                    "\n".join(_option_lines(cfg, f"env:{name}",
                                            "build_src_filter")))
                self.assertEqual(len(compiled), 1, compiled)
                self.assertEqual(mod.SKETCH[name].replace(os.sep, "/"),
                                 compiled[0])

    def test_unknown_env_fails_loud(self):
        with self.assertRaises(SystemExit):
            self._run("not_an_env")


class TestSketchConversion(unittest.TestCase):
    """teensy_pre.py: the converted sketch is identical across checkouts."""

    @staticmethod
    def _convert(project, temp_name):
        sketch = f"{Path(project).as_posix()}/targets/Phantasm/Phantasm.ino"
        os.makedirs(os.path.join(project, "targets", "Phantasm"))
        env = FakeEnv(PIOENV="phantasm", PROJECT_DIR=project, CONVERTED_SKETCH=(
            f'# 1 "C:\\\\Temp\\\\{temp_name}"\r\n'
            "#include <Arduino.h>\r\n"
            f'# 1 "{sketch}"\r\n'
            f'# 23 "{sketch}"\r\n'
            "void loop();\r\n"
            f'#line 30 "{sketch}"\r\n'
            "void loop() {}\r\n"))
        load_hook("teensy_pre.py", env=env)
        env.methods["ConvertInoToCpp"](env)
        with open(sketch + ".cpp", encoding="utf-8", newline="") as cpp:
            return cpp.read()

    def test_markers_name_the_sketch_by_relative_path_without_the_temp_file(self):
        with tempfile.TemporaryDirectory() as first, \
                tempfile.TemporaryDirectory() as second:
            converted = self._convert(first, "tmpa1b2")
            self.assertEqual(converted, self._convert(second, "tmpz9y8"))
        self.assertEqual(converted, (
            "#include <Arduino.h>\r\n"
            '# 1 "targets/Phantasm/Phantasm.ino"\r\n'
            '# 23 "targets/Phantasm/Phantasm.ino"\r\n'
            "void loop();\r\n"
            '#line 30 "targets/Phantasm/Phantasm.ino"\r\n'
            "void loop() {}\r\n"))


class TestCcacheLauncher(unittest.TestCase):
    """pio_ccache.py: compiles run through ccache only when it is on PATH."""

    def _run(self, path, targets=()):
        env = FakeEnv(ENV={"PATH": path}, LINKCOM="$LINK", ASCOM="$AS",
                      **COMPILE_COMMANDS)
        load_hook("pio_ccache.py", env=env, COMMAND_LINE_TARGETS=list(targets))
        return env

    def _with_ccache(self, targets=()):
        with tempfile.TemporaryDirectory() as bin_dir:
            for name in ("ccache", "ccache.exe"):
                tool = Path(bin_dir) / name
                tool.write_bytes(b"")
                tool.chmod(0o755)
            return self._run(bin_dir, targets)

    def test_compiles_are_prefixed_outside_the_signature(self):
        env = self._with_ccache()
        for com, command in COMPILE_COMMANDS.items():
            self.assertEqual(env[com], "$( ccache $) " + command)
        self.assertEqual((env["LINKCOM"], env["ASCOM"]), ("$LINK", "$AS"))

    def test_commands_are_unchanged_without_ccache(self):
        with tempfile.TemporaryDirectory() as empty:
            env = self._run(empty)
        for com, command in COMPILE_COMMANDS.items():
            self.assertEqual(env[com], command)

    def test_compile_database_keeps_plain_commands(self):
        env = self._with_ccache(targets=["compiledb"])
        for com, command in COMPILE_COMMANDS.items():
            self.assertEqual(env[com], command)


class TestLinkerMap(unittest.TestCase):
    """teensy_map.py: the per-symbol source the size investigations read."""

    def test_map_lands_in_this_envs_build_dir(self):
        env = FakeEnv(BUILD_DIR=os.path.join(".pio", "build", "phantasm"))
        load_hook("teensy_map.py", env=env)
        self.assertEqual(
            env["LINKFLAGS"],
            ["-Wl,-Map," + os.path.join(".pio", "build", "phantasm") + "/firmware.map"])


class TestNanoSpecs(unittest.TestCase):
    """teensy_nano.py: --specs=nano.specs must reach the LINK step, and appear
    exactly once per command (gcc fatals on a repeated spec)."""

    def test_reaches_compile_and_link_exactly_once(self):
        env = FakeEnv(CCFLAGS=["-O3"], CXXFLAGS=["-std=gnu++20"], LINKFLAGS=["-Os"])
        load_hook("teensy_nano.py", env=env)
        self.assertEqual(env["CCFLAGS"].count("--specs=nano.specs"), 1)
        self.assertEqual(env["LINKFLAGS"].count("--specs=nano.specs"), 1)

    def test_not_added_to_cxxflags(self):
        # A C++ command is `$CXXFLAGS $CCFLAGS`; a spec in both makes gcc fatal.
        env = FakeEnv(CCFLAGS=[], CXXFLAGS=[], LINKFLAGS=[])
        load_hook("teensy_nano.py", env=env)
        self.assertNotIn("--specs=nano.specs", env["CXXFLAGS"])


class TestIsystemDemotion(unittest.TestCase):
    """teensy_isystem.py: vendored include dirs must move from -I to -isystem."""

    def _run(self, cpppath, lib_builders=0, packages=None):
        projenv = FakeEnv(CPPPATH=list(cpppath))
        env = FakeEnv(
            CPPPATH=list(cpppath),
            PROJECT_PACKAGES_DIR=packages or _abs("home", "runner", ".platformio", "packages"),
            PROJECT_LIBDEPS_DIR=_abs("work", "Holosphere", ".pio", "libdeps"))
        env.lib_builders = [types.SimpleNamespace(env=FakeEnv())
                            for _ in range(lib_builders)]
        load_hook("teensy_isystem.py", env=env, projenv=projenv)
        return projenv, env

    def test_vendored_dirs_move_out_of_cpppath(self):
        for packages, include in [
                (_abs("home", "runner", ".platformio", "packages"), FRAMEWORK),
                (_abs("opt", "pio-core", "packages"), RELOCATED)]:
            with self.subTest(packages=packages):
                projenv, env = self._run(FIRST_PARTY + [LIBDEPS, include], packages=packages)
                for build_env in (projenv, env):
                    self.assertEqual(build_env["CPPPATH"], FIRST_PARTY)
                    self.assertEqual(build_env["CCFLAGS"],
                                     ["-isystem", LIBDEPS, "-isystem", include])

    def test_first_party_dirs_are_never_demoted(self):
        _, env = self._run(FIRST_PARTY + [FRAMEWORK])
        for path in FIRST_PARTY:
            self.assertIn(path, env["CPPPATH"])
            self.assertNotIn(path, env["CCFLAGS"])

    def test_library_builders_get_warnings_disabled(self):
        # -isystem cannot silence a warning in the body of the file being
        # compiled, only in its headers.
        _, env = self._run(FIRST_PARTY + [FRAMEWORK], lib_builders=3)
        for builder in env.GetLibBuilders():
            self.assertEqual(builder.env["CCFLAGS"], ["-w"])

    def test_demoting_nothing_fails_loud(self):
        with self.assertRaises(SystemExit) as cm:
            self._run(FIRST_PARTY)
        self.assertIn("demoted 0 third-party include dirs", str(cm.exception))

    def test_one_env_carrying_the_vendored_dirs_is_enough(self):
        projenv = FakeEnv(CPPPATH=list(FIRST_PARTY))
        env = FakeEnv(
            CPPPATH=FIRST_PARTY + [FRAMEWORK],
            PROJECT_PACKAGES_DIR=_abs("home", "runner", ".platformio", "packages"),
            PROJECT_LIBDEPS_DIR=_abs("work", "Holosphere", ".pio", "libdeps"))
        load_hook("teensy_isystem.py", env=env, projenv=projenv)
        self.assertEqual(projenv["CPPPATH"], FIRST_PARTY)
        self.assertEqual(env["CCFLAGS"], ["-isystem", FRAMEWORK])


class IncludeClassificationTests(unittest.TestCase):
    def test_checkout_packages_component_stays_first_party(self):
        teensy = "/pio/packages/teensy/core"
        fastled = "/work/packages/repo/.pio/libdeps/target/FastLED"
        env = FakeEnv(PROJECT_PACKAGES_DIR="/pio/packages",
                      PROJECT_LIBDEPS_DIR="/work/packages/repo/.pio/libdeps",
                      CPPPATH=[teensy, "/work/packages/repo/core",
                               "/pio/packages-other/core", fastled])
        load_hook("teensy_isystem.py", env=env, projenv=env)
        self.assertEqual(env["CPPPATH"],
                         ["/work/packages/repo/core", "/pio/packages-other/core"])
        self.assertEqual(env["CCFLAGS"], ["-isystem", os.path.normpath(teensy),
                                          "-isystem", os.path.normpath(fastled)])


if __name__ == "__main__":
    unittest.main()
