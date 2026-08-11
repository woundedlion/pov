#!/usr/bin/env python3
"""Host tests for the PlatformIO build hooks (platformio.ini extra_scripts).

The hooks only ever run inside a `pio run`, and most of them fail SILENTLY when
they stop working: teensy_isystem simply demotes nothing (and the warning
ratchet's first-party filter hides the resulting vendored-warning flood),
teensy_map stops emitting the map the size investigations read, teensy_nano
stops selecting the reduced libc without changing a single diagnostic. So each
hook is exec'd here against a fake SCons construction env and its effect on the
build asserted directly -- no ARM toolchain, no PlatformIO.

tools/teensy_gate_extra.py is covered by tools/teensy_gate_tests.

Run:  python -m unittest discover -s tools/teensy_hook_tests
"""

import configparser
import os
import re
import sys
import types
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TOOLS = REPO / "tools"
PIO_INI = REPO / "platformio.ini"

sys.path.insert(0, str(TOOLS))

import teensy_gate as tg   # noqa: E402

GATE_SCRIPT = "post:tools/teensy_gate_extra.py"

# Hooks every environment needs: without teensy_pre the sketch is never placed,
# without teensy_map the size investigations lose their per-symbol source, and
# without teensy_isystem the vendored warning flood re-enters the ratchet.
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


class FakeEnv(dict):
    """The subset of the SCons construction env the hooks touch."""

    def __init__(self, **kw):
        super().__init__(kw)
        self.methods = {}
        self.post_actions = []
        self.lib_builders = []

    def subst(self, s):
        return re.sub(r"\$(\w+)",
                      lambda m: str(self.get(m.group(1), m.group(0))), str(s))

    def Append(self, **kw):
        for key, val in kw.items():
            self.setdefault(key, []).extend(val)

    def Replace(self, **kw):
        self.update(kw)

    def AddMethod(self, fn, name):
        self.methods[name] = fn

    def AddPostAction(self, target, action):
        self.post_actions.append((target, action))

    def File(self, path):
        return path

    def GetLibBuilders(self):
        return self.lib_builders


def load_hook(name, **scons_globals):
    """Exec a hook the way PlatformIO does: SCons globals injected, not imported."""
    path = TOOLS / name
    mod = types.ModuleType(path.stem)
    mod.__file__ = str(path)
    mod.__dict__["Import"] = lambda *a, **k: None
    mod.__dict__.update(scons_globals)
    exec(compile(path.read_text(encoding="utf-8"), mod.__file__, "exec"), mod.__dict__)
    return mod


def _pio_config():
    cfg = configparser.RawConfigParser()
    cfg.read(PIO_INI, encoding="utf-8")
    return cfg


def _pio_envs():
    return [s.split(":", 1)[1] for s in _pio_config().sections()
            if s.startswith("env:")]


_INTERP = re.compile(r"\$\{([^.}]+)\.([^}]+)\}")


def _option_lines(cfg, section, option, depth=0):
    """One option's values as PlatformIO resolves them: [env] supplies the
    default for an [env:X] that omits the option (list options do NOT merge),
    and a ${section.option} line expands to that option's values."""
    if depth > 8:
        raise RecursionError(f"[{section}] {option} interpolates cyclically")
    if not cfg.has_option(section, option) and section.startswith("env:"):
        section = "env"
    out = []
    for line in cfg.get(section, option, fallback="").splitlines():
        line = line.strip()
        if not line:
            continue
        ref = _INTERP.fullmatch(line)
        if ref:
            out.extend(_option_lines(cfg, ref.group(1), ref.group(2), depth + 1))
        else:
            out.append(line)
    return out


class TestExtraScriptsExist(unittest.TestCase):
    """Every hook platformio.ini names must be a real file (a typo'd pre:/post:
    path is accepted by PlatformIO's parser and only surfaces at build time)."""

    def test_every_extra_script_resolves(self):
        cfg = _pio_config()
        seen = 0
        for section in cfg.sections():
            for line in cfg.get(section, "extra_scripts", fallback="").splitlines():
                line = line.strip()
                if not line or line.startswith("${"):
                    continue
                rel = line.split(":", 1)[1]
                self.assertTrue((REPO / rel).is_file(),
                                f"[{section}] extra_scripts names missing {rel}")
                seen += 1
        self.assertGreater(seen, 0, "parsed no extra_scripts from platformio.ini")


class TestBudgetedEnvsWireTheGate(unittest.TestCase):
    """The size/layout gate is a post hook, not a build step, and the envs that
    skip it re-type their extra_scripts block by hand (list options do not
    merge). A budgeted env that loses the line builds green with no ceiling
    enforced -- `pio run` and the teensy-size job both stay silent."""

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
    """PlatformIO cannot subtract from an inherited list, so the four envs that
    skip one hook re-type the whole block by hand. Adding a required hook to
    [env] therefore reaches only the envs that inherit it, and the four that do
    not build green without it."""

    def test_every_env_resolves_the_required_hooks(self):
        cfg = _pio_config()
        envs = _pio_envs()
        self.assertTrue(envs, "platformio.ini declares no environments")
        for name in envs:
            resolved = _option_lines(cfg, f"env:{name}", "extra_scripts")
            for script in REQUIRED_SCRIPTS:
                with self.subTest(env=name, script=script):
                    self.assertIn(script, resolved,
                                  f"env '{name}' does not wire {script}")

    def test_every_base_hook_but_the_gate_reaches_every_env(self):
        # Derived from [env], not from REQUIRED_SCRIPTS: a hook added there that
        # the re-typing envs never pick up is otherwise invisible, because the
        # hand-maintained tuple only has to be a subset of what [env] declares.
        # The gate hook is the one line those envs drop on purpose.
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

    def test_base_env_declares_every_required_hook(self):
        # The [env] block is where a new required hook is added; an entry missing
        # here would make the assertion above pass vacuously for the inheritors.
        declared = _option_lines(_pio_config(), "env", "extra_scripts")
        for script in REQUIRED_SCRIPTS:
            self.assertIn(script, declared, f"[env] does not declare {script}")

    def test_no_env_wires_a_hook_twice(self):
        # ${env.extra_scripts} plus a hand-listed copy runs the hook twice, which
        # double-appends its flags.
        cfg = _pio_config()
        for name in _pio_envs():
            resolved = _option_lines(cfg, f"env:{name}", "extra_scripts")
            self.assertEqual(len(resolved), len(set(resolved)),
                             f"env '{name}' lists a hook more than once: {resolved}")


class TestSketchSelection(unittest.TestCase):
    """teensy_pre.py: PlatformIO globs $PROJECT_SRC_DIR/*.ino, so the sketch is
    chosen here or setup()/loop() never link."""

    def _run(self, pioenv):
        env = FakeEnv(PIOENV=pioenv, PROJECT_DIR=str(REPO))
        return load_hook("teensy_pre.py", env=env), env

    def test_every_platformio_env_is_mapped(self):
        mod, _ = self._run("phantasm")
        for name in _pio_envs():
            self.assertIn(name, mod.SKETCH,
                          f"env '{name}' has no sketch mapping in teensy_pre.py")

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

    def test_unknown_env_fails_loud(self):
        with self.assertRaises(SystemExit):
            self._run("not_an_env")


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
        # A C++ command is `$CXXFLAGS $CCFLAGS`; adding it to both duplicates the
        # spec on one command line and gcc fatals ("spec 'link' already defined").
        env = FakeEnv(CCFLAGS=[], CXXFLAGS=[], LINKFLAGS=[])
        load_hook("teensy_nano.py", env=env)
        self.assertNotIn("--specs=nano.specs", env["CXXFLAGS"])


class TestIsystemDemotion(unittest.TestCase):
    """teensy_isystem.py: vendored include dirs must move from -I to -isystem."""

    def _run(self, cpppath, lib_builders=0):
        projenv = FakeEnv(CPPPATH=list(cpppath))
        env = FakeEnv(CPPPATH=list(cpppath))
        env.lib_builders = [types.SimpleNamespace(env=FakeEnv())
                            for _ in range(lib_builders)]
        load_hook("teensy_isystem.py", env=env, projenv=projenv)
        return projenv, env

    def test_vendored_dirs_move_out_of_cpppath(self):
        projenv, env = self._run(FIRST_PARTY + [LIBDEPS, FRAMEWORK, RELOCATED])
        for build_env in (projenv, env):
            self.assertEqual(build_env["CPPPATH"], FIRST_PARTY)
            self.assertEqual(build_env["CCFLAGS"],
                             ["-isystem", LIBDEPS,
                              "-isystem", FRAMEWORK,
                              "-isystem", RELOCATED])

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
        # The silent no-op this guard exists for: the markers stop matching and
        # the build quietly reverts to a vendored-warning flood.
        with self.assertRaises(SystemExit) as cm:
            self._run(FIRST_PARTY)
        self.assertIn("demoted 0 third-party include dirs", str(cm.exception))

    def test_one_env_carrying_the_vendored_dirs_is_enough(self):
        projenv = FakeEnv(CPPPATH=list(FIRST_PARTY))
        env = FakeEnv(CPPPATH=FIRST_PARTY + [FRAMEWORK])
        load_hook("teensy_isystem.py", env=env, projenv=projenv)
        self.assertEqual(projenv["CPPPATH"], FIRST_PARTY)
        self.assertEqual(env["CCFLAGS"], ["-isystem", FRAMEWORK])


if __name__ == "__main__":
    unittest.main()
