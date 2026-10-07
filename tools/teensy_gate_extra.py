"""PlatformIO post-build hook: run the Teensy 4 size/layout gate after link.

Glue for tools/teensy_gate.py: locates the built ELF and the ARM tools, feeds
their output to the gate, and exits non-zero on any violation (a post-action
fails `pio run` only by exiting non-zero or raising).
"""

import os
import subprocess
import sys

Import("env")  # noqa: F821  (SCons global injected by PlatformIO)

# SCons exec's this script without a __file__, so derive the tools dir from the
# project dir to import the toolchain-free gate logic.
_TOOLS_DIR = os.path.join(env.subst("$PROJECT_DIR"), "tools")
sys.path.insert(0, _TOOLS_DIR)
import teensy_gate  # noqa: E402

BUDGETS = os.path.join(_TOOLS_DIR, "teensy_budgets.json")
GATE_SOURCES = [BUDGETS,
                os.path.join(_TOOLS_DIR, "teensy_gate.py"),
                os.path.join(_TOOLS_DIR, "teensy_gate_extra.py"),
                os.path.join(_TOOLS_DIR, "phantasm.ld")]


def _tool(cc_path, suffix):
    """Derive a sibling ARM tool path from the C compiler (…-gcc -> …-<suffix>)."""
    base = cc_path
    for end in ("-gcc", "-gcc.exe"):
        if base.endswith(end):
            return base[: -len(end)] + "-" + suffix + (".exe" if end.endswith(".exe") else "")
    # Fall back to PATH lookup.
    return "arm-none-eabi-" + suffix


def _run(args, check=True):
    # teensy_size prints to STDERR (and exits non-zero on overflow, which we still
    # want to parse), so combine both streams and let the caller relax `check`.
    r = subprocess.run(args, capture_output=True, text=True, encoding="utf-8",
                       errors="replace", check=check)
    return r.stdout + r.stderr


_TEENSY_SIZE_NAMES = ("teensy_size", "teensy_size.exe")


def _teensy_size_candidates(env):
    """Paths to probe for teensy_size: installed tool package first, PATH last.

    teensy_size ships inside the PlatformIO `tool-teensy` package and is not on
    PATH.
    """
    roots = []
    pio_platform = getattr(env, "PioPlatform", None)
    if pio_platform is not None:
        root = pio_platform().get_package_dir("tool-teensy")
        if root:
            roots.append(root)
    packages = env.subst("$PROJECT_PACKAGES_DIR")
    if packages and not packages.startswith("$"):
        roots.append(os.path.join(packages, "tool-teensy"))
    core = os.environ.get("PLATFORMIO_CORE_DIR") or os.path.join(
        os.path.expanduser("~"), ".platformio")
    roots.append(os.path.join(core, "packages", "tool-teensy"))

    cands = [os.path.join(root, name) for root in roots
             for name in _TEENSY_SIZE_NAMES]
    cands.extend(_TEENSY_SIZE_NAMES)
    return list(dict.fromkeys(cands))


def _find_teensy_size(env):
    """Best-effort locate of teensy_size (ships with the Teensy platform tools).

    A candidate counts only if its --help output names teensy_size.
    """
    for cand in _teensy_size_candidates(env):
        try:
            r = subprocess.run([cand, "--help"], capture_output=True, text=True,
                               encoding="utf-8", errors="replace", check=False)
        except OSError:
            continue
        if "teensy_size" in (r.stdout + r.stderr).lower():
            return cand
    return None


def run_gate(source, target, env):
    elf = str(target[0])
    pioenv = env["PIOENV"]
    cc = env.subst("$CC")
    size_tool = _tool(cc, "size")
    readelf = _tool(cc, "readelf")

    try:
        teensy_size = _find_teensy_size(env)
        size_text = (_run([teensy_size, elf], check=False) if teensy_size
                     else _run([size_tool, "-A", "-x", elf]))
        syms_text = _run([readelf, "-sW", elf])
        secs_text = _run([readelf, "-SW", elf])
    except (OSError, subprocess.SubprocessError) as exc:
        print(f"::error::teensy-gate: a toolchain step failed before evaluation "
              f"({type(exc).__name__}: {exc}). This is a build/tooling break, "
              "not a size-budget violation.")
        sys.exit(2)

    lines, code = teensy_gate.run(
        pioenv, BUDGETS, size_text, from_teensy_size=teensy_size is not None,
        syms_text=syms_text, secs_text=secs_text, github=True)
    print("\n".join(lines))
    if code:
        sys.exit(code)


ELF = "$BUILD_DIR/${PROGNAME}.elf"

# A post-action runs only when its target relinks.
env.Depends(ELF, GATE_SOURCES)
env.AddPostAction(ELF, run_gate)
