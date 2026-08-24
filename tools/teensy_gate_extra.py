"""PlatformIO post-build hook: run the Teensy 4 size/layout gate after link.

Wired in via `extra_scripts = post:tools/teensy_gate_extra.py` (platformio.ini).
All decision logic lives in the toolchain-free, unit-tested tools/teensy_gate.py;
this file is glue: it locates the built ELF and the ARM tools, captures their
output, and fails `pio run` on any violation.

Why a post-ACTION that exits non-zero: a post-action that merely prints
does NOT fail `pio run`; only a non-zero exit / raised exception propagates. So
the gate raises on violation. Violations are emitted as GitHub `::error::`
annotations first so they render inline on the PR (ci.yml convention).
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
    r = subprocess.run(args, capture_output=True, text=True, check=check)
    return r.stdout + r.stderr


_TEENSY_SIZE_NAMES = ("teensy_size", "teensy_size.exe")


def _teensy_size_candidates(env):
    """Paths to probe for teensy_size: installed tool package first, PATH last.

    teensy_size ships inside the PlatformIO `tool-teensy` package and is not on
    PATH, so a bare-name probe alone never finds it and the gate silently drops
    to the uncalibrated `size -A` fallback.
    """
    roots = []
    pio_platform = getattr(env, "PioPlatform", None)
    if pio_platform is not None:
        try:
            root = pio_platform().get_package_dir("tool-teensy")
        except Exception:
            root = None
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

    Validates the probe output identifies itself as teensy_size, so an unrelated
    same-named binary on PATH (which would merely launch) is not accepted.
    """
    for cand in _teensy_size_candidates(env):
        try:
            r = subprocess.run([cand, "--help"], capture_output=True, text=True,
                               check=False)
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

    # Capturing the ARM-tool output is the one part of the gate that can fail for
    # reasons unrelated to firmware size: a missing/renamed tool (OSError), a
    # non-zero tool exit (CalledProcessError), or output the parser no longer
    # recognizes (empty regions). Attribute those to a DISTINCT ::error:: so a
    # toolchain/parser break never masquerades as a size-budget "region-missing"
    # violation or an opaque SCons traceback.
    try:
        # Region totals: prefer teensy_size (correct flash-LMA accounting),
        # fall back to `size -A` VMA bucketing (undercounts flash — see teensy_gate).
        teensy_size = _find_teensy_size(env)
        used_size_a_fallback = teensy_size is None
        if teensy_size:
            sizes = teensy_gate.parse_teensy_size(_run([teensy_size, elf], check=False))
        else:
            sizes = teensy_gate.fallback_sizes_from_size_a(
                _run([size_tool, "-A", "-x", elf]))
            print("::warning::teensy_size not found; using `size -A` fallback "
                  "(flash total undercounts; calibrate the flash ceiling against "
                  "teensy_size, not this).")

        symbols = teensy_gate.parse_readelf_symbols(_run([readelf, "-sW", elf]))
        sections = teensy_gate.parse_readelf_sections(_run([readelf, "-SW", elf]))
    except teensy_gate.TeensySizeFormatError as exc:
        print(f"::error::teensy-gate: invalid teensy_size output ({exc}). This "
              f"is a tooling/format error, not a size-budget violation.")
        sys.exit(2)
    except teensy_gate.SizeAFormatError as exc:
        print(f"::error::teensy-gate: invalid `size -A` output ({exc}). This is "
              f"a tooling/format error, not a size-budget violation.")
        sys.exit(2)
    except (OSError, subprocess.SubprocessError) as exc:
        print(f"::error::teensy-gate: a toolchain step failed before evaluation "
              f"({type(exc).__name__}: {exc}). This is a build/tooling break "
              f"(missing/renamed ARM tool or a non-zero tool exit), NOT a "
              f"size-budget violation — fix the toolchain, do not adjust budgets.")
        sys.exit(2)

    if not any(r in sizes for r in ("flash", "ram1", "ram2")):
        print("::error::teensy-gate: parsed no FLASH/RAM1/RAM2 regions from the "
              "size output. This is a toolchain/format break (the size tool's "
              "output shape changed), NOT a size-budget violation.")
        sys.exit(2)

    try:
        budgets = teensy_gate.load_budgets(BUDGETS)
    except teensy_gate.BudgetSchemaError as exc:
        print(f"::error::teensy-gate: invalid budgets schema in {BUDGETS} "
              f"({exc}). This is a budgets-file error (a key the gate never "
              f"reads disables its ceiling), NOT a size-budget violation.")
        sys.exit(2)
    if pioenv not in budgets:
        print(f"::error::teensy-gate: no budget for env '{pioenv}' in {BUDGETS}. "
              f"This is a budgets-file error (the env is unlisted or renamed), "
              f"NOT a size-budget violation.")
        sys.exit(2)

    # The fallback synthesizes region totals with no component breakdown, so a
    # per-component ceiling would evaluate to `component-missing` — a message
    # about a renamed field or a code-size regression, for a missing tool.
    if used_size_a_fallback and teensy_gate.declares_components(budgets[pioenv]):
        print(f"::error::teensy-gate: env '{pioenv}' declares per-component "
              f"ceilings, which the `size -A` fallback cannot measure. teensy_size "
              f"is unavailable, so no component figure was read: this is a "
              f"build/tooling break (install the Teensy platform tools), NOT a "
              f"size-budget violation — do not adjust budgets.")
        sys.exit(2)

    result = teensy_gate.evaluate(pioenv, budgets[pioenv], sizes, symbols, sections)
    if used_size_a_fallback:
        result.notes.insert(0, teensy_gate.UNCALIBRATED_NOTE)
    print(teensy_gate.render_report(result, github=True))
    if not result.passed:
        sys.exit(1)
    if used_size_a_fallback:
        # A bucketed guess must never read as a shipped verdict: CI accepts this
        # build on the gate's exit status, so an uncalibrated PASS exits non-zero.
        print("::error::teensy-gate: PASS is UNCALIBRATED - teensy_size was not "
              "found, so region totals come from `size -A` VMA bucketing. Install "
              "the Teensy platform tools (tool-teensy package) and re-run; do not "
              "record this as a gate PASS.")
        sys.exit(teensy_gate.EXIT_UNCALIBRATED_PASS)


ELF = "$BUILD_DIR/${PROGNAME}.elf"

# A post-action runs only when its target relinks.
env.Depends(ELF, GATE_SOURCES)
env.AddPostAction(ELF, run_gate)
