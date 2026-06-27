"""PlatformIO post-build hook: run the Teensy 4 size/layout gate after link.

Wired in via `extra_scripts = post:tools/teensy_gate_extra.py` (platformio.ini).
All decision logic lives in the toolchain-free, unit-tested tools/teensy_gate.py;
this file is glue: it locates the built ELF and the ARM tools, captures their
output, and fails `pio run` on any violation.

Why a post-ACTION that exits non-zero (spec §9): a post-action that merely prints
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


def _find_teensy_size():
    """Best-effort locate of teensy_size (ships with the Teensy platform tools).

    Validates the probe output identifies itself as teensy_size, so an unrelated
    same-named binary on PATH (which would merely launch) is not accepted.
    """
    for cand in ("teensy_size", "teensy_size.exe"):
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
        # Region totals: prefer teensy_size (correct flash-LMA accounting, §7.3),
        # fall back to `size -A` VMA bucketing (undercounts flash — see teensy_gate).
        teensy_size = _find_teensy_size()
        if teensy_size:
            sizes = teensy_gate.parse_teensy_size(_run([teensy_size, elf], check=False))
        else:
            totals = teensy_gate.region_totals_from_size_a(_run([size_tool, "-A", "-x", elf]))
            ram1 = totals.get("ITCM", 0) + totals.get("DTCM", 0)
            sizes = {
                "flash": {"used": totals.get("FLASH", 0), "free": 0},
                "ram1": {"used": ram1, "free": 0x80000 - ram1},
                "ram2": {"used": totals.get("OCRAM", 0),
                         "free": 0x80000 - totals.get("OCRAM", 0)},
            }
            print("::warning::teensy_size not found; using `size -A` fallback "
                  "(flash total undercounts; calibrate the flash ceiling against "
                  "teensy_size, not this).")

        symbols = teensy_gate.parse_readelf_symbols(_run([readelf, "-sW", elf]))
        sections = teensy_gate.parse_readelf_sections(_run([readelf, "-SW", elf]))
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

    budgets = teensy_gate.load_budgets(BUDGETS)
    if pioenv not in budgets:
        print(f"::error::no budget for env '{pioenv}' in {BUDGETS}")
        sys.exit(1)

    result = teensy_gate.evaluate(pioenv, budgets[pioenv], sizes, symbols, sections)
    print(teensy_gate.render_report(result, github=True))
    if not result.passed:
        sys.exit(1)


# Run after the .elf is linked.
env.AddPostAction("$BUILD_DIR/${PROGNAME}.elf", run_gate)
