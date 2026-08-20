#!/usr/bin/env python3
"""Single source for build values shared by CI, the pre-commit hook and just.

PINS are injected (--github-output / a named lookup) and must appear in no
build file literally. INLINE_PINS cannot be injected -- they name an action
input, an apt package, a pip requirement resolved before this script could run,
or an application installed by hand -- so --check pins them by consistency
instead: every spelling of them in INLINE_SCAN must derive from the value here.

SHARED_LITERALS are pinned the same way: strings two build files must spell
identically because neither can read the other.

ENGINE_RANGES name a manifest whose declared floor the pin must satisfy: the
manifest is read by a package manager, not by this script, so it cannot take
the pin -- --check asserts the two agree instead.

TEST_FILE_PIN reads the check_test_files.sh counts back out of the files that
spell them, so the copies cannot drift from each other or from the suite.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path


PINS = {
    "doxygen-awesome": "568f56cde6ac78b6dfcc14acd380b2e745c301ea",
    "emsdk": "5.0.0",
    # PyPI's rust-just, so the recipe runner is installed and held like ruff.
    "just": "1.52.0",
    "node": "24.13.0",
    "platformio": "6.1.19",
    "ruff": "0.14.4",
    # PyPI's shellcheck-py, whose version is the shellcheck release plus a
    # packaging suffix.
    "shellcheck": "0.11.0.1",
}

# Versions the build files must spell out literally: a setup-action input, an
# apt package name and a pip requirement are all resolved before this script
# could inject anything, and the interpreter pin bootstraps the script itself.
# They are single-sourced by CONSISTENCY instead of substitution -- --check
# asserts every occurrence equals the value here, so a partial bump fails.
INLINE_PINS = {
    "python": "3.12",
    "numpy": "2.4.3",
    "clang": "22",
    "clang-format": "22.1.8",
    "doxygen": "1.17.0",
    "doxygen-sha256":
        "75419ef4f446fc1c24ef12514b574e66e898ee6f527c6ae2ad84f91a905823c2",
    # apt.llvm.org's signing key, which authenticates every package that host
    # serves -- including the host's own copy of the key.
    "llvm-key-sha256":
        "8b2a587ffd672c4687e7581dad4b2f6c1bb2ad6b480cd9771ba2ff48e0b8c75d",
    # The major of the desktop KiCad install the fab gates shell out to. It is
    # installed from KiCad's own installer, so no build file can inject it.
    "kicad": "10",
}

ROOT = Path(__file__).resolve().parents[1]
CONSUMERS = {
    ROOT / ".github/workflows/ci.yml": (
        "build_pins.py --github-output",
        "build_pins.py platformio",
    ),
    ROOT / ".github/workflows/docs.yml": ("build_pins.py doxygen-awesome",),
    ROOT / "justfile": (
        "build_pins.py doxygen-awesome",
        "build_pins.py --check-tool doxygen",
        "build_pins.py --check-tool node",
        "build_pins.py --check-tool platformio",
        "build_pins.py --check-tool ruff",
        "build_pins.py --check-tool shellcheck",
    ),
}

# Files scanned for INLINE_PINS occurrences.
INLINE_SCAN = (
    ROOT / ".github/workflows/ci.yml",
    ROOT / ".github/workflows/docs.yml",
    ROOT / "justfile",
    ROOT / "scripts/generate_luts.py",
    ROOT / ".githooks/pre-commit",
    ROOT / "hardware/phantasm/gen/kicad_common.py",
    # Prose, but a contributor installs what it tells them to: a stale spelling
    # here sends them to a version the format gate then rejects.
    ROOT / "README.md",
    ROOT / "CONTRIBUTING.md",
    ROOT / "hardware/phantasm/README.md",
)

# (pattern, pin name, expected form of the pin value). The pattern's single
# capture group is the version as that spelling writes it; `form` maps the pin
# value to that spelling, so every occurrence is derived from one string.
INLINE_USES = (
    (r"python-version:\s*'([^']*)'", "python", lambda v: v),
    (r"\bnumpy==([\w.]+)", "numpy", lambda v: v),
    (r"\b(?:clang\+\+|clang|llvm)-(\d+)\b", "clang", lambda v: v),
    (r"\bllvm-\w+-(\d+)\b", "clang", lambda v: v),
    (r"\bclang-format==([\w.]+)", "clang-format", lambda v: v),
    (r"\bclang-format-(\d+)\b", "clang-format", lambda v: v.split(".")[0]),
    (r"\bclang-format (\d+)\b", "clang-format", lambda v: v.split(".")[0]),
    (r"EXPECTED_CLANG_FORMAT_MAJOR = (\d+)", "clang-format",
     lambda v: v.split(".")[0]),
    (r"HS_CLANG_FORMAT_MAJOR=(\d+)", "clang-format", lambda v: v.split(".")[0]),
    (r"Install Doxygen ([\w.]+) ", "doxygen", lambda v: v),
    (r"/Release_(\w+)/", "doxygen", lambda v: v.replace(".", "_")),
    (r"\bdoxygen-([\w.]+)\.linux", "doxygen", lambda v: v),
    (r"\bdoxygen-([\w.]+)/bin", "doxygen", lambda v: v),
    (r"([0-9a-f]{64})  doxygen\.tar\.gz", "doxygen-sha256", lambda v: v),
    (r"([0-9a-f]{64})  /tmp/llvm-snapshot\.gpg\.key", "llvm-key-sha256",
     lambda v: v),
    (r"KICAD_MAJOR = (\d+)", "kicad", lambda v: v),
    (r"\bKiCad (\d+)\b", "kicad", lambda v: v),
)

# The extensions the clang-format gate covers. ci.yml and the justfile select
# them as a git ls-files pathspec, the pre-commit hook as a grep alternation
# over staged paths; both spellings are derived from this one tuple, so an
# extension added here reaches all three copies or none.
FORMAT_EXTENSIONS = ("h", "hpp", "cpp", "cc", "inl")

# Strings that must read identically in several build files. The pre-commit hook
# is POSIX shell run per commit, ci.yml is workflow YAML and the justfile is a
# recipe list, so none can source a value from another; --check asserts every
# occurrence matches this one.
SHARED_LITERALS = {
    # Paths the clang-format gate skips: vendored sources and generated ones.
    "format-exclude": (
        r"(^|/)core/vendor/"
        r"|(^|/)core/color/color_luts\.h$"
        r"|(^|/)core/color/gamut_lut\.h$"
        r"|(^|/)core/spatial/reaction_graph\.cpp$"
        r"|(^|/)tests/mindsplatter_replay_corpus\.h$"
    ),
    "format-globs": " ".join(f"'*.{ext}'" for ext in FORMAT_EXTENSIONS),
    "format-extensions": r"\.(" + "|".join(FORMAT_EXTENSIONS) + r")$",
}

# (pattern, literal name, occurrences required across INLINE_SCAN). The
# pattern's single capture group is the literal as that file spells it; the
# count fails a copy that was dropped or another one added unnoticed. The
# pathspec pattern is anchored on the first glob so the shellcheck gate's own
# git ls-files pathspec is not swept in.
SHARED_LITERAL_USES = (
    (r"grep -vE '([^']*)'", "format-exclude", 3),
    (r"git ls-files -- ('\*\.h'(?: '\*\.\w+')*)", "format-globs", 2),
    (r"grep -E '([^']*)'", "format-extensions", 1),
)

# check_test_files.sh calls, whose count argument is exact. Every such glob is
# spelled once per file that runs the suite -- a CI step and usually a justfile
# recipe mirroring it -- and neither copy can read the other, so a suite that
# gains a file must move every pin at once.
TEST_FILE_PIN = re.compile(r"""check_test_files\.sh\s+(\d+)\s+['"]([^'"]+)['"]""")
TEST_FILE_PIN_SCAN = (
    ROOT / "justfile",
    ROOT / ".github/workflows/ci.yml",
)

# --check-tool targets: pin name -> (version command, how to install the pin,
# the form of the pin that command reports). `{pin}` in either is filled with
# the pin value. A pin absent here has nothing to interrogate: a git ref and
# two file digests name no program, and the emsdk and KiCad pins are checked
# where they are used (the WASM toolchain marker written beside the build,
# kicad_common.find_kicad_cli). The targets a recipe runs are invoked from
# it (CONSUMERS pins those call sites); clang, just, numpy and python are
# manual probes, no recipe spelling them -- the native build reaches clang
# through the CMake toolchain, which may be emsdk's.
CHECK_TOOLS = {
    "clang": (["clang-{pin}", "--version"], "apt install clang-{pin}",
              lambda v: v),
    "clang-format": (["clang-format", "--version"],
                     "pip install clang-format=={pin}", lambda v: v),
    "doxygen": (["doxygen", "--version"],
                "install Doxygen {pin} from doxygen.nl", lambda v: v),
    "just": (["just", "--version"], "pip install rust-just=={pin}",
             lambda v: v),
    "node": (["node", "--version"], "install Node {pin}", lambda v: v),
    # No console script; the pin is met by the module the interpreter imports.
    "numpy": ([sys.executable, "-c", "import numpy; print(numpy.__version__)"],
              "pip install numpy=={pin}", lambda v: v),
    "platformio": (["platformio", "--version"],
                   "pip install platformio=={pin}", lambda v: v),
    # The interpreter that matters is the one running this script, not
    # whatever `python` resolves to on PATH.
    "python": ([sys.executable, "--version"], "install Python {pin}",
               lambda v: v),
    "ruff": (["ruff", "--version"], "pip install ruff=={pin}", lambda v: v),
    # shellcheck-py's version is the shellcheck release plus a packaging
    # suffix, which the binary itself never reports.
    "shellcheck": (["shellcheck", "--version"],
                   "pip install shellcheck-py=={pin}",
                   lambda v: v.rsplit(".", 1)[0]),
}

# (manifest, JSON path to a `>=X` range, pin name). The range is the floor the
# tree's own tooling needs -- `node --test`'s glob expansion needs >= 22 -- and
# the pin is what every CI job installs, so a pin below the floor would run the
# suite on an interpreter the manifest already rejects.
ENGINE_RANGES = (
    (ROOT / "package.json", ("engines", "node"), "node"),
)


def duplicates_pin(text: str, name: str, value: str) -> bool:
    """Return whether a dependency context contains its literal pinned value."""
    aliases = (name.lower(), name.lower().replace("-", "_"))
    value_pattern = re.compile(
        rf"(?<![0-9A-Za-z.]){re.escape(value)}(?![0-9A-Za-z.])"
    )
    lines = text.splitlines()
    for index, line in enumerate(lines):
        code = line.split("#", 1)[0]
        if not value_pattern.search(code):
            continue
        context = "\n".join(lines[max(0, index - 2) : index + 1]).lower()
        if any(alias in context for alias in aliases):
            return True
    return False


def check_inline_pins() -> list[str]:
    """Return one error per occurrence of an INLINE_PINS value that disagrees
    with the pin, plus one per pin no scanned file spells at all."""
    errors: list[str] = []
    seen: dict[str, int] = {name: 0 for name in INLINE_PINS}
    for path in INLINE_SCAN:
        text = path.read_text(encoding="utf-8")
        for index, line in enumerate(text.splitlines(), 1):
            for pattern, name, form in INLINE_USES:
                want = form(INLINE_PINS[name])
                for found in re.findall(pattern, line):
                    seen[name] += 1
                    if found != want:
                        errors.append(
                            f"{path.relative_to(ROOT)}:{index}: {name} pinned to "
                            f"{want!r} in build_pins.py but written {found!r}"
                        )
    for name, count in seen.items():
        if count == 0:
            errors.append(f"{name} pin is spelled in no scanned file")
    return errors


def check_shared_literals() -> list[str]:
    """Return one error per occurrence of a SHARED_LITERALS value that disagrees
    with the literal, plus one per literal found the wrong number of times."""
    errors: list[str] = []
    for pattern, name, expected in SHARED_LITERAL_USES:
        want = SHARED_LITERALS[name]
        occurrences = 0
        for path in INLINE_SCAN:
            text = path.read_text(encoding="utf-8")
            for index, line in enumerate(text.splitlines(), 1):
                for found in re.findall(pattern, line):
                    occurrences += 1
                    if found != want:
                        errors.append(
                            f"{path.relative_to(ROOT)}:{index}: {name} differs "
                            f"from build_pins.py: {found!r}"
                        )
        if occurrences != expected:
            errors.append(
                f"{name} occurs {occurrences} time(s) in the scanned files, "
                f"expected {expected}"
            )
    return errors


def collect_test_file_pins() -> dict[str, list[tuple[Path, int, int]]]:
    """Return glob -> [(file, line, pinned count)] over TEST_FILE_PIN_SCAN."""
    pins: dict[str, list[tuple[Path, int, int]]] = {}
    for path in TEST_FILE_PIN_SCAN:
        text = path.read_text(encoding="utf-8")
        for index, line in enumerate(text.splitlines(), 1):
            for count, glob in TEST_FILE_PIN.findall(line):
                pins.setdefault(glob, []).append((path, index, int(count)))
    return pins


def check_test_file_pins() -> list[str]:
    """Return one error per check_test_files.sh glob whose copies disagree, or
    whose pin disagrees with the files the glob matches."""
    errors: list[str] = []
    for glob, uses in sorted(collect_test_file_pins().items()):
        counts = sorted({count for _, _, count in uses})
        found = sum(1 for path in ROOT.glob(glob) if path.is_file())
        for path, index, count in uses:
            where = f"{path.relative_to(ROOT)}:{index}"
            if len(counts) > 1:
                errors.append(
                    f"{where}: pins {count} for {glob!r}, which another copy "
                    f"pins {[c for c in counts if c != count]}")
            elif count != found:
                errors.append(
                    f"{where}: pins {count} for {glob!r}, which matches "
                    f"{found} file(s)")
    return errors


def _version_tuple(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.split("."))


def check_engine_ranges() -> list[str]:
    """Return one error per manifest whose declared floor the pin misses."""
    errors: list[str] = []
    for path, keys, name in ENGINE_RANGES:
        where = path.relative_to(ROOT)
        node = json.loads(path.read_text(encoding="utf-8"))
        for key in keys:
            node = node.get(key) if isinstance(node, dict) else None
        if not isinstance(node, str):
            errors.append(f"{where}: no {'.'.join(keys)} range for the "
                          f"{name} pin")
            continue
        found = re.fullmatch(r">=\s*([0-9]+(?:\.[0-9]+)*)", node.strip())
        if found is None:
            errors.append(
                f"{where}: {'.'.join(keys)} is {node!r}, which build_pins.py "
                f"cannot compare against the {name} pin (expected '>=X')")
            continue
        floor = _version_tuple(found.group(1))
        pinned = _version_tuple(PINS[name])
        width = max(len(floor), len(pinned))
        floor += (0,) * (width - len(floor))
        pinned += (0,) * (width - len(pinned))
        if pinned < floor:
            errors.append(
                f"{where}: {'.'.join(keys)} requires {node!r} but {name} is "
                f"pinned to {PINS[name]} in build_pins.py")
    return errors


def check_consumers() -> int:
    errors: list[str] = (check_inline_pins() + check_shared_literals()
                         + check_engine_ranges() + check_test_file_pins())
    for name in sorted(set(CHECK_TOOLS) - set(PINS | INLINE_PINS)):
        errors.append(f"CHECK_TOOLS names {name}, which is not a pin")
    for path, references in CONSUMERS.items():
        text = path.read_text(encoding="utf-8")
        for reference in references:
            if reference not in text:
                errors.append(f"{path.relative_to(ROOT)}: missing {reference!r}")
    workflow_paths = set((ROOT / ".github/workflows").glob("*.yml"))
    workflow_paths.update((ROOT / ".github/workflows").glob("*.yaml"))
    duplicate_paths = workflow_paths | set(CONSUMERS)
    for path in sorted(duplicate_paths):
        text = path.read_text(encoding="utf-8")
        for name, value in PINS.items():
            if duplicates_pin(text, name, value):
                errors.append(
                    f"{path.relative_to(ROOT)}: duplicates {name} pin {value}"
                )
    if errors:
        for error in errors:
            print(error)
        return 1
    print(f"build pins are single-sourced ({len(PINS)} injected, "
          f"{len(INLINE_PINS)} inline, {len(SHARED_LITERALS)} shared literal, "
          f"{len(ENGINE_RANGES)} engine range, "
          f"{len(collect_test_file_pins())} test-file glob)")
    return 0


def check_tool(name: str) -> int:
    """Fail unless the installed tool reports the version pinned here.

    A recipe that just invokes a linter runs whatever the developer installed,
    which gates the tree on a different rule set than CI's.

    The pin's precision is the comparison's: a major-only pin such as clang's
    is met by any release of that major, which is all the pin claims.
    """
    pin = (PINS | INLINE_PINS)[name]
    command, install, form = CHECK_TOOLS[name]
    command = [part.format(pin=pin) for part in command]
    want = form(pin).split(".")
    try:
        reported = subprocess.run(
            command, capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.SubprocessError):
        reported = ""
    found = re.search(r"\d[\w.]*", reported)
    if found is None or found.group().split(".")[: len(want)] != want:
        got = found.group() if found else "nothing runnable"
        print(f"{name} is pinned to {pin} but PATH has {got} "
              f"({install.format(pin=pin)})")
        return 1
    print(f"{name} {pin} matches the pin")
    return 0


def write_github_output() -> None:
    output = os.environ.get("GITHUB_OUTPUT")
    if not output:
        raise SystemExit("GITHUB_OUTPUT is not set")
    with open(output, "a", encoding="utf-8", newline="\n") as stream:
        for name, value in PINS.items():
            stream.write(f"{name.replace('-', '_')}={value}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", nargs="?", choices=sorted(PINS | INLINE_PINS))
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--check-tool", metavar="NAME",
                        choices=sorted(CHECK_TOOLS))
    parser.add_argument("--github-output", action="store_true")
    args = parser.parse_args()
    if args.check:
        return check_consumers()
    if args.check_tool:
        return check_tool(args.check_tool)
    if args.github_output:
        write_github_output()
        return 0
    if args.name:
        print((PINS | INLINE_PINS)[args.name])
        return 0
    parser.error(
        "specify a pin name, --check, --check-tool NAME, or --github-output")


if __name__ == "__main__":
    raise SystemExit(main())
