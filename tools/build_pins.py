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

The install-set check reads the same CMakeLists.txt rules that mirror files into
the daydream checkout and asserts each carries a line-ending pin, since those
bytes are digested there over LF.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import teensy_gate  # noqa: E402


PINS = {
    # PyPI's actionlint-py, whose version is the actionlint release plus a
    # packaging suffix.
    "actionlint": "1.7.12.24",
    # The daydream commit the tree-fence and cross-repo jobs check out. Those
    # jobs validate an exhaustive fence against that tree, so the checkout
    # tracks a committed pin: daydream's tip moves independently and would both
    # red unrelated engine changes and make a re-run of an unchanged commit flip
    # its verdict. daydream pins this repository the same way, in
    # holosphere_wasm.sha -- so the pair is circular and daydream lands first.
    "daydream": "c7ab8101386882e5849c4b7fcadf0ba9d0c39fad",
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
WORKFLOWS = ROOT / ".github/workflows"
ACTIONS = ROOT / ".github/actions"


def workflow_files() -> tuple[Path, ...]:
    """Every workflow and composite action, so a new one is scanned with no
    second edit here."""
    found = set(WORKFLOWS.glob("*.yml")) | set(WORKFLOWS.glob("*.yaml"))
    found |= set(ACTIONS.glob("*/action.yml"))
    found |= set(ACTIONS.glob("*/action.yaml"))
    return tuple(sorted(found))


CONSUMERS = {
    ROOT / ".github/workflows/ci.yml": (
        "build_pins.py --github-output",
        "python tools/build_pins.py --check",
    ),
    ROOT / ".github/actions/pinned-doxygen/action.yml": (
        "build_pins.py doxygen-awesome",
    ),
    ROOT / "justfile": (
        "build_pins.py doxygen-awesome",
        "build_pins.py --check-tool doxygen",
        "build_pins.py --check-tool node",
        "build_pins.py --check-tool platformio",
        "build_pins.py --check-tool ruff",
        "build_pins.py --check-tool shellcheck",
        "{{py}} tools/build_pins.py --check",
    ),
    ROOT / ".githooks/pre-commit": (
        '"$PYTHON_BIN" "$SNAPSHOT/tools/build_pins.py" --check',
    ),
}

# Files scanned for INLINE_PINS occurrences. Every workflow and composite action
# is covered, so a new one cannot spell a pin its own way and still pass --check.
INLINE_SCAN = (
    *workflow_files(),
    ROOT / "justfile",
    # The justfile's clang-format recipe body, where the gate's pathspec and
    # exclusion regex live.
    ROOT / "tools/clang_format_gate.sh",
    ROOT / "platformio.ini",
    ROOT / "CMakeLists.txt",
    ROOT / "scripts/generate_luts.py",
    ROOT / ".githooks/pre-commit",
    ROOT / "hardware/phantasm/gen/sexp.py",
    # Prose, but a contributor installs what it tells them to: a stale spelling
    # here sends them to a version the format gate then rejects.
    ROOT / "README.md",
    ROOT / "CONTRIBUTING.md",
    ROOT / "hardware/phantasm/README.md",
    # Both halves of each requirements pair: the hand-edited .in carries the
    # pin, the pip-compile'd .txt repeats it above the hashes.
    *(ROOT / "requirements" / f"{stem}{suffix}"
      for stem in ("actionlint", "clang-format", "just", "numpy",
                   "platformio", "ruff", "shellcheck")
      for suffix in (".in", ".txt")),
)

# (pattern, pin name, expected form of the pin value, occurrences required
# across INLINE_SCAN). The pattern's single capture group is the version as that
# spelling writes it; `form` maps the pin value to that spelling, so every
# occurrence is derived from one string. The count fails a site dropped or
# another added unnoticed: a bare "matches nothing anywhere" guard passes a
# re-spelling for as long as one sibling site keeps matching.
INLINE_USES = (
    (r"\bactionlint-py==([\w.]+)", "actionlint", lambda v: v, 2),
    # Quote-agnostic: setup-python's own README writes the input with double
    # quotes, which a single-quoted pattern reads as absent.
    (r"""python-version:\s*['"]?([^'"\s]+)['"]?""", "python", lambda v: v, 16),
    (r"\bnumpy==([\w.]+)", "numpy", lambda v: v, 2),
    (r"\b(?:clang\+\+|clang|llvm)-(\d+)\b", "clang", lambda v: v, 28),
    (r"\bllvm-\w+-(\d+)\b", "clang", lambda v: v, 7),
    (r"\bclang-format==([\w.]+)", "clang-format", lambda v: v, 4),
    (r"\bclang-format-(\d+)\b", "clang-format", lambda v: v.split(".")[0], 1),
    (r"\bclang-format (\d+)\b", "clang-format", lambda v: v.split(".")[0], 3),
    (r"\brust-just==([\w.]+)", "just", lambda v: v, 2),
    (r"\bplatformio==([\w.]+)", "platformio", lambda v: v, 2),
    (r"\bruff==([\w.]+)", "ruff", lambda v: v, 2),
    (r"\bshellcheck-py==([\w.]+)", "shellcheck", lambda v: v, 2),
    (r"EXPECTED_CLANG_FORMAT_MAJOR = (\d+)", "clang-format",
     lambda v: v.split(".")[0], 1),
    (r"HS_CLANG_FORMAT_MAJOR=(\d+)", "clang-format",
     lambda v: v.split(".")[0], 1),
    (r"Install Doxygen ([\w.]+) ", "doxygen", lambda v: v, 1),
    (r"/Release_(\w+)/", "doxygen", lambda v: v.replace(".", "_"), 1),
    (r"\bdoxygen-([\w.]+)\.linux", "doxygen", lambda v: v, 1),
    (r"\bdoxygen-([\w.]+)/bin", "doxygen", lambda v: v, 1),
    (r"([0-9a-f]{64})  doxygen\.tar\.gz", "doxygen-sha256", lambda v: v, 1),
    (r"([0-9a-f]{64})  /tmp/llvm-snapshot\.gpg\.key", "llvm-key-sha256",
     lambda v: v, 1),
    (r"KICAD_MAJOR = (\d+)", "kicad", lambda v: v, 1),
    (r"\bKiCad (\d+)\b", "kicad", lambda v: v, 6),
)

# The extensions the clang-format gate covers. ci.yml and the justfile select
# them as a git ls-files pathspec, the pre-commit hook as a grep alternation
# over staged paths; both spellings are derived from this one tuple, so an
# extension added here reaches all three copies or none.
FORMAT_EXTENSIONS = ("h", "hpp", "cpp", "cc", "inl")

# The float flags both shipping targets build with: the device firmware
# (platformio.ini) and the WASM modules (CMakeLists.txt). The CI leg runs the
# suite under the same pair and defines which test contract applies.
# -fno-finite-math-only must follow -ffast-math, which otherwise implies
# -ffinite-math-only and folds every std::isfinite() boundary predicate to
# constant true; a spelling that loses it builds green.
FLOAT_FLAGS = ("-ffast-math", "-fno-finite-math-only")
FAST_MATH_TEST_FLAGS = (*FLOAT_FLAGS, "-DHS_TEST_FAST_MATH=1")

# Strings that must read identically in several build files. The pre-commit hook
# is POSIX shell run per commit, ci.yml is workflow YAML, the justfile is a
# recipe list and the gate scripts are bash, so none can source a value from
# another; --check asserts every occurrence matches this one.
SHARED_LITERALS = {
    # Paths the clang-format gate skips: vendored sources and generated ones.
    "format-exclude": (
        r"(^|/)core/vendor/"
        r"|(^|/)core/color/color_luts\.h$"
        r"|(^|/)core/color/gamut_lut\.h$"
        r"|(^|/)core/color/srgb_decode_lut\.h$"
        r"|(^|/)core/color/triadic_palette_luts\.h$"
        r"|(^|/)core/mesh/relax_bakes_generated\.h$"
        r"|(^|/)core/spatial/reaction_graph\.cpp$"
        r"|(^|/)tests/mindsplatter_replay_corpus\.h$"
    ),
    "format-globs": " ".join(f"'*.{ext}'" for ext in FORMAT_EXTENSIONS),
    "format-extensions": r"\.(" + "|".join(FORMAT_EXTENSIONS) + r")$",
    # A flag list (platformio.ini's build_flags, ci.yml's matrix entry) and the
    # CMake quoted-argument form of the same pair.
    "float-flags": " ".join(FLOAT_FLAGS),
    "float-test-flags": " ".join(FAST_MATH_TEST_FLAGS),
    "float-flags-cmake": " ".join(f'"{flag}"' for flag in FLOAT_FLAGS),
    # The per-effect smoke window every CI leg drives. The justfile and
    # CONTRIBUTING repeat it because neither can read the workflow.
    "smoke-frames": "120",
}

# (pattern, literal name, occurrences required across INLINE_SCAN). The
# pattern's single capture group is the literal as that file spells it; the
# count fails a copy that was dropped or another one added unnoticed. The
# pathspec pattern is anchored on the first glob so the shellcheck gate's own
# git ls-files pathspec is not swept in.
SHARED_LITERAL_USES = (
    (r"grep -vE '([^']*)'", "format-exclude", 3),
    (r"git ls-files -- ('\*\.h'(?: '\*\.\w+')*)", "format-globs", 2),
    # Anchored on the extension alternation's own opening, so an unrelated
    # quoted `grep -E` elsewhere in a scanned file is not counted as a copy.
    (r"grep -E '(\\\.\([^']*)'", "format-extensions", 1),
    # Anchored on the start of a flag line (platformio.ini) and on the matrix key
    # that carries the pair plus its test-contract define (ci.yml), so prose
    # spellings are not swept in. Captures run to end of line: a partial edit
    # reads as a difference.
    (r"^\s+(-ffast-math\b.*)$", "float-flags", 1),
    (r"float_flags:\s+(-ffast-math\b.*)$", "float-test-flags", 1),
    # Every WASM target's compile and link line.
    (r'("-ffast-math" "[^"]+")', "float-flags-cmake", 4),
    # ci.yml declares the window once and aliases it into every other job; the
    # justfile spells it as a recipe parameter and CONTRIBUTING as prose.
    (r'HS_SMOKE_FRAMES(?:: &\w+ |="?)(\d+)', "smoke-frames", 3),
)

# --check-tool targets: pin name -> (version command, how to install the pin,
# the form of the pin that command reports). `{pin}` in either is filled with
# the pin value. A pin absent here has nothing to interrogate: two git refs and
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


# install(FILES|DIRECTORY) rules, matched to the first closing paren -- neither
# form nests one. The cross-repo install set is the subset destined for
# DAYDREAM_DIR; install(CODE) writes only files generated at install time.
INSTALL_RULE = re.compile(r"install\(\s*(FILES|DIRECTORY)\s(.*?)\)", re.DOTALL)
INSTALL_SOURCE = re.compile(r"\$\{CMAKE_CURRENT_SOURCE_DIR\}/([^\"\s]+)")
INSTALL_PATTERN = re.compile(r'PATTERN\s+"([^"]+)"')


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


def read_scanned(path: Path, errors: list[str]) -> str | None:
    """Text of a scanned file, or None with the failure recorded as an error.

    Every scanned path is named in a table here, so a rename or a deletion is a
    finding this gate reports; a traceback out of a hook is not a report.
    """
    try:
        return path.read_text(encoding="utf-8")
    except OSError as error:
        errors.append(f"{path.relative_to(ROOT)}: cannot be read ({error})")
        return None


def check_inline_pins() -> list[str]:
    """Return one error per occurrence of an INLINE_PINS value that disagrees
    with the pin, plus one per INLINE_USES spelling found the wrong number of
    times.

    The count is per pattern, not per pin: most pins have several spellings, so
    a per-pin total would let one spelling stop matching while a sibling absorbs
    its share.
    """
    errors: list[str] = []
    pin_values = {**PINS, **INLINE_PINS}
    seen: dict[str, int] = {pattern: 0 for pattern, _, _, _ in INLINE_USES}
    for path in INLINE_SCAN:
        text = read_scanned(path, errors)
        if text is None:
            continue
        for index, line in enumerate(text.splitlines(), 1):
            for pattern, name, form, _ in INLINE_USES:
                want = form(pin_values[name])
                for found in re.findall(pattern, line):
                    seen[pattern] += 1
                    if found != want:
                        errors.append(
                            f"{path.relative_to(ROOT)}:{index}: {name} pinned to "
                            f"{want!r} in build_pins.py but written {found!r}"
                        )
    for pattern, name, _, expected in INLINE_USES:
        if seen[pattern] != expected:
            errors.append(
                f"{name} spelling {pattern!r} occurs {seen[pattern]} time(s) in "
                f"the scanned files, expected {expected}"
            )
    return errors


def check_shared_literals() -> list[str]:
    """Return one error per occurrence of a SHARED_LITERALS value that disagrees
    with the literal, plus one per literal found the wrong number of times."""
    errors: list[str] = []
    for pattern, name, expected in SHARED_LITERAL_USES:
        want = SHARED_LITERALS[name]
        occurrences = 0
        for path in INLINE_SCAN:
            text = read_scanned(path, errors)
            if text is None:
                continue
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


def _version_tuple(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.split("."))


def check_engine_ranges() -> list[str]:
    """Return one error per manifest whose declared floor the pin misses."""
    errors: list[str] = []
    for path, keys, name in ENGINE_RANGES:
        where = path.relative_to(ROOT)
        text = read_scanned(path, errors)
        if text is None:
            continue
        try:
            node = json.loads(text)
        except json.JSONDecodeError as error:
            errors.append(f"{where}: is not valid JSON ({error})")
            continue
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


def check_flexram_geometry() -> list[str]:
    """Tie the budget and gate constants to the linker script geometry."""
    budgets_path = ROOT / "tools/teensy_budgets.json"
    try:
        budgets = teensy_gate.load_budgets(budgets_path)
    except (OSError, ValueError) as exc:
        # An unreadable file, a BudgetSchemaError, a JSON syntax error and an
        # unterminated block comment all land here; a traceback out of a hook is
        # not a report.
        return [f"tools/teensy_budgets.json: {exc}"]
    derived = budgets["phantasm"]["regions"]["ram1"][
        "components"]["code"]["max_banks_from_stack_floor"]
    bank_bytes = derived["bank_bytes"]
    total_banks = derived["total_banks"]
    errors: list[str] = []

    gate_text = read_scanned(ROOT / "tools/teensy_gate.py", errors)
    if gate_text is not None:
        gate_match = re.search(r"^FLEXRAM_BANK_BYTES = (0x[0-9a-fA-F]+|\d+)$",
                               gate_text, re.MULTILINE)
        if gate_match is None or int(gate_match.group(1), 0) != bank_bytes:
            errors.append("tools/teensy_gate.py: FlexRAM bank size differs "
                          "from teensy_budgets.json")

    shift = bank_bytes.bit_length() - 1
    linker = read_scanned(ROOT / "tools/phantasm.ld", errors)
    if linker is None:
        return errors
    linker_spellings = (
        f"+ 0x{bank_bytes - 1:X}) >> {shift}",
        f"(({total_banks} - _itcm_block_count) << {shift})",
    )
    for spelling in linker_spellings:
        if spelling not in linker:
            errors.append(
                f"tools/phantasm.ld: missing FlexRAM geometry derived from "
                f"teensy_budgets.json: {spelling!r}")
    return errors


def installed_sources() -> list[str]:
    """Return the repository files CMakeLists.txt installs into daydream.

    An install(DIRECTORY) rule contributes whatever its FILES_MATCHING patterns
    select right now, so a document dropped into patterns/ joins the set by
    existing rather than by being listed anywhere.
    """
    # An unreadable CMakeLists.txt yields no install source, which
    # check_install_eol reports rather than raising out of the hook.
    try:
        text = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    except OSError:
        return []
    sources: set[str] = set()
    for kind, body in INSTALL_RULE.findall(text):
        head, _, destination = body.partition("DESTINATION")
        if "DAYDREAM_DIR" not in destination:
            continue
        for source in INSTALL_SOURCE.findall(head):
            if kind == "FILES":
                sources.add(source)
                continue
            for pattern in INSTALL_PATTERN.findall(destination) or ["*"]:
                sources.update(
                    path.relative_to(ROOT).as_posix()
                    for path in (ROOT / source).rglob(pattern)
                    if path.is_file()
                )
    return sorted(sources)


def check_install_eol(paths: list[str]) -> list[str]:
    """Return one error per installed file carrying no line-ending pin.

    daydream's deploy gate byte-compares its copy of each installed file
    against this repository's, over LF bytes. Without an eol=lf (or binary)
    attribute the working copy holds whatever the host checked out, so an
    install run on Windows ships CRLF and forks the mirror -- which a Linux
    runner cannot reproduce.
    """
    if not paths:
        return ["CMakeLists.txt: no DAYDREAM_DIR install sources found"]
    # Bytes, not text=True: the text wrapper rewrites the input's newlines to
    # CRLF on Windows and git reads the CR as part of the path.
    query = subprocess.run(
        ["git", "-C", str(ROOT), "check-attr", "--stdin", "eol", "binary"],
        input="\n".join(paths).encode(), capture_output=True)
    if query.returncode != 0:
        return [f"git check-attr failed: {query.stderr.decode().strip()}"]
    attributes: dict[str, dict[str, str]] = {}
    for line in query.stdout.decode().splitlines():
        path, attribute, value = line.rsplit(": ", 2)
        attributes.setdefault(path, {})[attribute] = value
    errors: list[str] = []
    for path in paths:
        found = attributes.get(path, {})
        if found.get("eol") == "lf" or found.get("binary") == "set":
            continue
        errors.append(
            f"{path}: installed into daydream with no eol=lf pin "
            f"in .gitattributes"
        )
    return errors


def check_consumers() -> int:
    installed = installed_sources()
    errors: list[str] = (check_inline_pins() + check_shared_literals()
                         + check_engine_ranges() + check_flexram_geometry()
                         + check_install_eol(installed))
    for name in sorted(set(CHECK_TOOLS) - set(PINS | INLINE_PINS)):
        errors.append(f"CHECK_TOOLS names {name}, which is not a pin")
    for path, references in CONSUMERS.items():
        text = read_scanned(path, errors)
        if text is None:
            continue
        for reference in references:
            if reference not in text:
                errors.append(f"{path.relative_to(ROOT)}: missing {reference!r}")
    duplicate_paths = set(workflow_files()) | set(CONSUMERS)
    for path in sorted(duplicate_paths):
        text = read_scanned(path, errors)
        if text is None:
            continue
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
          f"{len(installed)} installed file)")
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
