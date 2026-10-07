#!/usr/bin/env python3
"""Single source for build values shared by CI, the pre-commit hook and just.

PINS are injected (--github-output / a named lookup) and must appear in no
build file literally. INLINE_PINS cannot be injected, so --check asserts every
spelling in INLINE_SCAN agrees; pip versions come from requirements/*.in.
SHARED_LITERALS are strings build files that cannot read one another must spell
identically. ENGINE_RANGES name a manifest floor the pin must satisfy. --check
also asserts every file CMakeLists.txt installs into daydream has an eol=lf pin.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import subprocess
import sys
from pathlib import Path

if sys.version_info < (3, 11):
    raise SystemExit("build_pins requires Python 3.11 or newer")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import teensy_gate  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]


def requirement_pin(name: str, package: str | None = None) -> str:
    """Read the single exact requirement from a tool's source manifest."""
    path = ROOT / "requirements" / f"{name}.in"
    package = package or name
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise SystemExit(f"{path.relative_to(ROOT)}: cannot be read ({error})")
    lines = [line.split("#", 1)[0].strip() for line in text.splitlines()]
    requirement = "\n".join(line for line in lines if line)
    found = re.fullmatch(rf"{re.escape(package)}==([0-9][\w.!+-]*)", requirement)
    if found is None:
        raise SystemExit(
            f"{path.relative_to(ROOT)}: expected one exact {package}==VERSION pin")
    return found.group(1)


PINS = {
    # PyPI's actionlint-py, whose version is the actionlint release plus a
    # packaging suffix.
    "actionlint": requirement_pin("actionlint", "actionlint-py"),
    # The daydream revision used to validate the exhaustive documentation tree.
    "daydream": "eb9d08e29e60a52bb2b5d9fc9162b4b1201d5302",
    "doxygen-awesome": "568f56cde6ac78b6dfcc14acd380b2e745c301ea",
    "emsdk": "5.0.0",
    # PyPI's rust-just, so the recipe runner is installed and held like ruff.
    "just": requirement_pin("just", "rust-just"),
    "node": "24.13.0",
    "platformio": requirement_pin("platformio"),
    "ruff": requirement_pin("ruff"),
    # PyPI's shellcheck-py, whose version is the shellcheck release plus a
    # packaging suffix.
    "shellcheck": requirement_pin("shellcheck", "shellcheck-py"),
}

# Versions the build files must spell out literally; --check asserts every
# occurrence equals the value here, so a partial bump fails.
INLINE_PINS = {
    "actionlint-release": PINS["actionlint"].rsplit(".", 1)[0],
    "actionlint-sha256":
        "8aca8db96f1b94770f1b0d72b6dddcb1ebb8123cb3712530b08cc387b349a3d8",
    "python": "3.12",
    "numpy": requirement_pin("numpy"),
    "clang": "22",
    "clang-format": requirement_pin("clang-format"),
    "cmake": requirement_pin("cmake"),
    "doxygen": "1.17.0",
    "doxygen-sha256":
        "75419ef4f446fc1c24ef12514b574e66e898ee6f527c6ae2ad84f91a905823c2",
    # apt.llvm.org's signing key.
    "llvm-key-sha256":
        "8b2a587ffd672c4687e7581dad4b2f6c1bb2ad6b480cd9771ba2ff48e0b8c75d",
    # Major of the hand-installed desktop KiCad the fab gates shell out to.
    "kicad": "10",
}

WORKFLOWS = ROOT / ".github/workflows"
ACTIONS = ROOT / ".github/actions"


def workflow_files() -> tuple[Path, ...]:
    """Every workflow and composite action under .github."""
    found = set(WORKFLOWS.glob("*.yml")) | set(WORKFLOWS.glob("*.yaml"))
    found |= set(ACTIONS.glob("*/action.yml"))
    found |= set(ACTIONS.glob("*/action.yaml"))
    return tuple(sorted(found))


CONSUMERS = {
    ROOT / ".github/workflows/ci.yml": (
        "build_pins.py --github-output",
        "python tools/build_pins.py --check",
        "python tools/build_pins.py --check-tool cmake",
    ),
    ROOT / ".github/actions/pinned-doxygen/action.yml": (
        "build_pins.py doxygen-awesome",
    ),
    ROOT / "justfile": (
        "build_pins.py doxygen-awesome",
        "build_pins.py --check-tool actionlint",
        "build_pins.py --check-tool just",
        "build_pins.py --check-tool clang-format",
        "build_pins.py --check-tool cmake",
        "build_pins.py --check-tool doxygen",
        "build_pins.py --check-tool node",
        "build_pins.py --check-tool numpy",
        "build_pins.py --check-tool platformio",
        "build_pins.py --check-tool ruff",
        "build_pins.py --check-tool shellcheck",
        '{{ python_command }} tools/build_pins.py --check',
    ),
    ROOT / ".githooks/pre-commit": (
        '"$PYTHON_BIN" "$SNAPSHOT/tools/build_pins.py" --check',
    ),
}

# Files scanned for INLINE_PINS occurrences.
INLINE_SCAN = (
    *workflow_files(),
    ROOT / "justfile",
    # The clang-format gate's pathspec and exclusion regex.
    ROOT / "tools/clang_format_gate.sh",
    ROOT / "tools/whitespace_gate.sh",
    ROOT / "tools/shellcheck_gate.sh",
    ROOT / "platformio.ini",
    ROOT / "CMakeLists.txt",
    ROOT / "scripts/generate_luts.py",
    ROOT / ".githooks/pre-commit",
    ROOT / "hardware/phantasm/gen/sexp.py",
    # Prose that contributors install from.
    ROOT / "README.md",
    ROOT / "CONTRIBUTING.md",
    ROOT / "hardware/phantasm/README.md",
    # Both halves of each requirements pair: the hand-edited .in carries the
    # pin, the pip-compile'd .txt repeats it above the hashes.
    *(ROOT / "requirements" / f"{stem}{suffix}"
      for stem in ("actionlint", "clang-format", "cmake", "just", "numpy",
                   "platformio", "ruff", "shellcheck")
      for suffix in (".in", ".txt")),
)

# Patterns and value transforms for each supported pin spelling.
INLINE_USES = (
    (r"/actionlint/releases/download/v([\d.]+)/", "actionlint-release", lambda v: v),
    (r"actionlint_([\d.]+)_linux_amd64", "actionlint-release", lambda v: v),
    (r"([0-9a-f]{64})  actionlint\.tar\.gz", "actionlint-sha256", lambda v: v),
    (r"\bactionlint-py==([\w.]+)", "actionlint", lambda v: v),
    # Quote-agnostic: setup-python's own README writes the input with double
    # quotes, which a single-quoted pattern reads as absent.
    (r"""python-version:\s*['"]?([^'"\s]+)['"]?""", "python", lambda v: v),
    (r"\bnumpy==([\w.]+)", "numpy", lambda v: v),
    (r"\b(?:clang\+\+|clang|llvm)-(\d+)\b", "clang", lambda v: v),
    (r"\bllvm-\w+-(\d+)\b", "clang", lambda v: v),
    (r"\bclang-format==([\w.]+)", "clang-format", lambda v: v),
    (r"\bcmake==([\w.]+)", "cmake", lambda v: v),
    (r"\bclang-format-(\d+)\b", "clang-format", lambda v: v.split(".")[0]),
    (r"\bclang-format (\d+)\b", "clang-format", lambda v: v.split(".")[0]),
    (r"\brust-just==([\w.]+)", "just", lambda v: v),
    (r"\bplatformio==([\w.]+)", "platformio", lambda v: v),
    (r"\bruff==([\w.]+)", "ruff", lambda v: v),
    (r"\bshellcheck-py==([\w.]+)", "shellcheck", lambda v: v),
    (r"EXPECTED_CLANG_FORMAT_MAJOR = (\d+)", "clang-format",
     lambda v: v.split(".")[0]),
    (r"HS_CLANG_FORMAT_MAJOR=(\d+)", "clang-format",
     lambda v: v.split(".")[0]),
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

# The shared clang-format gate selects these extensions through a git pathspec;
# the pre-commit hook matches the same set over staged paths.
FORMAT_EXTENSIONS = ("h", "hpp", "cpp", "cc", "inl", "ino")

# The float flags both shipping targets build with. -fno-finite-math-only must
# follow -ffast-math, which otherwise folds std::isfinite() to constant true.
FLOAT_FLAGS = ("-ffast-math", "-fno-finite-math-only")
FAST_MATH_TEST_FLAGS = (*FLOAT_FLAGS, "-DHS_TEST_FAST_MATH=1")

# Strings that must read identically in several build files that cannot source
# one another; --check asserts every occurrence matches this one.
SHARED_LITERALS = {
    "whitespace-rules": "blank-at-eol,blank-at-eof",
    "shell-selection": r"\.sh$|^\.githooks/",
    # Paths the clang-format gate skips: vendored FastNoiseLite and generated
    # tables.
    "format-exclude": (
        r"(^|/)core/vendor/FastNoiseLite\.h$"
        r"|(^|/)core/color/color_luts\.h$"
        r"|(^|/)core/color/gamut_lut\.h$"
        r"|(^|/)core/color/srgb_decode_lut\.h$"
        r"|(^|/)core/color/mindsplatter_palette_luts\.h$"
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
    # The per-effect smoke window every CI leg drives.
    "smoke-frames": "120",
    "smoke-stack-ceiling-debug": "6144",
}

# Patterns for literals shared across build tools.
SHARED_LITERAL_USES = (
    (r"core\.whitespace=([^\s]+)", "whitespace-rules"),
    (r"grep -E '(\\\.sh[^']*)'", "shell-selection"),
    (r"grep -vE '([^']*)'", "format-exclude"),
    (r"git ls-files -- ('\*\.h'(?: '\*\.\w+')*)", "format-globs"),
    # Anchored on the extension alternation's own opening, so an unrelated
    # quoted `grep -E` elsewhere in a scanned file is not counted as a copy.
    (r"grep -E '(\\\.\([^']*)'", "format-extensions"),
    # Anchored on a flag line's start (platformio.ini) and on ci.yml's matrix
    # key, so prose is not swept in; captures run to end of line so a partial
    # edit reads as a difference.
    (r"^\s+(-ffast-math\b.*)$", "float-flags"),
    (r"float_flags:\s+(-ffast-math\b.*)$", "float-test-flags"),
    # Every WASM target's compile and link line.
    (r'("-ffast-math" "[^"]+")', "float-flags-cmake"),
    # ci.yml anchor/alias, justfile recipe parameter, CONTRIBUTING prose.
    (r'HS_SMOKE_FRAMES(?:: &\w+ |="?)(\d+)', "smoke-frames"),
    (r'WASM_SMOKE_STACK_CEILING(?:: |=\"?)(\d+)', "smoke-stack-ceiling-debug"),
)

# --check-tool targets: pin name -> (version command, install hint, the form
# of the pin that command reports). `{pin}` is filled with the pin value.
CHECK_TOOLS = {
    "cmake": (["cmake", "--version"], "pip install cmake=={pin}", lambda v: v),
    "actionlint": (["actionlint", "-version"],
                   "pip install actionlint-py=={pin}",
                   lambda v: v.rsplit(".", 1)[0]),
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

# (manifest, JSON path to a `>=X` range, pin name); the pin must satisfy the
# range.
ENGINE_RANGES = (
    (ROOT / "package.json", ("engines", "node"), "node"),
)


# install(FILES|DIRECTORY) rules, matched to the first closing paren -- neither
# form nests one. The cross-repo install set is the subset destined for
# DAYDREAM_DIR; install(CODE) writes only files generated at install time.
INSTALL_RULE = re.compile(r"install\(\s*(FILES|DIRECTORY)\s(.*?)\)", re.DOTALL)
INSTALL_SOURCE = re.compile(r"\$\{CMAKE_CURRENT_SOURCE_DIR\}/([^\"\s]+)")
INSTALL_PATTERN = re.compile(
    r'\b(PATTERN|REGEX)\s+"([^"]+)"(?:\s+(EXCLUDE))?')


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
    """Read scanned text, or return None and record the read failure in errors."""
    try:
        return path.read_text(encoding="utf-8")
    except OSError as error:
        errors.append(f"{path.relative_to(ROOT)}: cannot be read ({error})")
        return None


INLINE_AUTHORITIES = {
    'cmake': ('requirements/cmake.in', 'requirements/cmake.txt'),
    'actionlint': (
        'requirements/actionlint.in',
        'requirements/actionlint.txt',
    ),
    'clang': (
        '.github/workflows/ci.yml',
        'README.md',
    ),
    'clang-format': (
        '.githooks/pre-commit',
        'CONTRIBUTING.md',
        'README.md',
        'requirements/clang-format.in',
        'requirements/clang-format.txt',
        'scripts/generate_luts.py',
        'tools/clang_format_gate.sh',
    ),
    'doxygen': (
        '.github/actions/pinned-doxygen/action.yml',
    ),
    'doxygen-sha256': (
        '.github/actions/pinned-doxygen/action.yml',
    ),
    'just': (
        'requirements/just.in',
        'requirements/just.txt',
    ),
    'kicad': (
        'README.md',
        'hardware/phantasm/README.md',
        'hardware/phantasm/gen/sexp.py',
    ),
    'llvm-key-sha256': (
        '.github/workflows/ci.yml',
    ),
    'numpy': (
        'requirements/numpy.in',
        'requirements/numpy.txt',
    ),
    'platformio': (
        'requirements/platformio.in',
        'requirements/platformio.txt',
    ),
    'python': (
        '.github/workflows/ci.yml',
        '.github/workflows/docs.yml',
    ),
    'ruff': (
        'requirements/ruff.in',
        'requirements/ruff.txt',
    ),
    'shellcheck': (
        'requirements/shellcheck.in',
        'requirements/shellcheck.txt',
    ),
}


LITERAL_AUTHORITIES = {
    'whitespace-rules': (
        '.githooks/pre-commit',
        'tools/whitespace_gate.sh',
    ),
    'smoke-stack-ceiling-debug': (
        '.github/workflows/ci.yml',
        'justfile',
    ),
    'shell-selection': (
        '.githooks/pre-commit',
        'tools/shellcheck_gate.sh',
    ),
    'float-flags-cmake': (
        'CMakeLists.txt',
    ),
    'format-exclude': (
        '.githooks/pre-commit',
        'tools/clang_format_gate.sh',
    ),
    'format-extensions': (
        '.githooks/pre-commit',
    ),
    'format-globs': (
        'tools/clang_format_gate.sh',
    ),
    'smoke-frames': (
        '.github/workflows/ci.yml',
        'CONTRIBUTING.md',
        'justfile',
    ),
}


def check_authorities(path, text, uses, authorities, errors):
    relative = path.relative_to(ROOT).as_posix()
    for name in {row[1] for row in uses}:
        if relative in authorities.get(name, ()) and not any(
                re.search(row[0], text) for row in uses if row[1] == name):
            errors.append(f"{relative}: missing {name} authority")


def check_inline_pins() -> list[str]:
    """Validate discovered pin spellings and report patterns with no matches."""
    errors: list[str] = []
    pin_values = {**PINS, **INLINE_PINS}
    seen: dict[str, int] = {pattern: 0 for pattern, _, _ in INLINE_USES}
    for path in INLINE_SCAN:
        text = read_scanned(path, errors)
        if text is None:
            continue
        check_authorities(path, text, INLINE_USES, INLINE_AUTHORITIES, errors)
        for index, line in enumerate(text.splitlines(), 1):
            for pattern, name, form in INLINE_USES:
                want = form(pin_values[name])
                for found in re.findall(pattern, line):
                    seen[pattern] += 1
                    if found != want:
                        errors.append(
                            f"{path.relative_to(ROOT)}:{index}: {name} pinned to "
                            f"{want!r} but written {found!r}"
                        )
    for pattern, name, _ in INLINE_USES:
        if not seen[pattern]:
            errors.append(
                f"{name} spelling {pattern!r} was not found in the scanned files"
            )
    return errors


def check_shared_literals() -> list[str]:
    """Validate shared literal values and report patterns with no matches."""
    errors: list[str] = []
    texts = {}
    for path in INLINE_SCAN:
        text = read_scanned(path, errors)
        if text is not None:
            texts[path] = text
    for pattern, name in SHARED_LITERAL_USES:
        want = SHARED_LITERALS[name]
        occurrences = 0
        for path, text in texts.items():
            check_authorities(path, text, ((pattern, name),),
                              LITERAL_AUTHORITIES, errors)
            for index, line in enumerate(text.splitlines(), 1):
                for found in re.findall(pattern, line):
                    occurrences += 1
                    if found != want:
                        errors.append(
                            f"{path.relative_to(ROOT)}:{index}: {name} differs "
                            f"from build_pins.py: {found!r}"
                        )
        if not occurrences:
            errors.append(
                f"{name} was not found in the scanned files"
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
        # BudgetSchemaError and JSON syntax errors are ValueErrors.
        return [f"tools/teensy_budgets.json: {exc}"]
    try:
        derived = budgets["phantasm"]["regions"]["ram1"][
            "components"]["code"]["max_banks_from_stack_floor"]
        bank_bytes = derived["bank_bytes"]
        total_banks = derived["total_banks"]
    except (KeyError, TypeError) as exc:
        return [f"tools/teensy_budgets.json: missing FlexRAM geometry: {exc}"]
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
    # Case-folded and whitespace-collapsed: a reflowed spelling is the same
    # geometry.
    def canonical(text: str) -> str:
        return re.sub(r"\s+", " ", text.lower())

    folded = canonical(linker)
    for spelling in linker_spellings:
        if canonical(spelling) not in folded:
            errors.append(
                f"tools/phantasm.ld: missing FlexRAM geometry derived from "
                f"teensy_budgets.json: {spelling!r}")
    return errors


def installed_sources() -> list[str]:
    """Return the repository files CMakeLists.txt installs into daydream.

    An install(DIRECTORY) rule contributes whatever its FILES_MATCHING patterns
    currently select.
    """
    # An unreadable CMakeLists.txt yields no install source.
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
            rules = INSTALL_PATTERN.findall(destination)
            inclusions = [(kind, pattern) for kind, pattern, excluded in rules
                          if not excluded]
            exclusions = [(kind, pattern) for kind, pattern, excluded in rules
                          if excluded]

            def matches(path: Path, rule: tuple[str, str]) -> bool:
                kind, pattern = rule
                spelling = path.as_posix()
                return (re.search(pattern, spelling) is not None if kind == "REGEX"
                        else fnmatch.fnmatchcase(spelling, f"*/{pattern}"))

            for path in (ROOT / source).rglob("*"):
                if (path.is_file()
                        and (not inclusions or any(matches(path, rule) for rule in inclusions))
                        and not any(matches(path, rule) for rule in exclusions)):
                    sources.add(path.relative_to(ROOT).as_posix())
    return sorted(sources)


def check_install_eol(paths: list[str]) -> list[str]:
    """Return one error per installed file carrying no line-ending pin.

    daydream's deploy gate byte-compares each installed file over LF bytes, so
    each needs an eol=lf (or binary) attribute.
    """
    if not paths:
        return ["CMakeLists.txt: no DAYDREAM_DIR install sources found"]
    # Bytes, not text=True: the text wrapper rewrites the input's newlines to
    # CRLF on Windows and git reads the CR as part of the path.
    try:
        query = subprocess.run(
            ["git", "-C", str(ROOT), "check-attr", "--stdin", "eol", "binary"],
            input="\n".join(paths).encode(), capture_output=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as error:
        return [f"git check-attr failed: {error}"]
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

    The pin's precision is the comparison's: a major-only pin such as clang's
    is met by any release of that major, which is all the pin claims.
    """
    pin = (PINS | INLINE_PINS)[name]
    command, install, form = CHECK_TOOLS[name]
    command = [part.format(pin=pin) for part in command]
    if name == "clang-format":
        command[0] = os.environ.get("CLANG_FORMAT", command[0])
    want = form(pin).split(".")
    try:
        reported = subprocess.run(
            command, capture_output=True, text=True, check=True, timeout=30
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
