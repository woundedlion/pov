#!/usr/bin/env python3
"""Fail when a Doxygen input symbol has no documentation comment.

The published Doxyfile sets EXTRACT_ALL, which silences doxygen's
undocumented-member warnings. This gate reruns doxygen with EXTRACT_ALL off,
drops the exempt symbol classes below, and reports each remaining source
location once (a template's member warns once per instantiation). Doxygen
skips global symbols of a file with no `@file` block, so such files fail too.
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

_WARNING_RE = re.compile(
    r"^(?P<path>.+?):(?P<line>\d+): warning: "
    r"(?:Member (?P<member>.+?) \((?P<kind>[\w ]+)\) of (?:\w+ )+(?P<scope>.+?)"
    r"|Compound (?P<compound>.+?)) is not documented\.$")
_ANY_WARNING_RE = re.compile(r"^(?P<path>.+?):(?P<line>\d+): warning: (?P<text>.*)$")

# Members every Pullback::Interp::Op operator model declares; the contract is
# documented once on the namespace.
OPERATOR_CONTRACT = frozenset({
    "ID", "NAME", "Input", "Output", "Params", "Prepared", "State",
    "TOPOLOGY", "FIELDS", "EDGE_DISTANCE_AVAILABLE", "ORACLE", "METRICS",
    "APPROXIMATE", "init", "prepare", "run", "advance", "validate_frame",
    "project",
})

# `Pullback::Spec` fields an effect's spec overrides; documented on the base.
SPEC_FIELDS = frozenset({
    "PROJECTION", "TRANSFER", "COVERAGE", "FIELD_COVERAGE", "HARMONY", "HUE",
    "BRIGHTNESS", "ANIMATED_PROJECTION", "SURFACE_PLACEMENT", "LensPolicy",
})

CONTAINER_ALIASES = frozenset({
    "value_type", "size_type", "difference_type", "reference",
    "const_reference", "pointer", "const_pointer", "iterator",
    "const_iterator",
})

_DEFAULTED_RE = re.compile(r"=\s*(?:default|delete)\s*$")


@dataclass(frozen=True)
class Undocumented:
    """One undocumented-symbol warning."""

    path: str  # Repository-relative POSIX path.
    line: int  # 1-based source line.
    name: str  # Member signature or compound name.
    kind: str  # Doxygen member kind, "compound", or "warning".
    scope: str  # Enclosing scope; empty for a compound.


def parse(text: str, root: Path) -> list[Undocumented]:
    """Return the source-file warnings in a doxygen warning log.

    A warning other than an undocumented symbol (a missing or unknown
    `@param`, say) has kind "warning" and its full text as the name.
    """
    prefix = root.resolve().as_posix().rstrip("/") + "/"
    found: list[Undocumented] = []
    for raw in text.splitlines():
        if raw[:1].isspace() and found and found[-1].kind == "warning":
            last = found.pop()
            found.append(Undocumented(last.path, last.line,
                                      f"{last.name} {raw.strip()}", "warning",
                                      ""))
            continue
        line = raw.strip().replace("\\", "/")
        match = _WARNING_RE.match(line)
        other = None if match else _ANY_WARNING_RE.match(line)
        if not match and not other:
            continue
        path = (match or other)["path"]
        if path.lower().startswith(prefix.lower()):
            path = path[len(prefix):]
        if other:
            found.append(Undocumented(path, int(other["line"]), other["text"],
                                      "warning", ""))
        elif match["compound"] is not None:
            found.append(Undocumented(path, int(match["line"]),
                                      match["compound"], "compound", ""))
        else:
            found.append(Undocumented(path, int(match["line"]), match["member"],
                                      match["kind"], match["scope"]))
    return found


def _base_name(signature: str) -> str:
    return re.split(r"[\s(<]", signature, maxsplit=1)[0]


_UNDEF_RE = re.compile(r"^\s*#\s*undef\s+(\w+)", re.M)


def undefined_macros(text: str) -> frozenset[str]:
    """Return the macros a source file `#undef`s."""
    return frozenset(_UNDEF_RE.findall(text))


def exempt(item: Undocumented, local_macros: frozenset[str] = frozenset()
           ) -> bool:
    """Whether the symbol belongs to a class that needs no comment of its own.

    `local_macros` names the macros the item's file `#undef`s; such a macro is
    file-local scaffolding (an X-macro helper, say).
    """
    name = _base_name(item.name)
    if item.kind == "macro definition" and name in local_macros:
        return True
    if item.kind == "function" and (_DEFAULTED_RE.search(item.name)
                                    or name.startswith("operator=")):
        return True
    if item.kind == "typedef" and name in CONTAINER_ALIASES:
        return True
    scope = item.scope
    if (scope.startswith("Pullback::Interp::Op::")
            and "::" not in scope[len("Pullback::Interp::Op::"):]
            and name in OPERATOR_CONTRACT):
        return True
    if scope.endswith("Spec") and scope != "Pullback::Spec" \
            and name in SPEC_FIELDS:
        return True
    # A trait specialization's `Type` result; the primary template is documented.
    if item.kind == "typedef" and name == "Type" and scope.endswith(">"):
        return True
    return False


def locations(items: list[Undocumented], root: Path | None = None
              ) -> list[Undocumented]:
    """Drop exempt items and keep one item per source location."""
    macros: dict[str, frozenset[str]] = {}

    def local_macros(path: str) -> frozenset[str]:
        if root is None:
            return frozenset()
        if path not in macros:
            source = root / path
            macros[path] = undefined_macros(source.read_text(
                encoding="utf-8", errors="replace")) if source.is_file() \
                else frozenset()
        return macros[path]

    seen: dict[tuple[str, int], Undocumented] = {}
    for item in items:
        if not exempt(item, local_macros(item.path)):
            seen.setdefault((item.path, item.line), item)
    return sorted(seen.values(), key=lambda i: (i.path, i.line))


def doxyfile_setting(text: str, key: str) -> list[str]:
    """Return the whitespace-separated values of a Doxyfile setting."""
    values: list[str] = []
    lines = iter(text.splitlines())
    for line in lines:
        name, sep, value = line.partition("=")
        if not sep or name.strip().rstrip("+").strip() != key:
            continue
        if not name.strip().endswith("+"):
            values = []
        while value.rstrip().endswith("\\"):
            values += value.rstrip()[:-1].split()
            value = next(lines, "")
        values += value.split()
    return values


def source_files(root: Path, paths: list[str]) -> list[str]:
    """Return the tracked Doxygen source inputs under `paths` or INPUT."""
    config = (root / "Doxyfile").read_text(encoding="utf-8")
    patterns = doxyfile_setting(config, "FILE_PATTERNS")
    excluded = [e.rstrip("/") for e in doxyfile_setting(config, "EXCLUDE")]
    inputs = paths or doxyfile_setting(config, "INPUT")
    listed = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z", "--", *inputs],
        check=True, capture_output=True).stdout.decode("utf-8").split("\0")
    return [
        path for path in listed
        if path and any(fnmatch.fnmatch(Path(path).name, p) for p in patterns)
        and not any(path == e or path.startswith(e + "/") for e in excluded)
    ]


_FILE_COMMAND_RE = re.compile(r"[@\\]file\b")


def missing_file_blocks(root: Path, files: list[str]) -> list[str]:
    """Return the files without an `@file` documentation block."""
    return [path for path in files if not _FILE_COMMAND_RE.search(
        (root / path).read_text(encoding="utf-8", errors="replace"))]


def run_doxygen(doxygen: str, root: Path, paths: list[str]) -> str:
    """Run doxygen over the Doxyfile inputs and return its warning log."""
    with tempfile.TemporaryDirectory() as scratch:
        log = Path(scratch) / "warnings.txt"
        overlay = [
            "EXTRACT_ALL = NO",
            "WARN_IF_UNDOCUMENTED = YES",
            "WARN_NO_PARAMDOC = YES",
            "WARN_IF_INCOMPLETE_DOC = YES",
            "WARN_AS_ERROR = NO",
            "GENERATE_HTML = NO",
            "GENERATE_LATEX = NO",
            "GENERATE_XML = NO",
            "HAVE_DOT = NO",
            f"OUTPUT_DIRECTORY = {Path(scratch).as_posix()}",
            f"WARN_LOGFILE = {log.as_posix()}",
        ]
        if paths:
            overlay.append("INPUT = " + " ".join(paths))
        config = (root / "Doxyfile").read_text(encoding="utf-8")
        config += "\n" + "\n".join(overlay) + "\n"
        subprocess.run([doxygen, "-"], input=config, text=True, cwd=root,
                       check=True, stdout=subprocess.DEVNULL)
        return log.read_text(encoding="utf-8", errors="replace")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="*",
                        help="limit the scan to these inputs (default: Doxyfile INPUT)")
    parser.add_argument("--doxygen", default=os.environ.get("DOXYGEN", "doxygen"),
                        help="doxygen executable (default: $DOXYGEN or doxygen)")
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    files = missing_file_blocks(args.root, source_files(args.root, args.paths))
    items = locations(parse(run_doxygen(args.doxygen, args.root, args.paths),
                            args.root), args.root)
    for path in files:
        print(f"{path}:1: missing @file block")
    for item in items:
        where = f"{item.path}:{item.line}"
        if item.kind == "warning":
            print(f"{where}: {item.name}")
            continue
        what = item.name if not item.scope else f"{item.scope}::{item.name}"
        print(f"{where}: undocumented {item.kind} {what}")
    if files or items:
        print(f"[doc-coverage] FAIL - {len(files)} file(s) without @file, "
              f"{len(items)} undocumented or incomplete symbol(s)", file=sys.stderr)
        return 1
    print("[doc-coverage] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
