#!/usr/bin/env python3
"""Check first-party C++ include spelling and layer dependency direction."""

from __future__ import annotations

import argparse
from pathlib import Path
import posixpath
import re
import subprocess
import sys

from doc_coverage import CPP_LEXEMES
from git_files import git_nul_list, staged_blobs

SOURCE_SUFFIXES = {".h", ".hpp", ".cpp", ".cc", ".inl", ".ino"}
LINKER_SCRIPTS = frozenset({"tools/phantasm.ld"})
CONFIG_HEADER = "core/vendor/FastNoiseLite_config.h"
INCLUDE = re.compile(r'^[ \t]*#[ \t]*include[ \t]*(?:"([^"\n]+)"|<([^>\n]+)>)', re.MULTILINE)
LINKER_INCLUDE = re.compile(r'^[ \t]*INCLUDE[ \t]+(?:"([^"\n]+)"|([^\s;]+))', re.MULTILINE)
FORBIDDEN_PREFIXES = {
    "core/": ("effects/", "targets/", "workbench/", "hardware/", "tests/", "tools/"),
    "core/pattern/": ("core/color/", "core/render/", "core/control/", "core/animation/"),
    "core/pattern/types.h": ("core/pattern/field.h", "core/pattern/sources.h",
                             "core/pattern/noise.h", "core/pattern/sphere_lut.h",
                             "core/vendor/FastNoiseLite.h", "core/math/noise_field.h"),
    "core/pattern/field.h": ("core/pattern/sources.h", "core/pattern/noise.h",
                             "core/pattern/sphere_lut.h", "core/vendor/FastNoiseLite.h"),
    "core/pattern/sources.h": ("core/pattern/noise.h", "core/vendor/FastNoiseLite.h",
                               "core/math/noise_field.h"),
    "core/pattern/sphere_lut.h": ("core/pattern/noise.h", "core/pattern/sources.h",
                                  "core/vendor/FastNoiseLite.h", "core/math/noise_field.h"),
    "core/color/modulation.h": ("core/pattern/field.h", "core/pattern/sources.h",
                                "core/pattern/noise.h", "core/pattern/sphere_lut.h",
                                "core/vendor/FastNoiseLite.h", "core/math/noise_field.h",
                                "core/render/", "core/control/"),
    "core/color/modulation_types.h": ("core/pattern/field.h", "core/pattern/sources.h",
                                      "core/pattern/noise.h", "core/pattern/sphere_lut.h",
                                      "core/vendor/FastNoiseLite.h", "core/math/noise_field.h",
                                      "core/render/", "core/control/", "core/color/"),
    "hardware/": ("effects/", "targets/", "workbench/", "tests/"),
    "effects/": ("targets/", "workbench/", "hardware/", "tests/"),
    "workbench/": ("effects/", "targets/"),
    "tools/pc_sampler.h": ("effects/", "targets/", "workbench/", "tests/"),
    "tools/profile_isr_window.h": ("effects/", "targets/", "workbench/", "tests/"),
    "tools/profile_clock.h": ("effects/", "targets/", "workbench/", "tests/"),
    "tools/phantasm.ld": ("effects/", "targets/", "workbench/", "tests/"),
}


def check_includes(root: Path, files: list[str],
                   sources: dict[str, bytes] | None = None) -> list[str]:
    """Return diagnostics for noncanonical includes of tracked headers.

    Includes resolve against `files`, the tracked paths. `sources` maps the
    paths to check to their content; by default every path in `files` is read
    from `root`.
    """
    known = set(files)
    diagnostics = []
    for name in files if sources is None else sorted(sources):
        if Path(name).suffix not in SOURCE_SUFFIXES and name not in LINKER_SCRIPTS:
            continue
        if name.startswith("core/vendor/") and name != CONFIG_HEADER:
            continue
        try:
            data = (root / name).read_bytes() if sources is None else sources[name]
            text = data.decode("utf-8")
        except (OSError, UnicodeError) as exc:
            diagnostics.append(f"{name}:1: cannot read source ({exc})")
            continue
        text = CPP_LEXEMES.sub(
            lambda m: re.sub(r"[^\n]", " ", m[0])
            if m["comment"] is not None or "\n" in m[0] else m[0], text)
        linker = name in LINKER_SCRIPTS
        for match in (LINKER_INCLUDE if linker else INCLUDE).finditer(text):
            include = match.group(1) or match.group(2)
            angle = not linker and match.group(2) is not None
            candidates = [include] if linker or angle else [
                posixpath.join(posixpath.dirname(name), include), include,
                "core/" + include, "effects/" + include,
                "hardware/" + include, "tests/" + include]
            canonical = next((resolved for resolved in map(posixpath.normpath, candidates)
                              if resolved in known), None)
            if canonical is None:
                continue
            line = text.count("\n", 0, match.start()) + 1
            if name != "core/engine/effects_legacy.h":
                for prefix, forbidden in FORBIDDEN_PREFIXES.items():
                    in_layer = name.startswith(prefix) if prefix.endswith("/") else name == prefix
                    if in_layer and canonical.startswith(forbidden):
                        diagnostics.append(f'{name}:{line}: forbidden include "{canonical}" from {prefix}')
            if include != canonical:
                diagnostics.append(f'{name}:{line}: include "{canonical}" '
                                   f'instead of "{include}"')
    return diagnostics


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--staged", action="store_true",
                        help="check the staged content of the paths read from stdin")
    args = parser.parse_args(argv)
    try:
        files = git_nul_list(args.root, "ls-files", "-z", timeout=None)
        sources = (staged_blobs(args.root, sys.stdin.read().splitlines())
                   if args.staged else None)
    except (OSError, subprocess.SubprocessError, UnicodeError) as error:
        print(f"[include-paths] tooling error: cannot read tracked sources ({error})",
              file=sys.stderr)
        return 2
    if not any(Path(name).suffix in SOURCE_SUFFIXES or name in LINKER_SCRIPTS for name in files):
        print("[include-paths] tooling error: no tracked C/C++ sources", file=sys.stderr)
        return 2
    diagnostics = check_includes(args.root, files, sources)
    if diagnostics:
        print("\n".join(diagnostics))
        return 1
    print("[include-paths] PASS - first-party includes use repository-root paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
