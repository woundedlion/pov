#!/usr/bin/env python3
"""Check repository-root spelling of first-party C++ includes."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys

from git_files import git_nul_list

SOURCE_SUFFIXES = {".h", ".hpp", ".cpp", ".cc", ".inl", ".ino"}
LINKER_SCRIPTS = frozenset({"tools/phantasm.ld"})
CONFIG_HEADER = "core/vendor/FastNoiseLite_config.h"
INCLUDE = re.compile(r'^[ \t]*#[ \t]*include[ \t]*"([^"\n]+)"', re.MULTILINE)
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
COMMENTS = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)


def check_includes(root: Path, files: list[str]) -> list[str]:
    """Return diagnostics for noncanonical includes of tracked headers."""
    root = root.resolve()
    known = set(files)
    diagnostics = []
    for name in files:
        if Path(name).suffix not in SOURCE_SUFFIXES and name not in LINKER_SCRIPTS:
            continue
        if name.startswith("core/vendor/") and name != CONFIG_HEADER:
            continue
        path = root / name
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as exc:
            diagnostics.append(f"{name}:1: cannot read source ({exc})")
            continue
        text = COMMENTS.sub(lambda match: "\n" * match.group().count("\n"), text)
        linker = name in LINKER_SCRIPTS
        for match in (LINKER_INCLUDE if linker else INCLUDE).finditer(text):
            include = match.group(1) or match.group(2)
            candidates = [root / include] if linker else [
                path.parent / include, root / include,
                root / "core" / include, root / "effects" / include,
                root / "hardware" / include, root / "tests" / include]
            target = next((candidate.resolve() for candidate in candidates
                           if candidate.is_file()), None)
            if target is None or not target.is_relative_to(root):
                continue
            canonical = target.relative_to(root).as_posix()
            line = text.count("\n", 0, match.start()) + 1
            if name != "core/engine/effects_legacy.h":
                for prefix, forbidden in FORBIDDEN_PREFIXES.items():
                    in_layer = name.startswith(prefix) if prefix.endswith("/") else name == prefix
                    if in_layer and canonical.startswith(forbidden):
                        diagnostics.append(f'{name}:{line}: forbidden include "{canonical}" from {prefix}')
            if canonical in known and include != canonical:
                diagnostics.append(f'{name}:{line}: include "{canonical}" '
                                   f'instead of "{include}"')
    return diagnostics


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    args = parser.parse_args(argv)
    try:
        files = git_nul_list(args.root, "ls-files", "-z", timeout=None)
    except (OSError, subprocess.SubprocessError, UnicodeError) as error:
        print(f"[include-paths] tooling error: cannot list tracked sources ({error})",
              file=sys.stderr)
        return 2
    if not any(Path(name).suffix in SOURCE_SUFFIXES or name in LINKER_SCRIPTS for name in files):
        print("[include-paths] tooling error: no tracked C/C++ sources", file=sys.stderr)
        return 2
    diagnostics = check_includes(args.root, files)
    if diagnostics:
        print("\n".join(diagnostics))
        return 1
    print("[include-paths] PASS - first-party includes use repository-root paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
