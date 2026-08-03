#!/usr/bin/env python3
"""Stage repo-relative images into the Doxygen HTML tree and verify references."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_HTML_ROOT = ROOT / "build" / "docs" / "html"

_IMG_SRC_RE = re.compile(
    r"<img\b[^>]*?\bsrc[ \t]*=[ \t]*[\"']([^\"']*)[\"']", re.IGNORECASE)


def references(html_root: Path) -> list[tuple[Path, str]]:
    """Return every (html file, <img> src) pair in the generated tree."""
    found: list[tuple[Path, str]] = []
    for path in sorted(html_root.rglob("*.html")):
        text = path.read_text(encoding="utf-8", errors="replace")
        found.extend((path, src) for src in _IMG_SRC_RE.findall(text))
    return found


def resolve(html_root: Path, repo_root: Path) -> tuple[list[str], int]:
    """Copy in missing repo-relative images; return (errors, staged count)."""
    errors: list[str] = []
    staged: set[Path] = set()
    for html_file, src in references(html_root):
        where = f"{html_file.relative_to(html_root).as_posix()}: {src!r}"
        parts = urlsplit(src)
        if parts.scheme or parts.netloc:
            continue
        if not parts.path:
            errors.append(f"{where} has no path component")
            continue
        target = (html_file.parent / unquote(parts.path)).resolve()
        if not target.is_relative_to(html_root):
            errors.append(f"{where} resolves outside the artifact")
            continue
        if target.is_file():
            continue
        source = repo_root / target.relative_to(html_root)
        if not source.is_file():
            errors.append(f"{where} is in neither the artifact nor the repository")
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        staged.add(target)
    return errors, len(staged)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "html_root", nargs="?", type=Path, default=DEFAULT_HTML_ROOT,
        help="generated Doxygen HTML directory (default: build/docs/html)")
    parser.add_argument(
        "--repo-root", type=Path, default=ROOT,
        help="repository root the unresolved references are staged from")
    args = parser.parse_args()

    html_root = args.html_root.resolve()
    if not html_root.is_dir():
        print(f"{html_root}: no generated HTML — build the docs first")
        return 1

    errors, staged = resolve(html_root, args.repo_root.resolve())
    if errors:
        for error in errors:
            print(f"::error::unresolvable docs image reference - {error}")
        return 1
    print(f"docs image references all resolve ({staged} staged from the repository)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
