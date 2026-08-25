#!/usr/bin/env python3
"""Verify documented image references; optionally stage them into Doxygen HTML."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_HTML_ROOT = ROOT / "build" / "docs" / "html"

_IMG_SRC_RE = re.compile(
    r"<img\b[^>]*?\bsrc[ \t]*=[ \t]*[\"']([^\"']*)[\"']", re.IGNORECASE)
_MARKDOWN_SUFFIXES = (".md", ".markdown")
_GIT_TIMEOUT_SECONDS = 30


def references(html_root: Path) -> list[tuple[Path, str]]:
    """Return every (html file, <img> src) pair in the generated tree."""
    found: list[tuple[Path, str]] = []
    for path in sorted(html_root.rglob("*.html")):
        text = path.read_text(encoding="utf-8", errors="replace")
        found.extend((path, src) for src in _IMG_SRC_RE.findall(text))
    return found


def markdown_references(
        repo_root: Path
) -> tuple[list[tuple[PurePosixPath, str]], list[str]]:
    """Return every (tracked Markdown file, <img> src) pair, and read errors."""
    command = ["git", "-c", f"safe.directory={repo_root.as_posix()}",
               "-C", str(repo_root), "ls-files", "-z"]
    found: list[tuple[PurePosixPath, str]] = []
    errors: list[str] = []
    try:
        result = subprocess.run(
            command, check=True, stdout=subprocess.PIPE,
            timeout=_GIT_TIMEOUT_SECONDS)
    except (OSError, subprocess.SubprocessError) as error:
        return found, [f"git ls-files failed: {error}"]
    for name in sorted(result.stdout.decode("utf-8").split("\0")):
        if not name or not name.casefold().endswith(_MARKDOWN_SUFFIXES):
            continue
        relative = PurePosixPath(name)
        # A tracked file can be absent from the working tree (an interrupted
        # checkout, a sparse one). That is an error to report, not a traceback.
        try:
            text = repo_root.joinpath(*relative.parts).read_text(
                encoding="utf-8", errors="replace")
        except OSError as error:
            errors.append(
                f"{relative.as_posix()}: tracked Markdown is unreadable: {error}")
            continue
        found.extend((relative, src) for src in _IMG_SRC_RE.findall(text))
    return found, errors


def verify(repo_root: Path) -> tuple[list[str], int]:
    """Resolve tracked Markdown references; return (errors, references checked)."""
    found, errors = markdown_references(repo_root)
    checked = 0
    for source, src in found:
        where = f"{source.as_posix()}: {src!r}"
        parts = urlsplit(src)
        if parts.scheme or parts.netloc:
            continue
        checked += 1
        if not parts.path:
            errors.append(f"{where} has no path component")
            continue
        decoded = unquote(parts.path)
        if decoded.startswith("/"):
            target = repo_root.joinpath(*PurePosixPath(decoded).parts[1:]).resolve()
        else:
            target = (repo_root.joinpath(*source.parent.parts) / decoded).resolve()
        if not target.is_relative_to(repo_root):
            errors.append(f"{where} resolves outside the repository")
        elif not target.is_file():
            errors.append(f"{where} names no file in the repository")
    return errors, checked


def stage(html_root: Path, repo_root: Path) -> tuple[list[str], int, int]:
    """Copy in missing repo-relative images.

    Returns (errors, staged count, references checked).
    """
    errors: list[str] = []
    staged: set[Path] = set()
    checked = 0
    for html_file, src in references(html_root):
        where = f"{html_file.relative_to(html_root).as_posix()}: {src!r}"
        parts = urlsplit(src)
        if parts.scheme or parts.netloc:
            continue
        checked += 1
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
    return errors, len(staged), checked


def report(errors: list[str]) -> None:
    for error in errors:
        print(f"::error::unresolvable docs image reference - {error}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", nargs="?", type=Path, const=DEFAULT_HTML_ROOT,
        metavar="HTML_ROOT",
        help="copy the referenced repository images into a generated Doxygen "
             "tree (default: build/docs/html) rather than verifying the "
             "tracked Markdown sources")
    parser.add_argument(
        "--repo-root", type=Path, default=ROOT,
        help="repository root the references are resolved against")
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()

    if args.stage is None:
        errors, checked = verify(repo_root)
        if errors:
            report(errors)
            return 1
        # No resolvable references means the checker was pointed somewhere
        # it cannot see the repository; passing would certify nothing.
        if not checked:
            print(f"[docs-images] tooling error: no repository-relative "
                  f"image references in tracked Markdown under "
                  f"{repo_root}", file=sys.stderr)
            return 2
        print(f"docs image references all resolve ({checked} checked)")
        return 0

    html_root = args.stage.resolve()
    if not html_root.is_dir():
        print(f"{html_root}: no generated HTML — build the docs first")
        return 1
    errors, staged, checked = stage(html_root, repo_root)
    if errors:
        report(errors)
        return 1
    # An artifact carrying no image reference is a gallery-less site;
    # publishing it would certify nothing.
    if not checked:
        print(f"[docs-images] tooling error: no repository-relative "
              f"image references under {html_root}", file=sys.stderr)
        return 2
    print(f"docs image references all resolve ({staged} staged from the repository)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
