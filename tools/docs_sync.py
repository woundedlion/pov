#!/usr/bin/env python3
"""Refresh documented repository maps and source-derived counts before validation."""

from __future__ import annotations

import fnmatch
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

import docs_check as dc


@dataclass
class Row:
    names: list[str]
    description: str
    continuation: list[str] = field(default_factory=list)
    children: list[Row] = field(default_factory=list)
    separators: list[str] = field(default_factory=list)


def parse_rows(body: list[str]) -> list[Row]:
    roots: list[Row] = []
    stack: list[Row] = []
    for line in body:
        match = dc.TREE_ROW_RE.match(line)
        if match is None:
            if stack and line.strip("│ "):
                stack[-1].continuation.append(line[len(stack) * 4:])
            elif stack:
                depth = max(0, (len(line.rstrip()) - 1) // 4)
                stack[min(depth, len(stack) - 1)].separators.append(line)
            continue
        depth = len(match["indent"]) // 4
        names = dc.tree_names(match["rest"])
        if not names or depth > len(stack):
            raise ValueError(f"invalid generated tree row: {line}")
        name_text = " / ".join(names)
        row = Row(names, match["rest"][len(name_text):])
        (roots if depth == 0 else stack[depth - 1].children).append(row)
        stack[depth:] = [row]
    return roots


def refresh_rows(rows: list[Row], entries: set[PurePosixPath],
                 allowed: tuple[str, ...], unmapped: tuple[str, ...],
                 parent: str = "", covered: set[str] | None = None) -> list[Row]:
    def exists(name: str) -> bool:
        candidate = f"{parent}/{name}".lstrip("/").rstrip("/")
        return (any(candidate == prefix.rstrip("/")
                    or candidate.startswith(prefix) for prefix in allowed)
                or any(fnmatch.fnmatchcase(path.as_posix(), candidate)
                       for path in entries))

    refreshed = []
    for row in rows:
        names = [name for name in row.names if exists(name)]
        for removed in set(row.names) - set(names):
            print(f"[docs-sync] removed absent tree path {parent}/{removed}".replace("path /", "path "))
        if not names:
            continue
        if names != row.names:
            width = len(" / ".join(row.names)) + len(row.description) - len(row.description.lstrip())
            row.description = " " * max(2, width - len(" / ".join(names))) + row.description.lstrip()
            row.names = names
        directory = f"{parent}/{names[0]}".strip("/")
        listed = set()
        if names[0].endswith("/"):
            def refresh_list(match: re.Match) -> str:
                stems = dc.tree_listed(match[0])
                if not stems:
                    return match[0]
                live = [stem for stem in stems if any(
                    path.parent.as_posix() == directory and path.stem == stem
                    for path in entries)]
                listed.update(path for stem in live
                              for path in dc.tree_stem_paths(directory, stem, entries))
                if live == stems:
                    return match[0]
                return "(" + ", ".join(live) + ")" if live else ""
            description = dc.TREE_LIST_RE.sub(refresh_list, "\n".join([row.description, *row.continuation]))
            row.description, *row.continuation = description.split("\n")
        if row.children or listed:
            row.children = refresh_rows(row.children, entries, allowed, unmapped, directory, listed)
        refreshed.append(row)

    drawn = {f"{parent}/{name}".strip("/") for row in refreshed for name in row.names}
    drawn.update(covered or ())
    drawn.update(ancestor.as_posix() for path in tuple(drawn)
                 for ancestor in PurePosixPath(path).parents if ancestor != PurePosixPath("."))
    for path in sorted(entries):
        if ("" if path.parent == PurePosixPath(".") else path.parent.as_posix()) != parent:
            continue
        candidate = path.as_posix()
        if (any(fnmatch.fnmatchcase(candidate, pattern) for pattern in drawn)
                or dc.tree_unmapped(candidate, unmapped)):
            continue
        raise ValueError(f"add a repository-map row with a role description for {candidate}")
    return refreshed


def render_rows(rows: list[Row], prefix: str = "") -> list[str]:
    lines = []
    for index, row in enumerate(rows):
        last = index == len(rows) - 1
        branch = "└── " if last else "├── "
        spine = prefix + ("    " if last else "│   ")
        lines.append(prefix + branch + " / ".join(row.names) + row.description)
        lines.extend(spine + line for line in row.continuation)
        lines.extend(render_rows(row.children, spine))
        lines.extend(row.separators)
    return lines


def sync_trees(text: str, entries: set[PurePosixPath],
               checkouts: dict[str, set[PurePosixPath]]) -> str:
    lines = text.splitlines(keepends=True)
    _, fences, issues = dc.visible_lines(PurePosixPath("README.md"), text)
    if issues:
        return text
    for fence in reversed(fences):
        directive = dc.tree_directive(fence.tag)
        if directive is None or not directive.exhaustive:
            continue
        target = checkouts.get(directive.checkout) if directive.checkout else entries
        if target is None:
            continue
        allowed = (dc.CHECKOUT_UNTRACKED_ALLOWED.get(directive.checkout, ())
                   if directive.checkout else dc.UNTRACKED_ALLOWED)
        unmapped = () if directive.checkout else dc.TREE_UNMAPPED
        rows = parse_rows([line for _, line in fence.body])
        replacement = render_rows(refresh_rows(rows, target, allowed, unmapped))
        if fence.body:
            lines[fence.body[0][0] - 1:fence.body[-1][0]] = [line + "\n" for line in replacement]
        else:
            lines[fence.start:fence.start] = [line + "\n" for line in replacement]
    return "".join(lines)


def replace_count(match: re.Match, counts: dict[str | int, int]) -> str:
    text = match[0]
    for group, count in sorted(counts.items(), key=lambda item: match.start(item[0]), reverse=True):
        if match[group] is not None:
            start, end = match.span(group)
            text = text[:start - match.start()] + str(count) + text[end - match.start():]
    return text


def source_counts(root: Path) -> dict[str, int]:
    header = root / dc.EFFECT_ROSTER_SOURCE
    playlist = root / dc.PHANTASM_PLAYLIST_SOURCE
    counts = {}
    if header.is_file() and playlist.is_file():
        source = header.read_text(encoding="utf-8")
        counts = {"HS_EFFECT_LIST": len(dc.effect_roster(source)),
                  "HS_PHANTASM_EFFECT_LIST": len(dc.phantasm_roster(playlist.read_text(encoding="utf-8"))),
                  "HS_SHADER_PRODUCT_GROUP": len(dc.shader_product_group(source))}
    return counts


def sync_text(relative: PurePosixPath, text: str, entries: set[PurePosixPath],
              checkouts: dict[str, set[PurePosixPath]], counts: dict[str, int]) -> str:
    if relative == PurePosixPath("README.md"):
        text = sync_trees(text, entries, checkouts)
        if counts:
            count = counts["HS_EFFECT_LIST"]
            headers = sum(entry.parent == dc.EFFECTS_DIR and entry.suffix == ".h" for entry in entries)
            text = dc.EFFECTS_ROW_RE.sub(lambda match: replace_count(match, {
                "headers": headers, "effects": count, "legacy_effects": count}), text)
            text = dc.EFFECTS_DIAGRAM_RE.sub(lambda match: replace_count(match, {"effects": count}), text)
    for document, pattern, macro, _ in dc.CARDINALITY_CLAIMS:
        if relative.as_posix() == document and macro in counts:
            text = pattern.sub(
                lambda match, expected=counts[macro]:
                    match[0] if dc.claimed_count(match[1]) == expected
                    else replace_count(match, {1: expected}),
                text)
    return text


def sync_repository(root: Path, checkout_roots: dict[str, Path],
                    revisions: dict[str, str]) -> None:
    markdown, entries = dc.tracked_entries(root)
    checkouts = {name: dc.tracked_entries(path, revisions.get(name))[1]
                 for name, path in checkout_roots.items()}
    counts = source_counts(root)
    for relative in markdown:
        path = root / relative
        before = path.read_text(encoding="utf-8-sig")
        after = sync_text(relative, before, entries, checkouts, counts)
        if after != before:
            path.write_text(after, encoding="utf-8", newline="\n")
            print(f"[docs-sync] updated {relative}")


def checkout_revisions(checkout_roots: dict[str, Path]) -> dict[str, str]:
    from build_pins import PINS
    return {name: PINS[name] for name in checkout_roots if name == "daydream"}


def discover_daydream(root: Path) -> Path | None:
    common = subprocess.check_output([
        "git", "-C", str(root), "rev-parse", "--path-format=absolute", "--git-common-dir"],
        text=True, timeout=30).strip()
    candidates = [Path(os.environ["DAYDREAM_CHECKOUT"])] if os.environ.get("DAYDREAM_CHECKOUT") else [
        root / "daydream", Path(common).parent.parent / "daydream"]
    revision = checkout_revisions({"daydream": root})["daydream"]
    for candidate in candidates:
        if candidate.is_dir():
            result = subprocess.run(["git", "-C", str(candidate), "cat-file", "-e", f"{revision}^{{tree}}"],
                                    capture_output=True, timeout=30)
            if result.returncode == 0:
                return candidate.resolve()
    print(f"[docs-sync] Daydream pin {revision} unavailable locally; leaving its map unchanged")
    return None
