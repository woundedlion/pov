#!/usr/bin/env python3
"""Validate fences, links, anchors, and path claims in tracked Markdown.

Structure only. A green run means every fence closes, every link, anchor and
backticked repo path resolves, and every tree fence matches the tracked tree --
not that the prose is true. A wrong number in a sentence, a renamed symbol in a
table, and any path written without backticks or a link are all invisible here.
"""

from __future__ import annotations

import argparse
import fnmatch
import functools
import posixpath
import re
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlsplit


@dataclass(frozen=True, order=True)
class Issue:
    path: str
    line: int
    message: str

    def __str__(self) -> str:
        return f"{self.path}:{self.line}: {self.message}"


_FENCE_RE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
_REFERENCE_DEFINITION_RE = re.compile(r"^ {0,3}\[([^]]+)\]:[ \t]*(.*)$")
_MARKDOWN_ESCAPE_RE = re.compile(r"\\([!\"#$%&'()*+,\-./:;<=>?@\[\\\]^_`{|}~])")

_ATX_HEADING_RE = re.compile(r"^ {0,3}#{1,6}(?:[ \t]+(.*?))?[ \t]*$")
_HTML_ANCHOR_RE = re.compile(
    r"<a\b[^>]*?\b(?:id|name)[ \t]*=[ \t]*[\"']([^\"']+)[\"']", re.IGNORECASE)
_EXPLICIT_HEADING_ID_RE = re.compile(r"[ \t]*\{#([^}\s]+)\}[ \t]*$")
_HTML_TAG_RE = re.compile(r"<[^>]*>")
_INLINE_LINK_TEXT_RE = re.compile(r"!?\[([^\]]*)\]\([^)]*\)|!?\[([^\]]*)\]\[[^\]]*\]")
_SLUG_DROP_RE = re.compile(r"[^\w\- ]", re.UNICODE)
# GitHub renders #L12 / #L12-L20 as a line range on a blob, not a document
# anchor. No structural check can tell whether the line still holds its subject,
# so a repo-relative link is not allowed to pin one: link the file and name the
# symbol in the link text, which a reader can verify and a rename breaks loudly.
_LINE_FRAGMENT_RE = re.compile(r"^L\d+(?:-L?\d+)?$", re.IGNORECASE)

# A backticked token is linted as a repo path only when it carries one of these
# suffixes; anything else in backticks is prose, an identifier, or a command.
_SOURCE_SUFFIXES = frozenset({
    ".c", ".cfg", ".cmake", ".cpp", ".css", ".h", ".hpp", ".html", ".ini",
    ".ino", ".js", ".json", ".kicad_mod", ".kicad_pcb", ".kicad_pro",
    ".kicad_sch", ".kicad_sym", ".ld", ".markdown", ".md", ".mjs", ".png",
    ".py", ".sh", ".svg", ".toml", ".ts", ".txt", ".wasm", ".wrl", ".yaml",
    ".yml",
})
# Optional trailing :line or :line-line, as the ledgers cite source spans.
_PATH_SPAN_RE = re.compile(r"^([A-Za-z0-9_.][\w.\-/]*)(?::\d+(?:-\d+)?)?$")

# Directory trees drawn inside a fence are path claims too, and the largest
# ones in the docs. The preceding directive says which repository a tree draws
# and how completely: `tree [<checkout>] [exhaustive]`. Bare `tree` is this
# repository, rooted at its root; `tree <checkout>` is a sibling repository,
# validated against the root --checkout supplies, and left unvalidated only
# when --skip-checkout names it. `exhaustive` adds the
# reverse direction — every tracked path under a drawn directory must have a
# row. An HTML comment carries the directive, so neither GitHub nor Doxygen
# renders it.
_DIRECTIVE_RE = re.compile(r"^ {0,3}<!--[ \t]*docs-check:[ \t]*(.*?)[ \t]*-->[ \t]*$")
_TREE_TAG = "tree"
_TREE_EXHAUSTIVE = "exhaustive"
_TREE_ROW_RE = re.compile(r"^(?P<indent>(?:│   |    )*)(?:├──|└──) +(?P<rest>\S.*)$")
_TREE_INDENT = 4
# A directory row may enumerate its children in prose rather than draw one row
# apiece, wrapping onto continuation lines that carry the spine but no branch.
# Within such a row every parenthesized group whose comma-separated parts are
# all bare stems names children of it, gated in both directions like a row.
_TREE_CONT_RE = re.compile(r"^[│ \t]*(?P<rest>\S.*)$")
_TREE_LIST_RE = re.compile(r"\(([^()]*)\)")
_TREE_STEM_RE = re.compile(r"[a-z0-9_]+")
# A drawn name: a path segment chain, optionally a directory's trailing slash,
# optionally a glob. Bare ellipsis rows elide a subtree and name nothing.
_TREE_NAME_RE = re.compile(r"^(?![.…]+/?$)[A-Za-z0-9_.*?][\w.\-*?]*(?:/[\w.\-*?]+)*/?$")

# Path prefixes the docs cite that this repository will never track: the files
# that live in the sibling daydream repository, where the same paths are real.
_SELF_REPO_HOSTS = frozenset({"github.com", "www.github.com"})
# The README is installed into the sibling daydream checkout, where relative
# paths break, so it cites this repository through absolute GitHub URLs. They
# name tracked paths and are validated like any repo-relative link.
_SELF_REPO_PATH_RE = re.compile(
    r"^/woundedlion/pov/(?:blob|tree|raw)/[^/]+/(.+)$")

# The `effects/` row summarizes its subtree instead of drawing it, so no tree
# gate sees the counts it states. Both are derivable: headers from the tracked
# tree, effects from the roster macro's cardinality.
_EFFECTS_TREE_ROW = "README.md"
_EFFECTS_ROW_RE = re.compile(
    r"\beffects/\s+(?P<headers>\d+) headers(?: covering (?P<effects>\d+) effects|: "
    r"one per effect \((?P<legacy_effects>\d+)\))")
# The architecture diagram restates the roster's cardinality in its own words,
# outside any tree fence and in a spelling the summary row's regex never
# matches, so it is derivable from the same source and gated the same way.
_EFFECTS_DIAGRAM_RE = re.compile(
    r"\beffects/\s+\((?P<effects>\d+) visual algorithms\)")
_EFFECTS_DIR = PurePosixPath("effects")
_EFFECT_ROSTER_SOURCE = PurePosixPath("core/engine/effects.h")
_EFFECT_ROSTER_DEFINE = "#define HS_EFFECT_LIST(X)"
# Shared X-row spelling with scripts/effect_roster.mjs and tools/profile_sweep.sh:
# whitespace inside the parens is tolerated so a reformat to `X( Foo )` cannot
# silently shrink one roster reader's count while the others keep the full set.
# Anchoring at the line start is what excludes a commented-out `// X(Foo)` or
# `/* X(Foo) */` row, matching the comment stripping effect_roster.mjs does.
_EFFECT_ROSTER_ENTRY_RE = re.compile(r"^\s*X\(\s*(\w+)\s*\)")

_UNTRACKED_ALLOWED = (
    ".github/workflows/deploy.yml",
)

_MARKDOWN_EXCLUDED = frozenset({PurePosixPath("docs/CODE_REVIEW.md")})

# Tracked paths an exhaustive tree deliberately leaves without a row: VCS
# metadata, the map's own document, and the test tree the map draws as one
# summary row plus its shared fixtures. A trailing slash covers a subtree.
_TREE_UNMAPPED = (
    ".gitattributes",
    ".gitignore",
    "README.md",
    "hardware/phantasm/.gitignore",
    "tests/",
)

# Rows a checkout's own repository gitignores, keyed by the directive's name.
_CHECKOUT_UNTRACKED_ALLOWED = {
    "daydream": ("node_modules/", "three.js/", "vendor/"),
}


def _untracked_allowance(candidate: str, used: set[str] | None) -> bool:
    """Reports whether an _UNTRACKED_ALLOWED prefix exempts candidate, recording it."""
    for prefix in _UNTRACKED_ALLOWED:
        if candidate.startswith(prefix):
            if used is not None:
                used.add(prefix)
            return True
    return False


def _stale_allowances(entries: set[PurePosixPath], used: set[str]) -> list[str]:
    """Names _UNTRACKED_ALLOWED entries the exemption no longer buys anything for."""
    stale = []
    for prefix in _UNTRACKED_ALLOWED:
        if PurePosixPath(prefix.rstrip("/")) in entries:
            stale.append(f"{prefix} (now tracked)")
        elif prefix not in used:
            stale.append(f"{prefix} (uncited)")
    return stale


def _normalize_label(label: str) -> str:
    return " ".join(label.split()).casefold()


def _is_escaped(text: str, index: int) -> bool:
    backslashes = 0
    index -= 1
    while index >= 0 and text[index] == "\\":
        backslashes += 1
        index -= 1
    return backslashes % 2 == 1


def _scan_code_spans(line: str) -> tuple[str, list[str]]:
    """Blanks every code span and returns the masked line plus the span bodies."""
    masked = list(line)
    spans = []
    index = 0
    while index < len(line):
        if line[index] != "`" or _is_escaped(line, index):
            index += 1
            continue
        end_run = index
        while end_run < len(line) and line[end_run] == "`":
            end_run += 1
        marker = line[index:end_run]
        close = line.find(marker, end_run)
        if close == -1:
            index = end_run
            continue
        spans.append(line[end_run:close].strip())
        for offset in range(index, close + len(marker)):
            masked[offset] = " "
        index = close + len(marker)
    return "".join(masked), spans


def _is_fence_close(line: str, marker: str) -> bool:
    stripped = line.lstrip(" ")
    if len(line) - len(stripped) > 3 or not stripped.startswith(marker[0]):
        return False
    run = len(stripped) - len(stripped.lstrip(marker[0]))
    return run >= len(marker) and not stripped[run:].strip()


@dataclass(frozen=True)
class VisibleLine:
    number: int
    raw: str
    masked: str
    spans: tuple[str, ...]


@dataclass(frozen=True)
class Fence:
    tag: str
    start: int
    body: tuple[tuple[int, str], ...]


def _visible_lines(path: PurePosixPath,
                   text: str) -> tuple[list[VisibleLine], list[Fence], list[Issue]]:
    visible = []
    fences = []
    issues = []
    fence: tuple[str, int, str] | None = None
    body: list[tuple[int, str]] = []
    pending = ""
    for line_number, line in enumerate(text.splitlines(), 1):
        if fence:
            if _is_fence_close(line, fence[0]):
                fences.append(Fence(fence[2], fence[1], tuple(body)))
                fence = None
                body = []
            else:
                body.append((line_number, line))
            continue
        match = _FENCE_RE.match(line)
        if match and not (match.group(1)[0] == "`" and "`" in match.group(2)):
            fence = (match.group(1), line_number, pending)
            pending = ""
            continue
        directive = _DIRECTIVE_RE.match(line)
        if directive:
            # A directive carries to the next fence across blank lines only.
            pending = directive.group(1)
        elif line.strip():
            pending = ""
        masked, spans = _scan_code_spans(line)
        visible.append(VisibleLine(line_number, line, masked, tuple(spans)))
    if fence:
        issues.append(Issue(path.as_posix(), fence[1], "unclosed fenced code block"))
    return visible, fences, issues


def _slug(heading: str) -> str:
    """Slugifies heading text the way GitHub derives its anchor ids."""
    text = _INLINE_LINK_TEXT_RE.sub(lambda m: m.group(1) or m.group(2) or "", heading)
    text = _HTML_TAG_RE.sub("", text).replace("`", "")
    text = _MARKDOWN_ESCAPE_RE.sub(r"\1", text)
    text = re.sub(r"[*~]", "", text).strip().casefold()
    return _SLUG_DROP_RE.sub("", text).replace(" ", "-")


def _anchors(visible: list[VisibleLine]) -> set[str]:
    """Collects every fragment a link may target: heading slugs and HTML ids."""
    found: set[str] = set()
    seen: dict[str, int] = {}
    for line in visible:
        found.update(name.casefold() for name in _HTML_ANCHOR_RE.findall(line.raw))
        match = _ATX_HEADING_RE.match(line.raw)
        if not match:
            continue
        heading = match.group(1) or ""
        explicit = _EXPLICIT_HEADING_ID_RE.search(heading)
        if explicit:
            found.add(explicit.group(1).casefold())
            heading = heading[:explicit.start()]
        heading = heading.rstrip("#").rstrip()
        slug = _slug(heading)
        if not slug:
            continue
        # Repeated headings get GitHub's -1, -2, ... disambiguating suffix.
        count = seen.get(slug, 0)
        seen[slug] = count + 1
        found.add(slug if count == 0 else f"{slug}-{count}")
    return found


def _find_closing_bracket(text: str, start: int) -> int:
    depth = 1
    for index in range(start + 1, len(text)):
        if _is_escaped(text, index):
            continue
        if text[index] == "[":
            depth += 1
        elif text[index] == "]":
            depth -= 1
            if depth == 0:
                return index
    return -1


def _find_closing_parenthesis(text: str, start: int) -> int:
    depth = 1
    in_angle = False
    for index in range(start + 1, len(text)):
        if _is_escaped(text, index):
            continue
        char = text[index]
        if char == "<":
            in_angle = True
        elif char == ">":
            in_angle = False
        elif not in_angle and char == "(":
            depth += 1
        elif not in_angle and char == ")":
            depth -= 1
            if depth == 0:
                return index
    return -1


def _destination(body: str) -> str | None:
    body = body.lstrip()
    if not body:
        return ""
    if body.startswith("<"):
        close = body.find(">", 1)
        return None if close == -1 else body[1:close]

    depth = 0
    end = 0
    while end < len(body):
        char = body[end]
        if char == "\\" and end + 1 < len(body):
            end += 2
            continue
        if char.isspace() and depth == 0:
            break
        if char == "(":
            depth += 1
        elif char == ")" and depth:
            depth -= 1
        end += 1
    return body[:end]


def _inline_links(line: str) -> tuple[list[str], list[str]]:
    destinations = []
    references = []
    index = 0
    while index < len(line):
        if line[index] != "[" or _is_escaped(line, index):
            index += 1
            continue
        close = _find_closing_bracket(line, index)
        if close == -1:
            index += 1
            continue
        label = line[index + 1:close]
        next_index = close + 1
        if next_index < len(line) and line[next_index] == "(":
            end = _find_closing_parenthesis(line, next_index)
            if end != -1:
                destination = _destination(line[next_index + 1:end])
                if destination is not None:
                    destinations.append(destination)
                index = end + 1
                continue
        elif next_index < len(line) and line[next_index] == "[":
            ref_close = _find_closing_bracket(line, next_index)
            if ref_close != -1:
                reference = line[next_index + 1:ref_close] or label
                references.append(_normalize_label(reference))
                index = ref_close + 1
                continue
        index = close + 1
    return destinations, references


def _resolved_target(source: PurePosixPath, target: str) -> PurePosixPath | None:
    target = _MARKDOWN_ESCAPE_RE.sub(r"\1", target.strip())
    if not target or target.startswith("#"):
        return None
    parsed = urlsplit(target)
    if parsed.scheme or parsed.netloc:
        if parsed.scheme.casefold() not in ("http", "https"):
            return None
        if parsed.netloc.casefold() not in _SELF_REPO_HOSTS:
            return None
        match = _SELF_REPO_PATH_RE.match(unquote(parsed.path))
        if not match:
            return None
        return PurePosixPath(posixpath.normpath(match.group(1)))
    decoded = unquote(parsed.path)
    if not decoded:
        return None
    if decoded.startswith("/"):
        resolved = posixpath.normpath(decoded.lstrip("/"))
    else:
        resolved = posixpath.normpath(posixpath.join(source.parent.as_posix(), decoded))
    return PurePosixPath(resolved)


def _fragment(target: str) -> str:
    parsed = urlsplit(_MARKDOWN_ESCAPE_RE.sub(r"\1", target.strip()))
    return unquote(parsed.fragment)


def _anchor_issue(source: PurePosixPath, line: int, target: str,
                  document: PurePosixPath,
                  anchors: dict[PurePosixPath, set[str]]) -> Issue | None:
    fragment = _fragment(target)
    # Only tracked Markdown has resolvable anchors; a blob line range is not one.
    if not fragment or document not in anchors:
        return None
    if _LINE_FRAGMENT_RE.match(fragment):
        return None
    known = anchors[document]
    if fragment.casefold() in known or _slug(fragment) in known:
        return None
    where = "this document" if document == source else document.as_posix()
    return Issue(source.as_posix(), line,
                 f"missing anchor {'#' + fragment!r} in {where}")


def _link_issue(source: PurePosixPath, line: int, target: str,
                entries: set[PurePosixPath],
                anchors: dict[PurePosixPath, set[str]]) -> Issue | None:
    resolved = _resolved_target(source, target)
    if resolved is None:
        cleaned = _MARKDOWN_ESCAPE_RE.sub(r"\1", target.strip())
        if cleaned.startswith("#"):
            return _anchor_issue(source, line, target, source, anchors)
        return None
    # A target above the root can never resolve on GitHub, whatever a local
    # checkout layout happens to put there.
    if resolved == PurePosixPath("..") or resolved.as_posix().startswith("../"):
        return Issue(source.as_posix(), line,
                     f"link target {target!r} escapes the repository root "
                     f"(resolved to {resolved.as_posix()!r})")
    if resolved not in entries:
        return Issue(source.as_posix(), line,
                     f"missing repo-relative link target {target!r} "
                     f"(resolved to {resolved.as_posix()!r})")
    if _LINE_FRAGMENT_RE.match(_fragment(target)):
        return Issue(source.as_posix(), line,
                     f"link target {target!r} pins a line number in "
                     f"{resolved.as_posix()!r}, which any edit above it moves; "
                     f"link the file and name the symbol in the link text")
    return _anchor_issue(source, line, target, resolved, anchors)


def _path_span_issue(source: PurePosixPath, line: int, span: str,
                     entries: set[PurePosixPath],
                     used: set[str] | None = None) -> Issue | None:
    """Reports a backticked repo path that no tracked file matches."""
    match = _PATH_SPAN_RE.match(span)
    if not match:
        return None
    candidate = match.group(1)
    if "/" not in candidate or PurePosixPath(candidate).suffix not in _SOURCE_SUFFIXES:
        return None
    if _untracked_allowance(candidate, used):
        return None
    # Scope to paths rooted at a real tracked directory, so references to build
    # output, external checkouts, and illustrative paths stay out of the gate.
    root = PurePosixPath(posixpath.normpath(candidate)).parts[0]
    if PurePosixPath(root) not in entries:
        return None
    for base in (PurePosixPath(""), source.parent):
        resolved = posixpath.normpath(posixpath.join(base.as_posix(), candidate))
        if PurePosixPath(resolved) in entries:
            return None
    return Issue(source.as_posix(), line, f"backticked path {candidate!r} does not exist")


def _tree_names(rest: str) -> list[str]:
    """Names one tree row draws: its leading token plus any ` / ` siblings."""
    tokens = rest.split()
    if not tokens or not _TREE_NAME_RE.match(tokens[0]):
        return []
    names = [tokens[0]]
    index = 1
    while (index + 1 < len(tokens) and tokens[index] == "/"
           and _TREE_NAME_RE.match(tokens[index + 1])):
        names.append(tokens[index + 1])
        index += 2
    return names


def _tree_listed(description: str) -> list[str]:
    """Child stems a directory row's prose enumerates, in the order drawn.

    A group holding anything but bare stems is prose, not a list, and names
    nothing; the row then elides its subtree as an undescribed one does.
    """
    stems = []
    for group in _TREE_LIST_RE.findall(description):
        parts = [part.strip() for part in group.split(",")]
        if parts and all(_TREE_STEM_RE.fullmatch(part) for part in parts):
            stems.extend(parts)
    return stems


def _tree_stem_paths(directory: str, stem: str,
                     entries: set[PurePosixPath]) -> list[str]:
    """Resolves a prose-named stem to the children it names, suffix and all.

    A stem no child carries resolves to itself, which the existence check then
    reports as the missing path the prose claims.
    """
    matches = sorted(entry.as_posix() for entry in entries
                     if entry.parent.as_posix() == directory
                     and entry.stem == stem)
    return matches or [posixpath.join(directory, stem)]


def _tree_entry_exists(candidate: str, entries: set[PurePosixPath],
                       allowed: Callable[[str], bool]) -> bool:
    if allowed(candidate):
        return True
    if "*" in candidate or "?" in candidate:
        return any(fnmatch.fnmatchcase(entry.as_posix(), candidate)
                   for entry in entries)
    return PurePosixPath(candidate) in entries


@dataclass(frozen=True)
class TreeDirective:
    checkout: str
    exhaustive: bool


def _tree_directive(tag: str) -> TreeDirective | None:
    """Parses `tree [<checkout>] [exhaustive]`; None when the tag is not one."""
    tokens = tag.split()
    if not tokens or tokens[0] != _TREE_TAG:
        return None
    rest = tokens[1:]
    exhaustive = bool(rest) and rest[-1] == _TREE_EXHAUSTIVE
    if exhaustive:
        rest = rest[:-1]
    if len(rest) > 1:
        return None
    return TreeDirective(rest[0] if rest else "", exhaustive)


def _checkout_allowance(candidate: str, prefixes: tuple[str, ...]) -> bool:
    return any(candidate == prefix.rstrip("/") or candidate.startswith(prefix)
               for prefix in prefixes)


def _tree_rows(source: PurePosixPath, fence: Fence,
               entries: set[PurePosixPath]) -> tuple[list[tuple[int, str]],
                                                     list[Issue]]:
    """Resolves a fence's rows to (line, repo-relative path) pairs."""
    rows = []
    issues = []
    stack: list[str] = []
    # The directory row whose description is still open, and its text so far.
    described: tuple[int, str] | None = None
    words: list[str] = []

    def close() -> None:
        nonlocal described
        if described is not None:
            line_number, directory = described
            rows.extend(
                (line_number, path)
                for stem in _tree_listed(" ".join(words))
                for path in _tree_stem_paths(directory, stem, entries))
        described = None
        words.clear()

    for line_number, line in fence.body:
        match = _TREE_ROW_RE.match(line)
        if not match:
            continuation = described and _TREE_CONT_RE.match(line)
            if continuation:
                words.append(continuation.group("rest"))
            else:
                close()
            continue
        close()
        depth = len(match.group("indent")) // _TREE_INDENT
        names = _tree_names(match.group("rest"))
        if not names:
            continue
        if depth > len(stack):
            issues.append(Issue(source.as_posix(), line_number,
                                f"tree row {names[0]!r} is indented past its parent"))
            continue
        parent = "/".join(stack[:depth])
        stack[depth:] = [names[0].rstrip("/")]
        rows.extend((line_number, posixpath.join(parent, name.rstrip("/")))
                    for name in names)
        if names[0].endswith("/"):
            described = (line_number,
                         posixpath.join(parent, names[0].rstrip("/")))
            words = [match.group("rest")]
    close()
    return rows, issues


def _tree_omissions(source: PurePosixPath, fence: Fence,
                    rows: list[tuple[int, str]],
                    entries: set[PurePosixPath]) -> list[Issue]:
    """Reports tracked paths under a drawn directory that no row names.

    A directory whose rows name none of its children is a summary row and its
    subtree is elided; one that names any child must name them all.
    """
    drawn = {""}
    for _, path in rows:
        drawn.add(path)
        drawn.update(parent.as_posix() for parent in PurePosixPath(path).parents
                     if parent != PurePosixPath("."))
    globs = tuple(path for path in drawn if "*" in path or "?" in path)

    def is_drawn(candidate: str) -> bool:
        return candidate in drawn or any(
            fnmatch.fnmatchcase(candidate, pattern) for pattern in globs)

    children: dict[str, set[str]] = {}
    for entry in entries:
        parent = entry.parent.as_posix()
        children.setdefault("" if parent == "." else parent, set()).add(
            entry.as_posix())

    omitted: set[str] = set()
    for directory in drawn:
        siblings = children.get(directory, ())
        if not any(is_drawn(child) for child in siblings):
            continue
        omitted.update(child for child in siblings
                       if not is_drawn(child) and not _tree_unmapped(child))
    return [Issue(source.as_posix(), fence.start,
                  f"tree omits tracked path {path!r}") for path in sorted(omitted)]


def _tree_unmapped(candidate: str) -> bool:
    return any(candidate == prefix or candidate.startswith(prefix)
               for prefix in _TREE_UNMAPPED)


def _tree_issues(source: PurePosixPath, fences: list[Fence],
                 entries: set[PurePosixPath],
                 used: set[str] | None = None,
                 checkouts: dict[str, set[PurePosixPath]] | None = None,
                 skipped: set[str] | None = None) -> list[Issue]:
    """Reports drawn tree rows that name a path the drawn repository lacks."""
    issues = []
    for fence in fences:
        if not fence.tag.split():
            continue
        directive = _tree_directive(fence.tag)
        if directive is None:
            issues.append(Issue(source.as_posix(), fence.start,
                                f"unknown docs-check directive {fence.tag!r}"))
            continue
        if directive.checkout:
            target = (checkouts or {}).get(directive.checkout)
            if target is None:
                if skipped is not None:
                    skipped.add(directive.checkout)
                continue
            prefixes = _CHECKOUT_UNTRACKED_ALLOWED.get(directive.checkout, ())
            allowed = functools.partial(_checkout_allowance, prefixes=prefixes)
        else:
            target = entries
            allowed = functools.partial(_untracked_allowance, used=used)
        rows, row_issues = _tree_rows(source, fence, target)
        issues.extend(row_issues)
        issues.extend(
            Issue(source.as_posix(), line_number,
                  f"tree path {candidate!r} does not exist")
            for line_number, candidate in rows
            if not _tree_entry_exists(candidate, target, allowed))
        if directive.exhaustive:
            issues.extend(_tree_omissions(source, fence, rows, target))
    return issues


def effect_roster(source: str) -> set[str]:
    """Names HS_EFFECT_LIST expands over, given the roster header's text.

    The list's one macro-expanded row, HS_SHADER_WORKBENCH_EFFECT(X), yields
    X(Shader) only under HS_ENABLE_SHADER_WORKBENCH, which neither the firmware
    playlist nor the gallery roster sets, so it is not among the names.
    """
    names: set[str] = set()
    inside = False
    for line in source.splitlines():
        if not inside:
            inside = line.startswith(_EFFECT_ROSTER_DEFINE)
            continue
        match = _EFFECT_ROSTER_ENTRY_RE.match(line)
        if match:
            names.add(match.group(1))
        if not line.rstrip().endswith("\\"):
            break
    return names


def effects_row_issues(text: str, entries: set[PurePosixPath],
                       roster: set[str] | None) -> list[Issue]:
    """Checks the summary row's header and effect counts against the tree.

    A count nobody derives drifts silently: the row elides its own subtree, so
    the exhaustive-tree gate never reaches it.
    """
    headers = sum(1 for entry in entries
                  if entry.parent == _EFFECTS_DIR and entry.suffix == ".h")
    issues = []
    matched = False
    diagram_matched = False
    for number, line in enumerate(text.splitlines(), 1):
        diagram = _EFFECTS_DIAGRAM_RE.search(line)
        if diagram:
            diagram_matched = True
            drawn = int(diagram.group("effects"))
            if roster is None:
                issues.append(Issue(
                    _EFFECTS_TREE_ROW, number,
                    f"architecture diagram claims {drawn} effects, but "
                    f"{_EFFECT_ROSTER_SOURCE} defines no HS_EFFECT_LIST to "
                    f"check it against"))
            elif drawn != len(roster):
                issues.append(Issue(
                    _EFFECTS_TREE_ROW, number,
                    f"architecture diagram claims {drawn} effects, "
                    f"HS_EFFECT_LIST names {len(roster)}"))
        match = _EFFECTS_ROW_RE.search(line)
        if not match:
            continue
        matched = True
        drawn_headers = int(match.group("headers"))
        if drawn_headers != headers:
            issues.append(Issue(
                _EFFECTS_TREE_ROW, number,
                f"effects/ row claims {drawn_headers} headers, "
                f"the tracked tree has {headers}"))
        drawn_effects = int(match.group("effects") or match.group("legacy_effects"))
        if roster is None:
            issues.append(Issue(
                _EFFECTS_TREE_ROW, number,
                f"effects/ row claims {drawn_effects} effects, but "
                f"{_EFFECT_ROSTER_SOURCE} defines no HS_EFFECT_LIST to "
                f"check it against"))
        elif drawn_effects != len(roster):
            issues.append(Issue(
                _EFFECTS_TREE_ROW, number,
                f"effects/ row claims {drawn_effects} effects, "
                f"HS_EFFECT_LIST names {len(roster)}"))
    if not matched:
        issues.append(Issue(
            _EFFECTS_TREE_ROW, 1,
            "no effects/ summary row, so its counts go unchecked"))
    if not diagram_matched:
        issues.append(Issue(
            _EFFECTS_TREE_ROW, 1,
            "no effects/ architecture-diagram row, so its effect count goes "
            "unchecked"))
    return issues


def check_text(source: PurePosixPath, text: str,
               entries: set[PurePosixPath],
               anchors: dict[PurePosixPath, set[str]] | None = None,
               used: set[str] | None = None,
               checkouts: dict[str, set[PurePosixPath]] | None = None,
               skipped: set[str] | None = None) -> list[Issue]:
    visible, fences, issues = _visible_lines(source, text)
    issues.extend(_tree_issues(source, fences, entries, used, checkouts, skipped))
    # The source's own anchors always resolve; cross-document ones need the
    # repository-wide map check_repository builds.
    anchors = dict(anchors or {})
    anchors[source] = _anchors(visible)

    definitions: dict[str, tuple[str, int]] = {}
    body_lines = []
    for line in visible:
        for span in line.spans:
            issue = _path_span_issue(source, line.number, span, entries, used)
            if issue:
                issues.append(issue)
        match = _REFERENCE_DEFINITION_RE.match(line.masked)
        if not match:
            body_lines.append((line.number, line.masked))
            continue
        target = _destination(match.group(2))
        if target is not None:
            definitions.setdefault(
                _normalize_label(match.group(1)), (target, line.number))

    for target, line_number in definitions.values():
        issue = _link_issue(source, line_number, target, entries, anchors)
        if issue:
            issues.append(issue)

    for line_number, line in body_lines:
        destinations, references = _inline_links(line)
        for target in destinations:
            issue = _link_issue(source, line_number, target, entries, anchors)
            if issue:
                issues.append(issue)
        for reference in references:
            if reference not in definitions:
                issues.append(Issue(source.as_posix(), line_number,
                                    f"undefined reference link [{reference}]"))
    return sorted(issues)


def _tracked_entries(root: Path) -> tuple[list[PurePosixPath], set[PurePosixPath]]:
    command = ["git", "-c", f"safe.directory={root.as_posix()}",
               "-C", str(root), "ls-files", "-z"]
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE)
    files = [PurePosixPath(name) for name in
             result.stdout.decode("utf-8").split("\0") if name]
    entries = set(files)
    for file_path in files:
        entries.update(parent for parent in file_path.parents
                       if parent != PurePosixPath("."))
    markdown = sorted(path for path in files
                      if path.suffix.casefold() in (".md", ".markdown")
                      and path not in _MARKDOWN_EXCLUDED)
    return markdown, entries


def check_repository(
        root: Path, checkout_roots: dict[str, Path] | None = None,
        skipped: set[str] | None = None
) -> tuple[list[PurePosixPath], list[Issue], list[str]]:
    markdown, entries = _tracked_entries(root)
    checkouts = {name: _tracked_entries(path)[1]
                 for name, path in (checkout_roots or {}).items()}
    issues = []
    used: set[str] = set()
    sources: dict[PurePosixPath, str] = {}
    for relative in markdown:
        path = root.joinpath(*relative.parts)
        try:
            sources[relative] = path.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as error:
            issues.append(Issue(relative.as_posix(), 1,
                                f"cannot read tracked Markdown as UTF-8: {error}"))

    # Anchors first: a link may point at a heading in any other tracked document.
    anchors = {relative: _anchors(_visible_lines(relative, text)[0])
               for relative, text in sources.items()}
    for relative, text in sources.items():
        issues.extend(check_text(relative, text, entries, anchors, used,
                                 checkouts, skipped))
    # Only this repository draws the row; a checkout without the roster header
    # is not the tree the claim is about.
    if _EFFECT_ROSTER_SOURCE in entries:
        try:
            roster = effect_roster(root.joinpath(
                *_EFFECT_ROSTER_SOURCE.parts).read_text(encoding="utf-8"))
        except (OSError, UnicodeError):
            roster = set()
        issues.extend(effects_row_issues(
            sources.get(PurePosixPath(_EFFECTS_TREE_ROW), ""), entries,
            roster or None))
    return markdown, sorted(issues), _stale_allowances(entries, used)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check tracked Markdown fences and repository links.")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--checkout", action="append", default=[], metavar="NAME=PATH",
        help="root of a sibling checkout a `tree <NAME>` fence draws")
    parser.add_argument(
        "--skip-checkout", action="append", default=[], metavar="NAME",
        help="accept `tree <NAME>` fences unvalidated; without it a fence "
             "naming a checkout given no --checkout root fails")
    args = parser.parse_args(argv)

    checkout_roots = {}
    for option in args.checkout:
        name, separator, path = option.partition("=")
        if not separator or not name or not path:
            parser.error(f"--checkout expects NAME=PATH, got {option!r}")
        checkout_root = Path(path).resolve()
        if not checkout_root.is_dir():
            parser.error(f"--checkout root for {name!r} is not a directory: "
                         f"{checkout_root}")
        checkout_roots[name] = checkout_root

    skipped: set[str] = set()
    try:
        markdown, issues, stale = check_repository(
            args.root.resolve(), checkout_roots, skipped)
    except (OSError, subprocess.SubprocessError, UnicodeError) as error:
        print(f"[docs-check] tooling error: {error}", file=sys.stderr)
        return 2
    # An unvalidated tree fence is not a pass: its rows went ungated, so every
    # one of them can be edited through a green run. Accepting that is allowed --
    # a developer without the sibling checkout must still be able to run the
    # checker -- but only by naming the checkout in --skip-checkout, and the
    # verdict line says so either way.
    unvalidated = skipped - set(args.skip_checkout)
    if skipped:
        print(f"::{'error' if unvalidated else 'warning'}::tree fences NOT "
              f"validated - no --checkout root for: {', '.join(sorted(skipped))}")
    # Warning only: dropping the last citation of an exempt path improves the
    # tree and must not red the build for whoever lands that commit.
    if markdown and stale:
        print("::warning::_UNTRACKED_ALLOWED in tools/docs_check.py is stale - "
              f"drop these entries: {', '.join(stale)}")
    if issues:
        for issue in issues:
            print(issue)
        print(f"[docs-check] FAIL - {len(issues)} issue(s)", file=sys.stderr)
        return 1
    # No tracked Markdown means the checker was pointed somewhere it cannot see
    # the repository (wrong --root, no git); passing would certify nothing.
    if not markdown:
        print(f"[docs-check] tooling error: no tracked Markdown under "
              f"{args.root.resolve()}", file=sys.stderr)
        return 2
    if unvalidated:
        print(f"[docs-check] FAIL - tree fences unvalidated: "
              f"{', '.join(sorted(unvalidated))} - pass --checkout NAME=PATH, or "
              f"--skip-checkout NAME to accept them unvalidated", file=sys.stderr)
        return 1
    note = (f"; tree fences NOT validated: {', '.join(sorted(skipped))}"
            if skipped else "")
    print(f"[docs-check] PASS - {len(markdown)} tracked Markdown file(s){note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
