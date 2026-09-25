#!/usr/bin/env python3
"""Validate fences, links, anchors, and path claims in tracked Markdown.

Structure only. A green run means every fence closes, every anchor resolves,
every link into this repository or a supplied sibling checkout resolves, every
recognized backticked path under a tracked repository root exists (bare
basenames and unknown first segments are not checked), every tree fence matches the tracked tree it
draws, the cardinalities CARDINALITY_CLAIMS names match their source
macros, and the composed-effect roster in docs/effects.md matches each
effect's PRESET_IDS and the product group -- not that the prose is true. A
link to any other host is never
visited, and a sibling checkout no --checkout root supplies leaves its fences
and links unvalidated, which the verdict line says. A renamed symbol in a
table, a number no claim names, and any path written without backticks or a
link are all invisible here.

The Doxyfile's PREDEFINED names must also appear in documented source.
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
TREE_ROW_RE = re.compile(r"^(?P<indent>(?:│   |    )*)(?:├──|└──) +(?P<rest>\S.*)$")
_TREE_INDENT = 4
# A directory row may enumerate its children in prose rather than draw one row
# apiece, wrapping onto continuation lines that carry the spine but no branch.
# Within such a row every parenthesized group whose comma-separated parts are
# all bare stems names children of it, gated in both directions like a row.
_TREE_CONT_RE = re.compile(r"^[│ \t]*(?P<rest>\S.*)$")
TREE_LIST_RE = re.compile(r"\(([^()]*)\)")
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
# A sibling repository's own GitHub URLs name paths that are real there. The
# `tree <NAME>` fences already resolve such a repository against a --checkout
# root, and the same root resolves these links.
_CHECKOUT_REPO_PATH_RE = re.compile(
    r"^/woundedlion/(?!pov/)(?P<checkout>[^/]+)/(?:blob|tree|raw)/[^/]+/"
    r"(?P<path>.+)$")

# The `effects/` row summarizes its subtree instead of drawing it, so no tree
# gate sees the counts it states. Both are derivable: headers from the tracked
# tree, effects from the roster macro's cardinality.
_EFFECTS_TREE_ROW = "README.md"
EFFECTS_ROW_RE = re.compile(
    r"\beffects/\s+(?P<headers>\d+) headers(?: covering (?P<effects>\d+) effects|: "
    r"one per effect \((?P<legacy_effects>\d+)\))")
# The architecture diagram restates the roster's cardinality in its own words,
# outside any tree fence and in a spelling the summary row's regex never
# matches, so it is derivable from the same source and gated the same way.
EFFECTS_DIAGRAM_RE = re.compile(
    r"\beffects/\s+\((?P<effects>\d+) visual algorithms\)")
EFFECTS_DIR = PurePosixPath("effects")
EFFECT_ROSTER_SOURCE = PurePosixPath("targets/effects.h")
_EFFECT_ROSTER_DEFINE = "#define HS_EFFECT_LIST(X)"
# Shared X-row spelling with scripts/effect_roster.mjs and tools/profile_sweep.sh:
# whitespace inside the parens is tolerated so a reformat to `X( Foo )` cannot
# silently shrink one roster reader's count while the others keep the full set.
# Comments are stripped before matching, as effect_roster.mjs does, so a
# commented-out row is not counted -- including one inside a block comment that
# spans lines, which no line-anchored pattern can exclude.
_COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)
_EFFECT_ROSTER_ENTRY_RE = re.compile(r"X\(\s*(\w+)\s*\)")

# The device playlist repeats the full roster's cardinality in README prose.
PHANTASM_PLAYLIST_SOURCE = PurePosixPath("targets/Phantasm/phantasm_playlist.h")
_PHANTASM_ROSTER_DEFINE = "#define HS_PHANTASM_EFFECT_LIST(X)"
# The shader promotion product group lives beside HS_EFFECT_LIST, and the ITCM
# ledger restates its cardinality in prose.
_PRODUCT_GROUP_DEFINE = "#define HS_SHADER_PRODUCT_GROUP(X)"
# Both macros spell a row `X(Name, seconds)`.
_NAMED_ROSTER_ENTRY_RE = re.compile(r"X\(\s*(\w+)\s*,")
_ITCM_LEDGER = "docs/ledgers/itcm_ledger.md"
# The effects reference spells the product group's cardinality in words and
# tables every composed effect with its preset count.
_EFFECTS_REFERENCE = "docs/effects.md"
_COMPOSED_ROSTER_ROW_RE = re.compile(
    r"^\|\s*`[^`|]+`\s*\|\s*`(?P<effect>\w+)`\s*\|\s*(?P<presets>\d+)\s*\|")
_PRESET_IDS_RE = re.compile(
    r"std::array<std::string_view,\s*(\d+)>\s*PRESET_IDS\b")
_NUMBER_WORDS = {
    word: value for value, word in enumerate((
        "zero", "one", "two", "three", "four", "five", "six", "seven",
        "eight", "nine", "ten", "eleven", "twelve", "thirteen", "fourteen",
        "fifteen", "sixteen", "seventeen", "eighteen", "nineteen", "twenty"))
} | {"thirty": 30, "forty": 40, "fifty": 50, "sixty": 60, "seventy": 70,
     "eighty": 80, "ninety": 90}
CARDINALITY_CLAIMS = (
    ("README.md",
     re.compile(r"compile-time roster and tests carry (\d+) firmware-capable"),
     "HS_EFFECT_LIST", "the firmware-capable roster size"),
    ("README.md", re.compile(r"\bcontains (\d+) effects\b"),
     "HS_PHANTASM_EFFECT_LIST", "the effects-reference playlist size"),
    ("README.md", re.compile(r"\b(\d+)-entry roster\b"),
     "HS_PHANTASM_EFFECT_LIST", "the frame-sync playlist size"),
    ("docs/specs/phantasm_frame_sync_spec.md",
     re.compile(r"\b(\d+)-entry roster\b"),
     "HS_PHANTASM_EFFECT_LIST", "the frame-sync playlist size"),
    (_ITCM_LEDGER, re.compile(r"\b(\d+) fixed-pipeline products\b"),
     "HS_SHADER_PRODUCT_GROUP", "the shader product-group size"),
    (_ITCM_LEDGER, re.compile(r"\b(\d+)-effect phantasm roster\b"),
     "HS_PHANTASM_EFFECT_LIST", "the ledger's phantasm roster size"),
    (_EFFECTS_REFERENCE,
     re.compile(r"\bThese ([\w-]+) effects form the product-only\b"),
     "HS_SHADER_PRODUCT_GROUP", "the composed roster size"),
    (_EFFECTS_REFERENCE,
     re.compile(r"\bthe ([\w-]+) promoted fixed descriptors\b"),
     "HS_SHADER_PRODUCT_GROUP", "the promoted-descriptor count"),
)

# Names are matched against the source file kinds Doxyfile documents.
_DOXYFILE = PurePosixPath("Doxyfile")
_GIT_TIMEOUT_SECONDS = 30
_DOXYGEN_SOURCE_SUFFIXES = frozenset({".h", ".cpp", ".ino"})
_DOXYFILE_TAG_RE = re.compile(r"^(\w+)[ \t]*\+?=[ \t]*(.*)$")
_DOXYFILE_ENTRY_RE = re.compile(r'"[^"]*"|\S+')
_PREDEFINED_NAME_RE = re.compile(r"^[A-Za-z_]\w*")
# DOXYGEN enables doc-only branches and need not appear in source.
_PREDEFINED_UNREFERENCED_ALLOWED = frozenset({"DOXYGEN"})

UNTRACKED_ALLOWED = (
    ".github/workflows/deploy.yml",
)
_UNTRACKED_LIST = "untracked-allowed"
_TREE_UNMAPPED_LIST = "tree-unmapped"

# Tracked paths an exhaustive tree deliberately leaves without a row: VCS
# metadata, the map's own document, and the test tree the map draws as one
# summary row plus its shared fixtures. A trailing slash covers a subtree.
TREE_UNMAPPED = (
    ".gitattributes",
    ".gitignore",
    "README.md",
    "hardware/phantasm/.gitignore",
    "tests/",
)

# Rows a checkout's own repository gitignores, keyed by the directive's name.
CHECKOUT_UNTRACKED_ALLOWED = {
    "daydream": ("node_modules/", "three.js/", "vendor/"),
}


# Prose also spells engine headers relative to core/ (`render/shading.h`), so
# every backticked candidate is resolved under this root as well.
_IMPLICIT_PATH_ROOT = PurePosixPath("core")


def _cited(allowlist: str, entry: str) -> str:
    """One citation token, namespaced by allowlist so the three cannot collide."""
    return f"{allowlist}:{entry}"


def _untracked_allowance(candidate: str, used: set[str] | None) -> bool:
    """Reports whether an UNTRACKED_ALLOWED prefix exempts candidate, recording it."""
    for prefix in UNTRACKED_ALLOWED:
        if candidate.startswith(prefix):
            if used is not None:
                used.add(_cited(_UNTRACKED_LIST, prefix))
            return True
    return False


def _stale_allowances(entries: set[PurePosixPath], used: set[str],
                      checkouts: dict[str, set[PurePosixPath]] | None = None
) -> list[str]:
    """Names allowlist entries the exemption no longer buys anything for.

    An UNTRACKED_ALLOWED prefix exempts a path the docs cite and this
    repository does not track, so tracking it makes the entry stale; a
    TREE_UNMAPPED prefix exempts a tracked path from needing a tree row, so
    untracking it does. A checkout allowance is judged only against a checkout
    a --checkout root supplied, since nothing else can say what it tracks.
    """
    stale = []
    for prefix in UNTRACKED_ALLOWED:
        if PurePosixPath(prefix.rstrip("/")) in entries:
            stale.append(f"{prefix} (now tracked)")
        elif _cited(_UNTRACKED_LIST, prefix) not in used:
            stale.append(f"{prefix} (uncited)")
    for prefix in TREE_UNMAPPED:
        entry = _cited(_TREE_UNMAPPED_LIST, prefix)
        if PurePosixPath(prefix.rstrip("/")) not in entries:
            stale.append(f"{entry} (untracked)")
        elif entry not in used:
            stale.append(f"{entry} (uncited)")
    for checkout, prefixes in CHECKOUT_UNTRACKED_ALLOWED.items():
        tracked = (checkouts or {}).get(checkout)
        if tracked is None:
            continue
        for prefix in prefixes:
            entry = _cited(checkout, prefix)
            if PurePosixPath(prefix.rstrip("/")) in tracked:
                stale.append(f"{entry} (now tracked)")
            elif entry not in used:
                stale.append(f"{entry} (uncited)")
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


def visible_lines(path: PurePosixPath,
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


def _checkout_link(target: str) -> tuple[str, PurePosixPath] | None:
    """Splits a sibling-repository GitHub URL into its checkout name and path."""
    parsed = urlsplit(_MARKDOWN_ESCAPE_RE.sub(r"\1", target.strip()))
    if parsed.scheme.casefold() not in ("http", "https"):
        return None
    if parsed.netloc.casefold() not in _SELF_REPO_HOSTS:
        return None
    match = _CHECKOUT_REPO_PATH_RE.match(unquote(parsed.path))
    if not match:
        return None
    return (match.group("checkout"),
            PurePosixPath(posixpath.normpath(match.group("path"))))


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
                anchors: dict[PurePosixPath, set[str]],
                checkouts: dict[str, set[PurePosixPath]] | None = None,
                skipped: set[str] | None = None) -> Issue | None:
    sibling = _checkout_link(target)
    if sibling is not None:
        name, path = sibling
        tracked = (checkouts or {}).get(name)
        if tracked is None:
            if skipped is not None:
                skipped.add(name)
            return None
        if path in tracked:
            return None
        return Issue(source.as_posix(), line,
                     f"missing {name} link target {target!r} "
                     f"(resolved to {path.as_posix()!r})")
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
    if PurePosixPath(candidate).suffix not in _SOURCE_SUFFIXES:
        return None
    if _untracked_allowance(candidate, used):
        return None
    if "/" not in candidate:
        if source.parent / candidate in entries:
            return None
        matches = sorted(path.as_posix() for path in entries
                         if path.name == candidate)
        if len(matches) > 1:
            return Issue(source.as_posix(), line,
                         f"ambiguous basename {candidate!r}: {', '.join(matches)}")
        return None
    in_scope = False
    for base in (PurePosixPath(""), source.parent, _IMPLICIT_PATH_ROOT):
        resolved = PurePosixPath(posixpath.normpath(
            posixpath.join(base.as_posix(), candidate)))
        # Scope each interpretation to its first directory under that base.
        scope = PurePosixPath(*resolved.parts[:len(base.parts) + 1])
        if scope not in entries:
            continue
        in_scope = True
        if resolved in entries:
            return None
    if not in_scope:
        return None
    return Issue(source.as_posix(), line, f"backticked path {candidate!r} does not exist")


def tree_names(rest: str) -> list[str]:
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


def tree_listed(description: str) -> list[str]:
    """Child stems a directory row's prose enumerates, in the order drawn.

    A group holding anything but bare stems is prose, not a list, and names
    nothing; the row then elides its subtree as an undescribed one does.
    """
    stems = []
    for group in TREE_LIST_RE.findall(description):
        parts = [part.strip() for part in group.split(",")]
        if parts and all(_TREE_STEM_RE.fullmatch(part) for part in parts):
            stems.extend(parts)
    return stems


def tree_stem_paths(directory: str, stem: str,
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


def tree_directive(tag: str) -> TreeDirective | None:
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


_REQUIRED_TREES = frozenset({
    (PurePosixPath("README.md"), TreeDirective("", True)),
    (PurePosixPath("README.md"), TreeDirective("daydream", True)),
})


def _checkout_allowance(candidate: str, checkout: str,
                        prefixes: tuple[str, ...],
                        used: set[str] | None = None) -> bool:
    for prefix in prefixes:
        if candidate == prefix.rstrip("/") or candidate.startswith(prefix):
            if used is not None:
                used.add(_cited(checkout, prefix))
            return True
    return False


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
                for stem in tree_listed(" ".join(words))
                for path in tree_stem_paths(directory, stem, entries))
        described = None
        words.clear()

    for line_number, line in fence.body:
        match = TREE_ROW_RE.match(line)
        if not match:
            continuation = described and _TREE_CONT_RE.match(line)
            if continuation:
                words.append(continuation.group("rest"))
            else:
                close()
            continue
        close()
        depth = len(match.group("indent")) // _TREE_INDENT
        names = tree_names(match.group("rest"))
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
                    entries: set[PurePosixPath],
                    unmapped: tuple[str, ...],
                    used: set[str] | None = None) -> list[Issue]:
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
                       if not is_drawn(child)
                       and not tree_unmapped(child, unmapped, used))
    return [Issue(source.as_posix(), fence.start,
                  f"tree omits tracked path {path!r}") for path in sorted(omitted)]


def tree_unmapped(candidate: str, prefixes: tuple[str, ...],
                   used: set[str] | None = None) -> bool:
    for prefix in prefixes:
        if candidate == prefix or candidate.startswith(prefix):
            if used is not None:
                used.add(_cited(_TREE_UNMAPPED_LIST, prefix))
            return True
    return False


def _tree_issues(source: PurePosixPath, fences: list[Fence],
                 entries: set[PurePosixPath],
                 used: set[str] | None = None,
                 checkouts: dict[str, set[PurePosixPath]] | None = None,
                 skipped: set[str] | None = None,
                 seen: set[tuple[PurePosixPath, TreeDirective]] | None = None
) -> list[Issue]:
    """Reports drawn tree rows that name a path the drawn repository lacks."""
    issues = []
    for fence in fences:
        if not fence.tag.split():
            continue
        directive = tree_directive(fence.tag)
        if directive is None:
            issues.append(Issue(source.as_posix(), fence.start,
                                f"unknown docs-check directive {fence.tag!r}"))
            continue
        if seen is not None:
            seen.add((source, directive))
        if directive.checkout:
            target = (checkouts or {}).get(directive.checkout)
            if target is None:
                if skipped is not None:
                    skipped.add(directive.checkout)
                continue
            prefixes = CHECKOUT_UNTRACKED_ALLOWED.get(directive.checkout, ())
            allowed = functools.partial(_checkout_allowance,
                                        checkout=directive.checkout,
                                        prefixes=prefixes, used=used)
            unmapped = ()
        else:
            target = entries
            allowed = functools.partial(_untracked_allowance, used=used)
            unmapped = TREE_UNMAPPED
        rows, row_issues = _tree_rows(source, fence, target)
        issues.extend(row_issues)
        issues.extend(
            Issue(source.as_posix(), line_number,
                  f"tree path {candidate!r} does not exist")
            for line_number, candidate in rows
            if not _tree_entry_exists(candidate, target, allowed))
        if directive.exhaustive:
            issues.extend(_tree_omissions(source, fence, rows, target,
                                          unmapped, used))
    return issues


def _macro_body(source: str, define: str) -> str:
    """The continued-line body of a macro definition, comments stripped.

    The body runs from the #define to the first line without a trailing
    backslash; the break is decided on the raw lines, so stripping cannot
    extend or truncate it.
    """
    body: list[str] = []
    inside = False
    for line in source.splitlines():
        if not inside:
            inside = line.startswith(define)
            if not inside:
                continue
            line = line[len(define):]
        body.append(line)
        if not line.rstrip().endswith("\\"):
            break
    return _COMMENT_RE.sub("", "\n".join(body))


def effect_roster(source: str) -> set[str]:
    """Names HS_EFFECT_LIST expands over, given the roster header's text.

    The list's one macro-expanded row, HS_SHADER_WORKBENCH_EFFECT(X), yields
    X(Shader) only under HS_ENABLE_SHADER_WORKBENCH, which neither the firmware
    playlist nor the gallery roster sets, so it is not among the names.
    """
    return set(_EFFECT_ROSTER_ENTRY_RE.findall(
        _macro_body(source, _EFFECT_ROSTER_DEFINE)))


def effects_row_issues(text: str, entries: set[PurePosixPath],
                       roster: set[str] | None) -> list[Issue]:
    """Checks the summary row's header and effect counts against the tree.

    A count nobody derives drifts silently: the row elides its own subtree, so
    the exhaustive-tree gate never reaches it.
    """
    headers = sum(1 for entry in entries
                  if entry.parent == EFFECTS_DIR and entry.suffix == ".h")
    issues = []
    matched = False
    diagram_matched = False
    for number, line in enumerate(text.splitlines(), 1):
        diagram = EFFECTS_DIAGRAM_RE.search(line)
        if diagram:
            diagram_matched = True
            drawn = int(diagram.group("effects"))
            if roster is None:
                issues.append(Issue(
                    _EFFECTS_TREE_ROW, number,
                    f"architecture diagram claims {drawn} effects, but "
                    f"{EFFECT_ROSTER_SOURCE} defines no HS_EFFECT_LIST to "
                    f"check it against"))
            elif drawn != len(roster):
                issues.append(Issue(
                    _EFFECTS_TREE_ROW, number,
                    f"architecture diagram claims {drawn} effects, "
                    f"HS_EFFECT_LIST names {len(roster)}"))
        match = EFFECTS_ROW_RE.search(line)
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
                f"{EFFECT_ROSTER_SOURCE} defines no HS_EFFECT_LIST to "
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


def doxyfile_predefined(text: str) -> list[tuple[int, str]]:
    """Macro names the Doxyfile's PREDEFINED tag defines, with their line numbers."""
    names: list[tuple[int, str]] = []
    inside = False
    for number, line in enumerate(text.splitlines(), 1):
        if inside:
            body = line
        else:
            match = _DOXYFILE_TAG_RE.match(line)
            if not match or match.group(1) != "PREDEFINED":
                continue
            body = match.group(2)
        stripped = body.rstrip()
        inside = stripped.endswith("\\")
        for entry in _DOXYFILE_ENTRY_RE.findall(stripped.rstrip("\\")):
            name = _PREDEFINED_NAME_RE.match(entry.strip('"'))
            if name:
                names.append((number, name.group()))
    return names


def unreferenced_predefined(root: Path, sources: list[PurePosixPath],
                            names: set[str]) -> set[str]:
    """Subset of names no documented source mentions."""
    pending = set(names) - _PREDEFINED_UNREFERENCED_ALLOWED
    if not pending:
        return pending
    pattern = re.compile(
        r"\b(" + "|".join(re.escape(name) for name in sorted(pending)) + r")\b")
    for relative in sources:
        try:
            text = root.joinpath(*relative.parts).read_text(
                encoding="utf-8", errors="replace")
        except OSError:
            continue
        pending.difference_update(
            match.group(1) for match in pattern.finditer(text))
        if not pending:
            break
    return pending


def doxyfile_predefined_issues(predefined: list[tuple[int, str]],
                               unreferenced: set[str]) -> list[Issue]:
    """Reports PREDEFINED names that guard nothing, and a tag that parsed empty."""
    path = _DOXYFILE.as_posix()
    if not predefined:
        return [Issue(path, 1, "no PREDEFINED tag, so its macro names go "
                               "unchecked")]
    return [Issue(path, number,
                  f"PREDEFINED defines {name}, which no documented source "
                  f"names")
            for number, name in predefined if name in unreferenced]


def phantasm_roster(source: str) -> set[str]:
    """Returns the effect names in HS_PHANTASM_EFFECT_LIST."""
    return set(_NAMED_ROSTER_ENTRY_RE.findall(
        _macro_body(source, _PHANTASM_ROSTER_DEFINE)))


def shader_product_group(source: str) -> set[str]:
    """Returns the effect names in HS_SHADER_PRODUCT_GROUP."""
    return set(_NAMED_ROSTER_ENTRY_RE.findall(
        _macro_body(source, _PRODUCT_GROUP_DEFINE)))


def claimed_count(token: str) -> int | None:
    """The count a prose claim spells, in digits or English number words."""
    digits = token.replace(",", "")
    if digits.isdigit():
        return int(digits)
    total = 0
    for part in token.casefold().split("-"):
        if part not in _NUMBER_WORDS:
            return None
        total += _NUMBER_WORDS[part]
    return total


def roster_claim_issues(sources: dict[PurePosixPath, str], roster: set[str],
                        playlist: set[str],
                        products: set[str]) -> list[Issue]:
    """Checks prose roster cardinalities against their source macros."""
    if not playlist:
        return [Issue(PHANTASM_PLAYLIST_SOURCE.as_posix(), 1,
                      "no HS_PHANTASM_EFFECT_LIST, so every playlist count "
                      "restated in prose goes unchecked")]
    if not products:
        return [Issue(EFFECT_ROSTER_SOURCE.as_posix(), 1,
                      "no HS_SHADER_PRODUCT_GROUP, so every product count "
                      "restated in prose goes unchecked")]
    counts = {"HS_EFFECT_LIST": len(roster),
              "HS_PHANTASM_EFFECT_LIST": len(playlist),
              "HS_SHADER_PRODUCT_GROUP": len(products)}
    issues = []
    for document, pattern, macro, subject in CARDINALITY_CLAIMS:
        expected = counts[macro]
        text = sources.get(PurePosixPath(document), "")
        matched = False
        for number, line in enumerate(text.splitlines(), 1):
            for match in pattern.finditer(line):
                matched = True
                claimed = claimed_count(match.group(1))
                if claimed is None:
                    issues.append(Issue(
                        document, number,
                        f"{subject} is stated as {match.group(1)!r}, "
                        "which is not a count"))
                elif claimed != expected:
                    issues.append(Issue(
                        document, number,
                        f"{subject} is stated as {claimed}, "
                        f"{macro} names {expected}"))
        if not matched:
            issues.append(Issue(document, 1,
                                f"no statement of {subject}, so it goes "
                                f"unchecked"))
    return issues


def composed_roster_issues(root: Path, effects_text: str,
                           entries: set[PurePosixPath],
                           products: set[str]) -> list[Issue]:
    """Checks the composed-effect roster table against each row's effect header
    (its PRESET_IDS size) and against HS_SHADER_PRODUCT_GROUP's membership."""
    issues = []
    rows: dict[str, tuple[int, int]] = {}
    for number, line in enumerate(effects_text.splitlines(), 1):
        match = _COMPOSED_ROSTER_ROW_RE.match(line)
        if not match:
            continue
        effect = match.group("effect")
        if effect in rows:
            issues.append(Issue(_EFFECTS_REFERENCE, number,
                                f"composed roster lists {effect} twice"))
            continue
        rows[effect] = (number, int(match.group("presets")))
    if not rows:
        return [Issue(_EFFECTS_REFERENCE, 1,
                      "no composed-effect roster table, so its preset counts "
                      "go unchecked")]
    for effect, (number, claimed) in rows.items():
        header = EFFECTS_DIR / f"{effect}.h"
        if header not in entries:
            issues.append(Issue(
                _EFFECTS_REFERENCE, number,
                f"composed roster names {header.as_posix()}, which is not "
                "tracked"))
            continue
        try:
            text = root.joinpath(*header.parts).read_text(encoding="utf-8")
        except (OSError, UnicodeError) as error:
            issues.append(Issue(header.as_posix(), 1,
                                f"cannot read as UTF-8: {error}"))
            continue
        declared = _PRESET_IDS_RE.search(_COMMENT_RE.sub("", text))
        if not declared:
            issues.append(Issue(
                _EFFECTS_REFERENCE, number,
                f"composed roster claims {claimed} presets for {effect}, "
                "which declares no PRESET_IDS"))
        elif int(declared.group(1)) != claimed:
            issues.append(Issue(
                _EFFECTS_REFERENCE, number,
                f"{effect} presets stated as {claimed}, PRESET_IDS holds "
                f"{declared.group(1)}"))
    for effect in sorted(products - rows.keys()):
        issues.append(Issue(
            _EFFECTS_REFERENCE, 1,
            f"composed roster omits {effect}, which HS_SHADER_PRODUCT_GROUP "
            "names"))
    for effect in sorted(rows.keys() - products):
        issues.append(Issue(
            _EFFECTS_REFERENCE, rows[effect][0],
            f"composed roster lists {effect}, which HS_SHADER_PRODUCT_GROUP "
            "does not name"))
    return issues


def check_text(source: PurePosixPath, text: str,
               entries: set[PurePosixPath],
               anchors: dict[PurePosixPath, set[str]] | None = None,
               used: set[str] | None = None,
               checkouts: dict[str, set[PurePosixPath]] | None = None,
               skipped: set[str] | None = None,
               seen_trees: set[tuple[PurePosixPath, TreeDirective]] | None = None
) -> list[Issue]:
    visible, fences, issues = visible_lines(source, text)
    issues.extend(_tree_issues(source, fences, entries, used, checkouts, skipped,
                               seen_trees))
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
        issue = _link_issue(source, line_number, target, entries, anchors,
                            checkouts, skipped)
        if issue:
            issues.append(issue)

    for line_number, line in body_lines:
        destinations, references = _inline_links(line)
        for target in destinations:
            issue = _link_issue(source, line_number, target, entries, anchors,
                                checkouts, skipped)
            if issue:
                issues.append(issue)
        for reference in references:
            if reference not in definitions:
                issues.append(Issue(source.as_posix(), line_number,
                                    f"undefined reference link [{reference}]"))
    return sorted(issues)


def tracked_entries(root: Path, revision: str | None = None
                     ) -> tuple[list[PurePosixPath], set[PurePosixPath]]:
    command = ["git", "-c", f"safe.directory={root.as_posix()}",
               "-C", str(root)]
    command.extend(["ls-tree", "-r", "--name-only", "-z", revision]
                   if revision else ["ls-files", "-z"])
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE,
                            timeout=_GIT_TIMEOUT_SECONDS)
    files = [PurePosixPath(name) for name in
             result.stdout.decode("utf-8").split("\0") if name]
    entries = set(files)
    for file_path in files:
        entries.update(parent for parent in file_path.parents
                       if parent != PurePosixPath("."))
    markdown = sorted(path for path in files
                      if path.suffix.casefold() in (".md", ".markdown"))
    return markdown, entries


def check_repository(
        root: Path, checkout_roots: dict[str, Path] | None = None,
        skipped: set[str] | None = None,
        checkout_revisions: dict[str, str] | None = None
) -> tuple[list[PurePosixPath], list[Issue], list[str]]:
    markdown, entries = tracked_entries(root)
    checkouts = {name: tracked_entries(path, (checkout_revisions or {}).get(name))[1]
                 for name, path in (checkout_roots or {}).items()}
    issues = []
    used: set[str] = set()
    sources: dict[PurePosixPath, str] = {}
    for relative in markdown:
        path = root.joinpath(*relative.parts)
        try:
            # utf-8-sig: a leading BOM would otherwise ride on the first line
            # and hide a directive or a fence from every check below.
            sources[relative] = path.read_text(encoding="utf-8-sig")
        except (OSError, UnicodeError) as error:
            issues.append(Issue(relative.as_posix(), 1,
                                f"cannot read tracked Markdown as UTF-8: {error}"))

    # Anchors first: a link may point at a heading in any other tracked document.
    anchors = {relative: _anchors(visible_lines(relative, text)[0])
               for relative, text in sources.items()}
    seen_trees: set[tuple[PurePosixPath, TreeDirective]] = set()
    for relative, text in sources.items():
        issues.extend(check_text(relative, text, entries, anchors, used,
                                 checkouts, skipped, seen_trees))
    if sources:
        for source, directive in sorted(
                _REQUIRED_TREES - seen_trees,
                key=lambda item: (item[0], item[1].checkout,
                                  item[1].exhaustive)):
            tag = " ".join(filter(None, (
                _TREE_TAG, directive.checkout,
                _TREE_EXHAUSTIVE if directive.exhaustive else "")))
            issues.append(Issue(
                source.as_posix(), 1,
                f"required docs-check directive {tag!r} is missing, so its "
                "tree goes unchecked"))
    effects_text = sources.get(PurePosixPath(_EFFECTS_TREE_ROW), "")
    roster_claimed = bool(
        EFFECTS_ROW_RE.search(effects_text)
        or EFFECTS_DIAGRAM_RE.search(effects_text)
        or any(macro == "HS_EFFECT_LIST"
               and pattern.search(sources.get(PurePosixPath(document), ""))
               for document, pattern, macro, _ in CARDINALITY_CLAIMS))
    if EFFECT_ROSTER_SOURCE in entries:
        try:
            header = root.joinpath(
                *EFFECT_ROSTER_SOURCE.parts).read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            header = ""
        roster = effect_roster(header)
        issues.extend(effects_row_issues(
            effects_text, entries, roster or None))
        try:
            playlist_header = root.joinpath(
                *PHANTASM_PLAYLIST_SOURCE.parts).read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            playlist_header = ""
        products = shader_product_group(header)
        issues.extend(roster_claim_issues(
            sources, roster, phantasm_roster(playlist_header), products))
        if products:
            issues.extend(composed_roster_issues(
                root, sources.get(PurePosixPath(_EFFECTS_REFERENCE), ""),
                entries, products))
    elif roster_claimed:
        issues.append(Issue(
            EFFECT_ROSTER_SOURCE.as_posix(), 1,
            "tracked documentation states effect-roster cardinalities, but "
            f"{EFFECT_ROSTER_SOURCE} is not tracked"))
    if _DOXYFILE in entries:
        try:
            doxyfile = root.joinpath(*_DOXYFILE.parts).read_text(encoding="utf-8")
        except (OSError, UnicodeError) as error:
            issues.append(Issue(_DOXYFILE.as_posix(), 1,
                                f"cannot read as UTF-8: {error}"))
        else:
            predefined = doxyfile_predefined(doxyfile)
            issues.extend(doxyfile_predefined_issues(
                predefined,
                unreferenced_predefined(
                    root,
                    sorted(entry for entry in entries
                           if entry.suffix in _DOXYGEN_SOURCE_SUFFIXES),
                    {name for _, name in predefined})))
    return markdown, sorted(issues), _stale_allowances(entries, used, checkouts)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check tracked Markdown fences and repository links.")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--sync", action="store_true",
                        help="refresh repository maps and derived counts before checking")
    parser.add_argument("--auto-checkout", action="store_true",
                        help="use locally available pinned Daydream; fail if unavailable (pass --checkout or --skip-checkout daydream)")
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
        revisions = {}
        if args.sync or args.auto_checkout or checkout_roots:
            import docs_sync
            if args.auto_checkout and "daydream" not in checkout_roots:
                checkout = docs_sync.discover_daydream(args.root.resolve())
                if checkout is None:
                    raise ValueError("pinned daydream checkout unavailable; supply --checkout daydream=PATH")
                else:
                    checkout_roots["daydream"] = checkout
            revisions = docs_sync.checkout_revisions(checkout_roots)
            if args.sync:
                docs_sync.sync_repository(args.root.resolve(), checkout_roots, revisions)
        markdown, issues, stale = check_repository(
            args.root.resolve(), checkout_roots, skipped, revisions)
    except (OSError, subprocess.SubprocessError, UnicodeError, ValueError) as error:
        print(f"[docs-check] tooling error: {error}", file=sys.stderr)
        return 2
    # An unvalidated tree fence or sibling link is not a pass: it went ungated,
    # so it can be edited through a green run. Accepting that is allowed --
    # a developer without the sibling checkout must still be able to run the
    # checker -- but only by naming the checkout in --skip-checkout, and the
    # verdict line says so either way.
    unvalidated = skipped - set(args.skip_checkout)
    # Warning only, like a stale UNTRACKED_ALLOWED entry: a name that skips
    # nothing is a stale exemption, and dropping it must not red a commit.
    unused_skips = sorted(set(args.skip_checkout) - skipped)
    if unused_skips:
        print(f"::warning::--skip-checkout names no unvalidated tree fence "
              f"or link: {', '.join(unused_skips)}")
    if skipped:
        print(f"::{'error' if unvalidated else 'warning'}::tree fences and "
              f"sibling links NOT validated - no --checkout root for: "
              f"{', '.join(sorted(skipped))}")
    # Warning only: dropping the last citation of an exempt path improves the
    # tree and must not red the build for whoever lands that commit.
    if markdown and stale:
        print("::warning::UNTRACKED_ALLOWED in tools/docs_check.py is stale - "
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
        print(f"[docs-check] FAIL - checkouts unvalidated: "
              f"{', '.join(sorted(unvalidated))} - pass --checkout NAME=PATH, or "
              f"--skip-checkout NAME to accept them unvalidated", file=sys.stderr)
        return 1
    note = (f"; checkouts NOT validated: {', '.join(sorted(skipped))}"
            if skipped else "")
    print(f"[docs-check] PASS - {len(markdown)} tracked Markdown file(s){note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
