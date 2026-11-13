# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Extract roster rows from continued X-macro bodies and the playlist .def."""

import re


_COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)
# A playlist row starts its line, so `//` lines never match.
_PLAYLIST_ROW_RE = re.compile(r"^[ \t]*(show|exclude)\([ \t]*(\w+)", re.MULTILINE)

#: The Phantasm playlist, relative to the repository root.
PHANTASM_PLAYLIST_SOURCE = "targets/Phantasm/phantasm_playlist.def"


def playlist_shows(source: str) -> list[str]:
    """Return the effect names of the playlist's show rows, in show order."""
    return [name for kind, name in _PLAYLIST_ROW_RE.findall(source)
            if kind == "show"]


def named_rows(body: str) -> list[str]:
    """Return names from X(Name, ...) roster rows."""
    return re.findall(r"\bX\(\s*(\w+)\s*,", body)


def macro_body(source: str, macro: str) -> str | None:
    """Return a comment-free, line-spliced body, or None if absent."""
    head = re.search(r"(?m)^[ \t]*#[ \t]*define[ \t]+" + re.escape(macro)
                     + r"[ \t]*\([ \t]*X[ \t]*\)", source)
    if not head:
        return None
    body = []
    for line in source[head.end():].splitlines():
        body.append(line)
        if not line.rstrip().endswith("\\"):
            break
    return _COMMENT_RE.sub("", re.sub(r"\\[ \t]*\n", " ", "\n".join(body)))
