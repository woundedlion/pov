# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Extract continued X-macro bodies for roster consumers."""

import re


_COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)


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
