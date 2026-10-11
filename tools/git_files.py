"""Git invocations shared by repository tools."""

from __future__ import annotations

import subprocess
from pathlib import Path

GIT_TIMEOUT_SECONDS = 30


def git_command(root: Path, *args: str) -> list[str]:
    """Return a git argv that runs `args` in `root`, trusting it as safe."""
    return ["git", "-c", f"safe.directory={root.as_posix()}", "-C", str(root), *args]


def git_nul_list(root: Path, *args: str,
                 timeout: float | None = GIT_TIMEOUT_SECONDS) -> list[str]:
    """Run a `-z` git listing in `root` and return its non-empty entries.

    Raises `subprocess.CalledProcessError` on a non-zero exit, and
    `UnicodeError` on output that is not UTF-8.
    """
    result = subprocess.run(git_command(root, *args), check=True,
                            stdout=subprocess.PIPE, timeout=timeout)
    return [name for name in result.stdout.decode("utf-8").split("\0") if name]
