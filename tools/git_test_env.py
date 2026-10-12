"""Isolated environments for Git fixture repositories."""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def isolated_env() -> dict[str, str]:
    """Return the process environment with no inherited repository or user git configuration."""
    local = subprocess.run(["git", "rev-parse", "--local-env-vars"],
                           check=True, capture_output=True, text=True).stdout.splitlines()
    env = {key: value for key, value in os.environ.items() if key not in local}
    env.update(GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_SYSTEM=os.devnull,
               GIT_AUTHOR_NAME="fixture", GIT_AUTHOR_EMAIL="fixture@example.invalid",
               GIT_COMMITTER_NAME="fixture", GIT_COMMITTER_EMAIL="fixture@example.invalid",
               HS_PYTHON=sys.executable)
    return env


def git_run(root: Path, *args: str, env: dict[str, str] | None = None,
            check: bool = True, text: bool = True) -> subprocess.CompletedProcess:
    """Run `git -C root args` with captured output, in `env` or an isolated environment."""
    return subprocess.run(["git", "-C", str(root), *args],
                          env=isolated_env() if env is None else env,
                          capture_output=True, text=text, check=check)


def require_tool(test_case, name: str) -> str:
    """Locate a tool, failing in CI and skipping locally when unavailable."""
    path = shutil.which(name)
    if path is None:
        if os.environ.get("CI"):
            test_case.fail(f"{name} is required in CI")
        test_case.skipTest(f"{name} is not installed")
    return path
