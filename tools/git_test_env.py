"""Isolated environments for Git fixture repositories."""

import os
import subprocess
import sys


def isolated_env() -> dict[str, str]:
    local = subprocess.run(["git", "rev-parse", "--local-env-vars"],
                           check=True, capture_output=True, text=True).stdout.splitlines()
    env = {key: value for key, value in os.environ.items() if key not in local}
    env.update(GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_SYSTEM=os.devnull,
               GIT_AUTHOR_NAME="fixture", GIT_AUTHOR_EMAIL="fixture@example.invalid",
               GIT_COMMITTER_NAME="fixture", GIT_COMMITTER_EMAIL="fixture@example.invalid",
               HS_PYTHON=sys.executable)
    return env
