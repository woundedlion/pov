#!/usr/bin/env bash
# Require lint coverage of every tracked Python source.
#
# A linter handed nothing exits 0, so a broken ruff.toml exclude list or
# .gitignore rule would leave the lint gate green over no files at all.
#
# usage: ruff_selection_guard.sh
set -eu

if [ "$#" -ne 0 ]; then
  echo "usage: $0" >&2
  exit 2
fi

tmp=$(mktemp)
trap 'rm -f -- "$tmp"' EXIT

# --show-files reports the selection; emptiness is the gate, not its status.
ruff check --no-cache --show-files . > "$tmp" || true

python - "$tmp" <<'PY'
import pathlib
import subprocess
import sys
expected = subprocess.check_output(
    ["git", "ls-files", "-z", "--", "*.py"], text=True).split("\0")
expected = {pathlib.Path(path).resolve() for path in expected if path}
selected = {pathlib.Path(path).resolve()
            for path in pathlib.Path(sys.argv[1]).read_text().splitlines() if path}
missing = sorted(expected - selected)
if not expected or missing:
    print("ruff omitted tracked sources:", *(str(path) for path in missing), sep="\n")
    sys.exit(1)
PY
