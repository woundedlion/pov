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

cd "$(git rev-parse --show-toplevel)"

tmp=$(mktemp)
trap 'rm -f -- "$tmp"' EXIT

# --show-files reports the selection; emptiness is the gate, not its status.
ruff check --no-cache --show-files . > "$tmp" || true

ruff_python=
for candidate in "${HS_PYTHON:-}" python3 python; do
  if "$candidate" --version >/dev/null 2>&1; then ruff_python="$candidate"; break; fi
done
if [ -z "$ruff_python" ]; then
  echo "no working Python found (set HS_PYTHON)" >&2
  exit 1
fi
"$ruff_python" - "$tmp" <<'PY'
import pathlib
import subprocess
import sys
expected = subprocess.check_output(
    ["git", "ls-files", "-z", "--", "*.py"], text=True).split("\0")
expected = {pathlib.Path(path).resolve() for path in expected if path}
selected = {pathlib.Path(path).resolve()
            for path in pathlib.Path(sys.argv[1]).read_text().splitlines() if path}
missing = sorted(expected - selected)
if not expected:
    print("ruff selection guard found no tracked Python sources")
    sys.exit(1)
if missing:
    print("ruff omitted tracked sources:", *(str(path) for path in missing), sep="\n")
    sys.exit(1)
PY
