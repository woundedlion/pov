#!/usr/bin/env bash
# Assert every tracked file whose .gitattributes entry declares an `eol` value
# actually carries those line endings, in the index blob and in the working copy.
#
# A working copy that diverged from its eol=lf blob -- checked out before the
# attribute was added, or rewritten by a CRLF-emitting tool -- is invisible to
# `git status`, so the shell, clang-format and mirror gates then read different
# bytes locally than they do on a fresh checkout.
#
# usage: eol_gate.sh [--fix-worktree]
set -eu

if [ "$#" -gt 1 ] || { [ "$#" -eq 1 ] && [ "$1" != --fix-worktree ]; }; then
  echo "usage: $0 [--fix-worktree]" >&2
  exit 2
fi

if [ "${1:-}" = --fix-worktree ]; then
  python=${HS_PYTHON:-python3}
  if [ -z "${HS_PYTHON:-}" ] && ! "$python" --version >/dev/null 2>&1; then
    python=python
  fi
  "$python" - <<'PY'
from pathlib import Path
import subprocess

records = subprocess.check_output(["git", "ls-files", "--eol", "-z"])
for record in records.split(b"\0"):
    if not record:
        continue
    info, raw_path = record.split(b"\t", 1)
    if b"eol=lf" not in info or not any(eol in info for eol in (b"w/crlf", b"w/mixed")):
        continue
    path = Path(raw_path.decode("utf-8"))
    data = path.read_bytes()
    path.write_bytes(data.replace(b"\r\n", b"\n"))
    print(f"normalized: {path}")
PY
fi

tmp=$(mktemp)
trap 'rm -f -- "$tmp"' EXIT

git ls-files --eol -z > "$tmp"

declared=0
diverged=0

# Record layout: "i/<kind>  w/<kind>  attr/<attributes>" then a TAB then the
# path. `none` is a file with no line ending at all and `-text` is binary
# content; neither contradicts a declared eol.
while IFS= read -r -d '' record; do
  info=${record%%$'\t'*}
  path=${record#*$'\t'}

  case $info in
    *'eol=lf'*) want=lf ;;
    *'eol=crlf'*) want=crlf ;;
    *) continue ;;
  esac
  declared=$((declared + 1))

  index=${info#i/}
  index=${index%% *}
  work=${info#*w/}
  work=${work%% *}

  for probe in "index:$index" "worktree:$work"; do
    where=${probe%%:*}
    have=${probe#*:}
    [ "$where" != worktree ] || [ -n "$have" ] || continue
    expected=$want
    [ "$where" != index ] || expected=lf
    case $have in
      "$expected" | none | -text) continue ;;
    esac
    echo "::error file=$path::$where line endings are $have, expected $expected (.gitattributes eol=$want)"
    diverged=$((diverged + 1))
  done
done < "$tmp"

if [ "$declared" -eq 0 ]; then
  echo "no tracked file declares an eol attribute -- the selection is broken" >&2
  exit 1
fi

if [ "$diverged" -ne 0 ]; then
  cat >&2 <<'REMEDY'

To normalize CRLF working copies, preserving edits: just normalize-eol
To restore a diverged index blob:     git add --renormalize -- <path>
REMEDY
  exit 1
fi

echo "eol clean: $declared files declare an eol attribute."
