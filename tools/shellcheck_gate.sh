#!/usr/bin/env bash
# Lint every tracked shell file: the gate scripts and the git hooks, where a
# defect of this class passes a gate silently.
#
# xargs handed an empty list runs nothing and exits 0, so the selection is
# asserted non-empty before it is linted.
#
# usage: shellcheck_gate.sh
set -euo pipefail

if [ "$#" -ne 0 ]; then
  echo "usage: $0" >&2
  exit 2
fi

tmp=$(mktemp)
trap 'rm -f -- "$tmp"' EXIT

git ls-files -- '*.sh' '.githooks/*' > "$tmp"

if [ ! -s "$tmp" ]; then
  echo "no shell files selected -- the shell path list is broken"
  exit 1
fi

xargs -d '\n' shellcheck -x < "$tmp"

while IFS= read -r action; do
  awk '
    /^[[:space:]]+run: \|$/ {
      match($0, /[^ ]/); indent = RSTART - 1; body_indent = 0; active = 1; next
    }
    active && /^[[:space:]]*$/ { print; next }
    active {
      match($0, /[^ ]/)
      if (RSTART - 1 > indent) {
        if (!body_indent) body_indent = RSTART
        print substr($0, body_indent); next
      }
    }
    { active = 0 }
  ' "$action" > "$tmp"
  if ! grep -q '[^[:space:]]' "$tmp"; then
    echo "no shell body extracted from $action" >&2
    exit 1
  fi
  shellcheck -s bash - < "$tmp"
done < <(git ls-files -- '.github/actions/*/action.yml' '.github/actions/*/action.yaml')
