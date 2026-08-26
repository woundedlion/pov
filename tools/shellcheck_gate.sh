#!/usr/bin/env bash
# Lint every tracked shell file: the gate scripts and the git hooks, where a
# defect of this class passes a gate silently.
#
# xargs handed an empty list runs nothing and exits 0, so the selection is
# asserted non-empty before it is linted.
#
# usage: shellcheck_gate.sh
set -eu

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

xargs shellcheck -x < "$tmp"
