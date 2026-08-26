#!/usr/bin/env bash
# Fail when ruff's file selection is empty.
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

if [ ! -s "$tmp" ]; then
  echo "ruff selected no files -- the ruff.toml exclude list or a .gitignore rule is broken"
  exit 1
fi
