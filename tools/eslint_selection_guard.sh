#!/usr/bin/env bash
# Fail when eslint's file selection is empty.
#
# A linter handed nothing exits 0, so a broken eslint.config.mjs ignores list
# would leave the lint gate green over no files at all. eslint reports one JSON
# object per file, so a report naming no filePath is an empty selection.
#
# usage: eslint_selection_guard.sh
set -eu

if [ "$#" -ne 0 ]; then
  echo "usage: $0" >&2
  exit 2
fi

tmp=$(mktemp)
trap 'rm -f -- "$tmp"' EXIT

# A lint failure still writes the report; the selection is what is gated here.
npx eslint . --format json > "$tmp" || true

if ! grep -q filePath "$tmp"; then
  echo "eslint selected no files -- the eslint.config.mjs ignores list is broken"
  exit 1
fi
