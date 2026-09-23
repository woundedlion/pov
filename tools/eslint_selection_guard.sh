#!/usr/bin/env bash
# Require lint coverage of every tracked JavaScript source.
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

node --input-type=commonjs - "$tmp" <<'JS'
const fs = require('node:fs');
const path = require('node:path');
const { execFileSync } = require('node:child_process');
const expected = execFileSync('git', ['ls-files', '-z', '--', '*.js', '*.mjs', '*.cjs'],
  { encoding: 'utf8' }).split('\0').filter(Boolean).map(name => path.resolve(name));
const selected = new Set(JSON.parse(fs.readFileSync(process.argv[2], 'utf8'))
  .filter(row => typeof row.filePath === 'string').map(row => path.resolve(row.filePath)));
const missing = expected.filter(name => !selected.has(name));
if (!expected.length || missing.length) {
  console.error('eslint omitted tracked sources:', ...missing);
  process.exit(1);
}
JS
