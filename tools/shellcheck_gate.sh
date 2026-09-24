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

scratch=$(mktemp -d)
tmp=$scratch/files
trap 'rm -rf -- "$scratch"' EXIT

git ls-files > "$scratch/tracked"
grep -E '\.sh$|^\.githooks/' "$scratch/tracked" > "$tmp" || [[ $? == 1 ]]

if [ ! -s "$tmp" ]; then
  echo "no shell files selected -- the shell path list is broken"
  exit 1
fi

xargs -d '\n' shellcheck -x < "$tmp"

while IFS= read -r action; do
  awk -v directory="$scratch" '
    function finish() {
      if (!has_run) return
      if (shell !~ /^(bash|sh|dash|ksh)$/ || body !~ /[^[:space:]]/) {
        print "unsupported shell or empty run body" > "/dev/stderr"; failed=1; exit 1
      }
      count++
      path=directory "/step-" count
      print body > path; close(path)
      print shell " " path
    }
    /^[[:space:]]*-[[:space:]]/ {
      finish(); shell=""; body=""; has_run=0; active=0
    }
    active {
      if ($0 ~ /^[[:space:]]*$/) { body=body "\n"; next }
      match($0, /[^ ]/)
      if (RSTART - 1 > indent) {
        if (!body_indent) body_indent=RSTART
        body=body substr($0, body_indent) "\n"; next
      }
      active=0
    }
    /^[[:space:]]+(-[[:space:]]+)?shell:/ {
      shell=$0; sub(/^[[:space:]]+(-[[:space:]]+)?shell:[[:space:]]*/, "", shell)
      next
    }
    /^[[:space:]]+(-[[:space:]]+)?run:/ {
      if ($0 !~ /run: \|[-+]?[[:space:]]*$/) {
        print "unsupported run scalar; use a literal | block" > "/dev/stderr"
        failed=1; exit 1
      }
      match($0, /run:/); indent=RSTART - 1; body_indent=0; active=1; has_run=1
    }
    END { if (!failed) { finish(); if (!count) exit 1 } }
  ' "$action" > "$tmp"
  while read -r shell body; do
    shellcheck -s "$shell" - < "$body"
  done < "$tmp"
done < <(git ls-files -- '.github/actions/*/action.yml' '.github/actions/*/action.yaml')
