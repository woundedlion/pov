#!/usr/bin/env bash
# Fail when a glob-discovered test suite is empty.
#
# Discovery should absorb additions without a maintenance edit, while an empty
# match must never be reported as a passing suite.
#
# usage: require_test_files.sh <glob>
set -eu

if [ "$#" -ne 1 ]; then
  echo "usage: $0 <glob>" >&2
  exit 2
fi
pattern=$1

files=()
while IFS= read -r file; do
  files+=("$file")
done < <({ compgen -G "$pattern" || true; } | sort)

if [ "${#files[@]}" -eq 0 ]; then
  echo "::error::no test files match '$pattern'"
  exit 1
fi

printf '%s: %d test file(s) discovered\n' "$pattern" "${#files[@]}"
printf '  %s\n' "${files[@]}"
