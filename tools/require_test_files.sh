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

# Scan conventional test names independently of the runner's extension/glob.
root=$(dirname "$pattern")
case "$pattern" in
  *.js|*.mjs|*.cjs|*.ts)
    family=javascript
    if git rev-parse --show-toplevel >/dev/null 2>&1; then
      root=.
    fi
    ;;
  *.py) family=python ;;
  *) family=other ;;
esac
unreachable=0
scan=$(mktemp)
trap 'rm -f "$scan"' EXIT
find "$root" -type d \( -name node_modules -o -name .git -o -name build \) -prune -o -type f -print0 >"$scan"
while IFS= read -r -d '' candidate; do
  candidate=${candidate#./}
  case "$family:$candidate" in
    javascript:*.test.js|javascript:*.test.mjs|javascript:*.test.cjs|javascript:*.test.ts|javascript:*.spec.js|javascript:*.spec.mjs|javascript:*.spec.cjs|javascript:*.spec.ts|python:*/test*.py) ;;
    *) continue ;;
  esac
  reached=0
  for file in "${files[@]}"; do
    if [ "$candidate" = "${file#./}" ]; then
      reached=1
      break
    fi
  done
  if [ "$reached" -eq 0 ]; then
    echo "::error::test file '$candidate' is unreachable from '$pattern'"
    unreachable=1
  fi
done <"$scan"
[ "$unreachable" -eq 0 ] || exit 1

printf '%s: %d test file(s) discovered\n' "$pattern" "${#files[@]}"
printf '  %s\n' "${files[@]}"
