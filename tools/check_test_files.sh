#!/usr/bin/env bash
# Assert a test-file glob matches exactly the pinned number of files.
#
# Every CI step that discovers its tests by pattern -- `unittest discover`,
# `node --test <glob>` -- runs fewer suites and stays green when a file is
# deleted, renamed out of the pattern, or the glob matches nothing at all. This
# is the wiring pin for that direction, shared by every such step so the guard
# is written once.
#
# The count is EXACT, not a floor: a floor degrades to "at least one file
# exists" as the suite grows, so adding a test file must move the pin here too.
#
# usage: check_test_files.sh <count> <glob>
set -eu

if [ "$#" -ne 2 ]; then
  echo "usage: $0 <count> <glob>" >&2
  exit 2
fi
want=$1
pattern=$2

case "$want" in
  '' | *[!0-9]*)
    echo "$0: <count> must be a number, got '$want'" >&2
    exit 2
    ;;
esac

# compgen expands the glob without leaving an unmatched pattern behind as a
# literal; it exits non-zero on no match, which is a legitimate (failing) count.
files=()
while IFS= read -r file; do
  files+=("$file")
done < <({ compgen -G "$pattern" || true; } | sort)

found=${#files[@]}
if [ "$found" -ne "$want" ]; then
  echo "::error::expected exactly $want file(s) matching '$pattern', found $found: ${files[*]-}"
  exit 1
fi
echo "$pattern: $found file(s), matching the pin."
