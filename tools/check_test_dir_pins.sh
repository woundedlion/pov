#!/usr/bin/env bash
# Assert every tools/*_tests/ suite directory is named by a check_test_files.sh
# line in the CI workflow.
#
# tools/check_test_files.sh pins the file count of a suite the workflow already
# names, which closes the direction where a test file is deleted or renamed out
# of a glob. It says nothing about a suite the workflow never names: a new
# tools/<x>_tests/ directory has no pin and no discover step, so it runs
# nowhere and every job stays green. This is the guard for that direction --
# the suite directories in the tree must all appear in a pin.
#
# usage: check_test_dir_pins.sh
set -eu

workflow=.github/workflows/ci.yml

if [ "$#" -ne 0 ]; then
  echo "usage: $0" >&2
  exit 2
fi

if [ ! -f "$workflow" ]; then
  echo "$0: no $workflow -- run from the repository root" >&2
  exit 2
fi

# compgen expands the glob without leaving an unmatched pattern behind as a
# literal, like check_test_files.sh; a plain file matching the glob is not a
# suite.
dirs=()
while IFS= read -r dir; do
  [ -d "$dir" ] && dirs+=("$dir")
done < <({ compgen -G 'tools/*_tests' || true; } | sort)

if [ "${#dirs[@]}" -eq 0 ]; then
  echo "::error::no tools/*_tests directories found -- the glob is broken"
  exit 2
fi

unpinned=()
for dir in "${dirs[@]}"; do
  if ! grep -qE "check_test_files\.sh .*${dir}/" "$workflow"; then
    unpinned+=("$dir")
  fi
done

if [ "${#unpinned[@]}" -ne 0 ]; then
  echo "::error::no check_test_files.sh pin in $workflow for: ${unpinned[*]} -- a suite the workflow never names runs nowhere. Add a pin and a discover step."
  exit 1
fi
echo "$workflow pins all ${#dirs[@]} tools/*_tests suite(s)."
