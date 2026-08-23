#!/usr/bin/env bash
# Assert every Python test-suite directory in the tree is named by a
# require_test_files.sh line in the CI workflow.
#
# tools/require_test_files.sh proves a discovered suite is non-empty. It says
# nothing about a suite the workflow never names: a new suite directory has no
# discovery step, so it runs nowhere and every job stays green. This is the
# guard for that direction -- every tracked suite directory must appear in a
# require call.
#
# The suites are read off the tracked test files rather than a glob of their
# parent directories, so one landing outside tools/*_tests -- as
# hardware/phantasm/gen/tests already does -- is swept without a list to widen.
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

mapfile -t dirs < <(git ls-files -- '*/test*.py' | sed 's:/[^/]*$::' | sort -u)

if [ "${#dirs[@]}" -eq 0 ]; then
  echo "::error::no test-suite directories found -- the sweep is broken"
  exit 2
fi

unpinned=()
for dir in "${dirs[@]}"; do
  if ! grep -qE "require_test_files\.sh .*${dir}/" "$workflow"; then
    unpinned+=("$dir")
  fi
done

if [ "${#unpinned[@]}" -ne 0 ]; then
  echo "::error::no require_test_files.sh call in $workflow for: ${unpinned[*]} -- a suite the workflow never names runs nowhere. Add a require call and a discover step."
  exit 1
fi
echo "$workflow discovers all ${#dirs[@]} test suite(s)."
