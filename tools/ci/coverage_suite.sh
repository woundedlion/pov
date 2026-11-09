#!/usr/bin/env bash
# Run every run_tests module of the instrumented build as its own process,
# HS_COVERAGE_JOBS at a time (default 4). LLVM_PROFILE_FILE must carry %p so
# each process writes its own profile.
#
# usage: coverage_suite.sh [run_tests path]
set -euo pipefail

runner="${1:-build/tests/tests/run_tests}"
logs=cov/logs
mkdir -p "$logs"
case "${LLVM_PROFILE_FILE:-}" in
  *%p*) ;;
  *) echo "coverage_suite: LLVM_PROFILE_FILE needs %p so processes do not share a profile" >&2; exit 2 ;;
esac

# The slowest modules start first so the last ones to finish are short.
heavy=(effects effects_smoke shapeshifter_oracle conway_soak plot_scan composed_effect pullback)
"$runner" --list | sed -n 's/^  //p' > "$logs/all"
{
  for module in "${heavy[@]}"; do grep -qx "$module" "$logs/all" && echo "$module"; done
  grep -vxF -f <(printf '%s\n' "${heavy[@]}") "$logs/all"
} > "$logs/order"

run_module() {
  if "$runner" "$1" > "$logs/$1.log" 2>&1; then
    echo "PASS $1"
  else
    echo "FAIL $1"
    tail -n 80 "$logs/$1.log"
    return 1
  fi
}
export -f run_module
export runner logs
# shellcheck disable=SC2016 # run_module expands $1 in the child shell
xargs -P "${HS_COVERAGE_JOBS:-4}" -n 1 bash -c 'run_module "$1"' _ < "$logs/order"
