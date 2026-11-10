#!/usr/bin/env bash
# Report which roster effects and presets the checked-out commit re-renders
# relative to a base commit, through render_diff.sh --report. Appends the
# report to $GITHUB_STEP_SUMMARY when set. Report only: always exits 0.
#
# usage: render_hash_report.sh [base-sha] [cmake configure args...]
#   [base-sha]  the commit to compare against; empty means HEAD's parent
set -uo pipefail

base="${1:-}"
[ $# -gt 0 ] && shift

summary() {
  printf '%s\n' "$@"
  [ -z "${GITHUB_STEP_SUMMARY:-}" ] || printf '%s\n' "$@" >> "$GITHUB_STEP_SUMMARY"
}

if [ -z "$base" ] || [ -z "${base//0/}" ]; then
  git fetch --quiet --no-tags --depth=2 origin "$(git rev-parse HEAD)" || true
  base="$(git rev-parse --verify --quiet 'HEAD^' || true)"
fi
if [ -n "$base" ] && ! git cat-file -e "$base^{commit}" 2>/dev/null; then
  git fetch --quiet --no-tags --depth=1 origin "$base" || true
fi
if [ -z "$base" ] || ! git cat-file -e "$base^{commit}" 2>/dev/null; then
  summary "## Render-hash report" "" "Base commit ${base:-(none)} is unavailable; no report."
  exit 0
fi

# Base and tip build in sibling trees; share compiler-cache entries across them.
export CCACHE_BASEDIR="${CCACHE_BASEDIR:-$(dirname "$PWD")}"
export CCACHE_NOHASHDIR=1
# Depth levers the test runner requires when CI is set.
export HS_SMOKE_FRAMES=120 HS_SKIPS_ARE_ERRORS=1 HS_BUFFER_FREE_WATCHDOG_US=30000000

start=$SECONDS
report="$(bash "$(dirname "$0")/../render_diff.sh" --report "$base" build/tests "$@" 2>&1)"
status=$?
if [ "$status" -ne 0 ]; then
  report="$(printf '%s\nrender_diff: exited with status %s; no report.' "$report" "$status")"
fi
echo "render_hash_report: $((SECONDS - start)) s"
lines="$(printf '%s\n' "$report" | grep '^render_diff: \|^  ' || true)"
headline="$(printf '%s\n' "$lines" | grep '^render_diff: ' | tail -1 | sed 's/^render_diff: //')"
headline="${headline^}"
changed="$(printf '%s\n' "$lines" | grep '^  ' | sed 's/^  /- /' || true)"
printf '%s\n' "$report"
summary "## Render-hash report" "" "Base \`$(git rev-parse --short "$base")\`: ${headline:-No output.}"
[ -z "$changed" ] || summary "" "$changed"
exit 0
