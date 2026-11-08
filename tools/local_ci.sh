#!/usr/bin/env bash
# Run every CI job this machine can run, stopping at the first failure. Each
# job calls the same commands as its .github/workflows/ci.yml job.
#
# usage: local_ci.sh                  every local job in this clean checkout
#        local_ci.sh <job>...         the named jobs in this clean checkout
#        local_ci.sh --commit <sha>...  every local job on each commit, in the
#                                     dedicated local-CI worktree (pre-push)
#        local_ci.sh --new-guards <tip> <remote-sha>
#                                     new guards in the pushed range need a
#                                     death case (pre-push)
#        local_ci.sh --list           the job table
#        local_ci.sh --check          verify the table names every ci.yml job
set -euo pipefail

# ci.yml job | the function that runs it here, or "skip: <reason>"
JOBS='
lint|job_lint
clang-format|job_clang_format
license-headers|job_license_headers
docs-markdown|job_docs_markdown
docs-images|job_docs_images
lut-provenance|job_lut_provenance
reaction-graph-provenance|job_reaction_graph_provenance
gamut-lut-provenance|job_gamut_lut_provenance
patterns-provenance|job_patterns_provenance
screenshot-gallery|job_screenshot_gallery
python-tests|job_python_tests
scripts-unit-tests|job_scripts_unit_tests
windows-tests|job_windows_tests
tests|job_tests
srgb-decode-provenance|job_srgb_decode_provenance
optimized-tests|job_optimized_tests
wasm|job_wasm
teensy-size|job_teensy_size
teensy-warnings|job_teensy_warnings
pcb-tests|job_pcb_tests
docs-doxygen|job_docs_doxygen
sanitizers|skip: the ASan and UBSan runtimes are Linux-only here
thread-sanitizer|skip: TSan is Linux-only here
code-coverage|skip: the emsdk LLVM has no Windows profile runtime
'

job_lint() { just lint; }
job_clang_format() { clang-format --version && bash tools/clang_format_gate.sh; }
job_license_headers() { python tools/license_check.py; }
# CI checks the docs against daydream at the revision build_pins.py pins.
job_docs_markdown() {
  local daydream="${DAYDREAM_DIR:-../daydream}" target="$PWD/build/local-ci/daydream" pin
  pin="$(python tools/build_pins.py daydream)"
  if [ -e "$target/.git" ]; then
    git -C "$target" checkout --quiet --detach --force "$pin"
  else
    git -C "$daydream" worktree prune
    git -C "$daydream" worktree add --quiet --detach "$target" "$pin"
  fi
  bash tools/ci/docs_markdown.sh "$target"
}
job_docs_images() { python tools/docs_images.py; }
job_lut_provenance() { bash tools/ci/lut_provenance.sh; }
job_reaction_graph_provenance() { bash tools/ci/reaction_graph_provenance.sh; }
job_gamut_lut_provenance() { python tools/gen_gamut_lut.py --check; }
job_patterns_provenance() { bash tools/ci/patterns_provenance.sh; }
job_screenshot_gallery() { bash tools/ci/screenshot_gallery.sh; }
job_python_tests() { bash tools/ci/python_tests.sh; }
job_scripts_unit_tests() { bash tools/ci/scripts_tests.sh; }
job_windows_tests() {
  bash tools/ci/native_suite.sh --exclude-generated
  bash tools/ci/render_goldens.sh
}
# The Linux job's suite runs under windows-tests; this is its remaining step.
job_tests() { bash tools/ci/p1_plot_scan.sh; }
job_srgb_decode_provenance() { bash tools/ci/srgb_decode_provenance.sh; }
# fastmath stays in CI: fast-math codegen differs between this toolchain and
# Linux clang, so its bit-exact results here do not predict CI's.
job_optimized_tests() {
  local flavor
  for flavor in ieee ndebug-death; do
    HS_CI_BUILD_DIR="build/optimized-$flavor" bash tools/ci/optimized.sh "$flavor"
    bash tools/ci/render_goldens.sh "build/optimized-$flavor"
  done
}
job_wasm() {
  # shellcheck disable=SC1091
  . "${EMSDK:-../emsdk}/emsdk_env.sh" >/dev/null
  rm -rf build/local-ci/daydream-install
  bash tools/ci/wasm.sh "$PWD/build/local-ci/daydream-install"
}
job_teensy_size() { bash tools/ci/teensy_size.sh "$PWD/build/local-ci/pio"; }
job_teensy_warnings() { bash tools/ci/teensy_warnings.sh; }
job_pcb_tests() {
  export KICAD_CONFIG_HOME="$PWD/build/local-ci/kicad-config"
  python tools/ci/kicad_config.py "$KICAD_CONFIG_HOME"
  python tools/run_python_tests.py --suite hardware/phantasm/gen/tests
}
job_docs_doxygen() {
  bash tools/ci/apply_doxygen_theme.sh
  local status=0
  bash tools/ci/docs_doxygen.sh || status=$?
  git checkout -- Doxyfile
  return "$status"
}

table() { printf '%s\n' "$JOBS" | sed '/^$/d'; }
runner_of() { table | awk -F'|' -v job="$1" '$1 == job { print $2 }'; }

check_table() {
  local workflow="${HS_CI_WORKFLOW:-.github/workflows/ci.yml}" missing extra
  # Job keys are the two-space-indented keys after jobs:.
  awk '
    /^jobs:[[:space:]]*(#.*)?$/ { injobs = 1; next }
    injobs && /^[^[:space:]]/ { injobs = 0 }
    injobs && /^  [A-Za-z0-9_-]+:[[:space:]]*(#.*)?$/ {
      sub(/^  /, ""); sub(/:[[:space:]]*(#.*)?$/, ""); print
    }
  ' "$workflow" | grep -vx ci-green | sort > "$scratch/workflow-jobs"
  [ -s "$scratch/workflow-jobs" ] || { echo "local_ci: no jobs parsed from $workflow" >&2; return 1; }
  table | cut -d'|' -f1 | sort > "$scratch/table-jobs"
  missing="$(comm -23 "$scratch/workflow-jobs" "$scratch/table-jobs")"
  extra="$(comm -13 "$scratch/workflow-jobs" "$scratch/table-jobs")"
  [ -z "$missing" ] || echo "local_ci: ci.yml jobs with no local runner or skip reason: $missing" >&2
  [ -z "$extra" ] || echo "local_ci: table jobs absent from ci.yml: $extra" >&2
  [ -z "$missing$extra" ]
}

run_jobs() {
  local logs="build/local-ci/logs" job runner start
  for job in "$@"; do
    [ -n "$(runner_of "$job")" ] || { echo "local_ci: unknown job $job" >&2; return 2; }
  done
  if [ -n "$(git status --porcelain)" ]; then
    echo "local_ci: the checkout has uncommitted changes; run it on a clean commit." >&2
    return 1
  fi
  check_table
  mkdir -p "$logs"
  table | grep '|skip:' | sed 's/|skip: */ skipped: /' || true
  for job in "$@"; do
    runner="$(runner_of "$job")"
    case "$runner" in skip:*) continue ;; esac
    start=$SECONDS
    # A separate process keeps errexit live inside the job function.
    if bash "$0" --run-one "$runner" > "$logs/$job.log" 2>&1; then
      echo "PASS $job ($((SECONDS - start)) s)"
    else
      echo "FAIL $job ($((SECONDS - start)) s); log: $PWD/$logs/$job.log"
      tail -n 60 "$logs/$job.log"
      return 1
    fi
  done
  echo "local CI passed at $(git rev-parse --short HEAD)."
}

# Prepares the dedicated worktree for each commit and runs the commit's own copy
# of this script there. Passed commits are recorded so a retried push is free.
run_commits() {
  local common main dir passed sha entry npm_stamp lock_blob
  common="$(git rev-parse --path-format=absolute --git-common-dir)"
  main="$(dirname "$common")"
  dir="${HS_LOCAL_CI_DIR:-$(dirname "$main")/temp/$(basename "$main")-local-ci}"
  lock="$common/hs-local-ci.lock"
  passed="$common/hs-local-ci-passed"
  until mkdir "$lock" 2>/dev/null; do
    if [ -f "$lock/pid" ] && ! kill -0 "$(cat "$lock/pid")" 2>/dev/null; then
      rm -rf "$lock"
      continue
    fi
    echo "local_ci: waiting for another local CI run ($lock)"
    sleep 30
  done
  echo $$ > "$lock/pid"
  trap 'rm -rf "$lock" "$scratch"' EXIT

  while IFS= read -r entry; do
    entry="$(cygpath -u "$entry" 2>/dev/null || printf '%s' "$entry")"
    case ":$PATH:" in *":$entry:"*) ;; *) PATH="$entry:$PATH" ;; esac
  done < <(git config --get-all hs.localCiPath || true)
  if [ -z "${DAYDREAM_DIR:-}" ] && [ -d "$(dirname "$main")/daydream" ]; then
    DAYDREAM_DIR="$(dirname "$main")/daydream"
    export DAYDREAM_DIR
  fi
  if [ -z "${EMSDK:-}" ] && [ -d "$(dirname "$main")/emsdk" ]; then
    EMSDK="$(dirname "$main")/emsdk"
    export EMSDK
  fi

  for sha in "$@"; do
    if [ -f "$passed" ] && grep -qx "$sha" "$passed"; then
      echo "local_ci: $sha already passed local CI."
      continue
    fi
    if [ -e "$dir/.git" ]; then
      git -C "$dir" checkout --quiet --detach --force "$sha"
      git -C "$dir" clean -ffdq
    else
      git worktree add --quiet --detach "$dir" "$sha"
    fi
    if [ ! -d "$dir/.pio/libdeps" ] && [ -d "$main/.pio/libdeps" ]; then
      mkdir -p "$dir/.pio"
      cp -r "$main/.pio/libdeps" "$dir/.pio/"
    fi
    npm_stamp="$dir/build/local-ci/package-lock.blob"
    lock_blob="$(git rev-parse "$sha:package-lock.json")"
    if [ "$(cat "$npm_stamp" 2>/dev/null)" != "$lock_blob" ]; then
      (cd "$dir" && npm ci --no-audit --no-fund)
      mkdir -p "$(dirname "$npm_stamp")"
      echo "$lock_blob" > "$npm_stamp"
    fi
    echo "local_ci: running local CI on $sha in $dir"
    (cd "$dir" && bash tools/local_ci.sh)
    echo "$sha" >> "$passed"
  done
}

# Fails when <tip> adds a guard with no death case over the remote's <base>; a
# new branch, or a base this clone lacks, is measured from origin/master.
check_new_guards() {
  local tip=$1 base=$2
  case $base in *[!0]*) ;; *) base= ;; esac
  if [ -z "$base" ] || ! git cat-file -e "$base^{commit}" 2>/dev/null; then
    base="$(git merge-base "$tip" origin/master 2>/dev/null)" || return 0
  fi
  python tests/check_death_pins.py . --new-guards "$base" "$tip"
}

scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT
cd "$(git rev-parse --show-toplevel)"

case "${1:-}" in
  --run-one) "$2" ;;
  --new-guards) check_new_guards "$2" "$3" ;;
  --list) table ;;
  --check) check_table ;;
  --commit) shift; run_commits "$@" ;;
  '') mapfile -t all < <(table | cut -d'|' -f1); run_jobs "${all[@]}" ;;
  *) run_jobs "$@" ;;
esac
