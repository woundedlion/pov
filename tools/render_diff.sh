#!/usr/bin/env bash
# Compare render hashes between a base commit and the working tree: proves a
# refactor or performance change left every rendered pixel unchanged, or names
# the effects and presets it changed.
#
# usage: render_diff.sh [--report] <base> [build-dir] [cmake configure args...]
#   --report     hash only the render_hashes module (every roster effect at
#                every preset), print the changed renders grouped per effect,
#                and always exit 0
#   <base>       commit to compare against, built in a dedicated worktree
#                (../temp/<repo>-render-diff, or HS_RENDER_DIFF_DIR)
#   [build-dir]  the working tree's test build (default build/tests),
#                configured with the tests preset when it has no cache yet
#   [cmake configure args...]  extra arguments for those configures
# Without --report: exits 0 when every hash matches, 1 when any differs, and
# 2 when the base prints no render hashes. Modules the base lacks are skipped.
set -euo pipefail

usage="usage: render_diff.sh [--report] <base> [build-dir] [cmake configure args...]"
report=0
if [ "${1:-}" = --report ]; then
  report=1
  shift
fi
base="${1:?$usage}"
shift
build_arg=build/tests
if [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; then
  build_arg="$1"
  shift
fi
configure_args=("$@")
root="$(git rev-parse --show-toplevel)"
tip_build="$(cd "$root" && realpath -m "$build_arg")"
common="$(git rev-parse --path-format=absolute --git-common-dir)"
main="$(dirname "$common")"
dir="${HS_RENDER_DIFF_DIR:-$(dirname "$main")/temp/$(basename "$main")-render-diff}"

if [ "$report" = 1 ]; then
  target=test_render_hashes
  modules=(render_hashes)
else
  target=run_tests
  modules=(render_hashes plot_parity hyper_lattice composed_effect shader_chain effects)
fi

sha="$(git rev-parse --verify "$base^{commit}")"
if [ "$report" = 1 ] && ! git cat-file -e "$sha:tests/test_render_hashes.h" 2>/dev/null; then
  echo "render_diff: base has no render hashes ($base predates the render_hashes module)."
  exit 0
fi

checkout_base() {
  if [ -e "$dir/.git" ]; then
    git -C "$dir" checkout --quiet --detach --force "$sha" &&
      git -C "$dir" clean -ffdq
  else
    git worktree add --quiet --detach "$dir" "$sha"
  fi
}

# Configures (once) and builds the hash target of the tree at $1 in $2.
build_tree() {
  [ -f "$2/CMakeCache.txt" ] ||
    (cd "$1" && cmake --preset tests -B "$2" -DHS_INSTALL_GIT_HOOKS=OFF \
      "${configure_args[@]}" >/dev/null)
  cmake --build "$2" --target "$target" >/dev/null
}

# Prints the render-hash lines of the build at $1 for every module the
# module list in $2 names.
run_hashes() {
  local module
  for module in "${modules[@]}"; do
    grep -qx "  $module" "$2" || continue
    HS_PRINT_RENDER_HASHES=1 "$1/tests/$target" "$module" | grep '^render-hash ' || true
  done
}

scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT

if [ "$report" = 1 ]; then
  for side in base tip; do
    if [ "$side" = base ]; then
      src="$dir" build="$dir/build/tests" label="$base"
      checkout_base || { echo "render_diff: checking out $base failed; no report."; exit 0; }
    else
      src="$root" build="$tip_build" label="the working tree"
    fi
    if ! build_tree "$src" "$build" ||
        ! "$build/tests/$target" --list > "$scratch/$side.modules" ||
        ! run_hashes "$build" "$scratch/$side.modules" > "$scratch/$side"; then
      echo "render_diff: building or running $label failed; no report."
      exit 0
    fi
    # "<effect> <preset> <hash>" rows.
    sed -E 's/^render-hash ([^ ]+) p([0-9]+)=([0-9a-f]+)$/\1 \2 \3/' "$scratch/$side" |
      sort -k1,1 -k2,2n > "$scratch/$side.rows"
  done
  if [ ! -s "$scratch/base.rows" ]; then
    echo "render_diff: base has no render hashes ($base printed none)."
    exit 0
  fi
  awk '
    NR == FNR { base[$1 " " $2] = $3; next }
    {
      key = $1 " " $2
      if (!(key in base)) print $1, $2, "new"
      else if (base[key] != $3) print $1, $2, "changed"
      delete base[key]
    }
    END { for (key in base) print key, "removed" }
  ' "$scratch/base.rows" "$scratch/tip.rows" | sort -k1,1 -k3,3 -k2,2n > "$scratch/delta"
  total="$(wc -l < "$scratch/tip.rows" | tr -d ' ')"
  if [ ! -s "$scratch/delta" ]; then
    echo "render_diff: no pixel changes across $total renders."
    exit 0
  fi
  echo "render_diff: $(wc -l < "$scratch/delta" | tr -d ' ') of $total renders differ in" \
    "$(cut -d' ' -f1 "$scratch/delta" | sort -u | wc -l | tr -d ' ') effect(s):"
  # One line per effect and kind, consecutive presets collapsed to ranges.
  awk '
    function flush(   out, i, lo, hi) {
      if (n == 0) return
      out = ""
      for (i = 1; i <= n; i++) {
        lo = presets[i]
        hi = lo
        while (i < n && presets[i + 1] == hi + 1) hi = presets[++i]
        out = out (out == "" ? "" : ", ") "p" lo (hi > lo ? "-p" hi : "")
      }
      print "  " effect " " out (kind == "changed" ? "" : " (" kind ")")
      n = 0
    }
    $1 != effect || $3 != kind { flush(); effect = $1; kind = $3 }
    { presets[++n] = $2 }
    END { flush() }
  ' "$scratch/delta"
  exit 0
fi

checkout_base
echo "render_diff: building $base ($sha) in $dir"
build_tree "$dir" "$dir/build/tests"
"$dir/build/tests/tests/run_tests" --list > "$scratch/modules"
run_hashes "$dir/build/tests" "$scratch/modules" > "$scratch/base"
echo "render_diff: building the working tree in $tip_build"
build_tree "$root" "$tip_build"
run_hashes "$tip_build" "$scratch/modules" > "$scratch/tip"
if [ ! -s "$scratch/base" ]; then
  echo "render_diff: $base prints no render hashes; pick a newer base." >&2
  exit 2
fi
if diff -u --label "$base" --label "working tree" "$scratch/base" "$scratch/tip"; then
  echo "render_diff: $(wc -l < "$scratch/tip") render(s) unchanged."
else
  echo "render_diff: rendered pixels changed." >&2
  exit 1
fi
