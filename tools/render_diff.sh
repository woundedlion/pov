#!/usr/bin/env bash
# Compare render hashes between a base commit and the working tree: proves a
# refactor or performance change left every rendered pixel unchanged.
#
# usage: render_diff.sh <base> [build-dir]
#   <base>       commit to compare against, built in a dedicated worktree
#                (../temp/<repo>-render-diff, or HS_RENDER_DIFF_DIR)
#   [build-dir]  the working tree's configured test build (default build/tests)
# Exits 0 when every hash matches, 1 when any differs.
set -euo pipefail

base="${1:?usage: render_diff.sh <base> [build-dir]}"
root="$(git rev-parse --show-toplevel)"
tip_build="$(cd "$root" && realpath -m "${2:-build/tests}")"
common="$(git rev-parse --path-format=absolute --git-common-dir)"
main="$(dirname "$common")"
dir="${HS_RENDER_DIFF_DIR:-$(dirname "$main")/temp/$(basename "$main")-render-diff}"
modules=(plot_parity hyper_lattice composed_effect shader_chain effects)

sha="$(git rev-parse --verify "$base^{commit}")"
if [ -e "$dir/.git" ]; then
  git -C "$dir" checkout --quiet --detach --force "$sha"
  git -C "$dir" clean -ffdq
else
  git worktree add --quiet --detach "$dir" "$sha"
fi
[ -f "$dir/build/tests/CMakeCache.txt" ] ||
  (cd "$dir" && cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF >/dev/null)

hashes() {
  cmake --build "$1" --target run_tests >/dev/null
  for module in "${modules[@]}"; do
    HS_PRINT_RENDER_HASHES=1 "$1/tests/run_tests" "$module" | grep '^render-hash ' || true
  done
}

scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT
echo "render_diff: building $base ($sha) in $dir"
hashes "$dir/build/tests" > "$scratch/base"
echo "render_diff: building the working tree in $tip_build"
hashes "$tip_build" > "$scratch/tip"
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
