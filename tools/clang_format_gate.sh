#!/usr/bin/env bash
# clang-format --dry-run --Werror over the whole tracked first-party C++ set.
#
# The pathspec and the exclusion regex (vendored and the seven generated
# sources) are
# pinned against the .githooks/pre-commit copy by
# tools/build_pins.py --check. xargs handed an empty list runs nothing and exits
# 0, so the selection is asserted non-empty before it is checked.
#
# usage: clang_format_gate.sh
set -eu

if [ "$#" -ne 0 ]; then
  echo "usage: $0" >&2
  exit 2
fi

CF="${CLANG_FORMAT:-clang-format}"
major="$("$CF" --version | sed -n 's/.*version \([0-9][0-9]*\).*/\1/p')"
if [ "$major" != 22 ]; then
  echo "clang-format 22 is required; found ${major:-unknown}" >&2
  exit 1
fi

tmp=$(mktemp)
trap 'rm -f -- "$tmp"' EXIT

git ls-files -- '*.h' '*.hpp' '*.cpp' '*.cc' '*.inl' '*.ino' | grep -vE '(^|/)core/vendor/FastNoiseLite\.h$|(^|/)core/color/color_luts\.h$|(^|/)core/color/gamut_lut\.h$|(^|/)core/color/srgb_decode_lut\.h$|(^|/)core/color/triadic_palette_luts\.h$|(^|/)core/mesh/relax_bakes_generated\.h$|(^|/)core/spatial/reaction_graph\.cpp$|(^|/)tests/mindsplatter_replay_corpus\.h$' > "$tmp" || true

if [ ! -s "$tmp" ]; then
  echo "no files selected -- the pathspec or exclusion regex is broken"
  exit 1
fi

xargs -d '\n' "$CF" --dry-run --Werror --style=file < "$tmp"
