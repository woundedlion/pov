#!/bin/bash
# profile_islamic_big.sh [seconds=24] [window=8] [extra build flags...]
#
# Repeats only IslamicStars' heaviest build chain so optimization captures reach
# the troublesome ambo raster in seconds instead of walking the whole roster.
set -euo pipefail

SECONDS_ARG=${1:-24}
WINDOW=${2:-8}
if [ "$#" -ge 1 ]; then shift; fi
if [ "$#" -ge 1 ]; then shift; fi

SUFFIX=""
if [ -n "${HS_PROFILE_DEEP:-}" ] && [ "${HS_PROFILE_DEEP:-}" != "0" ]; then
  SUFFIX="_deep"
fi

export HS_PROFILE_OUT="build/prof/islamicstars_big_ship${SUFFIX}.log"
export HS_PROFILE_EXPECT_SHAPE=truncatedIcosidodecahedron_truncate50d_ambo_dual
TREE=${HS_PROFILE_TREE:-$(git -C "$(dirname "$0")" rev-parse --show-toplevel)}
SHAPE=$(awk -v wanted="$HS_PROFILE_EXPECT_SHAPE" '
  /inline constexpr Entry islamic_registry\[\]/ { active=1; next }
  active && /^};/ { exit }
  active && /^[[:space:]]*\{"/ {
    name=$0; sub(/^[^"]*"/, "", name); sub(/".*/, "", name)
    if (name == wanted) print ordinal + 0
    ordinal++
  }
' "$TREE/core/mesh/solids.h")
[ -n "$SHAPE" ] || { echo "profile shape is absent: $HS_PROFILE_EXPECT_SHAPE" >&2; exit 1; }
bash "$(dirname "$0")/profile_one.sh" IslamicStars profile "$SECONDS_ARG" "$WINDOW" \
  "-D HS_ISLAMICSTARS_PROFILE_SHAPE=$SHAPE" \
  "-D HS_PROFILE_TRANS_SPEED=4" "$@"
