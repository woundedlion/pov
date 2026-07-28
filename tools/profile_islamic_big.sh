#!/bin/bash
# profile_islamic_big.sh [seconds=24] [window=8] [extra build flags...]
#
# Repeats only IslamicStars' largest build chain so optimization captures reach
# the troublesome ambo raster in seconds instead of walking the 24-shape roster.
set -eo pipefail

SECONDS_ARG=${1:-24}
WINDOW=${2:-8}
if [ "$#" -ge 1 ]; then shift; fi
if [ "$#" -ge 1 ]; then shift; fi

SUFFIX=""
if [ -n "$HS_PROFILE_DEEP" ] && [ "$HS_PROFILE_DEEP" != "0" ]; then
  SUFFIX="_deep"
fi

export HS_PROFILE_OUT="build/prof/islamicstars_big_ship${SUFFIX}.log"
bash "$(dirname "$0")/profile_one.sh" IslamicStars profile "$SECONDS_ARG" "$WINDOW" \
  "-D HS_ISLAMICSTARS_PROFILE_SHAPE=14" \
  "-D HS_PROFILE_TRANS_SPEED=4" "$@"
