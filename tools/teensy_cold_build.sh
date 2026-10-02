#!/usr/bin/env bash
# Cold `pio run -v` over every platformio.ini environment, teed to
# teensy_build.log for tools/teensy_warnings.py.
#
# A cached translation unit emits no warnings, so the object cache and
# .pio/build go first and the whole firmware tree recompiles; budget tens of
# minutes. pipefail so a failed build is not masked by tee's status.
#
# usage: teensy_cold_build.sh [log]
set -euo pipefail

if [ "$#" -gt 1 ]; then
  echo "usage: $0 [log]" >&2
  exit 2
fi
log=${1:-teensy_build.log}
case "$log" in
  /*|[A-Za-z]:[\\/]*) ;;
  *) log="$PWD/$log" ;;
esac
# shellcheck source-path=SCRIPTDIR source=device_lock.sh
. "$(dirname "$0")/device_lock.sh"
TREE=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$TREE"
cleanup() {
  [ -z "$TREE_TOKEN" ] || _hs_break_lock "$TREE_LOCK" "$TREE_TOKEN" || :
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
acquire_tree_lock

export LC_ALL=C
rm -rf .pio/build_cache .pio/build
pio run -v 2>&1 | tee "$log"
