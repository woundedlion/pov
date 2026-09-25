#!/bin/bash
set -euo pipefail
# shellcheck source-path=SCRIPTDIR source=device_lock.sh
. "$(dirname "$0")/device_lock.sh"
# shellcheck source-path=SCRIPTDIR source=teensy_flash.sh
. "$(dirname "$0")/teensy_flash.sh"
[ "$#" = 1 ] || { echo "usage: $0 <PlatformIO environment>" >&2; exit 2; }
ENV=$1
TREE=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$TREE"
cleanup() {
  hs_device_release
  [ -z "$TREE_TOKEN" ] || _hs_break_lock "$TREE_LOCK" "$TREE_TOKEN" || :
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
hs_device_acquire upload "$ENV" 900
acquire_tree_lock
pio run -e "$ENV"
hs_teensy_flash "$ENV"
