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
trap hs_device_release EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
hs_device_acquire upload "$ENV" 900
pio run -e "$ENV"
hs_teensy_flash "$ENV"
