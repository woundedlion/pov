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

rm -rf .pio/build_cache .pio/build
pio run -v 2>&1 | tee "$log"
