#!/usr/bin/env bash
# Cold-build every firmware target and fail on any first-party warning.
#
# usage: teensy_warnings.sh [teensy_warnings.py args...]
set -euo pipefail

bash tools/teensy_cold_build.sh teensy_build.log
python tools/teensy_warnings.py --build-log teensy_build.log "$@"
