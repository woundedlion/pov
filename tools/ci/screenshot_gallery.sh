#!/usr/bin/env bash
# Check the capture script parses and the gallery and profile archive cover
# every registered effect.
set -euo pipefail

node --check scripts/capture_screenshots.mjs
node scripts/check_screenshots.mjs
node scripts/check_profiles.mjs
