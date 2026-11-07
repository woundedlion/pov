#!/usr/bin/env bash
# Run every scripts/*.test.mjs suite, then the browser capture helpers.
#
# usage: scripts_tests.sh [playwright install args...]   e.g. --with-deps
set -euo pipefail

bash tools/require_test_files.sh 'scripts/*.test.mjs'
npm test
npm ci --no-audit --no-fund
npx playwright install "$@" chromium
HS_BROWSER_TESTS=1 node --test scripts/screenshot_resolution.test.mjs scripts/png_probe.test.mjs
