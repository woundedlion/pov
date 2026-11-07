#!/usr/bin/env bash
# Validate tracked Markdown against a daydream checkout, the build pins, and
# that no tracked path matches .gitignore.
#
# usage: docs_markdown.sh <daydream-checkout>
set -euo pipefail

python tools/docs_check.py --checkout "daydream=${1:?usage: docs_markdown.sh <daydream-checkout>}"
python tools/build_pins.py --check
tracked_ignored=$(git ls-files -i -c --exclude-standard)
if [ -n "$tracked_ignored" ]; then
  echo "Tracked paths matched by .gitignore:" >&2
  echo "$tracked_ignored" >&2
  exit 1
fi
