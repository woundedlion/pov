#!/usr/bin/env bash
# Check every tracked working-tree file for trailing whitespace and blank EOF lines.
# usage: whitespace_gate.sh
set -eu

[ "$#" -eq 0 ] || { echo "usage: $0" >&2; exit 2; }
git -c core.whitespace=blank-at-eol,blank-at-eof diff --check "$(git hash-object -t tree /dev/null)"
