#!/usr/bin/env bash
set -eu

git -c core.whitespace=blank-at-eol,blank-at-eof diff --check "$(git hash-object -t tree /dev/null)" HEAD
