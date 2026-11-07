#!/usr/bin/env bash
# Fetch the pinned doxygen-awesome theme into .doxygen-awesome and append the
# theme settings to Doxyfile.
set -euo pipefail

theme_sha="$(python tools/build_pins.py doxygen-awesome)"
if [ ! -d .doxygen-awesome/.git ]; then
  git init -q .doxygen-awesome
  git -C .doxygen-awesome remote add origin https://github.com/jothepro/doxygen-awesome-css.git
fi
git -C .doxygen-awesome fetch -q --depth 1 origin "$theme_sha"
git -C .doxygen-awesome checkout -q FETCH_HEAD
cat docs/doxygen-theme.cfg >> Doxyfile
