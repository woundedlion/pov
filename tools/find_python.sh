#!/bin/sh

# App Execution Alias stubs can resolve on PATH without running Python.
hs_find_python() {
  for hs_python_candidate in "${HS_PYTHON:-}" python3 python; do
    if "$hs_python_candidate" --version >/dev/null 2>&1; then
      printf '%s\n' "$hs_python_candidate"
      return 0
    fi
  done
  return 1
}
