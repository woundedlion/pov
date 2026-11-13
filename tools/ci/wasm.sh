#!/usr/bin/env bash
# Build the release, debug and strict-FP WASM configurations, gate first-party
# warnings, smoke-test the modules, and exercise the daydream install and
# provenance path. Needs an activated emsdk.
#
# usage: wasm.sh <daydream-install-dir>   a scratch directory the install fills
set -euo pipefail

dd="${1:?usage: wasm.sh <daydream-install-dir>}"
scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT

# DAYDREAM_DIR must exist at configure time and contain the install marker.
mkdir -p "$dd/src/app" "$dd/docs/screenshots" "$dd/generated/shader/patterns"
touch "$dd/src/app/daydream.js" "$dd/docs/screenshots/stale.png" \
  "$dd/docs/screenshots/unrelated.txt" "$dd/generated/shader/patterns/stale.shader.json"
mkdir -p build

cmake --preset wasm-release -DDAYDREAM_DIR="$dd"
cmake --build --preset wasm-release 2>&1 | tee build/wasm-release-build.log
cmake --build --preset wasm-release --target pullback_capture_wasm 2>&1 \
  | tee -a build/wasm-release-build.log

python tools/pullback_capture.py --manifest-dir tests/data/pullback \
  --operations-output "$scratch/ops.bin"
for resolution in 96x20 288x144; do
  node build/wasm-release/pullback_capture_wasm.cjs --resolution "$resolution" \
    --operations "$scratch/ops.bin" "$scratch/wasm-release-$resolution.bin"
done
python tools/pullback_capture.py --manifest-dir tests/data/pullback \
  --configuration wasm-release --backend-audit "$scratch"/wasm-release-*.bin

# Compile coverage for debug-only code and developer bindings.
cmake --preset wasm-debug -DDAYDREAM_DIR="$dd" -DHS_WASM_DEV_BINDINGS=ON
cmake --build --preset wasm-debug 2>&1 | tee build/wasm-debug-build.log

# Compile coverage for the strict-FP capture configuration.
cmake --preset wasm-strict-fp -DDAYDREAM_DIR="$dd"
cmake --build --preset wasm-strict-fp 2>&1 | tee build/wasm-strict-fp-build.log
cmake --build --preset wasm-strict-fp --target pullback_capture_wasm 2>&1 \
  | tee -a build/wasm-strict-fp-build.log

# Reject first-party warnings in all three configurations; sysroot warnings are excluded.
logs=(build/wasm-release-build.log build/wasm-debug-build.log build/wasm-strict-fp-build.log)
# grep exits 2 for a missing log, which would otherwise look like no warnings.
for log in "${logs[@]}"; do
  if [ ! -f "$log" ]; then
    echo "::error::$log missing — the warning capture broke, so the gate proves nothing"
    exit 1
  fi
done
if grep -nE '(^|[[:space:]/])(core|effects|workbench|hardware|targets|tools|tests)/[^:[:space:]]+:[0-9]+:([0-9]+:)?[[:space:]]*warning:' "${logs[@]}"; then
  echo "::error::first-party warning(s) in the WASM build — fix them (this target has no baseline)"
  exit 1
fi
if grep -nE '(^|[[:space:]])(<command line>(:[0-9]+){0,2}|em\+\+|emcc|wasm-ld|clang\+\+|clang)(\.exe)?:[[:space:]]*warning:' "${logs[@]}"; then
  echo "::error::fileless toolchain warning(s) in the WASM build — fix them (this target has no baseline)"
  exit 1
fi
echo "No first-party or fileless toolchain warnings in the release, debug or strict-FP WASM build."

# Debug frames have a 6 KB ceiling; optimized frames use 2 KB.
WASM_SMOKE_STACK_CEILING=6144 node scripts/wasm_smoke.mjs build/wasm-debug/holosphere_wasm.js
node scripts/wasm_smoke.mjs build/wasm-release/holosphere_wasm.js
# Engine payloads match targets/wasm/engine_types.d.ts and the emitted aliases;
# the debug module adds the developer bindings.
node scripts/engine_payload_types.mjs build/wasm-debug/holosphere_wasm.js
node scripts/engine_payload_types.mjs build/wasm-release/holosphere_wasm.js

# Verify the installed provenance set and both directory installs.
cmake --build --preset wasm-release-install
for f in holosphere_wasm.wasm holosphere_wasm.js holosphere_wasm.d.ts holosphere_bindings.d.ts engine_types.d.ts engine_contract.mjs segment_geometry.js segment_geometry.d.ts holosphere_wasm.sha holosphere_wasm.wasm.sha256 holosphere_wasm.toolchain pov_segment_map.json shader/composed_effect_roster.mjs shader/shader_workbench.mjs shader/sha256.mjs shader/engine_catalog.json; do
  if [ ! -f "$dd/generated/$f" ]; then
    echo "::error::install did not produce $f — a consumer in the daydream checkout depends on it"
    exit 1
  fi
done
if [ ! -f "$dd/README.md" ]; then
  echo "::error::install did not produce README.md"
  exit 1
fi
if ! cmp -s README.md "$dd/README.md"; then
  echo "::error::installed README.md differs from the engine source"
  exit 1
fi
if [ -e "$dd/docs/screenshots/stale.png" ]; then
  echo "::error::install retained a stale mirrored screenshot"
  exit 1
fi
if [ ! -e "$dd/docs/screenshots/unrelated.txt" ]; then
  echo "::error::install removed unrelated screenshot-directory content"
  exit 1
fi
diff -u \
  <(cd docs/screenshots && find . -type f -name '*.png' | sort) \
  <(cd "$dd/docs/screenshots" && find . -type f -name '*.png' | sort)
# Exercise installed glue after file(CONFIGURE) rewriting.
node scripts/wasm_smoke.mjs "$dd/generated/holosphere_wasm.js"
# A clean checkout must record a bare 40-hex SHA (no -dirty marker).
if grep -q -- '-dirty' "$dd/generated/holosphere_wasm.sha"; then
  echo "::error::engine .sha recorded as -dirty on a clean checkout"
  cat "$dd/generated/holosphere_wasm.sha"
  exit 1
fi
( cd "$dd/generated" && sha256sum -c holosphere_wasm.wasm.sha256 )
for workspace in "$PWD" "$(pwd -W 2>/dev/null || true)"; do
  if [ -n "$workspace" ] && grep -aqF "$workspace" "$dd/generated/holosphere_wasm.wasm"; then
    echo "::error::release WASM embeds its checkout path and cannot be reproduced in another workspace"
    exit 1
  fi
done
# Compare the toolchain marker against the compiler that built the artifact.
recorded="$(awk '$1 == "emsdk" { print $2 }' "$dd/generated/holosphere_wasm.toolchain")"
emcc="$EMSDK/upstream/emscripten/emcc"
[ -e "$emcc" ] || emcc="$emcc.bat"
built_with="$("$emcc" -dumpversion | tr -d '\r')"
if [ "$recorded" != "$built_with" ]; then
  echo "::error::holosphere_wasm.toolchain records emsdk '$recorded', but this job built with $built_with"
  cat "$dd/generated/holosphere_wasm.toolchain"
  exit 1
fi
# The shipped configuration is wasm-release with developer bindings off.
field() { awk -v k="$1" '$1 == k { print $2 }' "$dd/generated/holosphere_wasm.toolchain"; }
if [ "$(field build_type)" != "Release" ] || [ "$(field dev_bindings)" != "OFF" ]; then
  echo "::error::holosphere_wasm.toolchain records build_type '$(field build_type)' / dev_bindings '$(field dev_bindings)'; the release install must record Release / OFF"
  cat "$dd/generated/holosphere_wasm.toolchain"
  exit 1
fi
echo "Provenance set produced, artifact hashes verified, toolchain recorded:"
cat "$dd/generated/holosphere_wasm.toolchain"

# Mirrored shader documents match their sources.
cmp scripts/composed_effect_roster.mjs "$dd/generated/shader/composed_effect_roster.mjs"
( cd patterns && find . -type f \( -name '*.shader.json' -o -name catalog.json \) ! -name example.shader.json ) | sort > "$scratch/mirror-source.txt"
test ! -e "$dd/generated/shader/patterns/stale.shader.json"
( cd "$dd/generated/shader/patterns" && find . -maxdepth 1 -type f \( -name '*.shader.json' -o -name catalog.json \) ) | sort > "$scratch/mirror-installed.txt"
test -s "$scratch/mirror-source.txt"
test -s "$scratch/mirror-installed.txt"
diff -u "$scratch/mirror-source.txt" "$scratch/mirror-installed.txt"
