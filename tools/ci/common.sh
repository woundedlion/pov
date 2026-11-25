# shellcheck shell=bash
# Values and helpers shared by the CI job scripts and the CI workflow; sourced,
# not run.

# Generated-table pins: bit-exact generator comparisons, run by the coverage
# job's Debug build and excluded from every other suite run.
export GENERATED_TESTS='^unit_(color_luts|gamut_lut|tangent_field_constants)$'

# Paths whose change runs the PCB suite.
export PCB_INPUTS='^(hardware/phantasm/|hardware/pov_segmented\.h$|tools/ci/kicad_config\.py$|tools/run_python_tests\.py$|\.github/workflows/ci\.yml$)'

# Top-level job keys in a GitHub Actions workflow.
ci_job_keys() {
  awk '
    /^jobs:[[:space:]]*(#.*)?$/ { injobs = 1; next }
    injobs && /^[^[:space:]]/ { injobs = 0 }
    injobs && /^  [A-Za-z0-9_-]+:[[:space:]]*(#.*)?$/ {
      sub(/^  /, ""); sub(/:[[:space:]]*(#.*)?$/, ""); print
    }
  ' "$1"
}

# Fail unless a regenerated file matches its committed copy in full.
# usage: require_regenerated <committed> <generated> <generator> <regenerate-command>
require_regenerated() {
  local committed=$1 generated=$2 generator=$3 regenerate=$4
  local name=${committed##*/}
  local report="build/provenance/${name%.*}.diff"
  if ! diff -u "$committed" "$generated" > "$report"; then
    head -40 "$report"
    echo "::error::$committed is out of sync with $generator"
    echo "Regenerate with: $regenerate"
    exit 1
  fi
  echo "Regenerated output matches $committed in full."
}
