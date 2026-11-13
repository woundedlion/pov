# shellcheck shell=bash
# Values shared by the CI job scripts and the CI workflow; sourced, not run.

# Generated-table pins: bit-exact generator comparisons, run by the coverage
# job's Debug build and excluded from every other suite run.
export GENERATED_TESTS='^unit_(color_luts|gamut_lut|tangent_field_constants)$'

# Paths whose change runs the PCB suite.
export PCB_INPUTS='^(hardware/phantasm/|tools/ci/kicad_config\.py$|tools/run_python_tests\.py$|\.github/workflows/ci\.yml$)'

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
