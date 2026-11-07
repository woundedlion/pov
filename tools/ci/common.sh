# shellcheck shell=bash
# Values shared by the CI job scripts and the CI workflow; sourced, not run.

# Generated-table tests the sanitizer, optimized and Windows legs exclude.
export GENERATED_EXCLUDE='^unit_(color_luts|gamut_lut|mindsplatter_replay_form|tangent_field_constants)$'
