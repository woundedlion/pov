#!/usr/bin/env bash
# Regenerate the colour LUTs, diff the full header, and run the generator self-test.
set -euo pipefail
# shellcheck source=tools/ci/common.sh
. "$(dirname "$0")/common.sh"

"${CLANG_FORMAT:-clang-format}" --version
mkdir -p build/provenance
generated=build/provenance/color_luts.generated.h
python scripts/generate_luts.py -o "$generated"
require_regenerated core/color/color_luts.h "$generated" scripts/generate_luts.py \
  "python scripts/generate_luts.py -o core/color/color_luts.h"
python scripts/generate_luts.py --check
