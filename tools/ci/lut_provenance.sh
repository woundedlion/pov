#!/usr/bin/env bash
# Regenerate the colour LUTs, diff the full header, and run the generator self-test.
set -euo pipefail

"${CLANG_FORMAT:-clang-format}" --version
mkdir -p build/provenance
generated=build/provenance/color_luts.generated.h
python scripts/generate_luts.py -o "$generated"
if ! diff -u core/color/color_luts.h "$generated" > build/provenance/color_luts.diff; then
  head -40 build/provenance/color_luts.diff
  echo "::error::core/color/color_luts.h is out of sync with scripts/generate_luts.py"
  echo "Regenerate with: python scripts/generate_luts.py -o core/color/color_luts.h"
  exit 1
fi
echo "Regenerated header matches core/color/color_luts.h in full."
python scripts/generate_luts.py --check
