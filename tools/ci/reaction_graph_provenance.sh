#!/usr/bin/env bash
# Regenerate the reaction-graph table and diff it against the committed source.
set -euo pipefail
# shellcheck source=tools/ci/common.sh
. "$(dirname "$0")/common.sh"

mkdir -p build/provenance
generated=build/provenance/reaction_graph.generated.cpp
python scripts/generate_reaction_graph.py -o "$generated"
require_regenerated core/spatial/reaction_graph.cpp "$generated" \
  scripts/generate_reaction_graph.py "python scripts/generate_reaction_graph.py"
