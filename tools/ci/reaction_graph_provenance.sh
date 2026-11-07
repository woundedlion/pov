#!/usr/bin/env bash
# Regenerate the reaction-graph table and diff it against the committed source.
set -euo pipefail

mkdir -p build/provenance
generated=build/provenance/reaction_graph.generated.cpp
python scripts/generate_reaction_graph.py -o "$generated"
if ! diff -u core/spatial/reaction_graph.cpp "$generated" > build/provenance/reaction_graph.diff; then
  head -40 build/provenance/reaction_graph.diff
  echo "::error::core/spatial/reaction_graph.cpp is out of sync with scripts/generate_reaction_graph.py"
  echo "Regenerate with: python scripts/generate_reaction_graph.py"
  exit 1
fi
echo "Regenerated table matches core/spatial/reaction_graph.cpp in full."
