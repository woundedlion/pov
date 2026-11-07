#!/usr/bin/env bash
# Check promoted shader documents and composed presets match their generators.
set -euo pipefail

node scripts/generate_promoted_shader_documents.mjs --check
node scripts/generate_composed_presets.mjs --check
