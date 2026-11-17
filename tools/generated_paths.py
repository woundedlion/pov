"""Generated sources excluded from formatting and API comment checks."""

import re


GENERATED_PATHS = (
    "core/color/color_luts.h",
    "core/color/gamut_lut.h",
    "core/color/srgb_decode_lut.h",
    "effects/generated/",
    "core/mesh/relax_bakes_generated.h",
    "core/spatial/reaction_graph.cpp",
)


def generated_path_pattern(prefix: str) -> str:
    """Make a regex alternative per generated path, rooted by `prefix`."""
    return "|".join(prefix + re.escape(path) + ("" if path.endswith("/") else "$")
                    for path in GENERATED_PATHS)
