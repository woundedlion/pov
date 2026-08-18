# Pin the in-tree patches applied to core/vendor/FastNoiseLite.h against the
# patch record in core/vendor/FastNoiseLite_config.h. Dropping any of them
# compiles cleanly and fails silently: without the config include the
# FASTNOISELITE_ONLY_OPENSIMPLEX2 guards go inert and the full noise-type switch
# comes back (flash regression), and without HS_O3_FN the per-pixel sampler
# loses selective -O3 (perf regression). The upstream version is pinned too, so
# a bump has to re-apply the patches and re-record the version here.
# -D args: HEADER (FastNoiseLite.h), CONFIG (FastNoiseLite_config.h).

# Script mode inherits no policies from the project, so every policy would
# otherwise default to OLD. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

set(VENDORED_VERSION "1.1.1")

file(READ "${HEADER}" _header)
file(READ "${CONFIG}" _config)

set(_missing "")

if(NOT _header MATCHES "VERSION: ${VENDORED_VERSION}")
  list(APPEND _missing
    "FastNoiseLite.h no longer declares upstream VERSION: ${VENDORED_VERSION}")
endif()

if(NOT _header MATCHES "#include \"FastNoiseLite_config.h\"")
  list(APPEND _missing
    "FastNoiseLite.h no longer includes FastNoiseLite_config.h")
endif()

string(REGEX MATCHALL "#ifn?def FASTNOISELITE_ONLY_OPENSIMPLEX2"
  _guards "${_header}")
list(LENGTH _guards _guard_count)
if(NOT _guard_count EQUAL 4)
  list(APPEND _missing
    "FastNoiseLite.h has ${_guard_count} FASTNOISELITE_ONLY_OPENSIMPLEX2 guards, expected 4")
endif()

if(NOT _header MATCHES "HS_O3_FN float SingleOpenSimplex2\\(")
  list(APPEND _missing "SingleOpenSimplex2 has lost its HS_O3_FN")
endif()

if(NOT _header MATCHES "float GetNoiseSingle\\(" OR
   NOT _header MATCHES "float GetNoiseSingleTransformed\\(")
  list(APPEND _missing "FastNoiseLite has lost its raw-octave sampling paths")
endif()

if(NOT _header MATCHES "HS_HOT_FLASH_MEMBER void GetVectorNoiseSingle\\(" OR
   NOT _header MATCHES "HS_HOT_FLASH_MEMBER void[\r\n ]+SingleDomainWarpOpenSimplex2Gradient\\(")
  list(APPEND _missing "FastNoiseLite has lost its raw vector-noise path")
endif()

if(NOT _header MATCHES "void GetNoiseGradientSingle\\(" OR
   NOT _header MATCHES "void SingleOpenSimplex2Gradient\\(")
  list(APPEND _missing "FastNoiseLite has lost its analytic gradient path")
endif()

if(NOT _config MATCHES "#define FASTNOISELITE_ONLY_OPENSIMPLEX2")
  list(APPEND _missing
    "FastNoiseLite_config.h no longer defines FASTNOISELITE_ONLY_OPENSIMPLEX2")
endif()

if(NOT _config MATCHES "#include \"platform/platform\\.h\"")
  list(APPEND _missing
    "FastNoiseLite_config.h no longer includes platform/platform.h (HS_O3_FN)")
endif()

if(_missing)
  string(REPLACE ";" "\n  - " _report "${_missing}")
  message(FATAL_ERROR
    "Vendored FastNoiseLite patches are missing:\n  - ${_report}\n"
    "Re-apply the patches listed in core/vendor/FastNoiseLite_config.h and "
    "update VENDORED_VERSION in this script if the vendored version changed.")
endif()

message(STATUS
  "FastNoiseLite ${VENDORED_VERSION}: all 7 in-tree patches present")
