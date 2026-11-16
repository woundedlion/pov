/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Shared arenas and recipe timing access for Conway morph cases. */
#pragma once

#include "core/animation/recipe_build.h"
#include "core/memory.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/mesh.h"
#include "tests/conway_test_util.h"

#include <cstddef>
#include <cstdint>

namespace hs_test::conway_morph_tests {

inline uint8_t morph_target_buf[256 * 1024]; /**< Op output arena. */
inline uint8_t morph_temp_buf[256 * 1024];   /**< Op scratch arena. */
inline uint8_t morph_aux_buf[256 * 1024];    /**< Seed / second-result arena. */
inline uint8_t morph_persist_buf[64 * 1024]; /**< Persistent-seed arena. */
inline uint8_t morph_bank_buf[64 * 1024];    /**< Baked palette LUT arena. */

/**
 * @brief Builds a site's persistent seed through the auxiliary arena halves.
 * @tparam Site Site providing a seed builder.
 * @param site Site whose seed is built.
 * @param persist Arena owning the returned seed.
 * @return Seed mesh owned by persist.
 */
template <typename Site>
inline PolyMesh finalize_site_seed(const Site &site, Arena &persist) {
  constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
  Arena a(morph_aux_buf, HALF);
  Arena b(morph_aux_buf + HALF, HALF);
  return hs_test::finalize_site_seed(site, persist, a, b);
}

/** @brief Production recipe timings exposed to replay fixtures. */
struct RecipeLegLengths : Animation::RecipeBuild<RecipeLegLengths, 1, 1> {
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::HANKIN_LEG_FRAMES;
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::SWEEP_LEG_FRAMES;
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::RELAX_LEG_FRAMES;
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::clean_endpoint;
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::make_build_leg;
};

using ConwayGraph::T_EPS;

} // namespace hs_test::conway_morph_tests
