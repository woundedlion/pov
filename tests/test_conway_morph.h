/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Conway morph operators, OpLeg transitions, recipe chains, and walk policies.
 */
#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <vector>
#include "core/color/palettes.h"
#include "core/mesh/conway.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/hankin.h"
#include "core/mesh/recipe.h"
#include "core/mesh/solids.h"
#include "core/render/canvas.h"
#include "effects/HankinSolids.h"
#include "effects/IslamicStars.h"
#include "tests/conway_test_util.h"
#include "tests/mesh_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/test_conway.h" // check_euler_characteristic_two, face_type_histogram
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_morph_tests {

template <typename Mesh>
concept BorrowableSweepSeed = requires(Mesh &&mesh) {
  Animation::OpLeg::SweepSeed::borrow(std::forward<Mesh>(mesh));
};

static_assert(std::is_constructible_v<Animation::OpLeg::SweepSeed, PolyMesh &>);
static_assert(
    !std::is_constructible_v<Animation::OpLeg::SweepSeed, PolyMesh &&>);
static_assert(
    !std::is_constructible_v<Animation::OpLeg::SweepSeed, const PolyMesh &&>);
static_assert(BorrowableSweepSeed<PolyMesh &>);
static_assert(BorrowableSweepSeed<const PolyMesh &>);
static_assert(!BorrowableSweepSeed<PolyMesh>);
static_assert(!BorrowableSweepSeed<const PolyMesh>);

inline uint8_t morph_target_buf[256 * 1024]; /**< Op output arena. */
inline uint8_t morph_temp_buf[256 * 1024];   /**< Op scratch arena. */
inline uint8_t morph_aux_buf[256 * 1024];    /**< Seed / second-result arena. */
inline uint8_t morph_persist_buf[64 * 1024]; /**< Persistent-seed arena. */
inline uint8_t morph_bank_buf[64 * 1024];    /**< Baked palette LUT arena. */

/** @brief Production recipe timings exposed to replay fixtures. */
struct RecipeLegLengths : Animation::RecipeBuild<RecipeLegLengths, 1, 1> {
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::HANKIN_LEG_FRAMES;
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::SWEEP_LEG_FRAMES;
  using Animation::RecipeBuild<RecipeLegLengths, 1, 1>::RELAX_LEG_FRAMES;
};

using ConwayGraph::T_EPS;

#include "tests/conway_morph/seeds.h"
#include "tests/conway_morph/oracles.h"
#include "tests/conway_morph/settle.h"
#include "tests/conway_morph/bridge.h"
#include "tests/conway_morph/jitterbug.h"
#include "tests/conway_morph/clean_swap.h"
#include "tests/conway_morph/endpoints.h"
#include "tests/conway_morph/topology.h"
#include "tests/conway_morph/scratch_budget.h"
#include "tests/conway_morph/walk.h"
#include "tests/conway_morph/profile_tour.h"
#include "tests/conway_morph/ambo_hankin.h"
#include "tests/conway_morph/hankin.h"
#include "tests/conway_morph/hankin_stability.h"
#include "tests/conway_morph/recipe_steps.h"
#include "tests/conway_morph/medial.h"
#include "tests/conway_morph/recipe_replay.h"
#include "tests/conway_morph/reconcile.h"
/**
 * @brief Runs all OpLeg operator-level tests.
 * @return The module's failure count.
 */
inline int run_conway_morph_tests() {
  hs_test::ModuleFixture fixture("conway_morph");

  test_edge_endpoints_match_registry();
  test_edge_sweeps_hold_topology();
  test_edge_morph_frames_fit_scratch_budget();

  test_relax_is_vertex_order_identity();

  test_snub_tetrahedron_relax_converges_to_icosahedron();
  test_ambo_tetrahedron_is_regular_octahedron();

  test_jitterbug_icosa_point_is_regular();
  test_jitterbug_octa_end_covers_octahedron();
  test_jitterbug_sweep_holds_topology();

  test_truncate_near_half_merges_onto_ambo();
  test_ops_at_t_eps_primary_faces_match_seed();

  test_ambo_leg_on_hankin_seed_holds_topology();
  test_hankin_sweep_on_islamic_seeds_holds_topology();
  test_hankin_sweep_vertex_stability();
  test_opleg_hankin_sweep_smoke();

  test_truncate_leg_on_recipe_seeds_holds_topology();
  test_snub_leg_on_recipe_seeds_holds_topology();
  test_relax_leg_on_recipe_seeds_holds_topology();
  test_medial_dual_bridge_wellformed();
  test_opleg_medial_leg_smoke();
  test_opleg_dual_bridge_seam_correspondence();
  test_opleg_step_leg_smoke();
  test_opleg_step_paused_holds_frame();
  test_opleg_step_leg_overshooting_easing();
  test_opleg_gated_swap_smoke();
  test_opleg_edge_leg_crossfade();
  test_unsweepable_recipe_steps_are_gated();

  test_relax_source_hash_separates_bevel_inputs();
  test_recipe_chain_build_replay();
  test_reconcile_bijection_wellposed();

  test_walk_policy_coverage_and_balance();
  test_ordered_tour_full_coverage_and_wrap();

  return fixture.result();
}

} // namespace conway_morph_tests
} // namespace hs_test
