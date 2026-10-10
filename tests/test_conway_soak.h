/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Native soak of the OpLeg graph walk through the real HankinSolids frame
 * loop.
 */
#pragma once

#include <cstdint>
#include <cstdio>
#include <map>
#include <set>

#include "core/memory.h"
#include "core/mesh/conway_graph.h"
#include "effects/HankinSolids.h"
#include "tests/conway_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_soak_tests {

/** Soak render size. */
constexpr int SOAK_W = 96;
constexpr int SOAK_H = 20;

static_assert(!ConwayGraph::is_platonic(-1));
static_assert(!ConwayGraph::is_platonic(
    ConwayGraph::dual_platonic(ConwayGraph::TRUNCATED_TETRAHEDRON)));

/** Leg budget within which the seeded walk must have visited every node. */
constexpr int SOAK_LEG_BOUND = 192;

/** Extra legs run past full coverage so late-arriving leaf states also get
 * revisit (steady-state) checks. */
constexpr int SOAK_EXTRA_LEGS = 8;

/** Frame ceiling backstopping the leg bound. */
constexpr int SOAK_FRAME_CAP = (SOAK_LEG_BOUND + SOAK_EXTRA_LEGS) * 140;

/** Per-leg minimum lit pixels in a sampled frame. */
constexpr int SOAK_MIN_LIT_PIXELS = SOAK_W * SOAK_H * 3 / 4;

/** Summed-channel floor per sampled frame (~1% of an all-white frame). */
constexpr uint64_t SOAK_MIN_FRAME_ENERGY = 4000000ull;

/**
 * @brief Runs the full-graph soak: real frame loop, every node visited, no
 *        traps, steady-state persistent arena.
 */
inline void test_full_graph_walk_soak(uint32_t seed) {
  reset_globals();
  hs::random().seed(seed);
  HS_CONTEXT("walk seed", seed);
  std::printf("  [soak] walk seed=%u\n", seed);

  HankinSolids<SOAK_W, SOAK_H> fx;
  fx.init();

  bool visited[ConwayGraph::NUM_NODES] = {};
  int visited_count = 0;
  const auto mark = [&](int node) {
    if (node >= 0 && node < ConwayGraph::NUM_NODES && !visited[node]) {
      visited[node] = true;
      ++visited_count;
    }
  };

  // Arrival state keys. A sweep arrival's work is fixed by (node, held seed).
  // A pass-through arrival also starts the next leg in the same frame, whose
  // seed fix depends on the seed held over the completed leg.
  std::map<uint32_t, size_t> post_offset;
  std::set<uint32_t> transition_seen;

  int prev_node = HankinWalkProbe::node(fx);
  mark(prev_node);

  int legs = 0;
  int frames = 0;
  int legs_at_coverage = -1;
  int frames_at_coverage = -1;
  uint64_t render_energy = 0;
  // Dimmest sampled frame of the leg in flight; the sentinel distinguishes a
  // leg that got no sample from one that rendered a black frame.
  int leg_min_lit = SOAK_W * SOAK_H;
  uint64_t leg_min_energy = UINT64_MAX;
  size_t leg_hw = persistent_arena.get_high_water_mark();
  size_t leg_scratch_a_hw = scratch_arena_a.get_high_water_mark();
  size_t leg_scratch_b_hw = scratch_arena_b.get_high_water_mark();
  while (frames < SOAK_FRAME_CAP && legs < SOAK_LEG_BOUND + SOAK_EXTRA_LEGS) {
    const int leg_sid = HankinWalkProbe::seed_identity(fx);
    fx.draw_frame();
    fx.advance_display();
    ++frames;
    if (frames % 16 == 0) {
      int lit = 0;
      uint64_t frame_energy = 0;
      for (int y = 0; y < SOAK_H; ++y)
        for (int x = 0; x < SOAK_W; ++x) {
          const Pixel &p = fx.get_pixel(x, y);
          const uint64_t channels = static_cast<uint64_t>(p.r) + p.g + p.b;
          frame_energy += channels;
          if (channels != 0)
            ++lit;
        }
      render_energy += frame_energy;
      if (lit < leg_min_lit)
        leg_min_lit = lit;
      if (frame_energy < leg_min_energy)
        leg_min_energy = frame_energy;
    }

    const int node = HankinWalkProbe::node(fx);
    if (node == prev_node)
      continue;

    // Leg completion: the post-arrival persistent offset must reproduce
    // exactly on every revisit of the arrival key.
    ++legs;
    HS_EXPECT_NE(leg_min_energy, UINT64_MAX);
    if (leg_min_energy != UINT64_MAX) {
      HS_EXPECT_GE(leg_min_lit, SOAK_MIN_LIT_PIXELS);
      HS_EXPECT_GE(leg_min_energy, SOAK_MIN_FRAME_ENERGY);
    }
    leg_min_lit = SOAK_W * SOAK_H;
    leg_min_energy = UINT64_MAX;
    const int departed_node = prev_node;
    prev_node = node;
    mark(node);
    const size_t before_hw = leg_hw;
    const size_t before_scratch_a_hw = leg_scratch_a_hw;
    const size_t before_scratch_b_hw = leg_scratch_b_hw;
    leg_hw = persistent_arena.get_high_water_mark();
    leg_scratch_a_hw = scratch_arena_a.get_high_water_mark();
    leg_scratch_b_hw = scratch_arena_b.get_high_water_mark();

    const int sid = HankinWalkProbe::seed_identity(fx);
    const bool sid_ok = ConwayGraph::is_platonic(sid);
    HS_EXPECT_TRUE(sid_ok);
    if (!sid_ok || node < 0 || node >= ConwayGraph::NUM_NODES ||
        departed_node < 0 || departed_node >= ConwayGraph::NUM_NODES)
      continue;

    const bool swept = node == HankinWalkProbe::dest(fx);
    const uint32_t start =
        swept ? 0u
              : 1u + static_cast<uint32_t>(HankinWalkProbe::cur_edge(fx)) * 8u +
                    static_cast<uint32_t>(leg_sid);
    const uint32_t arrival_key = (start * 32u + node) * 8u + sid;

    // A repeated directed transition with the same arrival key cannot raise
    // arena peaks.
    const uint32_t transition_key = arrival_key * 32u + departed_node;
    if (transition_seen.count(transition_key)) {
      HS_EXPECT_EQ(leg_hw, before_hw);
      HS_EXPECT_EQ(leg_scratch_a_hw, before_scratch_a_hw);
      HS_EXPECT_EQ(leg_scratch_b_hw, before_scratch_b_hw);
    }
    transition_seen.insert(transition_key);

    const size_t off = persistent_arena.get_offset();
    const auto [it, first] = post_offset.emplace(arrival_key, off);
    if (!first) {
      if (off != it->second)
        std::printf("    [soak] persistent offset drift at '%s' (seed %d, "
                    "start %u): %zu -> %zu\n",
                    Solids::simple_registry[node].name, sid, start, it->second,
                    off);
      HS_EXPECT_EQ(off, it->second);
    }

    if (visited_count == ConwayGraph::NUM_NODES) {
      if (legs_at_coverage < 0) {
        legs_at_coverage = legs;
        frames_at_coverage = frames;
      }
      if (legs >= legs_at_coverage + SOAK_EXTRA_LEGS)
        break;
    }
  }

  HS_EXPECT_EQ(visited_count, ConwayGraph::NUM_NODES);
  HS_EXPECT_GT(legs_at_coverage, 0);
  HS_EXPECT_LE(legs_at_coverage, SOAK_LEG_BOUND);

  // The host persistent arena is over-provisioned; gate on the device figure.
  using Fx = HankinSolids<SOAK_W, SOAK_H>;
  HS_EXPECT_LE(persistent_arena.get_high_water_mark(),
               Fx::DEVICE_PERSISTENT_BYTES);

  std::printf(
      "  [soak] %d legs (%d frames) to full %d-node coverage; "
      "persistent hw=%zu/%zu scratch_a hw=%zu/%zu scratch_b hw=%zu/%zu "
      "sampled frame energy=%llu\n",
      legs_at_coverage, frames_at_coverage, ConwayGraph::NUM_NODES,
      persistent_arena.get_high_water_mark(), Fx::DEVICE_PERSISTENT_BYTES,
      scratch_arena_a.get_high_water_mark(), scratch_arena_a.get_capacity(),
      scratch_arena_b.get_high_water_mark(), scratch_arena_b.get_capacity(),
      static_cast<unsigned long long>(render_energy));
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs the OpLeg graph-walk soak.
 * @return The module's failure count.
 */
inline int run_conway_soak_tests() {
  hs_test::ModuleFixture fixture("conway_soak");
  for (uint32_t seed : {1337u, 42u})
    test_full_graph_walk_soak(seed);
  return fixture.result();
}

} // namespace conway_soak_tests
} // namespace hs_test
