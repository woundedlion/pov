/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Native soak for the ConwayMorph graph walk (docs/conway_morph_spec.md §7.8):
 * runs the real HankinSolids frame loop across enough legs to visit every
 * graph node, under native asserts. Surviving is most of the assertion — any
 * trap (seed reconciliation, per-frame hankin count guard, scratch overflow)
 * aborts the process. On top of that it pins:
 *   - full node coverage within a bounded leg count (deterministic RNG seed),
 *   - zero persistent-arena growth across leg compactions: the post-compaction
 *     offset is a pure function of (node, held seed), so every revisit must
 *     land on the byte-identical offset — steady state, not monotonic creep,
 *   - the render stays live (lit pixels sampled across the run).
 */
#pragma once

#include <cstdint>
#include <cstdio>

#include "core/engine/memory.h"
#include "core/mesh/conway_graph.h"
#include "effects/HankinSolids.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_soak_tests {

/**
 * @brief White-box accessor for HankinSolids' graph-walk state (befriended in
 *        effects/HankinSolids.h).
 * @details The soak needs the current node and held seed identity to pin
 *          coverage and the per-state post-compaction arena offset; neither is
 *          observable through the public effect surface.
 */
struct HankinWalkProbe {
  /**
   * @brief Current graph node (simple-registry index).
   */
  template <int W, int H> static int node(const HankinSolids<W, H> &fx) {
    return fx.node_;
  }
  /**
   * @brief Platonic solid the held seed mesh represents.
   */
  template <int W, int H>
  static int seed_identity(const HankinSolids<W, H> &fx) {
    return fx.seed_identity_;
  }
  /**
   * @brief In-flight leg's arrival data, or nullptr between legs.
   */
  template <int W, int H>
  static const Animation::ConwayMorph::Landing *
  pending_landing(const HankinSolids<W, H> &fx) {
    return fx.pending_landing_;
  }
  /**
   * @brief Face count of the current node's base mesh.
   */
  template <int W, int H>
  static size_t node_faces(const HankinSolids<W, H> &fx) {
    return fx.node_faces_;
  }
  /**
   * @brief Displayed palette per node base face (emission order).
   */
  template <int W, int H>
  static const uint8_t *node_face_palette(const HankinSolids<W, H> &fx) {
    return fx.node_face_palette_;
  }
  /**
   * @brief Live class-slot -> palette assignment the hankin cycle draws with.
   */
  template <int W, int H>
  static const std::array<int, HankinSolids<W, H>::NUM_PALETTES> &
  palette_idx(const HankinSolids<W, H> &fx) {
    return fx.palette_idx_;
  }
  /**
   * @brief Per-slot crossfade origin palette for the current hankin cycle.
   */
  template <int W, int H>
  static const std::array<int, HankinSolids<W, H>::NUM_PALETTES> &
  strap_from(const HankinSolids<W, H> &fx) {
    return fx.strap_from_;
  }
  /**
   * @brief Bitmask of slots crossfading over the current opening window.
   */
  template <int W, int H>
  static uint8_t strap_blend_mask(const HankinSolids<W, H> &fx) {
    return fx.strap_blend_mask_;
  }
  /**
   * @brief Sprite draws since the active hankin cycle's opening bookend.
   */
  template <int W, int H>
  static int hankin_cycle_frame(const HankinSolids<W, H> &fx) {
    return fx.hankin_cycle_frame_;
  }
  /**
   * @brief Strap-crossfade window length, in sprite frames.
   */
  template <int W, int H>
  static int strap_blend_frames(const HankinSolids<W, H> &) {
    return HankinSolids<W, H>::STRAP_BLEND_FRAMES;
  }
  /**
   * @brief The on-screen hankin mesh (topology carries the class slots).
   */
  template <int W, int H>
  static const MeshState &mesh(const HankinSolids<W, H> &fx) {
    return fx.mesh_;
  }
  /**
   * @brief Baked palette bank the effect shades with.
   */
  template <int W, int H>
  static const MeshPaletteBank &palette_bank(const HankinSolids<W, H> &fx) {
    return fx.palette_bank_;
  }
  /**
   * @brief Runs the production per-slot LUT resolution at a cycle frame.
   */
  template <int W, int H>
  static void resolve_slot_luts(
      HankinSolids<W, H> &fx, int cycle_frame,
      BakedPalette (&blended)[HankinSolids<W, H>::NUM_PALETTES],
      const BakedPalette *(&by_slot)[HankinSolids<W, H>::NUM_PALETTES],
      Arena &scratch) {
    fx.resolve_hankin_slot_luts(cycle_frame, blended, by_slot, scratch);
  }
};

/** Soak render size: small enough to keep the raster cheap, large enough that
 * every hankin face still covers pixels. */
constexpr int SOAK_W = 96;
constexpr int SOAK_H = 20;

/** Leg budget within which the seeded walk must have visited every node. */
constexpr int SOAK_LEG_BOUND = 96;

/** Extra legs run past full coverage so late-arriving leaf states also get
 * revisit (steady-state) checks. */
constexpr int SOAK_EXTRA_LEGS = 8;

/** Frame ceiling backstopping the leg bound (a leg is ~115-127 frames). */
constexpr int SOAK_FRAME_CAP = (SOAK_LEG_BOUND + SOAK_EXTRA_LEGS) * 140;

/**
 * @brief Runs the full-graph soak: real frame loop, every node visited, no
 *        traps, steady-state persistent arena.
 */
inline void test_full_graph_walk_soak() {
  reset_globals();

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

  // Post-compaction persistent offset per (node, held platonic seed); 0 =
  // state not yet seen.
  size_t post_offset[ConwayGraph::NUM_NODES][5] = {};

  int prev_node = HankinWalkProbe::node(fx);
  mark(prev_node);

  int legs = 0;
  int frames = 0;
  int legs_at_coverage = -1;
  uint64_t lit = 0;
  while (frames < SOAK_FRAME_CAP && legs < SOAK_LEG_BOUND) {
    fx.draw_frame();
    fx.advance_display();
    ++frames;
    if (frames % 16 == 0) {
      const Pixel &p = fx.get_pixel(SOAK_W / 2, SOAK_H / 2);
      lit += static_cast<uint64_t>(p.r) + p.g + p.b;
    }

    const int node = HankinWalkProbe::node(fx);
    if (node == prev_node)
      continue;

    // Leg completion: finish_morph_cycle compacted the persistent arena this
    // frame, so the offset now is the steady-state footprint of the arrived
    // (node, seed) pair and must reproduce exactly on every revisit.
    ++legs;
    prev_node = node;
    mark(node);

    const int sid = HankinWalkProbe::seed_identity(fx);
    HS_EXPECT_TRUE(ConwayGraph::is_platonic(sid));
    const size_t off = persistent_arena.get_offset();
    if (post_offset[node][sid] == 0) {
      post_offset[node][sid] = off;
    } else {
      if (off != post_offset[node][sid])
        std::printf("    [soak] persistent offset drift at '%s' (seed %d): "
                    "%zu -> %zu\n",
                    Solids::simple_registry[node].name, sid,
                    post_offset[node][sid], off);
      HS_EXPECT_EQ(off, post_offset[node][sid]);
    }

    if (visited_count == ConwayGraph::NUM_NODES) {
      if (legs_at_coverage < 0)
        legs_at_coverage = legs;
      if (legs >= legs_at_coverage + SOAK_EXTRA_LEGS)
        break;
    }
  }

  HS_EXPECT_EQ(visited_count, ConwayGraph::NUM_NODES);
  HS_EXPECT_GT(legs_at_coverage, 0);
  HS_EXPECT_LE(legs_at_coverage, SOAK_LEG_BOUND);
  HS_EXPECT_GT(lit, (uint64_t)0);

  std::printf(
      "  [soak] %d legs (%d frames) to full %d-node coverage; "
      "persistent hw=%zu scratch_a hw=%zu/%zu scratch_b hw=%zu/%zu "
      "center-sample sum=%llu\n",
      legs_at_coverage, frames, ConwayGraph::NUM_NODES,
      persistent_arena.get_high_water_mark(),
      scratch_arena_a.get_high_water_mark(), scratch_arena_a.get_capacity(),
      scratch_arena_b.get_high_water_mark(), scratch_arena_b.get_capacity(),
      static_cast<unsigned long long>(lit));
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs the ConwayMorph graph-walk soak.
 * @return The module's failure count.
 */
inline int run_conway_soak_tests() {
  hs_test::ModuleFixture fixture("conway_soak");
  test_full_graph_walk_soak();
  return fixture.result();
}

} // namespace conway_soak_tests
} // namespace hs_test
