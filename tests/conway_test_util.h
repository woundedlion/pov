/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/engine/memory.h"
#include "core/mesh/conway_graph.h"
#include "effects/HankinSolids.h"

namespace hs_test {
namespace conway_morph_tests {
/**
 * @brief Repartitions the global arena for the enclosing scope and restores the
 *        default split on the way out.
 * @details Declare it ahead of every Arena the scope carves out of the global
 * block, so the restore runs once those are gone. Without it the bespoke split
 * leaks into whatever case runs next.
 */
struct ScopedArenaSplit {
  /**
   * @brief Installs the split.
   * @param persistent Bytes for the persistent arena.
   * @param scratch_a Bytes for scratch arena A.
   * @param scratch_b Bytes for scratch arena B.
   */
  ScopedArenaSplit(size_t persistent, size_t scratch_a, size_t scratch_b) {
    configure_arenas(persistent, scratch_a, scratch_b);
  }
  ~ScopedArenaSplit() { configure_arenas_default(); }
  ScopedArenaSplit(const ScopedArenaSplit &) = delete;
  ScopedArenaSplit &operator=(const ScopedArenaSplit &) = delete;
};

/**
 * @brief Runs an edge's operator on a seed at one parameter point.
 * @param e Edge whose op kind is dispatched.
 * @param seed Seed mesh the op runs on.
 * @param target Arena receiving the output mesh.
 * @param temp Arena for the op's transient scratch.
 * @param t Operator parameter.
 * @param twist Snub twist (snub edges only).
 * @return The swept PolyMesh in `target`.
 */
inline PolyMesh run_edge_op(const ConwayGraph::EdgeSpec &e,
                            const PolyMesh &seed, Arena &target, Arena &temp,
                            float t, float twist) {
  switch (e.op) {
  case ConwayGraph::MorphOp::TRUNCATE:
    return MeshOps::truncate(seed, target, temp, t);
  case ConwayGraph::MorphOp::EXPAND:
    return MeshOps::expand(seed, target, temp, t);
  case ConwayGraph::MorphOp::SNUB:
    return MeshOps::snub(seed, target, temp, t, twist);
  case ConwayGraph::MorphOp::CHAMFER:
    return MeshOps::chamfer(seed, target, temp, t);
  }
  return PolyMesh{};
}

} // namespace conway_morph_tests
namespace conway_soak_tests {
/**
 * @brief White-box accessor for HankinSolids' graph-walk state (befriended in
 *        effects/HankinSolids.h).
 * @details The soak needs the current node and held seed identity to pin
 *          coverage and the per-state post-compaction arena offset; neither is
 *          observable through the public effect surface.
 */
struct HankinWalkProbe {
  /** @brief Builds a clean endpoint mesh from the held seed. */
  static PolyMesh node_mesh_at(const ConwayGraph::EdgeSpec &edge, bool to_end,
                               const PolyMesh &seed, Arena &a, Arena &b) {
    return HankinSolids<96, 20>::node_mesh_at(seed, edge, to_end, a, b);
  }

  /**
   * @brief Current graph node (simple-registry index).
   */
  template <int W, int H> static int node(const HankinSolids<W, H> &fx) {
    return fx.node;
  }
  /**
   * @brief Platonic solid the held seed mesh represents.
   */
  template <int W, int H>
  static int seed_identity(const HankinSolids<W, H> &fx) {
    return fx.seed_identity;
  }
  /**
   * @brief In-flight leg's arrival data, or nullptr between legs.
   */
  template <int W, int H>
  static const Animation::OpLeg::Landing *
  pending_landing(const HankinSolids<W, H> &fx) {
    return fx.pending_landing;
  }
  /**
   * @brief Face count of the current node's base mesh.
   */
  template <int W, int H>
  static size_t node_faces(const HankinSolids<W, H> &fx) {
    return fx.node_faces;
  }
  /**
   * @brief Displayed palette per node base face (emission order).
   */
  template <int W, int H>
  static const uint8_t *node_face_palette(const HankinSolids<W, H> &fx) {
    return fx.node_face_palette;
  }
  /**
   * @brief Per star face, the palette of the rosettes hosted inside it — the
   * ramp it cross-fades onto as it closes at the sweep's midpoint.
   */
  template <int W, int H>
  static const uint8_t *star_rim_palette(const HankinSolids<W, H> &fx) {
    return fx.star_rim_palette;
  }
  /**
   * @brief Live class-slot -> palette assignment the hankin cycle draws with.
   */
  template <int W, int H>
  static const std::array<int, HankinSolids<W, H>::NUM_PALETTES> &
  palette_idx(const HankinSolids<W, H> &fx) {
    return fx.palette_idx;
  }
  /**
   * @brief Per-slot crossfade origin palette for the current hankin cycle.
   */
  template <int W, int H>
  static const std::array<int, HankinSolids<W, H>::NUM_PALETTES> &
  strap_from(const HankinSolids<W, H> &fx) {
    return fx.strap_from;
  }
  /**
   * @brief Bitmask of slots crossfading over the current opening window.
   */
  template <int W, int H>
  static uint8_t strap_blend_mask(const HankinSolids<W, H> &fx) {
    return fx.strap_blend_mask;
  }
  /**
   * @brief Sprite draws since the active hankin cycle's opening bookend.
   */
  template <int W, int H>
  static int hankin_cycle_frame(const HankinSolids<W, H> &fx) {
    return fx.hankin_cycle_frame;
  }
  /**
   * @brief Strap-crossfade window length, in sprite frames.
   */
  template <int W, int H>
  static int strap_blend_frames(const HankinSolids<W, H> &) {
    return HankinSolids<W, H>::STRAP_BLEND_FRAMES;
  }
  /** @brief Shape blend weights at a cycle frame. */
  template <int W, int H>
  static auto shape_weights(const HankinSolids<W, H> &fx, int cycle_frame) {
    return fx.shape_weights(cycle_frame);
  }
  /** @brief Interlace sweep angle at a cycle frame. */
  template <int W, int H>
  static float sweep_angle(const HankinSolids<W, H> &fx, int cycle_frame) {
    return fx.sweep_wave()(math::ease_linear(static_cast<float>(cycle_frame) /
                                             fx.HANKIN_SWEEP_FRAMES));
  }
  /** @brief Interlace-angle sweep length, in sprite frames. */
  template <int W, int H> static int sweep_frames(const HankinSolids<W, H> &) {
    return HankinSolids<W, H>::HANKIN_SWEEP_FRAMES;
  }
  /**
   * @brief The on-screen hankin mesh (topology carries the class slots).
   */
  template <int W, int H>
  static const MeshState &mesh(const HankinSolids<W, H> &fx) {
    return fx.hankin_mesh;
  }
  /**
   * @brief Baked palette bank the effect shades with.
   */
  template <int W, int H>
  static const MeshPaletteBank &palette_bank(const HankinSolids<W, H> &fx) {
    return fx.palette_bank;
  }
  /**
   * @brief Runs the production per-slot LUT resolution at a cycle frame.
   */
  template <int W, int H>
  static void resolve_slot_luts(
      HankinSolids<W, H> &fx, int cycle_frame,
      BakedPalette (&blended)[HankinSolids<W, H>::NUM_PALETTES],
      const BakedPalette *(&star_by_slot)[HankinSolids<W, H>::NUM_PALETTES],
      const BakedPalette *(&strap_by_slot)[HankinSolids<W, H>::NUM_PALETTES],
      Arena &scratch) {
    fx.resolve_hankin_slot_luts(cycle_frame, blended, star_by_slot,
                                strap_by_slot, scratch);
  }
  /**
   * @brief Rebuilds the hankin mesh at `angle` and renders it through the
   * production draw path at the current (unadvanced) orientation, with the
   * given strap opening-fade weight.
   * @details Rendering without stepping the timeline keeps the camera fixed
   * across successive calls, so a diff isolates the mesh/fade change. Used to
   * pin the strap opening-fade: the newborn straps at a small angle must, when
   * faded, leave the angle-0 bookend nearly unchanged.
   */
  template <int W, int H>
  static void render_at_angle(HankinSolids<W, H> &fx, Canvas &canvas,
                              float angle, int cycle_frame, float strap_fade,
                              float close_blend, float terminal_fade,
                              float star_close, Arena &scratch) {
    MeshOps::update_hankin(fx.compiled_hankin, fx.hankin_mesh, persistent_arena,
                           angle);
    BakedPalette blended[HankinSolids<W, H>::NUM_PALETTES];
    const BakedPalette *star_by_slot[HankinSolids<W, H>::NUM_PALETTES];
    const BakedPalette *strap_by_slot[HankinSolids<W, H>::NUM_PALETTES];
    fx.resolve_hankin_slot_luts(cycle_frame, blended, star_by_slot,
                                strap_by_slot, scratch);
    fx.draw_mesh(canvas, fx.hankin_mesh, star_by_slot, strap_by_slot,
                 strap_fade, close_blend, terminal_fade, star_close);
  }
};

} // namespace conway_soak_tests
} // namespace hs_test
