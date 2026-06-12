/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

/**
 * @brief CRTP base for the spherical reaction-diffusion effects (BZ, Gray-Scott).
 * @tparam Derived Concrete effect supplying the system-specific render().
 * @tparam W Framebuffer width in pixels.
 * @tparam H Framebuffer height in pixels.
 * @details Both systems run on the same 7680-node Fibonacci-lattice K-NN graph,
 * share a Languid random-walk view orientation, build the cached node positions
 * once at init (the lattice is static), and interpolate with the same Wendland
 * C2 kernel. This base captures exactly that shared scaffolding. The physics
 * (Lotka-Volterra 3-species Q8 vs Gray-Scott 2-species Q16), state
 * representation, seeding, params, palette, and rendering are fundamentally
 * different and stay in the derived classes.
 *
 * Dispatch is static (CRTP): draw_frame() forwards to Derived::render() with no
 * virtual call, so the abstraction is zero-overhead. Derived classes befriend
 * this base so render() can stay private.
 */
template <typename Derived, int W, int H>
class ReactionDiffusionBase : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N; /**< Lattice node count. */
  static constexpr int RD_K = ReactionGraph::RD_K; /**< K-NN neighbors per node. */

  /**
   * @brief Constructs the base, disabling pixel persistence.
   * @details Each frame is fully repainted from lattice state, so the
   * framebuffer need not persist between frames.
   */
  ReactionDiffusionBase() : Effect(W, H) { persist_pixels = false; }

  /**
   * @brief Reports whether the effect composites over the framebuffer background.
   * @return Always true; both systems composite over the background.
   */
  bool show_bg() const override { return true; }

  /**
   * @brief Advances one animation frame and dispatches to the derived renderer.
   * @details Advances the orientation timeline, then statically dispatches to
   * Derived::render() for the system-specific physics and drawing.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    static_cast<Derived *>(this)->render(canvas);
  }

protected:
  // Wendland C2 compact kernel: w(d) = max(0, 1 - d²/R²)²
  static constexpr float D_AVG = ReactionGraph::D_AVG; /**< Mean inter-node spacing, sqrt(4π / RD_N). */
  static constexpr float KERNEL_R = 1.5f * D_AVG;      /**< Kernel support radius. */
  static constexpr float INV_R2 = 1.0f / (KERNEL_R * KERNEL_R); /**< Reciprocal of the squared support radius. */

  /**
   * @brief Squared Euclidean distance between two points (no sqrt).
   * @param a First point.
   * @param b Second point.
   * @return Squared distance between `a` and `b`.
   */
  static float dist2(const Vector &a, const Vector &b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }

  /**
   * @brief Refines a cubemap-LUT seed to its closest lattice node.
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Cached lattice node positions, indexed by node id.
   * @param center_node Seed node id from the cubemap LUT.
   * @return Node id of the genuine nearest node among the seed and its neighbors.
   * @details The CubemapLUT quantizes the query direction to a face cell, so its
   * seed can be a not-quite-nearest node; centering the kernel on the genuine
   * nearest node removes that quantization bias. Shared by both systems; without
   * it the interpolation inherits the cubemap grid artifact.
   */
  static int refine_nearest_node(const Vector &rv, const Vector *nodes,
                                 int center_node) {
    float best_d = dist2(rv, nodes[center_node]);
    int best_node = center_node;
    for_each_neighbor(center_node, [&](int ni) {
      float d = dist2(rv, nodes[ni]);
      if (d < best_d) {
        best_d = d;
        best_node = ni;
      }
    });
    return best_node;
  }

  /**
   * @brief Installs a Languid random-walk animation of the view orientation.
   * @details Shared by both systems; seeds OpenSimplex2 noise and adds the
   * random-walk animation to the timeline.
   */
  void init_orientation_animation() {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /**
   * @brief Builds the cached lattice node positions (RD_N entries) into `nodes`.
   * @param nodes Output array of at least RD_N positions to fill.
   * @details Called once at init; the Fibonacci lattice is static, not per-frame.
   */
  static void build_nodes(Vector *nodes) {
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);
  }

  /**
   * @brief Invokes `fn(ni)` for each valid K-NN neighbor index of `node`.
   * @tparam Fn Callable accepting a neighbor node id.
   * @param node Center node id whose neighbors are visited.
   * @param fn Callable invoked once per valid neighbor index.
   * @details The fundamental shared iteration over the lattice graph; the lambda
   * inlines, so this is identical codegen to a hand-written neighbor loop. The
   * discrete graph Laplacian of a field is `Σ_neighbors (field[ni] − center)`;
   * systems needing several fields per node fuse them into one walk here rather
   * than calling a single-field helper per field, which would re-read the
   * neighbor list once per species.
   */
  template <typename Fn> static void for_each_neighbor(int node, Fn &&fn) {
    for (int k = 0; k < RD_K; ++k) {
      int ni = ReactionGraph::neighbors[node][k];
      if (ni >= 0)
        fn(ni);
    }
  }

  /**
   * @brief Wendland C2 kernel walk over `center` and its neighbors.
   * @tparam OnWeight Callable accepting (node_index, weight).
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Cached lattice node positions, indexed by node id.
   * @param center Center node id of the kernel.
   * @param on_weight Callable invoked as `on_weight(node_index, weight)` for
   * every node inside the support radius.
   * @details The caller accumulates whatever fields it needs and applies the
   * total-weight guard, so this stays agnostic to species count and fixed-point.
   */
  template <typename OnWeight>
  static void kernel_accumulate(const Vector &rv, const Vector *nodes,
                                int center, OnWeight &&on_weight) {
    auto visit = [&](int i) {
      float u = 1.0f - dist2(rv, nodes[i]) * INV_R2;
      if (u > 0)
        on_weight(i, u * u);
    };
    visit(center);
    for_each_neighbor(center, visit);
  }

  Orientation<> orientation;        /**< Current view orientation on the sphere. */
  FastNoiseLite noise;              /**< Noise source driving the orientation walk. */
  ReactionGraph::CubemapLUT cube_lut; /**< Cubemap LUT for fast nearest-node seeding. */
  Timeline timeline;                /**< Animation timeline advancing the orientation. */
};
