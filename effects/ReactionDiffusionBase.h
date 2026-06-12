/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

/**
 * @brief CRTP base for the spherical reaction-diffusion effects (BZ, Gray-Scott).
 *
 * Both systems run on the same 7680-node Fibonacci-lattice K-NN graph, share a
 * Languid random-walk view orientation, build the cached node positions once at
 * init (the lattice is static), and interpolate with the same Wendland C2
 * kernel. This base captures exactly that shared scaffolding. The physics
 * (Lotka-Volterra 3-species Q8 vs
 * Gray-Scott 2-species Q16), state representation, seeding, params, palette, and
 * rendering are fundamentally different and stay in the derived classes.
 *
 * Dispatch is static (CRTP): draw_frame() forwards to Derived::render() with no
 * virtual call, so the abstraction is zero-overhead. Derived classes befriend
 * this base so render() can stay private.
 */
template <typename Derived, int W, int H>
class ReactionDiffusionBase : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N;
  static constexpr int RD_K = ReactionGraph::RD_K;

  /** Each frame is fully repainted from lattice state, so the framebuffer
   *  need not persist between frames. */
  ReactionDiffusionBase() : Effect(W, H) { persist_pixels = false; }

  /** Both systems composite over the framebuffer background. */
  bool show_bg() const override { return true; }

  /** Advance the orientation timeline, then statically dispatch to
   *  Derived::render() for the system-specific physics and drawing. */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    static_cast<Derived *>(this)->render(canvas);
  }

protected:
  // Wendland C2 compact kernel: w(d) = max(0, 1 - d²/R²)²
  static constexpr float D_AVG = ReactionGraph::D_AVG; // sqrt(4π / RD_N)
  static constexpr float KERNEL_R = 1.5f * D_AVG;
  static constexpr float INV_R2 = 1.0f / (KERNEL_R * KERNEL_R);

  /** Squared Euclidean distance between two points (no sqrt). */
  static float dist2(const Vector &a, const Vector &b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }

  /** Refine a cubemap-LUT seed to its closest lattice node by checking the
   *  seed's direct neighbors. The CubemapLUT quantizes the query direction to a
   *  face cell, so its seed can be a not-quite-nearest node; centering the kernel
   *  on the genuine nearest node removes that quantization bias. Shared by both
   *  systems — without it the interpolation inherits the cubemap grid artifact. */
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

  /** Languid random-walk of the view orientation (shared by both systems). */
  void init_orientation_animation() {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /** Build the cached lattice node positions (RD_N entries) into `nodes`.
   *  Called once at init — the Fibonacci lattice is static, not per-frame. */
  static void build_nodes(Vector *nodes) {
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);
  }

  /** Invoke `fn(ni)` for each valid K-NN neighbor index of `node`.
   *  The fundamental shared iteration over the lattice graph; the lambda
   *  inlines, so this is identical codegen to a hand-written neighbor loop.
   *  The discrete graph Laplacian of a field is `Σ_neighbors (field[ni] −
   *  center)`; systems needing several fields per node fuse them into one walk
   *  here rather than calling a single-field helper per field, which would
   *  re-read the neighbor list once per species. */
  template <typename Fn> static void for_each_neighbor(int node, Fn &&fn) {
    for (int k = 0; k < RD_K; ++k) {
      int ni = ReactionGraph::neighbors[node][k];
      if (ni >= 0)
        fn(ni);
    }
  }

  /** Wendland C2 kernel walk over `center` and its neighbors. Calls
   *  `on_weight(node_index, weight)` for every node inside the support radius.
   *  The caller accumulates whatever fields it needs and applies the total-
   *  weight guard, so this stays agnostic to species count and fixed-point. */
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

  Orientation<> orientation;
  FastNoiseLite noise;
  ReactionGraph::CubemapLUT cube_lut;
  Timeline timeline;
};
