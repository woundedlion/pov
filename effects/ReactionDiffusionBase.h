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
 * Languid random-walk view orientation, rebuild the cached node positions each
 * frame, and interpolate with the same Wendland C2 kernel. This base captures
 * exactly that shared scaffolding. The physics (Lotka-Volterra 3-species Q8 vs
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
  static constexpr int H_VIRT = H + hs::H_OFFSET;

  ReactionDiffusionBase() : Effect(W, H) { persist_pixels = false; }

  bool show_bg() const override { return true; }

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

  static float dist2(const Vector &a, const Vector &b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }

  /** Languid random-walk of the view orientation (shared by both systems). */
  void init_orientation_animation() {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /** Rebuild the cached lattice node positions (RD_N entries) into `nodes`. */
  static void build_nodes(Vector *nodes) {
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);
  }

  Orientation<W> orientation;
  FastNoiseLite noise;
  ReactionGraph::CubemapLUT cube_lut;
  Timeline<W> timeline;
};
