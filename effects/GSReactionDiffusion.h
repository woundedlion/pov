/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <algorithm>
#include <cstring>
#include "core/effects_engine.h"
#include "effects/ReactionDiffusionBase.h"

/**
 * @brief Gray-Scott reaction-diffusion on a Fibonacci lattice sphere.
 *
 * Two species (A, B) evolve via Gray-Scott dynamics (A·B² autocatalysis with
 * feed/kill) on the shared 7680-node lattice, producing spots/stripes/mazes.
 * State is Q16 (uint16_t) for the cubic reaction-term precision. Shared
 * lattice/orientation/kernel scaffolding lives in ReactionDiffusionBase.
 */
template <int W, int H>
class GSReactionDiffusion
    : public ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::build_nodes;
  using Base::cube_lut;
  using Base::dist2;
  using Base::for_each_neighbor;
  using Base::init_orientation_animation;
  using Base::INV_R2;
  using Base::kernel_accumulate;
  using Base::orientation;
  using Base::RD_K;
  using Base::RD_N;
  using Base::registerParam;

public:
  FLASHMEM GSReactionDiffusion() = default;

  void init() override {
    // 170KB holds the 48KB Cubemap LUT + 30KB State + 90KB node positions. The
    // node array (7680 × Vector) is the fixed Fibonacci lattice — queries are
    // un-oriented onto it (not the reverse), so it does not depend on the
    // per-frame view orientation and is built ONCE here instead of every frame
    // (mirrors BZReactionDiffusion). Peak footprint is unchanged: the array was
    // already allocated per-frame in the scratch arena, so it just moves from
    // transient scratch to persistent and stops being recomputed each frame.
    configure_arenas(170 * 1024, GLOBAL_ARENA_SIZE - 170 * 1024, 0);

    registerParam("Feed", &params.feed, 0.0f, 0.1f);
    registerParam("Kill", &params.k, 0.0f, 0.1f);
    registerParam("dA", &params.dA, 0.0f, 1.0f);
    registerParam("dB", &params.dB, 0.0f, 1.0f);
    registerParam("Speed", &params.dt, 0.1f, 3.0f);

    state.A = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    state.B = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    // A starts at 1.0 (65535), B starts at 0.0
    for (int i = 0; i < RD_N; i++) {
      state.A[i] = 65535;
      state.B[i] = 0;
    }

    cube_lut.build(persistent_arena);
    nodes = static_cast<Vector *>(
        persistent_arena.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    build_nodes(nodes);
    seed_clusters();
    init_orientation_animation();
  }

private:
  // Q16: 16-bit fixed point needed for GS cubic reaction term precision
  static constexpr float Q16_SCALE = 65535.0f;
  static constexpr float Q16_INV = 1.0f / Q16_SCALE;
  static inline float from_q16(uint16_t v) { return v * Q16_INV; }
  static inline uint16_t to_q16(float v) {
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 1.0f) * Q16_SCALE);
  }

  void seed_clusters() {
    for (int i = 0; i < 30; i++) {
      int idx = hs::rand_int(0, RD_N - 1);
      state.B[idx] = 65535;
      for (int nb : ReactionGraph::neighbors[idx])
        if (nb >= 0)
          state.B[nb] = 65535;
    }
  }

  // Gray-Scott: dA/dt = dA·∇²A - A·B² + feed·(1-A)
  //             dB/dt = dB·∇²B + A·B² - (k+feed)·B
  void step_physics(uint16_t *nA, uint16_t *nB) {
    for (int i = 0; i < RD_N; i++) {
      float a = from_q16(state.A[i]);
      float b = from_q16(state.B[i]);

      // Both species' Laplacians share one neighbor walk (fused on purpose:
      // two single-field graph_laplacian() calls would double the lattice reads).
      float lA = 0, lB = 0;
      for_each_neighbor(i, [&](int ni) {
        lA += from_q16(state.A[ni]) - a;
        lB += from_q16(state.B[ni]) - b;
      });

      float abb = a * b * b;
      nA[i] = to_q16(a + (params.dA * lA - abb + params.feed * (1.0f - a)) *
                             params.dt);
      nB[i] = to_q16(b + (params.dB * lB + abb - (params.k + params.feed) * b) *
                             params.dt);
    }

    std::swap(state.A, nA);
    std::swap(state.B, nB);
  }

  float interpolate_b(const Vector &p, int nearest, const Vector *nodes) const {
    float tw = 0, wb = 0;
    kernel_accumulate(p, nodes, nearest, [&](int i, float w) {
      wb += from_q16(state.B[i]) * w;
      tw += w;
    });
    // Empty kernel (no node within the support radius): guard the division. A
    // 0/0 NaN would slip past the b < 0.05 cull in the shader (NaN compares
    // false) and silently poison palette.get(). Mirrors BZ's sample_kernel.
    if (tw <= 0.0001f)
      return 0.0f;
    return wb / tw;
  }

  void render(Canvas &canvas) {
    ScratchScope _frame(scratch_arena_a);
    uint16_t *sA = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    uint16_t *sB = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));

    for (int k = 0; k < 16; k++)
      step_physics(sA, sB);

    // `nodes` is the fixed lattice, built once in init() — not rebuilt here.
    auto shader = [&](const Vector &v) -> Color4 {
      Vector rv = orientation.unorient(v);
      int nearest = cube_lut.lookup(rv);
      float b = interpolate_b(rv, nearest, nodes);

      if (b < 0.05f)
        return Color4(Pixel(0, 0, 0), 0.0f);

      float t = hs::clamp((b - 0.1f) * 4.0f, 0.0f, 1.0f);
      return palette.get(t);
    };

    Scan::Shader::draw<W, H, 4>(canvas, shader);
  }

  struct {
    uint16_t *A = nullptr, *B = nullptr;
  } state;

  // Fixed Fibonacci-lattice node positions, built once in init() (see BZ).
  Vector *nodes = nullptr;

  GenerativePalette palette{
      GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY,
      BrightnessProfile::ASCENDING, SaturationProfile::VIBRANT};

  struct Params {
    float feed = 0.04f;
    float k = 0.06f;
    float dA = 0.02f;
    float dB = 0.01f;
    float dt = 2.5f;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(GSReactionDiffusion)
