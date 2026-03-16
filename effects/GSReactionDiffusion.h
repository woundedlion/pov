/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <algorithm>
#include <cstring>
#include "../effects_engine.h"
#include "../reaction_graph.h"
#include "../scan.h"

template <int W, int H> class GSReactionDiffusion : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N;
  static constexpr int RD_K = ReactionGraph::RD_K;
  static constexpr int H_VIRT = H + hs::H_OFFSET;

  FLASHMEM GSReactionDiffusion() : Effect(W, H) { persist_pixels = false; }

  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - 120 * 1024, 120 * 1024, 0);

    registerParam("Feed", &params.feed, 0.0f, 0.1f);
    registerParam("Kill", &params.k, 0.0f, 0.1f);
    registerParam("dA", &params.dA, 0.0f, 1.0f);
    registerParam("dB", &params.dB, 0.0f, 1.0f);
    registerParam("Speed", &params.dt, 0.1f, 2.0f);

    state.A = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    state.B = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    // A starts at 1.0 (65535), B starts at 0.0
    for (int i = 0; i < RD_N; i++) { state.A[i] = 65535; state.B[i] = 0; }

    pixel_to_node = static_cast<uint16_t *>(
        persistent_arena.allocate(W * H * sizeof(uint16_t), alignof(uint16_t)));
    build_pixel_to_node_lut();
    seed_clusters();

    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    render(canvas);
  }

private:
  // Q16: 16-bit fixed point needed for GS cubic reaction term precision
  static constexpr float Q16_SCALE = 65535.0f;
  static constexpr float Q16_INV = 1.0f / Q16_SCALE;
  static inline float from_q16(uint16_t v) { return v * Q16_INV; }
  static inline uint16_t to_q16(float v) {
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 1.0f) * Q16_SCALE);
  }

  static float dist2(const Vector &a, const Vector &b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }

  // Greedy walk on K-NN graph from Y-seeded estimate
  int find_nearest_node(const Vector &p) const {
    int cur = hs::clamp(
        static_cast<int>((1.0f - p.y) * 0.5f * (RD_N - 1) + 0.5f), 0,
        RD_N - 1);
    float best_d = dist2(p, ReactionGraph::node(cur));

    for (int iter = 0; iter < 64; ++iter) {
      bool improved = false;
      for (int k = 0; k < RD_K; ++k) {
        int ni = ReactionGraph::neighbors[cur][k];
        if (ni < 0) continue;
        float d = dist2(p, ReactionGraph::node(ni));
        if (d < best_d) { best_d = d; cur = ni; improved = true; }
      }
      if (!improved) break;
    }
    return cur;
  }

  void build_pixel_to_node_lut() {
    for (int y = 0; y < H; ++y) {
      float phi = (static_cast<float>(y) * PI_F) / (H_VIRT - 1);
      float sp = sinf(phi), cp = cosf(phi);
      for (int x = 0; x < W; ++x) {
        float theta = (static_cast<float>(x) * 2.0f * PI_F) / W;
        Vector p(sp * cosf(theta), cp, sp * sinf(theta));
        pixel_to_node[y * W + x] = static_cast<uint16_t>(find_nearest_node(p));
      }
    }
  }

  void seed_clusters() {
    for (int i = 0; i < 30; i++) {
      int idx = hs::rand_int(0, RD_N - 1);
      state.B[idx] = 65535;
      for (int nb : ReactionGraph::neighbors[idx])
        if (nb >= 0) state.B[nb] = 65535;
    }
  }

  // Gray-Scott: dA/dt = dA·∇²A - A·B² + feed·(1-A)
  //             dB/dt = dB·∇²B + A·B² - (k+feed)·B
  void step_physics(uint16_t *nA, uint16_t *nB) {
    for (int i = 0; i < RD_N; i++) {
      float a = from_q16(state.A[i]);
      float b = from_q16(state.B[i]);

      float lA = 0, lB = 0;
      for (int k = 0; k < RD_K; k++) {
        int ni = ReactionGraph::neighbors[i][k];
        if (ni < 0) continue;
        lA += from_q16(state.A[ni]) - a;
        lB += from_q16(state.B[ni]) - b;
      }

      float abb = a * b * b;
      nA[i] = to_q16(a + (params.dA * lA - abb + params.feed * (1.0f - a)) * params.dt);
      nB[i] = to_q16(b + (params.dB * lB + abb - (params.k + params.feed) * b) * params.dt);
    }

    std::swap(state.A, nA);
    std::swap(state.B, nB);
  }

  // Wendland C2 compact kernel: w(d) = max(0, 1 - d²/R²)²
  static constexpr float D_AVG = 0.04044f; // sqrt(4π / RD_N)
  static constexpr float KERNEL_R = 1.5f * D_AVG;
  static constexpr float INV_R2 = 1.0f / (KERNEL_R * KERNEL_R);

  float interpolate_b(const Vector &p, int nearest,
                      const Vector *nodes) const {
    float tw = 0, wb = 0;
    auto acc = [&](int i) {
      float u = 1.0f - dist2(p, nodes[i]) * INV_R2;
      if (u <= 0) return;
      float w = u * u;
      wb += from_q16(state.B[i]) * w;
      tw += w;
    };

    acc(nearest);
    for (int k = 0; k < RD_K; ++k) {
      int ni = ReactionGraph::neighbors[nearest][k];
      if (ni >= 0) acc(ni);
    }
    return wb / tw;
  }

  int lookup_nearest(const Vector &rv) const {
    Spherical s(rv);
    int x = static_cast<int>((s.theta * W) / (2.0f * PI_F) + 0.5f) % W;
    int y = hs::clamp(
        static_cast<int>((s.phi * (H_VIRT - 1)) / PI_F + 0.5f), 0, H - 1);
    if (x < 0) x += W;
    return pixel_to_node[y * W + x];
  }

  void render(Canvas &canvas) {
    ScratchScope _frame(scratch_arena_a);
    uint16_t *sA = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    uint16_t *sB = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));

    for (int k = 0; k < 16; k++)
      step_physics(sA, sB);

    Vector *nodes = static_cast<Vector *>(
        scratch_arena_a.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);

    auto shader = [&](const Vector &v) -> Color4 {
      Vector rv = orientation.unorient(v);
      int nearest = lookup_nearest(rv);
      float b = interpolate_b(rv, nearest, nodes);

      if (b < 0.05f) return Color4(Pixel(0, 0, 0), 0.0f);

      float t = hs::clamp((b - 0.1f) * 4.0f, 0.0f, 1.0f);
      return palette.get(t);
    };

    Scan::Shader::draw<W, H, 4>(canvas, shader);
  }

  struct { uint16_t *A = nullptr, *B = nullptr; } state;

  GenerativePalette palette{GradientShape::STRAIGHT,
                            HarmonyType::SPLIT_COMPLEMENTARY,
                            BrightnessProfile::ASCENDING,
                            SaturationProfile::VIBRANT};
  Orientation<W> orientation;
  FastNoiseLite noise;
  uint16_t *pixel_to_node = nullptr;
  Timeline<W> timeline;

  struct Params {
    float feed = 0.04f;
    float k = 0.06f;
    float dA = 0.02f;
    float dB = 0.01f;
    float dt = 2.5f;
  } params;
};

#include "../effect_registry.h"
REGISTER_EFFECT(GSReactionDiffusion)
