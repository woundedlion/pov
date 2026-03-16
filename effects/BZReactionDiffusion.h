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

/**
 * @brief Belousov-Zhabotinsky reaction-diffusion on a Fibonacci lattice sphere.
 *
 * Three competing chemical species (A, B, C) evolve via Lotka-Volterra dynamics
 * on a 7680-node Fibonacci lattice with K=6 nearest neighbors. The cyclic
 * competition (A→B→C→A) creates self-sustaining spiral waves that persist
 * indefinitely. State is stored as Q8 (uint8_t). Per-pixel rendering uses
 * Wendland C2 kernel interpolation for smooth Voronoi boundaries.
 *
 * Memory budget (persistent arena):
 *   - State:  3 arrays × 7680 × 1B = 23,040 B
 *   - LUT:    288 × 144 × 2B       = 82,944 B
 *   - Total:  105,984 B (103 KB)
 *
 * Scratch arena (per frame):
 *   - Physics:  3 × 7680 × 1B      = 23,040 B
 *   - Node XYZ: 7680 × 12B         = 92,160 B
 *   - Total:    115,200 B (112 KB)
 */
template <int W, int H> class BZReactionDiffusion : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N;
  static constexpr int RD_K = ReactionGraph::RD_K;
  static constexpr int H_VIRT = H + hs::H_OFFSET;

  FLASHMEM BZReactionDiffusion() : Effect(W, H) { persist_pixels = false; }

  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - 140 * 1024, 140 * 1024, 0);

    registerParam("Alpha", &params.alpha, 0.0f, 4.0f);
    registerParam("Diff", &params.D, 0.001f, 0.1f);
    registerParam("Speed", &params.dt, 0.0f, 1.0f);

    state.A = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    state.B = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    state.C = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    memset(state.A, 0, RD_N);
    memset(state.B, 0, RD_N);
    memset(state.C, 0, RD_N);

    pixel_to_node = static_cast<uint16_t *>(
        persistent_arena.allocate(W * H * sizeof(uint16_t), alignof(uint16_t)));
    build_pixel_to_node_lut();
    seed_spiral_nuclei();

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
  static inline float from_q8(uint8_t v) { return v * (1.0f / 255.0f); }
  static inline uint8_t to_q8(float v) {
    return static_cast<uint8_t>(hs::clamp(v, 0.0f, 1.0f) * 255.0f);
  }

  static float dist2(const Vector &a, const Vector &b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }


  // Greedy walk on K-NN graph from Y-seeded estimate
  int find_nearest_node(const Vector &p) const {
    int cur = hs::clamp(
        static_cast<int>((1.0f - p.y) * 0.5f * (RD_N - 1) + 0.5f), 0, RD_N - 1);
    float best_d = dist2(p, ReactionGraph::node(cur));

    for (int iter = 0; iter < 64; ++iter) {
      bool improved = false;
      for (int k = 0; k < RD_K; ++k) {
        int ni = ReactionGraph::neighbors[cur][k];
        if (ni < 0)
          continue;
        float d = dist2(p, ReactionGraph::node(ni));
        if (d < best_d) {
          best_d = d;
          cur = ni;
          improved = true;
        }
      }
      if (!improved)
        break;
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

  // 3 clusters per species ensures all 3 are always present
  void seed_spiral_nuclei() {
    uint8_t *species[] = {state.A, state.B, state.C};
    for (int s = 0; s < 3; s++) {
      for (int k = 0; k < 3; k++) {
        int center = hs::rand_int(0, RD_N - 1);
        species[s][center] = 255;
        for (int nb : ReactionGraph::neighbors[center])
          if (nb >= 0)
            species[s][nb] = 255;
      }
    }
  }

  // Lotka-Volterra reaction + graph Laplacian diffusion
  void step_physics(uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    for (int i = 0; i < RD_N; i++) {
      float a = from_q8(state.A[i]);
      float b = from_q8(state.B[i]);
      float c = from_q8(state.C[i]);

      float lA = 0, lB = 0, lC = 0;
      for (int k = 0; k < RD_K; k++) {
        int ni = ReactionGraph::neighbors[i][k];
        if (ni < 0)
          continue;
        lA += from_q8(state.A[ni]) - a;
        lB += from_q8(state.B[ni]) - b;
        lC += from_q8(state.C[ni]) - c;
      }

      nA[i] = to_q8(a + (params.D * lA + a * (1 - a - params.alpha * c)) *
                            params.dt);
      nB[i] = to_q8(b + (params.D * lB + b * (1 - b - params.alpha * a)) *
                            params.dt);
      nC[i] = to_q8(c + (params.D * lC + c * (1 - c - params.alpha * b)) *
                            params.dt);
    }

    // Stochastic perturbation prevents convergence on closed manifold
    for (int p = 0; p < 8; p++) {
      int idx = hs::rand_int(0, RD_N - 1);
      int s = hs::rand_int(0, 2);
      uint8_t *t = (s == 0) ? nA : (s == 1) ? nB : nC;
      t[idx] =
          static_cast<uint8_t>(std::min(static_cast<int>(t[idx]) + 3, 255));
    }

    std::swap(state.A, nA);
    std::swap(state.B, nB);
    std::swap(state.C, nC);
  }

  // Wendland C2 compact kernel: w(d) = max(0, 1 - d²/R²)²
  static constexpr float D_AVG = 0.04044f; // sqrt(4π / RD_N)
  static constexpr float KERNEL_R = 1.5f * D_AVG;
  static constexpr float INV_R2 = 1.0f / (KERNEL_R * KERNEL_R);

  struct InterpolatedState {
    float a, b, c;
  };

  InterpolatedState interpolate_at(const Vector &p, int nearest,
                                   const Vector *nodes) const {
    float tw = 0, wa = 0, wb = 0, wc = 0;
    auto acc = [&](int i) {
      float u = 1.0f - dist2(p, nodes[i]) * INV_R2;
      if (u <= 0)
        return;
      float w = u * u;
      wa += from_q8(state.A[i]) * w;
      wb += from_q8(state.B[i]) * w;
      wc += from_q8(state.C[i]) * w;
      tw += w;
    };

    acc(nearest);
    for (int k = 0; k < RD_K; ++k) {
      int ni = ReactionGraph::neighbors[nearest][k];
      if (ni >= 0)
        acc(ni);
    }
    float inv = 1.0f / tw;
    return {wa * inv, wb * inv, wc * inv};
  }

  static Pixel blend_species(float a, float b, float c, const Color4 &ca,
                             const Color4 &cb, const Color4 &cc) {
    float r = ca.color.r * a, g = ca.color.g * a, bl = ca.color.b * a;

    float ib = 1.0f - b;
    r = r * ib + cb.color.r * b;
    g = g * ib + cb.color.g * b;
    bl = bl * ib + cb.color.b * b;

    float ic = 1.0f - c;
    r = r * ic + cc.color.r * c;
    g = g * ic + cc.color.g * c;
    bl = bl * ic + cc.color.b * c;

    return Pixel(static_cast<uint16_t>(hs::clamp(r, 0.0f, 65535.0f)),
                 static_cast<uint16_t>(hs::clamp(g, 0.0f, 65535.0f)),
                 static_cast<uint16_t>(hs::clamp(bl, 0.0f, 65535.0f)));
  }

  int lookup_nearest(const Vector &rv) const {
    Spherical s(rv);
    int x = static_cast<int>((s.theta * W) / (2.0f * PI_F) + 0.5f) % W;
    int y = hs::clamp(static_cast<int>((s.phi * (H_VIRT - 1)) / PI_F + 0.5f), 0,
                      H - 1);
    if (x < 0)
      x += W;
    return pixel_to_node[y * W + x];
  }

  void render(Canvas &canvas) {
    ScratchScope _frame(scratch_arena_a);
    uint8_t *sA = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *sB = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *sC = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));

    for (int k = 0; k < 2; k++)
      step_physics(sA, sB, sC);

    Vector *nodes = static_cast<Vector *>(
        scratch_arena_a.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);

    Color4 ca = palette.get(0.0f);
    Color4 cb = palette.get(0.5f);
    Color4 cc = palette.get(1.0f);

    auto shader = [&](const Vector &v) -> Color4 {
      Vector rv = orientation.unorient(v);
      int nearest = lookup_nearest(rv);
      auto [a, b, c] = interpolate_at(rv, nearest, nodes);
      if (a + b + c < 0.01f)
        return Color4(Pixel(0, 0, 0), 0.0f);
      return Color4(blend_species(a, b, c, ca, cb, cc), 1.0f);
    };

    Scan::Shader::draw<W, H, 4>(canvas, shader);
  }

  struct {
    uint8_t *A = nullptr, *B = nullptr, *C = nullptr;
  } state;

  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::DESCENDING,
                            SaturationProfile::VIBRANT, 42};
  Orientation<W> orientation;
  FastNoiseLite noise;
  uint16_t *pixel_to_node = nullptr;
  Timeline<W> timeline;

  struct Params {
    float alpha = 3.0f;
    float D = 0.05f;
    float dt = 0.35f;
  } params;
};

#include "../effect_registry.h"
REGISTER_EFFECT(BZReactionDiffusion)
