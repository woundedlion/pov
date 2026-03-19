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

template <int W, int H> class BZReactionDiffusion : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N;
  static constexpr int RD_K = ReactionGraph::RD_K;
  static constexpr int H_VIRT = H + hs::H_OFFSET;

  // 64x64 per face gives incredible accuracy for 3840 nodes
  static constexpr int CUBE_RES = 64;

  FLASHMEM BZReactionDiffusion() : Effect(W, H) { persist_pixels = false; }

  void init() override {
    // 75KB easily holds the 48KB Cubemap LUT + 23KB State
    configure_arenas(75 * 1024, GLOBAL_ARENA_SIZE - 75 * 1024, 0);

    registerParam("Alpha", &params.alpha, 0.0f, 4.0f);
    registerParam("Diff", &params.D, 0.001f, 0.1f);
    registerParam("Speed", &params.dt, 0.0f, 1.0f);

    state.A = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    state.B = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    state.C = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    memset(state.A, 0, RD_N);
    memset(state.B, 0, RD_N);
    memset(state.C, 0, RD_N);

    // Allocate and build the fast Cubemap LUT
    cube_lut = static_cast<uint16_t *>(persistent_arena.allocate(
        6 * CUBE_RES * CUBE_RES * sizeof(uint16_t), alignof(uint16_t)));
    build_cubemap_lut();

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

  // Used only once at startup
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

  void build_cubemap_lut() {
    for (int face = 0; face < 6; ++face) {
      for (int y = 0; y < CUBE_RES; ++y) {
        for (int x = 0; x < CUBE_RES; ++x) {
          float u = (x + 0.5f) / CUBE_RES * 2.0f - 1.0f;
          float v = (y + 0.5f) / CUBE_RES * 2.0f - 1.0f;
          Vector dir;
          if (face == 0)
            dir = Vector(1.0f, v, -u);       // +X
          else if (face == 1)
            dir = Vector(-1.0f, v, u);        // -X
          else if (face == 2)
            dir = Vector(u, 1.0f, -v);        // +Y
          else if (face == 3)
            dir = Vector(u, -1.0f, v);        // -Y
          else if (face == 4)
            dir = Vector(u, v, 1.0f);         // +Z
          else
            dir = Vector(-u, v, -1.0f);       // -Z

          float len =
              sqrtf(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
          dir.x /= len;
          dir.y /= len;
          dir.z /= len;
          cube_lut[(face * CUBE_RES + y) * CUBE_RES + x] =
              find_nearest_node(dir);
        }
      }
    }
  }

  // ZERO Trigonometry. Runs in ~10 cycles!
  int lookup_cubemap(const Vector &p) const {
    float ax = fabsf(p.x), ay = fabsf(p.y), az = fabsf(p.z);
    int face = 0;
    float u = 0, v = 0;

    if (ax >= ay && ax >= az) {
      float inv = 1.0f / ax;
      if (p.x >= 0) {
        face = 0; u = -p.z * inv; v = p.y * inv;
      } else {
        face = 1; u = p.z * inv; v = p.y * inv;
      }
    } else if (ay >= ax && ay >= az) {
      float inv = 1.0f / ay;
      if (p.y >= 0) {
        face = 2; u = p.x * inv; v = -p.z * inv;
      } else {
        face = 3; u = p.x * inv; v = p.z * inv;
      }
    } else {
      float inv = 1.0f / az;
      if (p.z >= 0) {
        face = 4; u = p.x * inv; v = p.y * inv;
      } else {
        face = 5; u = -p.x * inv; v = p.y * inv;
      }
    }

    int ui = hs::clamp(static_cast<int>((u + 1.0f) * 0.5f * CUBE_RES), 0,
                       CUBE_RES - 1);
    int vi = hs::clamp(static_cast<int>((v + 1.0f) * 0.5f * CUBE_RES), 0,
                       CUBE_RES - 1);
    return cube_lut[(face * CUBE_RES + vi) * CUBE_RES + ui];
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

    // 1. VERTEX SHADER: O(1) predictable cubemap lookup (Runs 1× per pixel)
    auto vertex_shader = [&](Fragment &frag) {
      Vector rv = orientation.unorient(frag.pos);
      frag.v0 = static_cast<float>(lookup_cubemap(rv));
    };

    // 2. FRAGMENT SHADER: Fully inlined kernel (Runs 4× per pixel for SSAA)
    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      int center_node = static_cast<int>(frag.v0);
      Vector rv = orientation.unorient(v);

      // Pass A: Find the absolute best node out of the center's 7-node cluster
      float best_d = dist2(rv, nodes[center_node]);
      int best_node = center_node;

      for (int k = 0; k < RD_K; ++k) {
        int ni = ReactionGraph::neighbors[center_node][k];
        if (ni >= 0) {
          float d = dist2(rv, nodes[ni]);
          if (d < best_d) {
            best_d = d;
            best_node = ni;
          }
        }
      }

      // Pass B: Calculate the Wendland Kernel around the true best node
      float tw = 0, wa = 0, wb = 0, wc = 0;

      // Reuse best_d from Pass A (saves one dist2 call)
      float u0 = 1.0f - best_d * INV_R2;
      if (u0 > 0) {
        float w = u0 * u0;
        wa += state.A[best_node] * w;
        wb += state.B[best_node] * w;
        wc += state.C[best_node] * w;
        tw += w;
      }

      // Loop neighbors of best node
      for (int k = 0; k < RD_K; ++k) {
        int ni = ReactionGraph::neighbors[best_node][k];
        if (ni >= 0) {
          float d = dist2(rv, nodes[ni]);
          float u = 1.0f - d * INV_R2;
          if (u > 0) {
            float w = u * u;
            wa += state.A[ni] * w;
            wb += state.B[ni] * w;
            wc += state.C[ni] * w;
            tw += w;
          }
        }
      }

      // Final output with deferred Q8 div
      if (tw <= 0.0001f) {
        frag.color = Color4(Pixel(0, 0, 0), 0.0f);
      } else {
        float inv = (1.0f / 255.0f) / tw;
        frag.color = Color4(
            blend_species(wa * inv, wb * inv, wc * inv, ca, cb, cc), 1.0f);
      }
    };

    Scan::Shader::draw<W, H, 4>(canvas, fragment_shader, vertex_shader);
  }

  struct {
    uint8_t *A = nullptr, *B = nullptr, *C = nullptr;
  } state;

  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::DESCENDING,
                            SaturationProfile::VIBRANT, 42};
  Orientation<W> orientation;
  FastNoiseLite noise;
  uint16_t *cube_lut = nullptr;
  Timeline<W> timeline;

  struct Params {
    float alpha = 3.0f;
    float D = 0.05f;
    float dt = 0.35f;
  } params;
};

#include "../effect_registry.h"
REGISTER_EFFECT(BZReactionDiffusion)
