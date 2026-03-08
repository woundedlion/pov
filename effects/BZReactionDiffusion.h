/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <algorithm>
#include "../effects_engine.h"
#include "../generators.h"
#include "../reaction_graph.h"

template <int W, int H> class BZReactionDiffusion : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N;
  static constexpr int RD_K = ReactionGraph::RD_K;

  FLASHMEM BZReactionDiffusion()
      : Effect(W, H), filters(Filter::World::Orient<W>(orientation),
                              Filter::Screen::AntiAlias<W, H>()) {}

  void init() override {
    registerParam("Alpha", &params.alpha, 0.0f, 2.0f);
    registerParam("Diff", &params.D, 0.001f, 0.1f);
    registerParam("Speed", &params.dt, 0.0f, 1.0f);
    registerParam("GlobalAlpha", &params.global_alpha, 0.0f, 1.0f);

    // Compute nodes from formula (neighbors are PROGMEM)
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);
    timeline
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, PI_F / 2, 600,
                                       ease_mid, true))
        .add(0, Animation::PeriodicTimer(
                    96, [this](Canvas &c) { this->spawn(); }, true));
    spawn();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  /**
   * @brief Holds the chemical state and parameters for a single
   * Belousov-Zhabotinsky reaction instance.
   * @details Simulates a 3-species cyclic competition (Rock-Paper-Scissors).
   */
  struct BZReactionContext {
    std::array<float, RD_N> A;
    std::array<float, RD_N> B;
    std::array<float, RD_N> C;
    std::array<float, RD_N> nextA;
    std::array<float, RD_N> nextB;
    std::array<float, RD_N> nextC;

    // Simulation Parameters Removed from Context (Use Global Params)

    GenerativePalette palette;

    BZReactionContext()
        : palette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                  BrightnessProfile::DESCENDING, SaturationProfile::VIBRANT) {
      reset();
    }

    void reset() {
      A.fill(0.0f);
      B.fill(0.0f);
      C.fill(0.0f);
      palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                                  BrightnessProfile::DESCENDING,
                                  SaturationProfile::VIBRANT);
    }
  };

  struct BZGridGenerator : public IGenerator<std::array<Vector, RD_N>> {
    std::array<Vector, RD_N> generate(Arena &geom,
                                      MemoryCtx &ctx) const override {
      std::array<Vector, RD_N> out_nodes;
      const float phi = PI_F * (3.0f - sqrtf(5.0f));
      for (int i = 0; i < RD_N; i++) {
        float y = 1.0f - (static_cast<float>(i) / (RD_N - 1)) * 2.0f;
        float radius = sqrtf(1.0f - y * y);
        float theta = phi * i;
        out_nodes[i] = Vector(cosf(theta) * radius, y, sinf(theta) * radius);
      }
      return out_nodes;
    }
  };

  void spawn() {
    contexts.push_back(BZReactionContext());
    BZReactionContext &ctx = contexts.back();
    seed(ctx);

    timeline.add(0, Animation::Sprite(
                        [this, &ctx](Canvas &c, float opacity) {
                          this->render_reaction(c, ctx, opacity);
                        },
                        192, 32, ease_mid, 32, ease_mid));
  }

  void seed(BZReactionContext &ctx) {
    // Sparse Seeding for Spirals (Droplets)
    for (int k = 0; k < 50; k++) {
      int center = hs::rand_int(0, RD_N - 1);
      float r = hs::rand_f();

      // Determine which species to seed
      float *target = (r < 0.33f)   ? ctx.A.data()
                      : (r < 0.66f) ? ctx.B.data()
                                    : ctx.C.data();

      target[center] = 1.0f;
      for (int neighbor : ReactionGraph::neighbors[center]) {
        if (neighbor >= 0)
          target[neighbor] = 1.0f;
      }
    }
  }

  void render_reaction(Canvas &canvas, BZReactionContext &ctx, float opacity) {
    // Simulate Physics
    // JS Parity: 2 steps per frame (was 12)
    for (int k = 0; k < 2; k++) {
      update_physics(ctx);
    }

    // Pre-fetch palette colors
    Color4 ca = ctx.palette.get(0.0f);
    Color4 cb = ctx.palette.get(0.5f);
    Color4 cc = ctx.palette.get(1.0f);

    // Draw
    for (int i = 0; i < RD_N; i++) {
      float a = ctx.A[i];
      float b = ctx.B[i];
      float c = ctx.C[i];
      float sum = a + b + c;

      if (sum > 0.01f) {
        // Blending logic: Start black, lerp towards A, then B, then C
        // Step 1: Base Black -> A
        float r = ca.color.r * a;
        float g = ca.color.g * a;
        float b_ch = ca.color.b * a;

        // Step 2: -> B
        float inv_b = 1.0f - b;
        r = r * inv_b + cb.color.r * b;
        g = g * inv_b + cb.color.g * b;
        b_ch = b_ch * inv_b + cb.color.b * b;

        // Step 3: -> C
        float inv_c = 1.0f - c;
        r = r * inv_c + cc.color.r * c;
        g = g * inv_c + cc.color.g * c;
        b_ch = b_ch * inv_c + cc.color.b * c;

        Pixel color(static_cast<uint16_t>(hs::clamp(r, 0.0f, 65535.0f)),
                    static_cast<uint16_t>(hs::clamp(g, 0.0f, 65535.0f)),
                    static_cast<uint16_t>(hs::clamp(b_ch, 0.0f, 65535.0f)));

        Color4 c_final(color);
        c_final.alpha *= params.global_alpha * opacity;
        auto shader = [c_final](const Vector &p, Fragment &f) {
          f.color = c_final;
        };
        Plot::Point::draw(filters, canvas, Fragment(nodes[i]), shader);
      }
    }
  }

  void update_physics(BZReactionContext &ctx) {
    for (int i = 0; i < RD_N; i++) {
      float a = ctx.A[i];
      float b = ctx.B[i];
      float c = ctx.C[i];

      float lapA = 0.0f;
      float lapB = 0.0f;
      float lapC = 0.0f;

      for (int k = 0; k < RD_K; k++) {
        int n_idx = ReactionGraph::neighbors[i][k];
        if (n_idx < 0)
          continue;

        lapA += (ctx.A[n_idx] - a);
        lapB += (ctx.B[n_idx] - b);
        lapC += (ctx.C[n_idx] - c);
      }

      // 3-Species Cyclic Model
      float da = a * (1.0f - a - params.alpha * c);
      float db = b * (1.0f - b - params.alpha * a);
      float dc = c * (1.0f - c - params.alpha * b);

      ctx.nextA[i] =
          hs::clamp(a + (params.D * lapA + da) * params.dt, 0.0f, 1.0f);
      ctx.nextB[i] =
          hs::clamp(b + (params.D * lapB + db) * params.dt, 0.0f, 1.0f);
      ctx.nextC[i] =
          hs::clamp(c + (params.D * lapC + dc) * params.dt, 0.0f, 1.0f);
    }

    // Swap buffers
    std::swap(ctx.A, ctx.nextA);
    std::swap(ctx.B, ctx.nextB);
    std::swap(ctx.C, ctx.nextC);
  }

  std::array<Vector, RD_N> nodes;
  Orientation<W> orientation;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  StaticCircularBuffer<BZReactionContext, 4> contexts;
  Timeline<W> timeline;

  struct Params {
    float alpha = 1.6f;
    float D = 0.03f;
    float dt = 0.2f;
    float global_alpha = 0.3f;
  } params;
};