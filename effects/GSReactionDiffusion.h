/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../reaction_graph.h"
#include <algorithm>
#include <array>

template <int W, int H> class GSReactionDiffusion : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N;
  static constexpr int RD_K = ReactionGraph::RD_K;

  FLASHMEM GSReactionDiffusion()
      : Effect(W, H), filters(Filter::World::Orient<W>(orientation),
                              Filter::Screen::AntiAlias<W, H>()) {}

  void init() override {
    registerParam("Feed", &params.feed, 0.0f, 0.1f);
    registerParam("Kill (k)", &params.k, 0.0f, 0.1f);
    registerParam("dA", &params.dA, 0.0f, 1.0f);
    registerParam("dB", &params.dB, 0.0f, 1.0f);
    registerParam("dt", &params.dt, 0.1f, 2.0f);
    registerParam("Global Alpha", &params.global_alpha, 0.0f, 1.0f);

    // Compute nodes from formula (neighbors are PROGMEM)
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);
    timeline
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, PI_F / 2, 64,
                                       ease_mid, true))
        .add(0, Animation::PeriodicTimer(
                    96, [this](Canvas &c) { this->spawn(); }, true));

    // Seed and spawn
    spawn();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  /**
   * @brief Holds the chemical state and parameters for a single Gray-Scott
   * reaction instance.
   */
  struct GSReactionContext {
    std::array<float, RD_N> A;
    std::array<float, RD_N> B;
    std::array<float, RD_N> nextA;
    std::array<float, RD_N> nextB;

    // Simulation params moved to main Effect class

    GenerativePalette palette;

    GSReactionContext()
        : palette(GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY,
                  BrightnessProfile::ASCENDING, SaturationProfile::VIBRANT) {
      reset();
    }

    void reset() {
      A.fill(1.0f);
      B.fill(0.0f);
      // Regenerate palette for variety on every spawn
      palette = GenerativePalette(
          GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY,
          BrightnessProfile::ASCENDING, SaturationProfile::VIBRANT);
    }
  };

  void spawn() {
    contexts.push_back(GSReactionContext());
    GSReactionContext &ctx = contexts.back();
    seed(ctx);

    timeline.add(0, Animation::Sprite(
                        [this, &ctx](Canvas &c, float opacity) {
                          this->render_reaction(c, ctx, opacity);
                        },
                        192, 32, ease_mid, 32, ease_mid));
  }

  void seed(GSReactionContext &ctx) {
    for (int i = 0; i < 5; i++) {
      int idx = hs::rand_int(0, RD_N);
      ctx.B[idx] = 1.0f;
      for (int neighbor : ReactionGraph::neighbors[idx]) {
        if (neighbor >= 0)
          ctx.B[neighbor] = 1.0f;
      }
    }
  }

  void render_reaction(Canvas &canvas, GSReactionContext &ctx, float opacity) {
    // Simulate Physics (8 steps per frame)
    for (int k = 0; k < 8; k++) {
      update_physics(ctx);
    }

    // Draw
    for (int i = 0; i < RD_N; i++) {
      float b = ctx.B[i];
      if (b > 0.1f) {
        float t = hs::clamp((b - 0.15f) * 4.0f, 0.0f, 1.0f);
        Color4 c = ctx.palette.get(t);
        c.alpha *= opacity * params.global_alpha;
        auto shader = [c](const Vector &p, Fragment &f) { f.color = c; };
        Plot::Point::draw(filters, canvas, Fragment(nodes[i]), shader);
      }
    }
  }

  void update_physics(GSReactionContext &ctx) {
    for (int i = 0; i < RD_N; i++) {
      float a = ctx.A[i];
      float b = ctx.B[i];

      // Laplacian (Graph-based)
      float lapA = 0.0f;
      float lapB = 0.0f;

      for (int k = 0; k < RD_K; k++) {
        int n_idx = ReactionGraph::neighbors[i][k];
        if (n_idx < 0)
          continue;

        lapA += (ctx.A[n_idx] - a);
        lapB += (ctx.B[n_idx] - b);
      }

      // Reaction-Diffusion logic
      float reaction = a * b * b;
      float feed_term = params.feed * (1.0f - a);
      float kill_term = (params.k + params.feed) * b;

      ctx.nextA[i] =
          hs::clamp(a + (params.dA * lapA - reaction + feed_term) * params.dt,
                    0.0f, 1.0f);
      ctx.nextB[i] =
          hs::clamp(b + (params.dB * lapB + reaction - kill_term) * params.dt,
                    0.0f, 1.0f);
    }

    // Swap buffers
    std::swap(ctx.A, ctx.nextA);
    std::swap(ctx.B, ctx.nextB);
  }

  std::array<Vector, RD_N> nodes;
  Orientation<W> orientation;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  StaticCircularBuffer<GSReactionContext, 4> contexts;
  Timeline<W> timeline;

  struct Params {
    float feed = 0.0545f;
    float k = 0.062f;
    float dA = 0.15f;
    float dB = 0.075f;
    float dt = 1.0f;
    float global_alpha = 0.3f;
  } params;
};

#include "../effect_registry.h"
REGISTER_EFFECT(GSReactionDiffusion)
