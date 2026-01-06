/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include "../color.h"
#include "../geometry.h"

template <int W>
class GSReactionDiffusion : public Effect {
public:

  static constexpr int RD_N = 1024; // Number of nodes in the graph
  static constexpr int RD_K = 6; // Number of neighbors per node

  GSReactionDiffusion() :
    Effect(W),
    filters(FilterOrient<W>(orientation), FilterAntiAlias<W>())
  {
    persist_pixels = false;
    build_graph();
    timeline
      .add(0, Rotation<W>(orientation, Y_AXIS, PI_F / 2, 64, ease_mid, true))
      .add(0, PeriodicTimer(96, [this](Canvas& c) { this->spawn(); }, true));
    spawn();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:

  /**
   * @brief Holds the chemical state and parameters for a single Gray-Scott reaction instance.
   */
  struct GSReactionContext {
    std::array<float, GSReactionDiffusion::RD_N> A;
    std::array<float, RD_N> B;
    std::array<float, RD_N> nextA;
    std::array<float, RD_N> nextB;

    // Simulation Parameters
    float feed = 0.0545f;
    float k = 0.062f;
    float dA = 0.15f;
    float dB = 0.075f;
    float dt = 1.0f;

    GenerativePalette palette;

    GSReactionContext() :
      palette(GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::ASCENDING, SaturationProfile::VIBRANT)
    {
      reset();
    }

    void reset() {
      A.fill(1.0f);
      B.fill(0.0f);
      // Regenerate palette for variety on every spawn
      palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::ASCENDING, SaturationProfile::VIBRANT);
    }
  };

  void build_graph() {
    // Generate Nodes (Fibonacci Lattice)
    const float phi = PI_F * (3.0f - sqrtf(5.0f));

    for (int i = 0; i < RD_N; i++) {
      float y = 1.0f - (static_cast<float>(i) / (RD_N - 1)) * 2.0f;
      float radius = sqrtf(1.0f - y * y);
      float theta = phi * i;

      nodes[i] = Vector(
        cosf(theta) * radius,
        y,
        sinf(theta) * radius
      );
    }

    // Build Neighbors (Brute force K-NN is acceptable for run-once setup)
    for (int i = 0; i < RD_N; i++) {
      const Vector& p1 = nodes[i];

      // Track K nearest neighbors: {distance_squared, index}
      std::array<std::pair<float, int>, RD_K> best;
      best.fill({ std::numeric_limits<float>::max(), -1 });

      for (int j = 0; j < RD_N; j++) {
        if (i == j) continue;

        float d2 = distance_squared(p1, nodes[j]);

        // Insertion sort into fixed size array
        if (d2 < best[RD_K - 1].first) {
          int pos = RD_K - 1;
          while (pos > 0 && d2 < best[pos - 1].first) {
            best[pos] = best[pos - 1];
            pos--;
          }
          best[pos] = { d2, j };
        }
      }

      // Store indices
      for (int k = 0; k < RD_K; k++) {
        neighbors[i][k] = best[k].second;
      }
    }
  }

  void spawn() {
    contexts.push_back(GSReactionContext());
    GSReactionContext& ctx = contexts.back();
    seed(ctx);

    timeline.add(0, Sprite(
      [this, &ctx](Canvas& c, float opacity) {
        this->render_reaction(c, ctx, opacity);
      },
      192, 32, ease_mid, 32, ease_mid
    ));
  }

  void seed(GSReactionContext& ctx) {
    for (int i = 0; i < 5; i++) {
      int idx = hs::rand_int(0, RD_N);
      ctx.B[idx] = 1.0f;
      for (int neighbor : neighbors[idx]) {
        if (neighbor >= 0) ctx.B[neighbor] = 1.0f;
      }
    }
  }

  void render_reaction(Canvas& canvas, GSReactionContext& ctx, float opacity) {
    // Simulate Physics (12 steps per frame for stability)
    for (int k = 0; k < 12; k++) {
      update_physics(ctx);
    }

    // Draw
    for (int i = 0; i < RD_N; i++) {
      float b = ctx.B[i];
      if (b > 0.1f) {
        float t = std::clamp((b - 0.15f) * 4.0f, 0.0f, 1.0f);
        Color4 c = ctx.palette.get(t);
        filters.plot(canvas, nodes[i], c.color, 0, global_alpha * opacity * c.alpha);
      }
    }
  }

  void update_physics(GSReactionContext& ctx) {
    for (int i = 0; i < RD_N; i++) {
      float a = ctx.A[i];
      float b = ctx.B[i];

      // Laplacian (Graph-based)
      float lapA = 0.0f;
      float lapB = 0.0f;

      // Unrolled loop for fixed K=6
      for (int k = 0; k < RD_K; k++) {
        int n_idx = neighbors[i][k];
        if (n_idx == -1) continue;

        // Weight is 1.0 for all neighbors in this uniform graph
        lapA += (ctx.A[n_idx] - a);
        lapB += (ctx.B[n_idx] - b);
      }

      // Reaction-Diffusion logic
      float reaction = a * b * b;
      float feed_term = ctx.feed * (1.0f - a);
      float kill_term = (ctx.k + ctx.feed) * b;

      ctx.nextA[i] = std::clamp(a + (ctx.dA * lapA - reaction + feed_term) * ctx.dt, 0.0f, 1.0f);
      ctx.nextB[i] = std::clamp(b + (ctx.dB * lapB + reaction - kill_term) * ctx.dt, 0.0f, 1.0f);
    }

    // Swap buffers
    std::swap(ctx.A, ctx.nextA);
    std::swap(ctx.B, ctx.nextB);
  }

  std::array<Vector, RD_N> nodes;
  std::array<std::array<int, RD_K>, RD_N> neighbors;
  Orientation orientation;

  Pipeline<W,
    FilterOrient<W>,
    FilterAntiAlias<W>> filters;

  StaticCircularBuffer<GSReactionContext, 4> contexts;
  Timeline timeline;
  float global_alpha = 0.3f;
};