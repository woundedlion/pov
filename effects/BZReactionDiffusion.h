/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include "../color.h"
#include "../geometry.h"

template <int W>
class BZReactionDiffusion : public Effect {
public:

  static constexpr int RD_H = 20;   // Height of the grid
  static constexpr int RD_N = W * RD_H * 2; // Number of nodes in the graph
  static constexpr int RD_K = 6;    // Number of neighbors per node

  BZReactionDiffusion() :
    Effect(W),
    filters(Filter::World::Orient<W>(orientation), Filter::Screen::AntiAlias<W>())
  {
    persist_pixels = false;
    build_graph();
    timeline
      .add(0, Rotation<W>(orientation, Y_AXIS, PI_F / 2, 600, ease_mid, true))
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
   * @brief Holds the chemical state and parameters for a single Belousov-Zhabotinsky reaction instance.
   * @details Simulates a 3-species cyclic competition (Rock-Paper-Scissors).
   */
  struct BZReactionContext {
    std::array<float, RD_N> A;
    std::array<float, RD_N> B;
    std::array<float, RD_N> C;
    std::array<float, RD_N> nextA;
    std::array<float, RD_N> nextB;
    std::array<float, RD_N> nextC;

    // Simulation Parameters
    float alpha = 1.6f; // Predation rate
    float D = 0.03f;    // Diffusion rate
    float dt = 0.2f;    // Time step

    GenerativePalette palette;

    BZReactionContext() :
      palette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::DESCENDING, SaturationProfile::VIBRANT)
    {
      reset();
    }

    void reset() {
      A.fill(0.0f);
      B.fill(0.0f);
      C.fill(0.0f);
      palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::DESCENDING, SaturationProfile::VIBRANT);
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

    // Build Neighbors using Spatial Hashing (Grid Optimization)
    static constexpr int GRID_SIZE = 20;
    static constexpr float CELL_SIZE = 2.0f / GRID_SIZE; // Domain [-1, 1], size 2.0
    
    // Temporary grid structure
    std::vector<std::vector<int>> grid(GRID_SIZE * GRID_SIZE * GRID_SIZE);

    auto get_grid_idx = [&](const Vector& p) {
      int gx = std::clamp(static_cast<int>((p.i + 1.0f) / CELL_SIZE), 0, GRID_SIZE - 1);
      int gy = std::clamp(static_cast<int>((p.j + 1.0f) / CELL_SIZE), 0, GRID_SIZE - 1);
      int gz = std::clamp(static_cast<int>((p.k + 1.0f) / CELL_SIZE), 0, GRID_SIZE - 1);
      return gx + gy * GRID_SIZE + gz * GRID_SIZE * GRID_SIZE;
    };

    // Bin points
    for (int i = 0; i < RD_N; i++) {
      grid[get_grid_idx(nodes[i])].push_back(i);
    }

    // Neighbor search
    for (int i = 0; i < RD_N; i++) {
      const Vector& p1 = nodes[i];
      
      std::array<std::pair<float, int>, RD_K> best;
      best.fill({ std::numeric_limits<float>::max(), -1 });

      int gx = std::clamp(static_cast<int>((p1.i + 1.0f) / CELL_SIZE), 0, GRID_SIZE - 1);
      int gy = std::clamp(static_cast<int>((p1.j + 1.0f) / CELL_SIZE), 0, GRID_SIZE - 1);
      int gz = std::clamp(static_cast<int>((p1.k + 1.0f) / CELL_SIZE), 0, GRID_SIZE - 1);

      for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
          for (int z = -1; z <= 1; z++) {
            int nx = gx + x;
            int ny = gy + y;
            int nz = gz + z;

             if (nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE && nz >= 0 && nz < GRID_SIZE) {
               int cell_idx = nx + ny * GRID_SIZE + nz * GRID_SIZE * GRID_SIZE;
               const auto& cell_points = grid[cell_idx];

               for (int j : cell_points) {
                 if (i == j) continue;

                 float d2 = distance_squared(p1, nodes[j]);

                 if (d2 < best[RD_K - 1].first) {
                   int pos = RD_K - 1;
                   while (pos > 0 && d2 < best[pos - 1].first) {
                     best[pos] = best[pos - 1];
                     pos--;
                   }
                   best[pos] = { d2, j };
                 }
               }
            }
          }
        }
      }

      for (int k = 0; k < RD_K; k++) {
        neighbors[i][k] = best[k].second;
      }
    }
  }

  void spawn() {
    contexts.push_back(BZReactionContext());
    BZReactionContext& ctx = contexts.back();
    seed(ctx);

    timeline.add(0, Sprite(
      [this, &ctx](Canvas& c, float opacity) {
        this->render_reaction(c, ctx, opacity);
      },
      192, 32, ease_mid, 32, ease_mid
    ));
  }

  void seed(BZReactionContext& ctx) {
    // Sparse Seeding for Spirals (Droplets)
    for (int k = 0; k < 50; k++) {
      int center = hs::rand_int(0, RD_N - 1);
      float r = hs::rand_f();

      // Determine which species to seed
      float* target = (r < 0.33f) ? ctx.A.data() : (r < 0.66f) ? ctx.B.data() : ctx.C.data();

      target[center] = 1.0f;
      for (int neighbor : neighbors[center]) {
        if (neighbor >= 0) target[neighbor] = 1.0f;
      }
    }
  }

  void render_reaction(Canvas& canvas, BZReactionContext& ctx, float opacity) {
    // Simulate Physics (12 steps per frame)
    for (int k = 0; k < 12; k++) {
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

        Pixel color(
          static_cast<uint8_t>(std::clamp(r, 0.0f, 255.0f)),
          static_cast<uint8_t>(std::clamp(g, 0.0f, 255.0f)),
          static_cast<uint8_t>(std::clamp(b_ch, 0.0f, 255.0f))
        );

        Color4 c_final(color);
        c_final.alpha *= global_alpha * opacity;
						auto shader = [c_final](const Vector& p, const Fragment& f) -> Fragment {
                            Fragment out = f;
                            out.color = c_final;
                            return out;
                        };
                        Plot::Point::draw(filters, canvas, nodes[i], shader);
      }
    }
  }

  void update_physics(BZReactionContext& ctx) {
    for (int i = 0; i < RD_N; i++) {
      float a = ctx.A[i];
      float b = ctx.B[i];
      float c = ctx.C[i];

      float lapA = 0.0f;
      float lapB = 0.0f;
      float lapC = 0.0f;

      for (int k = 0; k < RD_K; k++) {
        int n_idx = neighbors[i][k];
        if (n_idx == -1) continue;

        lapA += (ctx.A[n_idx] - a);
        lapB += (ctx.B[n_idx] - b);
        lapC += (ctx.C[n_idx] - c);
      }

      // 3-Species Cyclic Model
      float da = a * (1.0f - a - ctx.alpha * c);
      float db = b * (1.0f - b - ctx.alpha * a);
      float dc = c * (1.0f - c - ctx.alpha * b);

      ctx.nextA[i] = std::clamp(a + (ctx.D * lapA + da) * ctx.dt, 0.0f, 1.0f);
      ctx.nextB[i] = std::clamp(b + (ctx.D * lapB + db) * ctx.dt, 0.0f, 1.0f);
      ctx.nextC[i] = std::clamp(c + (ctx.D * lapC + dc) * ctx.dt, 0.0f, 1.0f);
    }

    // Swap buffers
    std::swap(ctx.A, ctx.nextA);
    std::swap(ctx.B, ctx.nextB);
    std::swap(ctx.C, ctx.nextC);
  }

  std::array<Vector, RD_N> nodes;
  std::array<std::array<int, RD_K>, RD_N> neighbors;
  Orientation orientation;

  Pipeline<W,
    Filter::World::Orient<W>,
    Filter::Screen::AntiAlias<W>> filters;

  StaticCircularBuffer<BZReactionContext, 4> contexts;
  Timeline timeline;

  float global_alpha = 0.3f;
};