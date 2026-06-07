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
 * @brief Belousov-Zhabotinsky reaction-diffusion on a Fibonacci lattice sphere.
 *
 * Three competing chemical species (A, B, C) evolve via Lotka-Volterra dynamics
 * on a 7680-node Fibonacci lattice with K=6 nearest neighbors. The cyclic
 * competition (A→B→C→A) creates self-sustaining spiral waves that persist
 * indefinitely. State is stored as Q8 (uint8_t). Per-pixel rendering uses
 * Wendland C2 kernel interpolation for smooth Voronoi boundaries.
 *
 * Shared lattice/orientation/kernel scaffolding lives in ReactionDiffusionBase.
 *
 * Memory budget (persistent arena):
 *   - State:  3 arrays × 7680 × 1B = 23,040 B
 *   - Total:  23,040 B (22 KB)
 *
 * Scratch arena (per frame):
 *   - Physics:  3 × 7680 × 1B      = 23,040 B
 *   - Node XYZ: 7680 × 12B         = 92,160 B
 *   - Total:    115,200 B (112 KB)
 */
template <int W, int H>
class BZReactionDiffusion
    : public ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::build_nodes;
  using Base::cube_lut;
  using Base::dist2;
  using Base::for_each_neighbor;
  using Base::graph_laplacian;
  using Base::init_orientation_animation;
  using Base::INV_R2;
  using Base::kernel_accumulate;
  using Base::orientation;
  using Base::RD_K;
  using Base::RD_N;
  using Base::registerParam;

public:
  FLASHMEM BZReactionDiffusion() = default;

  void init() override {
    // 75KB easily holds the 48KB Cubemap LUT + 23KB State
    configure_arenas(75 * 1024, GLOBAL_ARENA_SIZE - 75 * 1024, 0);

    registerParam("Alpha", &params.alpha, 0.0f, 4.0f);
    registerParam("Diff", &params.D, 0.001f, 0.1f);
    registerParam("Speed", &params.dt, 0.0f, 1.0f);

    allocate_state();
    cube_lut.build(persistent_arena);
    seed_spiral_nuclei();
    init_orientation_animation();
  }

private:
  // ---------------------------------------------------------------------------
  // Q8 fixed-point helpers
  // ---------------------------------------------------------------------------

  static inline float from_q8(uint8_t v) { return v * (1.0f / 255.0f); }
  static inline uint8_t to_q8(float v) {
    return static_cast<uint8_t>(hs::clamp(v, 0.0f, 1.0f) * 255.0f);
  }

  // ---------------------------------------------------------------------------
  // Initialization helpers
  // ---------------------------------------------------------------------------

  void allocate_state() {
    state.A = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    state.B = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    state.C = static_cast<uint8_t *>(persistent_arena.allocate(RD_N, 1));
    memset(state.A, 0, RD_N);
    memset(state.B, 0, RD_N);
    memset(state.C, 0, RD_N);
  }

  // ---------------------------------------------------------------------------
  // Seeding
  // ---------------------------------------------------------------------------

  /** Seed 3 clusters per species to ensure all 3 are always present. */
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

  // ---------------------------------------------------------------------------
  // Physics: Lotka-Volterra reaction + graph Laplacian diffusion
  // ---------------------------------------------------------------------------

  /** Advance one species: diffusion + Lotka-Volterra competition step. */
  uint8_t advance_species(float conc, float predator, float laplacian) const {
    return to_q8(conc + (params.D * laplacian +
                         conc * (1 - conc - params.alpha * predator)) *
                            params.dt);
  }

  /** Apply stochastic perturbation to prevent convergence on closed manifold.
   */
  static void perturb_state(uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    for (int p = 0; p < 8; p++) {
      int idx = hs::rand_int(0, RD_N - 1);
      int s = hs::rand_int(0, 2);
      uint8_t *t = (s == 0) ? nA : (s == 1) ? nB : nC;
      t[idx] =
          static_cast<uint8_t>(std::min(static_cast<int>(t[idx]) + 3, 255));
    }
  }

  /** Full physics step: reaction-diffusion + perturbation + swap. */
  void step_physics(uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    auto from = [](uint8_t v) { return from_q8(v); };
    for (int i = 0; i < RD_N; i++) {
      float a = from_q8(state.A[i]);
      float b = from_q8(state.B[i]);
      float c = from_q8(state.C[i]);

      nA[i] = advance_species(a, c, graph_laplacian(state.A, i, a, from));
      nB[i] = advance_species(b, a, graph_laplacian(state.B, i, b, from));
      nC[i] = advance_species(c, b, graph_laplacian(state.C, i, c, from));
    }

    perturb_state(nA, nB, nC);

    std::swap(state.A, nA);
    std::swap(state.B, nB);
    std::swap(state.C, nC);
  }

  // ---------------------------------------------------------------------------
  // Rendering: Wendland C2 kernel interpolation
  // ---------------------------------------------------------------------------

  /** Blend three species concentrations into a single pixel via palette. */
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

  /** Find the closest node by walking from a cubemap seed through neighbors. */
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

  /** Sample the kernel-interpolated color at a world-space point. */
  Color4 sample_kernel(const Vector &rv, const Vector *nodes, int best_node,
                       const Color4 &ca, const Color4 &cb,
                       const Color4 &cc) const {
    float tw = 0, wa = 0, wb = 0, wc = 0;
    kernel_accumulate(rv, nodes, best_node, [&](int i, float w) {
      wa += state.A[i] * w;
      wb += state.B[i] * w;
      wc += state.C[i] * w;
      tw += w;
    });

    if (tw <= 0.0001f)
      return Color4(Pixel(0, 0, 0), 0.0f);

    float inv = (1.0f / 255.0f) / tw;
    return Color4(blend_species(wa * inv, wb * inv, wc * inv, ca, cb, cc),
                  1.0f);
  }

  /** Allocate scratch, run physics steps, then rasterize with 4× SSAA. */
  void render(Canvas &canvas) {
    ScratchScope _frame(scratch_arena_a);

    uint8_t *sA = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *sB = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *sC = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));

    for (int k = 0; k < 2; k++)
      step_physics(sA, sB, sC);

    Vector *nodes = static_cast<Vector *>(
        scratch_arena_a.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    build_nodes(nodes);

    Color4 ca = palette.get(0.0f);
    Color4 cb = palette.get(0.5f);
    Color4 cc = palette.get(1.0f);

    auto vertex_shader = [&](Fragment &frag) {
      Vector rv = orientation.unorient(frag.pos);
      frag.v0 = static_cast<float>(cube_lut.lookup(rv));
    };

    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      int center_node = static_cast<int>(frag.v0);
      Vector rv = orientation.unorient(v);
      int best_node = refine_nearest_node(rv, nodes, center_node);
      frag.color = sample_kernel(rv, nodes, best_node, ca, cb, cc);
    };

    Scan::Shader::draw<W, H, 4>(canvas, fragment_shader, vertex_shader);
  }

  // ---------------------------------------------------------------------------
  // Member data
  // ---------------------------------------------------------------------------

  struct {
    uint8_t *A = nullptr, *B = nullptr, *C = nullptr;
  } state;

  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::DESCENDING,
                            SaturationProfile::VIBRANT, 42};

  struct Params {
    float alpha = 3.0f;
    float D = 0.05f;
    float dt = 0.35f;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(BZReactionDiffusion)
