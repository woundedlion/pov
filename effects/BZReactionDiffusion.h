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
 *   - Cubemap LUT:                 ~49,152 B
 *   - State:  3 arrays × 7680 × 1B = 23,040 B
 *   - Node XYZ: 7680 × 12B         = 92,160 B  (fixed lattice, built once)
 *
 * Scratch arena (per frame):
 *   - Physics:  3 × 7680 × 1B      = 23,040 B
 *   - Total:    23,040 B (22 KB)
 */
template <int W, int H>
class BZReactionDiffusion
    : public ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::build_nodes;
  using Base::cube_lut;
  using Base::for_each_neighbor;
  using Base::init_orientation_animation;
  using Base::kernel_accumulate;
  using Base::orientation;
  using Base::RD_N;
  using Base::refine_nearest_node;
  using Base::registerParam;

public:
  FLASHMEM BZReactionDiffusion() = default;

  void init() override {
    // 165KB holds the 48KB Cubemap LUT + 23KB State + 90KB node positions. The
    // node array (7680 × Vector) is the fixed Fibonacci lattice — it does not
    // depend on the per-frame view orientation (queries are un-oriented onto it,
    // not the reverse), so it is built ONCE here instead of every frame. This is
    // not a net RAM increase: the array was already allocated per-frame in the
    // scratch arena, so the peak footprint is unchanged — it just moves from
    // transient scratch to persistent and stops being recomputed each frame.
    configure_arenas(165 * 1024, GLOBAL_ARENA_SIZE - 165 * 1024, 0);

    // "Compete", not "Alpha": this is the Lotka-Volterra predation coefficient,
    // not an opacity. Naming it "Alpha" collided with the engine-wide
    // Alpha-as-opacity convention every other effect follows.
    registerParam("Compete", &params.alpha, 0.0f, 4.0f);
    registerParam("Diff", &params.D, 0.001f, 0.1f);
    registerParam("Speed", &params.dt, 0.0f, 1.0f);

    allocate_state();
    cube_lut.build(persistent_arena);
    nodes = static_cast<Vector *>(
        persistent_arena.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    build_nodes(nodes);
    seed_spiral_nuclei();
    init_orientation_animation();
  }

private:
  // ---------------------------------------------------------------------------
  // Q8 fixed-point helpers
  // ---------------------------------------------------------------------------

  static inline float from_q8(uint8_t v) { return v * (1.0f / 255.0f); }
  static inline uint8_t to_q8(float v) {
    // Round to nearest (+0.5f), not truncate: truncation loses every sub-LSB
    // positive update while negatives still decrement, biasing the dynamics
    // downward. clamp bounds the product to [0, 255], so +0.5f tops out at
    // 255.5 -> 255 with no overflow.
    return static_cast<uint8_t>(hs::clamp(v, 0.0f, 1.0f) * 255.0f + 0.5f);
  }

  // Simulation tuning.
  static constexpr int CLUSTERS_PER_SPECIES = 3; // seed blobs per species at init
  static constexpr int STEPS_PER_FRAME = 2;      // physics substeps per render
  static constexpr int NUM_PERTURBATIONS = 8;    // random nudges per physics step
  static constexpr int PERTURB_AMOUNT = 3;       // Q8 magnitude of each nudge

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

  /** Seed CLUSTERS_PER_SPECIES clusters per species to ensure all 3 are always
   *  present. */
  void seed_spiral_nuclei() {
    uint8_t *species[] = {state.A, state.B, state.C};
    for (int s = 0; s < 3; s++) {
      for (int k = 0; k < CLUSTERS_PER_SPECIES; k++) {
        int center = hs::rand_int(0, RD_N);
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
    for (int p = 0; p < NUM_PERTURBATIONS; p++) {
      int idx = hs::rand_int(0, RD_N);
      int s = hs::rand_int(0, 3);  // half-open [0,3): all three species, incl. C
      uint8_t *t = (s == 0) ? nA : (s == 1) ? nB : nC;
      t[idx] = static_cast<uint8_t>(
          std::min(static_cast<int>(t[idx]) + PERTURB_AMOUNT, 255));
    }
  }

  /** Full physics step: reaction-diffusion + perturbation.
   *
   * Pure double-buffered (Jacobi): reads the current buffers, writes the next
   * ones. The caller owns the ping-pong so the result can be landed back in the
   * persistent state regardless of substep parity (see render()). */
  void step_physics(const uint8_t *cA, const uint8_t *cB, const uint8_t *cC,
                    uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    for (int i = 0; i < RD_N; i++) {
      float a = from_q8(cA[i]);
      float b = from_q8(cB[i]);
      float c = from_q8(cC[i]);

      // All three species' Laplacians share one neighbor walk (fused on
      // purpose, per for_each_neighbor's note: three single-field walks would
      // triple the lattice reads on this dominant loop). Each Σ runs over the
      // same neighbor order, so the result is identical to the split form.
      float lA = 0, lB = 0, lC = 0;
      for_each_neighbor(i, [&](int ni) {
        lA += from_q8(cA[ni]) - a;
        lB += from_q8(cB[ni]) - b;
        lC += from_q8(cC[ni]) - c;
      });

      nA[i] = advance_species(a, c, lA);
      nB[i] = advance_species(b, a, lB);
      nC[i] = advance_species(c, b, lC);
    }

    perturb_state(nA, nB, nC);
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

    // Ping-pong between the persistent state and the scratch buffers. cur always
    // names the latest generation; nxt is the write target for the next step.
    uint8_t *curA = state.A, *curB = state.B, *curC = state.C;
    uint8_t *nxtA = sA, *nxtB = sB, *nxtC = sC;
    for (int k = 0; k < STEPS_PER_FRAME; k++) {
      step_physics(curA, curB, curC, nxtA, nxtB, nxtC);
      std::swap(curA, nxtA);
      std::swap(curB, nxtB);
      std::swap(curC, nxtC);
    }

    // Land the final generation back in the persistent buffers so it survives
    // the ScratchScope pop. An even substep count already leaves cur == state
    // (no copy); an odd count leaves it in scratch and is copied back — making
    // the effect correct for any STEPS_PER_FRAME, not just even values.
    if (curA != state.A) {
      std::memcpy(state.A, curA, RD_N);
      std::memcpy(state.B, curB, RD_N);
      std::memcpy(state.C, curC, RD_N);
    }

    // `nodes` is the fixed lattice, built once in init() — not rebuilt here.
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

  // Fixed Fibonacci-lattice node positions, built once in init() (persistent).
  Vector *nodes = nullptr;

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
