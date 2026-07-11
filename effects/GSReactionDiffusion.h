/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <algorithm>
#include <cstring>
#include <utility>
#include "core/engine/engine.h"
#include "effects/ReactionDiffusionBase.h"

// Unit-test accessor (tests/test_effects.h) reaching the private Q16
// conversions and one Gray-Scott substep that the smoke/determinism harness
// cannot pin.
namespace hs_test {
namespace effects_tests {
struct GSWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Builds per-node two-ring "renderable" flags for the B field.
 * @param b Per-node B concentrations, Q16.
 * @param hot1 Scratch: per-node flag, set when any of {node, neighbors}
 *        reaches the threshold.
 * @param hot2 Output: per-node flag, set when any node within two hops
 *        reaches the threshold.
 * @param count Node count.
 * @param threshold Q16 render floor.
 * @details A kernel sample is a convex average over the refined stencil and
 * the refined center is at most one hop from the seed, so a seed whose
 * two-ring sits entirely below the floor cannot produce a renderable sample —
 * culling on !hot2[seed] is exact, not approximate. Non-template so HS_COLD
 * reliably keeps the once-per-frame pass off ITCM.
 */
[[maybe_unused]] HS_COLD static void
fill_hot_flags(const uint16_t *b, uint8_t *hot1, uint8_t *hot2, int count,
               uint16_t threshold) {
  for (int i = 0; i < count; ++i) {
    bool hot = b[i] >= threshold;
    for (int k = 0; k < ReactionGraph::RD_K && !hot; ++k)
      hot = b[ReactionGraph::neighbors[i][k]] >= threshold;
    hot1[i] = hot;
  }
  for (int i = 0; i < count; ++i) {
    bool hot = hot1[i];
    for (int k = 0; k < ReactionGraph::RD_K && !hot; ++k)
      hot = hot1[ReactionGraph::neighbors[i][k]];
    hot2[i] = hot;
  }
}

/**
 * @brief Gray-Scott reaction-diffusion on a Fibonacci lattice sphere.
 *
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 *
 * @details
 * Two species (A, B) evolve via Gray-Scott dynamics (A·B² autocatalysis with
 * feed/kill) on the shared 7680-node lattice, producing spots/stripes/mazes.
 * Persistent state is Q16 (uint16_t) for the cubic reaction-term precision;
 * substeps integrate in float and quantize back once per frame. Shared
 * lattice/orientation/kernel scaffolding lives in ReactionDiffusionBase.
 *
 * Memory budget (persistent arena, configured 174 KB):
 *   - Cubemap LUT:                  6 × 64² × 2B = 49,152 B
 *   - State:   2 arrays × 7680 × 2B (Q16)        = 30,720 B
 *   - Node XYZ: 7680 × 12B                       = 92,160 B  (fixed lattice, built once)
 *   - Palette LUT: 256 × 12B                     =  3,072 B  (the extra tenant vs BZ)
 *   - Total:                                       175,104 B (171 KB)
 *
 * Scratch arena (per frame, disjoint phases):
 *   - Physics: float ping-pong 4 × 7680 × 4B                             = 122,880 B
 *   - Raster:  oriented lattice 7680 × 12B + cull flags 2 × 7680 × 1B    = 107,520 B
 */
template <int W, int H>
class GSReactionDiffusion
    : public ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::cube_lut;
  using Base::for_each_neighbor;
  using Base::init_lattice;
  using Base::nodes;
  using Base::orientation;
  using Base::RD_K;
  using Base::RD_N;
  using Base::refine_and_accumulate;
  using Base::register_param;
  using Base::seed_face_lut;

public:
  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  FLASHMEM GSReactionDiffusion() = default;

  /**
   * @brief One-time setup: arenas, GUI params, A/B state, cubemap LUT, lattice.
   * @details Carves the persistent arena, registers the GUI params, allocates
   * and seeds the A/B state, and builds the cubemap LUT and lattice nodes once.
   */
  void init() override {
    constexpr size_t CUBE_LUT_BYTES = 6u * ReactionGraph::CubemapLUT::RES *
                                     ReactionGraph::CubemapLUT::RES *
                                     sizeof(uint16_t); // cube_lut.build
    constexpr size_t STATE_BYTES = 2u * RD_N * sizeof(uint16_t); // A + B, Q16
    constexpr size_t NODE_BYTES = RD_N * sizeof(Vector);         // build_nodes
    constexpr size_t PALETTE_BYTES =
        BakedPalette::LUT_SIZE * sizeof(Color4); // palette.bake
    constexpr size_t PERSISTENT_BYTES = 174 * 1024;
    static_assert(CUBE_LUT_BYTES + STATE_BYTES + NODE_BYTES + PALETTE_BYTES <=
                      PERSISTENT_BYTES,
                  "GS persistent arena too small for LUT + state + nodes + "
                  "palette");
    // render()'s scratch peaks at the larger of the physics phase (4 float
    // ping-pong buffers) and the raster phase (the oriented lattice + cull
    // flags); the two run under disjoint scopes.
    constexpr size_t PHYSICS_SCRATCH_BYTES = 4u * RD_N * sizeof(float);
    constexpr size_t RASTER_SCRATCH_BYTES =
        RD_N * sizeof(Vector) + 2u * RD_N * sizeof(uint8_t);
    constexpr size_t SCRATCH_BYTES =
        PHYSICS_SCRATCH_BYTES > RASTER_SCRATCH_BYTES ? PHYSICS_SCRATCH_BYTES
                                                     : RASTER_SCRATCH_BYTES;
    static_assert(SCRATCH_BYTES <= GLOBAL_ARENA_SIZE - PERSISTENT_BYTES,
                  "GS scratch arena too small for render()'s phase peak");
    configure_arenas(PERSISTENT_BYTES, GLOBAL_ARENA_SIZE - PERSISTENT_BYTES, 0);

    register_param("Feed", &params.feed, 0.0f, 0.1f);
    register_param("Kill", &params.k, 0.0f, 0.1f);
    // dA/dB cap at 0.05: explicit Euler is stable only while dt·D·λmax ≤ 2. The
    // graph Laplacian on a degree-RD_K lattice has |λ|max ≤ 2·RD_K (= 12 at
    // RD_K=6), giving 3·0.05·12 = 1.8 ≤ 2 at Speed top.
    register_param("dA", &params.d_a, 0.0f, 0.05f);
    register_param("dB", &params.d_b, 0.0f, 0.05f);
    register_param("Speed", &params.dt, 0.1f, 3.0f);

    state.A = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    state.B = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    for (int i = 0; i < RD_N; i++) {
      state.A[i] = 65535;
      state.B[i] = 0;
    }

    palette.bake(persistent_arena,
                 GenerativePalette{GradientShape::STRAIGHT,
                                   HarmonyType::SPLIT_COMPLEMENTARY,
                                   BrightnessProfile::ASCENDING,
                                   SaturationProfile::VIBRANT});

    cube_lut.build(persistent_arena);
    init_lattice();
    seed_clusters();
  }

private:
  // Test seam: lets unit tests reach the Q16 helpers, step_physics, and params
  // without exposing them to production callers.
  friend struct ::hs_test::effects_tests::GSWhiteBox;

  /**
   * @brief Q16 full-scale factor: maps 1.0 to 65535 for the 16-bit fixed-point
   * state needed by the GS cubic reaction term.
   */
  static constexpr float Q16_SCALE = 65535.0f;
  static constexpr float Q16_INV = 1.0f / Q16_SCALE; /**< Reciprocal of
                                                          Q16_SCALE. */

  static constexpr int NUM_SEED_CLUSTERS = 30; /**< Initial B blobs seeded at
                                                    init. */
  /**
   * @brief Physics substeps advanced per rendered frame.
   * @details 16 substeps/frame advance the slow GS morphogenesis at a visible
   * rate; each stays within the explicit-Euler stability bound (see "Speed").
   */
  static constexpr int STEPS_PER_FRAME = 16;
  /**
   * @brief Lower bound of the B render band: below this, pixels are transparent;
   * [B_COLOR_FLOOR, B_COLOR_FLOOR + 1/B_COLOR_SCALE] maps to the full palette
   * range [0,1].
   */
  static constexpr float B_COLOR_FLOOR = 0.1f;
  static constexpr float B_COLOR_SCALE = 4.0f; /**< Slope mapping B above the
                                                    floor into palette t. */
  /**
   * @brief Cull threshold; coincides with the color floor so there is no band
   * between the two: a pixel is either fully transparent or on the gradient.
   * @details A cull below the floor would map b in [cull, floor) to t==0,
   * rendering an opaque flat plateau of the lowest palette color.
   */
  static constexpr float B_CULL_THRESHOLD = B_COLOR_FLOOR;
  /**
   * @brief Converts a Q16 fixed-point value to a float in [0, 1].
   * @param v Q16 value, where 65535 represents 1.0.
   * @return Concentration as a float in [0, 1].
   */
  static inline float from_q16(uint16_t v) { return v * Q16_INV; }
  /**
   * @brief Converts a float concentration to a clamped, rounded Q16 value.
   * @param v Concentration as a float (clamped to [0, 1]).
   * @return Q16 value in [0, 65535].
   * @details Rounds to nearest (+0.5f); truncating would bias the dynamics down
   * by dropping sub-LSB positive updates. clamp bounds the input so 65535.5 ->
   * 65535 with no overflow.
   */
  static inline uint16_t to_q16(float v) {
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 1.0f) * Q16_SCALE + 0.5f);
  }

  /**
   * @brief Seeds NUM_SEED_CLUSTERS fully-saturated B blobs as nucleation sites.
   * @details Drops blobs of fully-saturated B (a random node plus its immediate
   * neighbors) so the otherwise-uniform A=1/B=0 field has nucleation sites for
   * the GS instability to grow from.
   */
  HS_COLD_MEMBER void seed_clusters() {
    for (int i = 0; i < NUM_SEED_CLUSTERS; i++) {
      int idx = hs::rand_int(0, RD_N);
      state.B[idx] = 65535;
      for_each_neighbor(idx, [&](int nb) { state.B[nb] = 65535; });
    }
  }

  /**
   * @brief Advances one Gray-Scott substep into the next buffers (Jacobi).
   * @param c_a Current A field (read-only), float in [0, 1] per node.
   * @param c_b Current B field (read-only), float in [0, 1] per node.
   * @param n_a Next A field (write target), float in [0, 1] per node.
   * @param n_b Next B field (write target), float in [0, 1] per node.
   * @details Gray-Scott: dA/dt = dA·∇²A - A·B² + feed·(1-A);
   * dB/dt = dB·∇²B + A·B² - (k+feed)·B. Double-buffered Jacobi: reads current
   * buffers, writes next; the caller owns the ping-pong (see render()). The
   * [0, 1] clamp saturates explicit-Euler overshoot past the stability bound
   * (see "Speed"); substeps stay in float so the Q16 state quantizes once per
   * frame, not once per substep.
   */
  void step_physics(const float *c_a, const float *c_b, float *n_a,
                    float *n_b) {
    for (int i = 0; i < RD_N; i++) {
      float a = c_a[i];
      float b = c_b[i];

      float sum_a = 0, sum_b = 0;
      for_each_neighbor(i, [&](int ni) {
        sum_a += c_a[ni];
        sum_b += c_b[ni];
      });
      float l_a = sum_a - RD_K * a;
      float l_b = sum_b - RD_K * b;

      float abb = a * b * b;
      n_a[i] = hs::clamp(
          a + (params.d_a * l_a - abb + params.feed * (1.0f - a)) * params.dt,
          0.0f, 1.0f);
      n_b[i] = hs::clamp(
          b + (params.d_b * l_b + abb - (params.k + params.feed) * b) *
                  params.dt,
          0.0f, 1.0f);
    }
  }

  /**
   * @brief Kernel-weighted sample of the B concentration at a point.
   * @param p Query point on the sphere.
   * @param seed Seed node id from the cubemap LUT, refined to the true nearest
   * inside the fused stencil walk.
   * @param nodes Node positions in the same frame as `p`.
   * @return Support-radius weighted average of B in [0, 1]; 0 if no node is
   * within the support radius.
   */
  float interpolate_b(const Vector &p, int seed, const Vector *nodes) const {
    float tw = 0, wb = 0;
    refine_and_accumulate(p, nodes, seed, [&](int i, float w) {
      wb += from_q16(state.B[i]) * w;
      tw += w;
    });
    // Guard the division: returning 0 (not a 0/0 NaN) keeps the value cullable
    // by render()'s `b < B_CULL_THRESHOLD` test downstream — a NaN compares false
    // there, slips through, and poisons palette.get().
    if (tw <= Base::KERNEL_MIN_TOTAL_WEIGHT)
      return 0.0f;
    return wb / tw;
  }

  /**
   * @brief Advances the sim STEPS_PER_FRAME substeps and rasterizes the B field.
   * @param canvas Destination canvas to draw the sphere into.
   * @details Rasterizes the B field onto the sphere via the orientation-aware
   * SSAA shader pipeline after advancing the simulation.
   */
  void render(Canvas &canvas) {
    ScratchScope frame_guard(scratch_arena_a);
    {
      // The substeps ping-pong in float; the Q16 state converts in once and
      // quantizes back once per frame, not once per substep.
      ScratchScope physics_guard(scratch_arena_a);
      float *cur_a = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *cur_b = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *nxt_a = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *nxt_b = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));

      for (int i = 0; i < RD_N; i++) {
        cur_a[i] = from_q16(state.A[i]);
        cur_b[i] = from_q16(state.B[i]);
      }
      for (int k = 0; k < STEPS_PER_FRAME; k++) {
        step_physics(cur_a, cur_b, nxt_a, nxt_b);
        std::swap(cur_a, nxt_a);
        std::swap(cur_b, nxt_b);
      }
      for (int i = 0; i < RD_N; i++) {
        state.A[i] = to_q16(cur_a[i]);
        state.B[i] = to_q16(cur_b[i]);
      }
    }

    // Physics scratch is popped; the raster phase reuses the arena for the
    // oriented lattice so the kernel walks stay in world space, plus the
    // two-ring cull flags.
    Vector *world_nodes = static_cast<Vector *>(
        scratch_arena_a.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    orient_nodes(nodes, world_nodes, RD_N, orientation.get());
    uint8_t *hot1 = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *hot2 = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    fill_hot_flags(state.B, hot1, hot2, RD_N, to_q16(B_CULL_THRESHOLD));

    // Seed the cubemap lookup once per pixel center; a seed whose two-ring
    // sits below the render floor is culled for the whole pixel (v0 = -1).
    auto vertex_shader = [&](Fragment &frag) {
      seed_face_lut(frag);
      if (!hot2[static_cast<int>(frag.v0)])
        frag.v0 = -1.0f;
    };

    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      int seed = static_cast<int>(frag.v0);
      if (seed < 0) {
        frag.color = Color4(Pixel(0, 0, 0), 0.0f);
        return;
      }
      float b = interpolate_b(v, seed, world_nodes);

      if (b < B_CULL_THRESHOLD) {
        frag.color = Color4(Pixel(0, 0, 0), 0.0f);
        return;
      }

      float t = hs::clamp((b - B_COLOR_FLOOR) * B_COLOR_SCALE, 0.0f, 1.0f);
      frag.color = palette.get(t);
    };

    Scan::Shader::draw<W, H, 4>(canvas, fragment_shader, vertex_shader);
  }

  /**
   * @brief Persistent Q16 state buffers for the two species.
   */
  struct {
    uint16_t *A = nullptr, *B = nullptr; /**< Per-node A/B concentrations, Q16. */
  } state;

  /** @brief 16-bit LUT baked from the generative palette mapping B to RGB. */
  BakedPalette palette;

  /**
   * @brief GUI-tunable Gray-Scott parameters.
   */
  struct Params {
    float feed = 0.04f; /**< Feed rate of A. */
    float k = 0.06f;    /**< Kill rate of B. */
    float d_a = 0.02f;   /**< Diffusion coefficient of A. */
    float d_b = 0.01f;   /**< Diffusion coefficient of B. */
    float dt = 2.5f;    /**< Integration timestep (Speed slider). */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(GSReactionDiffusion)
