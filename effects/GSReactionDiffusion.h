/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <algorithm>
#include <cstring>
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
 * @brief Gray-Scott reaction-diffusion on a Fibonacci lattice sphere.
 *
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 *
 * @details
 * Two species (A, B) evolve via Gray-Scott dynamics (A·B² autocatalysis with
 * feed/kill) on the shared 7680-node lattice, producing spots/stripes/mazes.
 * State is Q16 (uint16_t) for the cubic reaction-term precision. Shared
 * lattice/orientation/kernel scaffolding lives in ReactionDiffusionBase.
 *
 * Memory budget (persistent arena, configured 174 KB):
 *   - Cubemap LUT:                  6 × 64² × 2B = 49,152 B
 *   - State:   2 arrays × 7680 × 2B (Q16)        = 30,720 B
 *   - Node XYZ: 7680 × 12B                       = 92,160 B  (fixed lattice, built once)
 *   - Palette LUT: 256 × 12B                     =  3,072 B  (the extra tenant vs BZ)
 *   - Total:                                       175,104 B (171 KB)
 *
 * Scratch arena (per frame, STEPS_PER_FRAME = 16):
 *   - Physics ping-pong: 2 × 7680 × 2B   = 30,720 B
 *   - Float pre-convert: 2 × 7680 × 4B   = 61,440 B
 *   - Total:    92,160 B (90 KB)
 */
template <int W, int H>
class GSReactionDiffusion
    : public ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::advance_substeps;
  using Base::cube_lut;
  using Base::for_each_neighbor;
  using Base::init_lattice;
  using Base::kernel_accumulate;
  using Base::nodes;
  using Base::orientation;
  using Base::RD_N;
  using Base::refine_nearest_node;
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
    // render() carves 2 uint16 + 2 float species buffers from scratch_a.
    constexpr size_t SCRATCH_BYTES =
        2u * RD_N * sizeof(uint16_t) + 2u * RD_N * sizeof(float);
    static_assert(SCRATCH_BYTES <= GLOBAL_ARENA_SIZE - PERSISTENT_BYTES,
                  "GS scratch arena too small for render()'s species buffers");
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
   * @details GS morphogenesis is slow per step, so 16 substeps/frame are needed
   * for the pattern to grow at a visible rate. Each substep stays inside the
   * explicit-Euler stability bound (see the "Speed" param), so the higher count
   * is pure throughput, not a stability risk.
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
   * @details Rounds to nearest (+0.5f), not truncate: plain truncation drops
   * every sub-LSB positive update while still applying negative ones, biasing
   * the RD dynamics downward into a diffusion dead-zone. clamp keeps the product
   * in [0, 65535], so +0.5f tops out at 65535.5 -> 65535 with no overflow.
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
  void seed_clusters() {
    for (int i = 0; i < NUM_SEED_CLUSTERS; i++) {
      int idx = hs::rand_int(0, RD_N);
      state.B[idx] = 65535;
      for_each_neighbor(idx, [&](int nb) { state.B[nb] = 65535; });
    }
  }

  /**
   * @brief Advances one Gray-Scott substep into the next buffers (Jacobi).
   * @param c_a Current A field (read-only), Q16 per node.
   * @param c_b Current B field (read-only), Q16 per node.
   * @param n_a Next A field (write target), Q16 per node.
   * @param n_b Next B field (write target), Q16 per node.
   * @param f_a Float scratch (RD_N) for the current A generation.
   * @param f_b Float scratch (RD_N) for the current B generation.
   * @details Gray-Scott: dA/dt = dA·∇²A - A·B² + feed·(1-A);
   * dB/dt = dB·∇²B + A·B² - (k+feed)·B. Pure double-buffered (Jacobi) step:
   * reads the current buffers, writes the next ones. The caller owns the
   * ping-pong so the result can be landed back in the persistent state
   * regardless of substep parity (see render()). The current generation is
   * pre-converted into f_a/f_b once, so the neighbor loop reads floats instead of
   * reconverting each node's Q16 value on every neighbor visit (~6-7x per node).
   */
  void step_physics(const uint16_t *c_a, const uint16_t *c_b, uint16_t *n_a,
                    uint16_t *n_b, float *f_a, float *f_b) {
    for (int i = 0; i < RD_N; i++) {
      f_a[i] = from_q16(c_a[i]);
      f_b[i] = from_q16(c_b[i]);
    }
    for (int i = 0; i < RD_N; i++) {
      float a = f_a[i];
      float b = f_b[i];

      float l_a = 0, l_b = 0;
      for_each_neighbor(i, [&](int ni) {
        l_a += f_a[ni] - a;
        l_b += f_b[ni] - b;
      });

      float abb = a * b * b;
      n_a[i] = to_q16(a + (params.d_a * l_a - abb + params.feed * (1.0f - a)) *
                             params.dt);
      n_b[i] = to_q16(b + (params.d_b * l_b + abb - (params.k + params.feed) * b) *
                             params.dt);
    }
  }

  /**
   * @brief Kernel-weighted sample of the B concentration at a point.
   * @param p Query point on the sphere (unoriented lattice space).
   * @param nearest Precomputed index of the nearest lattice node to seed from.
   * @param nodes Fixed lattice node-position array.
   * @return Support-radius weighted average of B in [0, 1]; 0 if no node is
   * within the support radius.
   */
  float interpolate_b(const Vector &p, int nearest, const Vector *nodes) const {
    float tw = 0, wb = 0;
    kernel_accumulate(p, nodes, nearest, [&](int i, float w) {
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
    uint16_t *s_a = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    uint16_t *s_b = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    float *f_a = static_cast<float *>(
        scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
    float *f_b = static_cast<float *>(
        scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));

    advance_substeps(STEPS_PER_FRAME,
                     std::array<uint16_t *, 2>{state.A, state.B},
                     std::array<uint16_t *, 2>{s_a, s_b},
                     [&](auto &cur, auto &nxt) {
                       step_physics(cur[0], cur[1], nxt[0], nxt[1], f_a, f_b);
                     });

    // Seed the cubemap lookup once per pixel center; it only feeds
    // refine_nearest_node, which re-finds the true nearest per sub-sample.
    auto vertex_shader = [this](Fragment &frag) { seed_face_lut(frag); };

    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      Vector rv = orientation.unorient(v);
      // The CubemapLUT seed is quantized to a face cell; refine to the true
      // nearest node.
      int nearest =
          refine_nearest_node(rv, nodes, static_cast<int>(frag.v0));
      float b = interpolate_b(rv, nearest, nodes);

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
