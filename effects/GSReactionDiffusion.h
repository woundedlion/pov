/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <algorithm>
#include <cstring>
#include "core/engine.h"
#include "effects/ReactionDiffusionBase.h"

// Forward declaration of the unit-test accessor (tests/test_effects.h) that
// reaches the private Q16 conversions and one Gray-Scott substep. The
// smoke/determinism harness cannot see a fixed-point round-trip error or a
// numerically-unstable-but-deterministic blow-up, so the dynamics are pinned
// directly through this seam.
namespace hs_test {
namespace effects_tests {
struct GSWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Gray-Scott reaction-diffusion on a Fibonacci lattice sphere.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Two species (A, B) evolve via Gray-Scott dynamics (A·B²
 * autocatalysis with feed/kill) on the shared 7680-node lattice, producing
 * spots/stripes/mazes. State is Q16 (uint16_t) for the cubic reaction-term
 * precision. Shared lattice/orientation/kernel scaffolding lives in
 * ReactionDiffusionBase.
 */
template <int W, int H>
class GSReactionDiffusion
    : public ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H>;
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
    // 170KB holds the 48KB Cubemap LUT + 30KB State + 90KB node positions. The
    // node array (7680 × Vector) is the fixed Fibonacci lattice — queries are
    // un-oriented onto it (not the reverse), so it does not depend on the
    // per-frame view orientation and is built ONCE here instead of every frame.
    //
    // Derive the requirement from sizeof so a future type-size change (the Q16
    // state element, Vector, RD_N, or the LUT resolution) trips a compile-time
    // static_assert here rather than the runtime init trap on an undersize.
    constexpr size_t kCubeLutBytes = 6u * ReactionGraph::CubemapLUT::RES *
                                     ReactionGraph::CubemapLUT::RES *
                                     sizeof(uint16_t); // cube_lut.build
    constexpr size_t kStateBytes = 2u * RD_N * sizeof(uint16_t); // A + B, Q16
    constexpr size_t kNodeBytes = RD_N * sizeof(Vector);         // build_nodes
    constexpr size_t kPersistentBytes = 170 * 1024;
    static_assert(kCubeLutBytes + kStateBytes + kNodeBytes <= kPersistentBytes,
                  "GS persistent arena too small for LUT + state + nodes");
    configure_arenas(kPersistentBytes, GLOBAL_ARENA_SIZE - kPersistentBytes, 0);

    registerParam("Feed", &params.feed, 0.0f, 0.1f);
    registerParam("Kill", &params.k, 0.0f, 0.1f);
    // dA/dB cap at 0.05: this is explicit Euler, stable only while
    // dt·D·λmax ≤ 2. The graph Laplacian's |λ|max ≤ 2·6 = 12 (6-NN lattice), so
    // at the Speed slider's top (dt = 3.0) the joint worst case is 3·0.05·12 =
    // 1.8 ≤ 2 — stable across the whole GUI. This bound covers only the diffusion
    // term; the full forward-Euler stability region also depends on the reaction
    // Jacobian (the A·B² autocatalysis and the feed/kill terms), which is assumed
    // sub-dominant within the registered Feed/Kill/dt ranges below. If any of
    // those slider ranges are ever widened, the reaction Jacobian must be
    // re-checked — it is not bounded by the diffusion argument above.
    //
    // Backstop: the Euler bound above is a *visual-quality* guarantee, not the
    // hard stability backstop. to_q16() clamps every written state to [0, 1]
    // (Q16) each step, so even a range-widening or reaction-Jacobian overshoot
    // that violates the Euler bound cannot make the simulation numerically blow
    // up — it can only degrade visually (oscillation / saturation artifacts).
    registerParam("dA", &params.dA, 0.0f, 0.05f);
    registerParam("dB", &params.dB, 0.0f, 0.05f);
    registerParam("Speed", &params.dt, 0.1f, 3.0f);

    state.A = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    state.B = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    // A starts at 1.0 (65535), B starts at 0.0
    for (int i = 0; i < RD_N; i++) {
      state.A[i] = 65535;
      state.B[i] = 0;
    }

    cube_lut.build(persistent_arena);
    nodes = static_cast<Vector *>(
        persistent_arena.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    build_nodes(nodes);
    seed_clusters();
    init_orientation_animation();
  }

private:
  // Test seam: lets the unit tests reach to_q16/from_q16, step_physics, and
  // params without exposing them to production callers (zero runtime cost).
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
   * @details GS morphogenesis is slow per step: the feed/kill reaction (feed
   * ~0.04, kill ~0.06) and the small diffusion rates move the field only
   * fractionally each substep, so 16 substeps/frame are needed for the pattern to
   * grow and drift at a visible rate. BZ's oscillatory dynamics move far faster
   * per step and read visibly in 2. Every substep stays inside the explicit-Euler
   * stability bound documented at the "Speed" (dt) param, so the higher count is
   * pure evolution throughput, not a stability risk.
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
      for (int nb : ReactionGraph::neighbors[idx])
        if (nb >= 0)
          state.B[nb] = 65535;
    }
  }

  /**
   * @brief Advances one Gray-Scott substep into the next buffers (Jacobi).
   * @param cA Current A field (read-only), Q16 per node.
   * @param cB Current B field (read-only), Q16 per node.
   * @param nA Next A field (write target), Q16 per node.
   * @param nB Next B field (write target), Q16 per node.
   * @details Gray-Scott: dA/dt = dA·∇²A - A·B² + feed·(1-A);
   * dB/dt = dB·∇²B + A·B² - (k+feed)·B. Pure double-buffered (Jacobi) step:
   * reads the current buffers, writes the next ones. The caller owns the
   * ping-pong so the result can be landed back in the persistent state
   * regardless of substep parity (see render()).
   */
  void step_physics(const uint16_t *cA, const uint16_t *cB, uint16_t *nA,
                    uint16_t *nB) {
    for (int i = 0; i < RD_N; i++) {
      float a = from_q16(cA[i]);
      float b = from_q16(cB[i]);

      // Both species' Laplacians share one neighbor walk (fused on purpose:
      // two single-field neighbor walks would double the lattice reads).
      float lA = 0, lB = 0;
      for_each_neighbor(i, [&](int ni) {
        lA += from_q16(cA[ni]) - a;
        lB += from_q16(cB[ni]) - b;
      });

      float abb = a * b * b;
      nA[i] = to_q16(a + (params.dA * lA - abb + params.feed * (1.0f - a)) *
                             params.dt);
      nB[i] = to_q16(b + (params.dB * lB + abb - (params.k + params.feed) * b) *
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
    // Empty kernel (no node within the support radius): guard the division. A
    // 0/0 NaN would slip past the b < B_CULL_THRESHOLD cull in the shader (NaN
    // compares false) and silently poison palette.get().
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
    ScratchScope _frame(scratch_arena_a);
    uint16_t *sA = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    uint16_t *sB = static_cast<uint16_t *>(
        scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));

    // Ping-pong between the persistent state and the scratch buffers. cur always
    // names the latest generation; nxt is the write target for the next step.
    uint16_t *curA = state.A, *curB = state.B;
    uint16_t *nxtA = sA, *nxtB = sB;
    for (int k = 0; k < STEPS_PER_FRAME; k++) {
      step_physics(curA, curB, nxtA, nxtB);
      std::swap(curA, nxtA);
      std::swap(curB, nxtB);
    }

    // Land the final generation back in the persistent buffers so it survives
    // the ScratchScope pop. An even substep count already leaves cur == state
    // (no copy); an odd count leaves it in scratch and is copied back — correct
    // for any STEPS_PER_FRAME parity.
    if (curA != state.A) {
      std::memcpy(state.A, curA, RD_N * sizeof(uint16_t));
      std::memcpy(state.B, curB, RD_N * sizeof(uint16_t));
    }

    // `nodes` is the fixed lattice, built once in init(). Hoist the per-pixel
    // cubemap lookup into the vertex shader (run once at the pixel center) and
    // keep the per-sub-sample refine + interpolation in the fragment shader,
    // sparing ~3 redundant CubemapLUT lookups per pixel under 4× SSAA. The
    // pixel-center seed only feeds refine_nearest_node, which re-finds the
    // genuine nearest from the sub-sample position, so the rendered result is
    // identical to looking the cubemap up per sub-sample.
    auto vertex_shader = [&](Fragment &frag) {
      Vector rv = orientation.unorient(frag.pos);
      frag.v0 = static_cast<float>(cube_lut.lookup(rv));
    };

    // Per sub-sample: refine the seed, sample B, and map to a palette color or
    // cull to transparent below B_CULL_THRESHOLD.
    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      Vector rv = orientation.unorient(v);
      // Refine the cubemap seed to the genuine nearest node before centering the
      // kernel — the CubemapLUT quantizes `rv` to a face cell and can return a
      // not-quite-nearest seed; refine_nearest_node removes that quantization
      // bias.
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

  Vector *nodes = nullptr; /**< Fixed Fibonacci-lattice node positions, built
                                once in init(). */

  /** @brief Color palette mapping the B gradient to RGB. */
  GenerativePalette palette{
      GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY,
      BrightnessProfile::ASCENDING, SaturationProfile::VIBRANT};

  /**
   * @brief GUI-tunable Gray-Scott parameters.
   */
  struct Params {
    float feed = 0.04f; /**< Feed rate of A. */
    float k = 0.06f;    /**< Kill rate of B. */
    float dA = 0.02f;   /**< Diffusion coefficient of A. */
    float dB = 0.01f;   /**< Diffusion coefficient of B. */
    float dt = 2.5f;    /**< Integration timestep (Speed slider). */
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(GSReactionDiffusion)
