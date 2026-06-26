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

// Unit-test accessor (tests/test_effects.h) reaching the private Q8
// conversions, advance_species, perturb_state, and one physics substep, which
// the smoke/determinism harness cannot pin.
namespace hs_test {
namespace effects_tests {
struct BZWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Belousov-Zhabotinsky reaction-diffusion on a Fibonacci lattice sphere.
 *
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 *
 * @details
 * Three competing chemical species (A, B, C) evolve via Lotka-Volterra dynamics
 * on a 7680-node Fibonacci lattice with K=6 nearest neighbors. The cyclic
 * competition (A→B→C→A) creates self-sustaining spiral waves that persist
 * indefinitely. State is stored as Q8 (uint8_t). Per-pixel rendering uses
 * Wendland C2 kernel interpolation for smooth cell boundaries between lattice
 * nodes.
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
  using Base::cube_lut;
  using Base::for_each_neighbor;
  using Base::init_lattice;
  using Base::kernel_accumulate;
  using Base::land_back;
  using Base::nodes;
  using Base::orientation;
  using Base::RD_N;
  using Base::refine_nearest_node;
  using Base::registerParam;
  using Base::seed_face_lut;

public:
  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  FLASHMEM BZReactionDiffusion() = default;

  /**
   * @brief Carves the arenas, registers tunable params, builds the lattice,
   *        and seeds the initial spiral nuclei.
   */
  void init() override {
    constexpr size_t kCubeLutBytes = 6u * ReactionGraph::CubemapLUT::RES *
                                     ReactionGraph::CubemapLUT::RES *
                                     sizeof(uint16_t); // cube_lut.build
    constexpr size_t kStateBytes = 3u * RD_N * sizeof(uint8_t); // allocate_state
    constexpr size_t kNodeBytes = RD_N * sizeof(Vector);        // build_nodes
    constexpr size_t kPersistentBytes = 165 * 1024;
    static_assert(kCubeLutBytes + kStateBytes + kNodeBytes <= kPersistentBytes,
                  "BZ persistent arena too small for LUT + state + nodes");
    configure_arenas(kPersistentBytes, GLOBAL_ARENA_SIZE - kPersistentBytes, 0);

    // "Compete" is the Lotka-Volterra predation coefficient, not an opacity.
    registerParam("Compete", &params.alpha, 0.0f, 4.0f);
    // Explicit Euler is stable only while dt·D·λmax ≤ 2 (|λ|max ≤ 12 on the 6-NN
    // lattice), bounding these Diff/Speed tops.
    registerParam("Diff", &params.D, 0.001f, 0.1f);
    registerParam("Speed", &params.dt, 0.0f, 1.0f);

    allocate_state();
    cube_lut.build(persistent_arena);
    init_lattice();
    seed_spiral_nuclei();

    // Fixed-seed palette: the three species colors never change after init.
    color_a = palette.get(0.0f);
    color_b = palette.get(0.5f);
    color_c = palette.get(1.0f);
  }

private:
  // Test seam: lets the unit tests reach the private Q8 helpers, physics, and
  // params without exposing them to production callers.
  friend struct ::hs_test::effects_tests::BZWhiteBox;

  // ---------------------------------------------------------------------------
  // Q8 fixed-point helpers
  // ---------------------------------------------------------------------------

  /**
   * @brief Q8 full-scale factor: maps the [0, 255] byte state to [0.0, 1.0].
   */
  static constexpr float Q8_SCALE = 255.0f;
  static constexpr float Q8_INV = 1.0f / Q8_SCALE; /**< Reciprocal of Q8_SCALE. */

  /**
   * @brief Concentration-sum floor below which a location is treated as empty.
   * @details The kernel-blended species concentrations are in [0, 1]; when their
   * sum falls below this, every species has decayed to ~0 there, so blend_species
   * has no hue to mix and sample_kernel culls to transparent. Distinct from
   * KERNEL_MIN_TOTAL_WEIGHT (which guards the Wendland weight sum, not the
   * concentrations): a kernel can carry full weight yet still average to ~0 if
   * all three species are absent at the covered nodes.
   */
  static constexpr float SPECIES_EMPTY_EPS = 1e-6f;

  /**
   * @brief Converts a Q8 fixed-point byte to a normalized float.
   * @param v Q8 value in [0, 255].
   * @return Concentration in [0.0, 1.0].
   */
  static inline float from_q8(uint8_t v) { return v * Q8_INV; }
  /**
   * @brief Converts a normalized concentration to a Q8 fixed-point byte.
   * @param v Concentration; clamped to [0.0, 1.0] before scaling.
   * @return Q8 value in [0, 255], rounded to nearest.
   * @details Rounds to nearest (+0.5f) rather than truncating: truncation loses
   *          every sub-LSB positive update while negatives still decrement,
   *          biasing the dynamics downward. clamp bounds the product to
   *          [0, 255], so +0.5f tops out at 255.5 -> 255 with no overflow.
   */
  static inline uint8_t to_q8(float v) {
    return static_cast<uint8_t>(hs::clamp(v, 0.0f, 1.0f) * Q8_SCALE + 0.5f);
  }

  // Simulation tuning.
  static constexpr int CLUSTERS_PER_SPECIES = 3; // seed blobs per species at init
  static constexpr int STEPS_PER_FRAME = 2;      // physics substeps per render
  static constexpr int NUM_PERTURBATIONS = 8;    // random nudges per physics step
  static constexpr int PERTURB_AMOUNT = 3;       // Q8 magnitude of each nudge

  // ---------------------------------------------------------------------------
  // Initialization helpers
  // ---------------------------------------------------------------------------

  /**
   * @brief Allocates the three persistent Q8 species arrays and zeroes them.
   */
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

  /**
   * @brief Seeds CLUSTERS_PER_SPECIES clusters per species at random nodes.
   * @details Saturates each cluster center and its neighbors to Q8 255,
   *          ensuring all three species are present so the cyclic competition
   *          can sustain spiral waves.
   */
  void seed_spiral_nuclei() {
    uint8_t *species[] = {state.A, state.B, state.C};
    for (int s = 0; s < 3; s++) {
      for (int k = 0; k < CLUSTERS_PER_SPECIES; k++) {
        int center = hs::rand_int(0, RD_N);
        species[s][center] = 255;
        for_each_neighbor(center, [&](int nb) { species[s][nb] = 255; });
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Physics: Lotka-Volterra reaction + graph Laplacian diffusion
  // ---------------------------------------------------------------------------

  /**
   * @brief Advances one species by a single diffusion + competition step.
   * @param conc Current concentration of this species in [0.0, 1.0].
   * @param predator Concentration of the species that preys on this one,
   *                 in [0.0, 1.0].
   * @param laplacian Graph-Laplacian of this species over its neighbors.
   * @return Updated concentration as a Q8 byte in [0, 255].
   * @details Combines Fickian diffusion (D * laplacian) with Lotka-Volterra
   *          competition (conc * (1 - conc - alpha * predator)), scaled by the
   *          timestep dt.
   */
  uint8_t advance_species(float conc, float predator, float laplacian) const {
    return to_q8(conc + (params.D * laplacian +
                         conc * (1 - conc - params.alpha * predator)) *
                            params.dt);
  }

  /**
   * @brief Applies stochastic perturbations to prevent convergence.
   * @param nA Species A buffer to nudge (Q8, modified in place).
   * @param nB Species B buffer to nudge (Q8, modified in place).
   * @param nC Species C buffer to nudge (Q8, modified in place).
   * @details Nudges NUM_PERTURBATIONS random nodes by PERTURB_AMOUNT (Q8),
   *          saturating at 255, to keep the dynamics from settling on the
   *          closed manifold.
   * @note Draws from the global deterministic RNG, so this advances the shared
   *       stream by exactly 2*NUM_PERTURBATIONS draws (idx + species each
   *       iteration). Keeps sim-vs-device byte-determinism but couples this
   *       effect's stream position to the substep count: retuning the draw count
   *       is a global-determinism change.
   */
  static void perturb_state(uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    for (int p = 0; p < NUM_PERTURBATIONS; p++) {
      int idx = hs::rand_int(0, RD_N);
      int s = hs::rand_int(0, 3);
      uint8_t *t = (s == 0) ? nA : (s == 1) ? nB : nC;
      t[idx] = static_cast<uint8_t>(
          std::min(static_cast<int>(t[idx]) + PERTURB_AMOUNT, 255));
    }
  }

  /**
   * @brief Runs one full physics step: reaction-diffusion plus perturbation.
   * @param cA Current species A buffer (Q8, read-only).
   * @param cB Current species B buffer (Q8, read-only).
   * @param cC Current species C buffer (Q8, read-only).
   * @param nA Next species A buffer (Q8, written).
   * @param nB Next species B buffer (Q8, written).
   * @param nC Next species C buffer (Q8, written).
   * @details Pure double-buffered (Jacobi): reads the current buffers, writes
   *          the next ones. The caller owns the ping-pong so the result can be
   *          landed back in the persistent state regardless of substep parity
   *          (see render()).
   */
  void step_physics(const uint8_t *cA, const uint8_t *cB, const uint8_t *cC,
                    uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    for (int i = 0; i < RD_N; i++) {
      float a = from_q8(cA[i]);
      float b = from_q8(cB[i]);
      float c = from_q8(cC[i]);

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

  /**
   * @brief Blends three species concentrations into a single pixel.
   * @param a Species A concentration in [0.0, 1.0].
   * @param b Species B concentration in [0.0, 1.0].
   * @param c Species C concentration in [0.0, 1.0].
   * @param ca Palette color for species A.
   * @param cb Palette color for species B.
   * @param cc Palette color for species C.
   * @return Composited RGB pixel; the convex blend keeps each 16-bit channel in [0, 65535] without clamping.
   * @details Concentration-weighted average of the three species colors:
   *          (ca·a + cb·b + cc·c) / (a + b + c). Each species contributes in
   *          proportion to its local concentration with no order bias.
   *          Normalizing by total concentration makes the output a pure mix of
   *          the palette colors (hue-driven), not concentration-dimmed.
   */
  static Pixel blend_species(float a, float b, float c, const Color4 &ca,
                             const Color4 &cb, const Color4 &cc) {
    float wsum = a + b + c;
    if (wsum < SPECIES_EMPTY_EPS)
      return Pixel(0, 0, 0);
    float inv = 1.0f / wsum;
    float r = (ca.color.r * a + cb.color.r * b + cc.color.r * c) * inv;
    float g = (ca.color.g * a + cb.color.g * b + cc.color.g * c) * inv;
    float bl = (ca.color.b * a + cb.color.b * b + cc.color.b * c) * inv;

    return Pixel(static_cast<uint16_t>(r), static_cast<uint16_t>(g),
                 static_cast<uint16_t>(bl));
  }

  /**
   * @brief Samples the kernel-interpolated color at a world-space point.
   * @param rv World-space query direction (un-oriented onto the lattice).
   * @param nodes Fixed Fibonacci-lattice node positions.
   * @param best_node Index of the nearest lattice node to rv.
   * @param ca Palette color for species A.
   * @param cb Palette color for species B.
   * @param cc Palette color for species C.
   * @return Opaque blended Color4, or transparent black if no kernel weight
   *         accumulates.
   * @details Accumulates Wendland C2 kernel weights over the neighborhood of
   *          best_node, normalizes by total weight, and blends the species.
   */
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

    if (tw <= Base::KERNEL_MIN_TOTAL_WEIGHT)
      return Color4(Pixel(0, 0, 0), 0.0f);

    float inv = Q8_INV / tw;
    float a = wa * inv, b = wb * inv, c = wc * inv;
    float total = a + b + c;
    if (total < SPECIES_EMPTY_EPS)
      return Color4(Pixel(0, 0, 0), 0.0f);
    return Color4(blend_species(a, b, c, ca, cb, cc),
                  hs::clamp(total, 0.0f, 1.0f));
  }

  /**
   * @brief Runs `steps` ping-ponged physics substeps and lands the final
   *        generation back in the persistent state.
   * @param steps Number of physics substeps to advance.
   * @param sA Scratch buffer A for the ping-pong (RD_N bytes).
   * @param sB Scratch buffer B for the ping-pong (RD_N bytes).
   * @param sC Scratch buffer C for the ping-pong (RD_N bytes).
   * @details Pure double-buffered: ping-pongs between the persistent state and
   *          the scratch buffers. An odd `steps` count leaves the final
   *          generation in scratch, so it is copied back into the persistent
   *          buffers before the caller's ScratchScope pops.
   */
  void advance_substeps(int steps, uint8_t *sA, uint8_t *sB, uint8_t *sC) {
    uint8_t *curA = state.A, *curB = state.B, *curC = state.C;
    uint8_t *nxtA = sA, *nxtB = sB, *nxtC = sC;
    for (int k = 0; k < steps; k++) {
      step_physics(curA, curB, curC, nxtA, nxtB, nxtC);
      std::swap(curA, nxtA);
      std::swap(curB, nxtB);
      std::swap(curC, nxtC);
    }

    land_back(std::array<uint8_t *, 3>{state.A, state.B, state.C},
              std::array<uint8_t *, 3>{curA, curB, curC}, RD_N);
  }

  /**
   * @brief Allocates scratch, advances the simulation, then rasterizes a frame.
   * @param canvas Destination canvas to draw into.
   * @details Runs STEPS_PER_FRAME ping-ponged physics substeps between the
   *          persistent state and scratch buffers, lands the final generation
   *          back in the persistent buffers, then rasterizes with 4x SSAA using
   *          a cubemap-LUT vertex shader and a kernel-sampling fragment shader.
   */
  void render(Canvas &canvas) {
    ScratchScope frame_guard(scratch_arena_a);

    uint8_t *sA = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *sB = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *sC = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));

    advance_substeps(STEPS_PER_FRAME, sA, sB, sC);

    const Color4 &ca = color_a;
    const Color4 &cb = color_b;
    const Color4 &cc = color_c;

    auto vertex_shader = [this](Fragment &frag) { seed_face_lut(frag); };

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

  /**
   * @brief Persistent Q8 concentration buffers for the three species.
   * @details A, B, and C each point to RD_N bytes; values are Q8 in [0, 255].
   */
  struct {
    uint8_t *A = nullptr, *B = nullptr, *C = nullptr;
  } state;

  /** @brief Triadic generative palette mapping species concentration to color. */
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::DESCENDING,
                            SaturationProfile::VIBRANT, 42};

  /** @brief Per-species palette colors, evaluated once in init(). */
  Color4 color_a, color_b, color_c;

  /** @brief User-tunable Lotka-Volterra and diffusion parameters. */
  struct Params {
    float alpha = 3.0f; /**< Competition (predation) coefficient. */
    float D = 0.05f;    /**< Diffusion coefficient. */
    float dt = 0.35f;   /**< Integration timestep per substep. */
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(BZReactionDiffusion)
