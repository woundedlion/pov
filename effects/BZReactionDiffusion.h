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
  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  FLASHMEM BZReactionDiffusion() = default;

  /**
   * @brief Carves the arenas, registers tunable params, builds the lattice,
   *        and seeds the initial spiral nuclei.
   */
  void init() override {
    // 165KB persistent: 48KB Cubemap LUT + 23KB State + 90KB node positions. The
    // node array (7680 × Vector) is the fixed Fibonacci lattice and does not
    // depend on the per-frame view orientation (queries are un-oriented onto it,
    // not the reverse), so it lives in the persistent arena and is built once.
    configure_arenas(165 * 1024, GLOBAL_ARENA_SIZE - 165 * 1024, 0);

    // "Compete" is the Lotka-Volterra predation coefficient, not an opacity;
    // the name avoids clashing with the engine-wide Alpha-as-opacity convention.
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

  /**
   * @brief Converts a Q8 fixed-point byte to a normalized float.
   * @param v Q8 value in [0, 255].
   * @return Concentration in [0.0, 1.0].
   */
  static inline float from_q8(uint8_t v) { return v * (1.0f / 255.0f); }
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
        for (int nb : ReactionGraph::neighbors[center])
          if (nb >= 0)
            species[s][nb] = 255;
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

  /**
   * @brief Blends three species concentrations into a single pixel.
   * @param a Species A concentration in [0.0, 1.0].
   * @param b Species B concentration in [0.0, 1.0].
   * @param c Species C concentration in [0.0, 1.0].
   * @param ca Palette color for species A.
   * @param cb Palette color for species B.
   * @param cc Palette color for species C.
   * @return Composited RGB pixel (16-bit channels clamped to [0, 65535]).
   * @details Concentration-weighted average of the three species colors:
   *          (ca·a + cb·b + cc·c) / (a + b + c). Each species contributes in
   *          proportion to its local concentration with no order bias — unlike
   *          a series alpha-over composite, which let whichever species was
   *          composited last (C) dominate every region where species overlapped.
   *          The normalization by total concentration makes the output a pure
   *          mix of the palette colors (hue-driven), not concentration-dimmed.
   */
  static Pixel blend_species(float a, float b, float c, const Color4 &ca,
                             const Color4 &cb, const Color4 &cc) {
    float wsum = a + b + c;
    if (wsum < 1e-6f)
      return Pixel(0, 0, 0);
    float inv = 1.0f / wsum;
    float r = (ca.color.r * a + cb.color.r * b + cc.color.r * c) * inv;
    float g = (ca.color.g * a + cb.color.g * b + cc.color.g * c) * inv;
    float bl = (ca.color.b * a + cb.color.b * b + cc.color.b * c) * inv;

    return Pixel(static_cast<uint16_t>(hs::clamp(r, 0.0f, 65535.0f)),
                 static_cast<uint16_t>(hs::clamp(g, 0.0f, 65535.0f)),
                 static_cast<uint16_t>(hs::clamp(bl, 0.0f, 65535.0f)));
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

    if (tw <= 0.0001f)
      return Color4(Pixel(0, 0, 0), 0.0f);

    float inv = (1.0f / 255.0f) / tw;
    float a = wa * inv, b = wb * inv, c = wc * inv;
    // Cull a location where every species has decayed to ~0 to transparent,
    // matching GSReactionDiffusion's B_CULL_THRESHOLD cull. Such a point is
    // semantically empty, so it must not paint opaque black over a (possibly
    // non-black) background. This is the same emptiness test blend_species
    // applies internally; deciding it here lets the alpha follow suit.
    if (a + b + c < 1e-6f)
      return Color4(Pixel(0, 0, 0), 0.0f);
    return Color4(blend_species(a, b, c, ca, cb, cc), 1.0f);
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
    // the ScratchScope pop. An even substep count leaves cur == state (no copy);
    // an odd count leaves the result in scratch, so copy it back. Correct for any
    // STEPS_PER_FRAME parity.
    if (curA != state.A) {
      std::memcpy(state.A, curA, RD_N);
      std::memcpy(state.B, curB, RD_N);
      std::memcpy(state.C, curC, RD_N);
    }

    // `nodes` is the fixed lattice, built once in init().
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

  /**
   * @brief Persistent Q8 concentration buffers for the three species.
   * @details A, B, and C each point to RD_N bytes; values are Q8 in [0, 255].
   */
  struct {
    uint8_t *A = nullptr, *B = nullptr, *C = nullptr;
  } state;

  /**
   * @brief Fixed Fibonacci-lattice node positions, built once in init().
   * @details Lives in the persistent arena; independent of the per-frame view
   *          orientation.
   */
  Vector *nodes = nullptr;

  /** @brief Triadic generative palette mapping species concentration to color. */
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::DESCENDING,
                            SaturationProfile::VIBRANT, 42};

  /** @brief User-tunable Lotka-Volterra and diffusion parameters. */
  struct Params {
    float alpha = 3.0f; /**< Competition (predation) coefficient. */
    float D = 0.05f;    /**< Diffusion coefficient. */
    float dt = 0.35f;   /**< Integration timestep per substep. */
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(BZReactionDiffusion)
