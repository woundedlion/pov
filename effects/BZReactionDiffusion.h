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
 * Scratch arena (per frame, disjoint phases):
 *   - Physics: ping-pong 3 × 7680 × 1B + float pre-convert 3 × 7680 × 4B = 115,200 B
 *   - Raster:  oriented lattice 7680 × 12B                               =  92,160 B
 */
template <int W, int H>
class BZReactionDiffusion
    : public ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::advance_substeps;
  using Base::cube_lut;
  using Base::dist2;
  using Base::for_each_neighbor;
  using Base::init_lattice;
  using Base::INV_R2;
  using Base::nodes;
  using Base::orientation;
  using Base::RD_K;
  using Base::RD_N;
  using Base::refine_center;
  using Base::register_param;
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
    constexpr size_t CUBE_LUT_BYTES = 6u * ReactionGraph::CubemapLUT::RES *
                                     ReactionGraph::CubemapLUT::RES *
                                     sizeof(uint16_t); // cube_lut.build
    constexpr size_t STATE_BYTES = 3u * RD_N * sizeof(uint8_t); // allocate_state
    // NODE_BYTES bounds both the resident node array and the equal-size transient
    // lattice cube_lut.build() carves and rewinds before init_lattice() allocates
    // the resident one; the build() peak (state + LUT + transient) is the max.
    constexpr size_t NODE_BYTES = RD_N * sizeof(Vector);        // build_nodes
    constexpr size_t PERSISTENT_BYTES = 165 * 1024;
    static_assert(CUBE_LUT_BYTES + STATE_BYTES + NODE_BYTES <= PERSISTENT_BYTES,
                  "BZ persistent arena too small for LUT + state + build peak");
    // render()'s scratch peaks at the larger of the physics phase (3 uint8 +
    // 3 float species buffers) and the raster phase (the oriented lattice);
    // the two run under disjoint scopes.
    constexpr size_t PHYSICS_SCRATCH_BYTES =
        3u * RD_N * sizeof(uint8_t) + 3u * RD_N * sizeof(float);
    constexpr size_t RASTER_SCRATCH_BYTES = RD_N * sizeof(Vector);
    constexpr size_t SCRATCH_BYTES =
        PHYSICS_SCRATCH_BYTES > RASTER_SCRATCH_BYTES ? PHYSICS_SCRATCH_BYTES
                                                     : RASTER_SCRATCH_BYTES;
    static_assert(SCRATCH_BYTES <= GLOBAL_ARENA_SIZE - PERSISTENT_BYTES,
                  "BZ scratch arena too small for render()'s phase peak");
    configure_arenas(PERSISTENT_BYTES, GLOBAL_ARENA_SIZE - PERSISTENT_BYTES, 0);

    // Lotka-Volterra predation coefficient; bounded only by to_q8's [0,1] clamp
    // (not the diffusion stability bound below), so a high value saturates.
    register_param("Compete", &params.alpha, 0.0f, 4.0f);
    // Explicit Euler is stable only while dt·D·λmax ≤ 2. The graph Laplacian on a
    // degree-RD_K lattice has |λ|max ≤ 2·RD_K (= 12 at RD_K=6), bounding these
    // Diff/Speed tops.
    register_param("Diff", &params.D, 0.001f, 0.1f);
    register_param("Speed", &params.dt, 0.0f, 1.0f);

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
   * @details Distinct from KERNEL_MIN_TOTAL_WEIGHT: that guards the Wendland
   * weight sum, this the blended concentration (a full-weight kernel can still
   * average to ~0 if all species are absent), and sample_kernel culls below it.
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
   * @details Rounds to nearest (+0.5f); truncating would bias the dynamics down
   *          by dropping sub-LSB positive updates. clamp bounds the input so
   *          255.5 -> 255 with no overflow.
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
  HS_COLD_MEMBER void seed_spiral_nuclei() {
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
  HS_O3_FN uint8_t advance_species(float conc, float predator,
                                   float laplacian) const {
    return to_q8(conc + (params.D * laplacian +
                         conc * (1 - conc - params.alpha * predator)) *
                            params.dt);
  }

  /**
   * @brief Applies stochastic perturbations to prevent convergence.
   * @param n_a Species A buffer to nudge (Q8, modified in place).
   * @param n_b Species B buffer to nudge (Q8, modified in place).
   * @param n_c Species C buffer to nudge (Q8, modified in place).
   * @details Nudges NUM_PERTURBATIONS random nodes by PERTURB_AMOUNT (Q8),
   *          saturating at 255, to keep the dynamics from settling on the
   *          closed manifold.
   * @note Draws from the global deterministic RNG (2*NUM_PERTURBATIONS draws per
   *       call), so retuning the draw count is a global-determinism change.
   */
  static void perturb_state(uint8_t *n_a, uint8_t *n_b, uint8_t *n_c) {
    for (int p = 0; p < NUM_PERTURBATIONS; p++) {
      int idx = hs::rand_int(0, RD_N);
      int s = hs::rand_int(0, 3);
      uint8_t *t = (s == 0) ? n_a : (s == 1) ? n_b : n_c;
      t[idx] = static_cast<uint8_t>(
          std::min(static_cast<int>(t[idx]) + PERTURB_AMOUNT, 255));
    }
  }

  /**
   * @brief Runs one full physics step: reaction-diffusion plus perturbation.
   * @param c_a Current species A buffer (Q8, read-only).
   * @param c_b Current species B buffer (Q8, read-only).
   * @param c_c Current species C buffer (Q8, read-only).
   * @param n_a Next species A buffer (Q8, written).
   * @param n_b Next species B buffer (Q8, written).
   * @param n_c Next species C buffer (Q8, written).
   * @param f_a Float scratch (RD_N) for the current A generation.
   * @param f_b Float scratch (RD_N) for the current B generation.
   * @param f_c Float scratch (RD_N) for the current C generation.
   * @details Double-buffered Jacobi: reads current buffers, writes next; the
   *          caller owns the ping-pong (see render()). Pre-converts the current
   *          generation into f_a/f_b/f_c once so the neighbor loop reads floats.
   */
  HS_O3_FN void step_physics(const uint8_t *c_a, const uint8_t *c_b,
                             const uint8_t *c_c, uint8_t *n_a, uint8_t *n_b,
                             uint8_t *n_c, float *f_a, float *f_b, float *f_c) {
    for (int i = 0; i < RD_N; i++) {
      f_a[i] = from_q8(c_a[i]);
      f_b[i] = from_q8(c_b[i]);
      f_c[i] = from_q8(c_c[i]);
    }
    for (int i = 0; i < RD_N; i++) {
      float a = f_a[i];
      float b = f_b[i];
      float c = f_c[i];

      float sum_a = 0, sum_b = 0, sum_c = 0;
      for_each_neighbor(i, [&](int ni) {
        sum_a += f_a[ni];
        sum_b += f_b[ni];
        sum_c += f_c[ni];
      });
      float l_a = sum_a - RD_K * a;
      float l_b = sum_b - RD_K * b;
      float l_c = sum_c - RD_K * c;

      n_a[i] = advance_species(a, c, l_a);
      n_b[i] = advance_species(b, a, l_b);
      n_c[i] = advance_species(c, b, l_c);
    }

    perturb_state(n_a, n_b, n_c);
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
   * @return Composited RGB pixel. Requires non-negative @p a, @p b, @p c: the
   *         blend then normalizes by their sum into a convex combination, which
   *         keeps each 16-bit channel in [0, 65535] so the cast needs no clamp.
   * @pre a + b + c >= SPECIES_EMPTY_EPS (the sole caller returns transparent
   *      below that floor), so the reciprocal is finite.
   * @details Concentration-weighted average: (ca·a + cb·b + cc·c) / (a + b + c),
   *          a hue-driven mix that is not concentration-dimmed.
   */
  HS_O3_FN static Pixel blend_species(float a, float b, float c,
                                      const Color4 &ca, const Color4 &cb,
                                      const Color4 &cc) {
    float inv = 1.0f / (a + b + c);
    float r = (ca.color.r * a + cb.color.r * b + cc.color.r * c) * inv;
    float g = (ca.color.g * a + cb.color.g * b + cc.color.g * c) * inv;
    float bl = (ca.color.b * a + cb.color.b * b + cc.color.b * c) * inv;

    return Pixel(static_cast<uint16_t>(r + 0.5f), static_cast<uint16_t>(g + 0.5f),
                 static_cast<uint16_t>(bl + 0.5f));
  }

  /**
   * @brief Turns accumulated kernel weights into a blended pixel color.
   * @param tw Total Wendland weight over the stencil.
   * @param wa Weighted sum of species A (Q8) over the stencil.
   * @param wb Weighted sum of species B (Q8) over the stencil.
   * @param wc Weighted sum of species C (Q8) over the stencil.
   * @param ca Palette color for species A.
   * @param cb Palette color for species B.
   * @param cc Palette color for species C.
   * @return Opaque blended Color4, or transparent black if no kernel weight
   *         accumulates.
   * @details The tail shared by every sub-sample: normalizes by total weight,
   *          culls empty walks, and blends the three species.
   */
  HS_O3_FN static Color4 finalize_sample(float tw, float wa, float wb, float wc,
                                         const Color4 &ca, const Color4 &cb,
                                         const Color4 &cc) {
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
   * @brief Hoisted 4× SSAA body for one pixel: refine once, re-weight per sample.
   * @tparam Grid Scan::Shader::SsaaGrid type supplying the sub-pixel offsets.
   * @param seed Cubemap-LUT seed node id for this pixel (from the vertex shader).
   * @param center_rv World-space direction at the pixel center.
   * @param world_nodes Oriented lattice node positions.
   * @param grid Row's SSAA sub-pixel grid.
   * @param x Pixel column.
   * @param ca Palette color for species A.
   * @param cb Palette color for species B.
   * @param cc Palette color for species C.
   * @return The finished, alpha-premultiplied pixel.
   * @details The four ±0.25 px sub-samples share one interpolation stencil (the
   * nearest node and its neighbors), refined and gathered once at the pixel
   * center; only the Wendland weights vary per sub-sample. A sub-sample
   * straddling a Voronoi boundary reuses the center's stencil rather than its
   * own — the ±0.25 px offset keeps that difference below one node spacing.
   * The whole body carries HS_O3_FN: an -Os loop around -O3 leaf calls forfeits
   * most of the codegen win.
   */
  template <typename Grid>
  HS_O3_FN Pixel shade_pixel(int seed, const Vector &center_rv,
                             const Vector *world_nodes, const Grid &grid, int x,
                             const Color4 &ca, const Color4 &cb,
                             const Color4 &cc) const {
    int center = refine_center(center_rv, world_nodes, seed);
    Vector spos[RD_K + 1];
    uint8_t sa[RD_K + 1], sb[RD_K + 1], sc[RD_K + 1];
    spos[0] = world_nodes[center];
    sa[0] = state.A[center];
    sb[0] = state.B[center];
    sc[0] = state.C[center];
    int k = 1;
    for_each_neighbor(center, [&](int ni) {
      spos[k] = world_nodes[ni];
      sa[k] = state.A[ni];
      sb[k] = state.B[ni];
      sc[k] = state.C[ni];
      ++k;
    });

    constexpr float inv_samples = 1.0f / 4.0f;
    Pixel accum(0, 0, 0);
    for (int i = 0; i < 4; ++i) {
      Vector v = grid.at(x, i);
      float tw = 0, wa = 0, wb = 0, wc = 0;
      for (int j = 0; j < RD_K + 1; ++j) {
        float u = 1.0f - dist2(v, spos[j]) * INV_R2;
        if (u > 0) {
          float w = u * u;
          wa += sa[j] * w;
          wb += sb[j] * w;
          wc += sc[j] * w;
          tw += w;
        }
      }
      Color4 c = finalize_sample(tw, wa, wb, wc, ca, cb, cc);
      accum += c.color * (c.alpha * inv_samples);
    }
    return accum;
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
    HS_PROFILE(bz_render);
    {
      ScratchScope physics_guard(scratch_arena_a);
      HS_PROFILE(bz_physics);
      uint8_t *s_a = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
      uint8_t *s_b = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
      uint8_t *s_c = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
      float *f_a = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *f_b = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *f_c = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));

      advance_substeps(STEPS_PER_FRAME,
                       std::array<uint8_t *, 3>{state.A, state.B, state.C},
                       std::array<uint8_t *, 3>{s_a, s_b, s_c},
                       [&](auto &cur, auto &nxt) {
                         step_physics(cur[0], cur[1], cur[2],
                                      nxt[0], nxt[1], nxt[2], f_a, f_b, f_c);
                       });
    }

    // Physics scratch is popped; the raster phase reuses the arena for the
    // oriented lattice so the kernel walks stay in world space.
    Vector *world_nodes = static_cast<Vector *>(
        scratch_arena_a.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    {
      HS_PROFILE(bz_orient);
      orient_nodes(nodes, world_nodes, RD_N, orientation.get());
    }

    const Color4 &ca = color_a;
    const Color4 &cb = color_b;
    const Color4 &cc = color_c;

    auto vertex_shader = [this](Fragment &frag) { seed_face_lut(frag); };

    auto pixel_shader = [&](Fragment &frag, const auto &grid, int x) -> Pixel {
      return shade_pixel(static_cast<int>(frag.v0), frag.pos, world_nodes, grid,
                         x, ca, cb, cc);
    };

    {
      HS_PROFILE(bz_raster);
      Scan::Shader::draw_grid<W, H>(canvas, vertex_shader, pixel_shader);
    }
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

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(BZReactionDiffusion)
