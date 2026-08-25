/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file BZReactionDiffusion.h
 * @brief Belousov-Zhabotinsky reaction-diffusion on a Fibonacci lattice
 *        sphere.
 */

#include <algorithm>
#include <cstring>
#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "effects/ReactionDiffusionBase.h"

// Unit-test accessor (tests/test_effects.h) reaching the private Q16
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
 * indefinitely. State is stored as Q16 (uint16_t): at the low end of the Diff
 * slider the diffusion term moves a node by only ~3e-4 of full scale per
 * substep, and a store coarser than that rounds the spatial coupling away
 * entirely, leaving uncoupled per-node ODEs. Per-pixel rendering uses Wendland
 * C2 kernel interpolation for smooth cell boundaries between lattice nodes.
 *
 * Shared lattice/orientation/kernel scaffolding lives in ReactionDiffusionBase.
 *
 * Memory budget (persistent arena, configured 184 KB):
 *   - Cubemap LUT:                  6 × 64² × 2B = 49,152 B
 *   - State:   3 arrays × 7680 × 2B (Q16)        = 46,080 B
 *   - Node XYZ: 7680 × 12B                       = 92,160 B  (fixed lattice, built once)
 *   - Total:                                       187,392 B (183 KB)
 *
 * Scratch arena (per frame, disjoint phases):
 *   - Physics: float generation mirror 3 × 7680 × 4B = 92,160 B
 *   - Raster:  oriented lattice 7680 × 12B           = 92,160 B
 */
template <int W, int H>
class BZReactionDiffusion
    : public ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<BZReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::cube_lut;
  using Base::dist2;
  using Base::for_each_neighbor;
  using Base::from_q16;
  using Base::init_lattice;
  using Base::orient_lattice;
  using Base::Q16_SCALE;
  using Base::RD_K;
  using Base::RD_N;
  using Base::refine_center;
  using Base::refine_render_center;
  using Base::register_param;
  using Base::seed_blobs;
  using Base::seed_face_lut;
  using Base::to_q16;
  using Base::with_wendland_weight;

public:
  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  BZReactionDiffusion() = default;

  /**
   * @brief Carves the arenas, registers tunable params, builds the lattice,
   *        and seeds the initial spiral nuclei.
   */
  void init() override {
    constexpr size_t PERSISTENT_BYTES = 184 * 1024;
    // render()'s scratch peaks at the larger of the physics phase (the 3 float
    // generation mirrors) and the raster phase (the oriented lattice); the two
    // run under disjoint scopes.
    constexpr size_t PHYSICS_SCRATCH_BYTES = 3u * RD_N * sizeof(float);
    constexpr size_t RASTER_SCRATCH_BYTES = RD_N * sizeof(Vector);
    constexpr size_t SCRATCH_BYTES =
        PHYSICS_SCRATCH_BYTES > RASTER_SCRATCH_BYTES ? PHYSICS_SCRATCH_BYTES
                                                     : RASTER_SCRATCH_BYTES;
    Base::template configure_rd_arenas<uint16_t, 3, PERSISTENT_BYTES, 0,
                                       SCRATCH_BYTES>();

    // Lotka-Volterra predation coefficient; bounded only by to_q16's [0,1] clamp
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
    color_a = FloatRgb(palette.get(0.0f).color);
    color_b = FloatRgb(palette.get(0.5f).color);
    color_c = FloatRgb(palette.get(1.0f).color);
  }

private:
  // Test seam: lets the unit tests reach the private Q16 helpers, physics, and
  // params without exposing them to production callers.
  friend struct ::hs_test::effects_tests::BZWhiteBox;

  struct FloatRgb {
    float r = 0, g = 0, b = 0;

    FloatRgb() = default;
    explicit FloatRgb(const Pixel &color)
        : r(color.r), g(color.g), b(color.b) {}
  };

  /**
   * @brief Concentration-sum floor below which a location is treated as empty.
   * @details Distinct from KERNEL_MIN_TOTAL_WEIGHT: that guards the Wendland
   * weight sum, this the blended concentration; a full-weight kernel can still
   * average to ~0 if all species are absent.
   */
  static constexpr float SPECIES_EMPTY_EPS = 1e-6f;

  // Simulation tuning.
  static constexpr int CLUSTERS_PER_SPECIES =
      3;                                      // seed blobs per species at init
  static constexpr int STEPS_PER_FRAME = 2;   // physics substeps per render
  static constexpr int NUM_PERTURBATIONS = 8; // random nudges per physics step
  static constexpr int PERTURB_AMOUNT = 771;  // Q16 nudge, 1.18% of full scale

  // ---------------------------------------------------------------------------
  // Initialization helpers
  // ---------------------------------------------------------------------------

  /**
   * @brief Allocates the three persistent Q16 species arrays and zeroes them.
   */
  void allocate_state() {
    constexpr size_t BYTES = RD_N * sizeof(uint16_t);
    state.A = static_cast<uint16_t *>(
        persistent_arena.allocate(BYTES, alignof(uint16_t)));
    state.B = static_cast<uint16_t *>(
        persistent_arena.allocate(BYTES, alignof(uint16_t)));
    state.C = static_cast<uint16_t *>(
        persistent_arena.allocate(BYTES, alignof(uint16_t)));
    memset(state.A, 0, BYTES);
    memset(state.B, 0, BYTES);
    memset(state.C, 0, BYTES);
  }

  // ---------------------------------------------------------------------------
  // Seeding
  // ---------------------------------------------------------------------------

  /**
   * @brief Seeds CLUSTERS_PER_SPECIES saturated blobs per species.
   * @details Puts all three species on the sphere so the cyclic competition has
   *          something to sustain spiral waves from. Seeds A, then B, then C:
   *          the order fixes the RNG stream position for everything downstream.
   */
  HS_COLD_MEMBER void seed_spiral_nuclei() {
    uint16_t *species[] = {state.A, state.B, state.C};
    for (uint16_t *field : species)
      seed_blobs(field, CLUSTERS_PER_SPECIES);
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
   * @return Updated concentration as a Q16 sample in [0, 65535].
   * @details Combines Fickian diffusion (D * laplacian) with Lotka-Volterra
   *          competition (conc * (1 - conc - alpha * predator)), scaled by the
   *          timestep dt.
   */
  HS_O3_FN uint16_t advance_species(float conc, float predator,
                                    float laplacian) const {
    return to_q16(conc + (params.D * laplacian +
                          conc * (1 - conc - params.alpha * predator)) *
                             params.dt);
  }

  /**
   * @brief Applies stochastic perturbations to prevent convergence.
   * @param n_a Species A buffer to nudge (Q16, modified in place).
   * @param n_b Species B buffer to nudge (Q16, modified in place).
   * @param n_c Species C buffer to nudge (Q16, modified in place).
   * @details Nudges NUM_PERTURBATIONS random nodes by PERTURB_AMOUNT (Q16)
   *          scaled by the timestep, saturating at 65535, to keep the dynamics
   *          from settling on the closed manifold. The scaling matches the
   *          reaction-diffusion terms: at dt = 0 the physics is frozen and the
   *          nudge is 0, so a stopped sphere cannot drift to saturation.
   * @note Draws from the global deterministic RNG (2*NUM_PERTURBATIONS draws per
   *       call) whatever the timestep, so retuning the draw count is a
   *       global-determinism change.
   */
  void perturb_state(uint16_t *n_a, uint16_t *n_b, uint16_t *n_c) const {
    const int amount = static_cast<int>(PERTURB_AMOUNT * params.dt);
    for (int p = 0; p < NUM_PERTURBATIONS; p++) {
      int idx = hs::rand_int(0, RD_N);
      int s = hs::rand_int(0, 3);
      uint16_t *t = (s == 0) ? n_a : (s == 1) ? n_b : n_c;
      t[idx] = static_cast<uint16_t>(
          std::min(static_cast<int>(t[idx]) + amount, 65535));
    }
  }

  /**
   * @brief Runs one full physics step: reaction-diffusion plus perturbation.
   * @param s_a Species A state (Q16, advanced in place).
   * @param s_b Species B state (Q16, advanced in place).
   * @param s_c Species C state (Q16, advanced in place).
   * @param f_a Float mirror (RD_N) of the current A generation.
   * @param f_b Float mirror (RD_N) of the current B generation.
   * @param f_c Float mirror (RD_N) of the current C generation.
   * @details Jacobi: the whole current generation is mirrored into
   *          f_a/f_b/f_c first, and the update loop reads only that mirror, so
   *          writing the new generation over the state in place cannot feed a
   *          half-updated neighbor back into the Laplacian. The mirror also
   *          converts each node's Q16 sample once instead of once per neighbor
   *          visit.
   */
  HS_O3_FN void step_physics(uint16_t *s_a, uint16_t *s_b, uint16_t *s_c,
                             float *f_a, float *f_b, float *f_c) {
    for (int i = 0; i < RD_N; i++) {
      f_a[i] = from_q16(s_a[i]);
      f_b[i] = from_q16(s_b[i]);
      f_c[i] = from_q16(s_c[i]);
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

      s_a[i] = advance_species(a, c, l_a);
      s_b[i] = advance_species(b, a, l_b);
      s_c[i] = advance_species(c, b, l_c);
    }

    perturb_state(s_a, s_b, s_c);
  }

  // ---------------------------------------------------------------------------
  // Rendering: Wendland C2 kernel interpolation
  // ---------------------------------------------------------------------------

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
                             const FloatRgb &ca, const FloatRgb &cb,
                             const FloatRgb &cc) const {
    int center = refine_render_center(center_rv, world_nodes, seed);
    Vector spos[RD_K + 1];
    uint16_t sa[RD_K + 1], sb[RD_K + 1], sc[RD_K + 1];
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

    constexpr float INV_SAMPLES = 1.0f / Grid::SAMPLES;
    float mix_a = 0, mix_b = 0, mix_c = 0;
    for (int i = 0; i < Grid::SAMPLES; ++i) {
      Vector v = grid.at(x, i);
      float tw = 0, wa = 0, wb = 0, wc = 0;
      // always_inline on the accumulator: without it GCC spends +32 B of ITCM
      // on this stencil walk.
      for (int j = 0; j < RD_K + 1; ++j)
        with_wendland_weight(dist2(v, spos[j]),
                             [&](float w) __attribute__((always_inline)) {
                               wa += sa[j] * w;
                               wb += sb[j] * w;
                               wc += sc[j] * w;
                               tw += w;
                             });
      float species_sum = wa + wb + wc;
      if (tw <= Base::KERNEL_MIN_TOTAL_WEIGHT ||
          species_sum < SPECIES_EMPTY_EPS * Q16_SCALE * tw)
        continue;
      // Normalize species mass, floored by the weighted Q16 unit concentration.
      const float sample_normalizer = std::max(species_sum, Q16_SCALE * tw);
      float scale = INV_SAMPLES / sample_normalizer;
      mix_a += wa * scale;
      mix_b += wb * scale;
      mix_c += wc * scale;
    }
    float accum_r = ca.r * mix_a + cb.r * mix_b + cc.r * mix_c;
    float accum_g = ca.g * mix_a + cb.g * mix_b + cc.g * mix_c;
    float accum_b = ca.b * mix_a + cb.b * mix_b + cc.b * mix_c;
    return Pixel(
        static_cast<uint16_t>(hs::clamp(accum_r + 0.5f, 0.0f, 65535.0f)),
        static_cast<uint16_t>(hs::clamp(accum_g + 0.5f, 0.0f, 65535.0f)),
        static_cast<uint16_t>(hs::clamp(accum_b + 0.5f, 0.0f, 65535.0f)));
  }

  /**
   * @brief Allocates scratch, advances the simulation, then rasterizes a frame.
   * @param canvas Destination canvas to draw into.
   * @details Advances the persistent state STEPS_PER_FRAME substeps in place
   *          over one set of float generation mirrors, then rasterizes with 4x
   *          SSAA using a cubemap-LUT vertex shader and a kernel-sampling
   *          fragment shader.
   */
  void render(Canvas &canvas) {
    ScratchScope frame_guard(scratch_arena_a);
    HS_PROFILE(bz_render);
    {
      ScratchScope physics_guard(scratch_arena_a);
      HS_PROFILE(bz_physics);
      float *f_a = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *f_b = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *f_c = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));

      for (int k = 0; k < STEPS_PER_FRAME; ++k)
        step_physics(state.A, state.B, state.C, f_a, f_b, f_c);
    }

    // Physics scratch is popped; the raster phase reuses the arena for the
    // oriented lattice so the kernel walks stay in world space.
    Vector *world_nodes;
    {
      HS_PROFILE(bz_orient);
      world_nodes = orient_lattice();
    }

    const FloatRgb &ca = color_a;
    const FloatRgb &cb = color_b;
    const FloatRgb &cc = color_c;

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
   * @brief Persistent Q16 concentration buffers for the three species.
   * @details A, B, and C each point to RD_N samples in [0, 65535].
   */
  struct {
    uint16_t *A = nullptr, *B = nullptr, *C = nullptr;
  } state;

  /** @brief Legacy species colors expressed as an explicit generative recipe. */
  GenerativePalette palette{EffectPaletteRecipes::bz_reaction_diffusion()};

  /** @brief Per-species palette channels converted once in init(). */
  FloatRgb color_a, color_b, color_c;

  /** @brief User-tunable Lotka-Volterra and diffusion parameters. */
  struct Params {
    float alpha = 3.0f; /**< Competition (predation) coefficient. */
    float D = 0.05f;    /**< Diffusion coefficient. */
    float dt = 0.35f;   /**< Integration timestep per substep. */
  } params;
};

#include "core/control/registry.h"
REGISTER_EFFECT(BZReactionDiffusion)
