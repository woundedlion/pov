/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file GSReactionDiffusion.h
 * @brief Gray-Scott reaction-diffusion on a Fibonacci lattice sphere.
 */

#include <array>
#include <cmath>
#include <utility>
#include "core/color/effect_palette_recipes.h"
#include "core/color/noise_shimmer_palette.h"
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
 * Persistent state is Q16 (uint16_t) for the cubic reaction-term precision;
 * substeps integrate in float and quantize back once per frame. Shared
 * lattice/orientation/kernel scaffolding lives in ReactionDiffusionBase.
 *
 * A reaction runs until its field has all but stopped moving, then dissolves off
 * the sphere and reseeds at new sites, so each cycle grows a different form from
 * the same constants. Editing the constants dissolves the current field too.
 *
 * Each node carries the two strongest seed palettes and their mixing weight.
 * Palette pigment follows B diffusion and autocatalysis; rendering blends RGB
 * in linear light. Hue and shimmer share one sphere-domain noise field.
 */
template <int W, int H>
class GSReactionDiffusion
    : public ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::accumulate_stencil;
  using Base::cube_lut;
  using Base::for_each_neighbor;
  using Base::from_q16;
  using Base::gather_stencil;
  using Base::init_lattice;
  using Base::orient_lattice;
  using Base::Q16_INV;
  using Base::RD_K;
  using Base::RD_N;
  using Base::refine_and_accumulate;
  using Base::refine_render_center;
  using Base::register_param;
  using Base::rasterize_lattice;
  using Base::to_q16;

public:
  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  GSReactionDiffusion() = default;

  /**
   * @brief One-time setup: arenas, GUI params, A/B state, cubemap LUT, lattice.
   * @details Carves the persistent arena, registers the GUI params, allocates
   * and seeds the A/B state, and builds the cubemap LUT and lattice nodes once.
   */
  void init() override {
    constexpr size_t PALETTE_BYTES =
        NUM_SEED_CLUSTERS * PALETTE_SIZE * sizeof(Pixel) + alignof(Pixel) +
        HueNoiseLutView::SIZE * sizeof(int8_t);
    constexpr size_t PERSISTENT_BYTES = 192 * 1024 + 32;
    constexpr size_t PHYSICS_SCRATCH_BYTES =
        2u * RD_N * sizeof(float) + RD_N * sizeof(uint16_t);
    constexpr size_t RASTER_SCRATCH_BYTES =
        RD_N * sizeof(math::Vector) + 2u * RD_N * sizeof(uint8_t);
    constexpr size_t SCRATCH_BYTES =
        std::max(PHYSICS_SCRATCH_BYTES, RASTER_SCRATCH_BYTES);
    Base::template configure_rd_arenas<uint16_t, 3, PERSISTENT_BYTES,
                                       PALETTE_BYTES, SCRATCH_BYTES, 3>();

    register_param("Feed", &params.feed, 0.0f, 0.1f);
    register_param("Kill", &params.k, 0.0f, 0.1f);
    // The six-step temporal block scales the Euler timestep by 5/3, which at
    // the joint Speed/diffusion maximum exceeds the linear diffusion bound;
    // step_physics' per-substep [0,1] clamp bounds that corner.
    register_param("dA", &params.d_a, 0.0f, 0.05f);
    register_param("dB", &params.d_b, 0.0f, 0.05f);
    register_param("Speed", &params.dt, 0.1f, 3.0f);
    register_param("Noise Speed", &params.noise_speed, -0.001f, 0.001f);
    register_param("Noise Scale", &params.noise_scale, 1.0f / 64.0f, 8.0f);
    register_param("Hue Shift", &params.hue_shift, -4.0f, 4.0f);
    register_param("Shimmer", &params.shimmer, 0.0f, 1.0f);

    state.A = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    state.B = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    state.pigment = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    palettes = static_cast<Pixel *>(persistent_arena.allocate(
        NUM_SEED_CLUSTERS * PALETTE_SIZE * sizeof(Pixel), alignof(Pixel)));
    color_noise_lut = static_cast<int8_t *>(
        persistent_arena.allocate(HueNoiseLutView::SIZE, alignof(int8_t)));
    color_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    color_noise.SetSeed(6047);
    color_noise.SetFrequency(1.0f);
    refresh_color_noise();

    validate_physics_neighbors();
    init_lattice();
    seed_reaction();
    reaction_edited(); // latch the defaults; frame 1 is not an edit
  }

private:
  // Test seam: lets unit tests reach the Q16 helpers, step_physics, and params
  // without exposing them to production callers.
  friend struct ::hs_test::effects_tests::GSWhiteBox;

  /**
   * @brief Fully-saturated B blobs seeded per reaction.
   * @details Nucleation sites for the GS instability to grow from; without them
   * the uniform A=1/B=0 field never moves.
   */
  static constexpr int NUM_SEED_CLUSTERS = 30;
  static constexpr int PALETTE_SIZE = 32;
  static_assert(NUM_SEED_CLUSTERS <= 32);
  /** @brief Substep budget used to calibrate the stabilization threshold. */
  static constexpr int BASELINE_STEPS_PER_FRAME = 16;
  /** @brief Rendered frames the dissolve takes to convert every node back to
   * rest; 6.4 s at the 16 fps cadence. */
  static constexpr int DISSOLVE_FRAMES = 103;
  /**
   * @brief Rendered frames a fresh reaction runs before the stabilization
   * detector arms; 9.6 s at the 16 fps cadence.
   * @details A young field has only NUM_SEED_CLUSTERS active sites, so its mean
   * |dB| sits under MEAN_DB_STABLE and would read as stalled at birth.
   */
  static constexpr int MIN_GROW_FRAMES = 154;
  /** @brief Consecutive sub-floor frames that count as stabilized. */
  static constexpr int STABLE_HOLD_FRAMES = 39;
  /** @brief Fraction of the dissolve spent fading a node before conversion. */
  static constexpr float DISSOLVE_FADE_FRACTION = 1.0f / 16.0f;
  /**
   * @brief Mean per-node |dB| per frame below which the field counts as
   * settled, at DEFAULT_DT and BASELINE_STEPS_PER_FRAME; the detector rescales
   * it by params.dt / DEFAULT_DT and by EVOLUTION_STEPS_PER_FRAME /
   * BASELINE_STEPS_PER_FRAME. This is a calibration heuristic: Q16
   * quantization makes the low-Speed response nonlinear.
   * @details Loose relative to the 1.1e-6..4.0e-6 Q16 chatter a converged field
   * floors at; fires at ~222 baseline frames.
   */
  static constexpr float MEAN_DB_STABLE = 2.0e-4f;
  /** @brief Speed the stabilization floor is calibrated at. */
  static constexpr float DEFAULT_DT = 2.5f;
  /** @brief Base-dt substep equivalents advanced per rendered frame. */
  static constexpr int EVOLUTION_STEPS_PER_FRAME = 10;
  /**
   * @brief Euler integrations performed per rendered frame.
   * @details Six 5/3-sized integrations cover the same simulated interval as
   * ten base-dt integrations.
   */
  static constexpr int STEPS_PER_FRAME = 6;
  static constexpr float STEP_DT_SCALE =
      static_cast<float>(EVOLUTION_STEPS_PER_FRAME) / STEPS_PER_FRAME;
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

  HS_COLD_MEMBER static GenerativePalette make_palette() {
    return GenerativePalette{EffectPaletteRecipes::gs_reaction_diffusion(
        EffectPaletteRecipes::random_base_turns())};
  }

  // Two 5-bit palette indices and a 6-bit weight for the first palette.
  static void add_pigment(float *weights, uint16_t pigment, float mass) {
    float first = static_cast<float>(pigment >> 10) * (1.0f / 63.0f);
    weights[pigment & 31] += mass * first;
    weights[(pigment >> 5) & 31] += mass * (1.0f - first);
  }

  static uint16_t pack_pigment(const float *weights) {
    int first = 0, second = 1;
    if (weights[second] > weights[first])
      std::swap(first, second);
    for (int i = 2; i < NUM_SEED_CLUSTERS; ++i) {
      if (weights[i] > weights[first]) {
        second = first;
        first = i;
      } else if (weights[i] > weights[second]) {
        second = i;
      }
    }
    float mass = weights[first] + weights[second];
    int mix = mass > 0.0f
                  ? static_cast<int>(63.0f * weights[first] / mass + 0.5f)
                  : 63;
    return static_cast<uint16_t>(first | (second << 5) | (mix << 10));
  }

  HS_COLD_MEMBER void seed_reaction() {
    for (int i = 0; i < RD_N; ++i) {
      state.A[i] = 65535;
      state.B[i] = 0;
      state.pigment[i] = 63u << 10;
    }
    for (int seed = 0; seed < NUM_SEED_CLUSTERS; ++seed) {
      auto palette = make_palette();
      for (int j = 0; j < PALETTE_SIZE; ++j)
        palettes[seed * PALETTE_SIZE + j] =
            palette.get(static_cast<float>(j) / (PALETTE_SIZE - 1)).color;
      int center = hs::rand_int(0, RD_N);
      auto plant = [&](int node) {
        float weights[NUM_SEED_CLUSTERS] = {};
        add_pigment(weights, state.pigment[node], from_q16(state.B[node]));
        weights[seed] += 1.0f;
        state.pigment[node] = pack_pigment(weights);
        state.B[node] = 65535;
      };
      plant(center);
      for_each_neighbor(center, plant);
    }
  }

  void step_pigment(const float *a, const float *b, uint16_t *next) {
    const float DT = params.dt * STEP_DT_SCALE;
    const float DIFFUSION = params.d_b * DT;
    for (int i = 0; i < RD_N; ++i) {
      float weights[NUM_SEED_CLUSTERS] = {};
      float retained = std::max(0.0f, b[i] * (1.0f - RD_K * DIFFUSION -
                                              (params.k + params.feed) * DT) +
                                          a[i] * b[i] * b[i] * DT);
      add_pigment(weights, state.pigment[i], retained);
      Base::template for_each_neighbor<true>(i, [&](int nb) {
        add_pigment(weights, state.pigment[nb], b[nb] * DIFFUSION);
      });
      next[i] = pack_pigment(weights);
    }
    std::copy_n(next, RD_N, state.pigment);
  }

  Pixel palette_color(int seed, float t) const {
    return lut_sample_pixel(palettes + seed * PALETTE_SIZE, PALETTE_SIZE,
                            t * (PALETTE_SIZE - 1));
  }

  struct SeedPalette {
    const Pixel *colors;

    Color4 get(float t) const {
      return Color4(
          lut_sample_pixel(colors, PALETTE_SIZE, t * (PALETTE_SIZE - 1)), 1.0f);
    }
  };

  Pixel modified_palette_color(int seed, float t, float shift,
                               float lightness) const {
    if (shift == 0.0f && lightness == 0.0f)
      return palette_color(seed, t);
    const Color4 SOURCE(palette_color(seed, t), 1.0f);
    const Color4 SHIFTED =
        shift == 0.0f ? SOURCE : hue_rotate_lut_gamut(SOURCE, shift);
    return NoiseShimmerPalette<SeedPalette>::lift_color(SHIFTED, lightness)
        .color;
  }

  void refresh_color_noise() {
    color_noise_cache.refresh(std::span<int8_t, HueNoiseLutView::SIZE>(
                                  color_noise_lut, HueNoiseLutView::SIZE),
                              color_noise, params.noise_scale,
                              color_noise_phase);
  }

  void advance_color_noise() {
    color_noise_phase = math::wrap_t(color_noise_phase + params.noise_speed);
    if (params.hue_shift != 0.0f || params.shimmer != 0.0f)
      refresh_color_noise();
  }

  float sample_color_noise(const math::Vector &direction) const {
    return sample_hue_noise_lut({color_noise_lut, true}, direction);
  }

  template <typename Sample>
  static Pixel mix_pigment(uint16_t pigment, Sample &&sample) {
    Pixel first = sample(pigment & 31);
    if ((pigment >> 10) == 63)
      return first;
    Pixel second = sample((pigment >> 5) & 31);
    return second.lerp16(
        first, static_cast<uint16_t>(((pigment >> 10) * 65535u + 31u) / 63u));
  }

  Pixel pigment_color(uint16_t pigment, float t) const {
    return mix_pigment(pigment,
                       [&](int seed) { return palette_color(seed, t); });
  }

  /**
   * @brief Fades nodes approaching the dissolve frontier, then holds them at
   * rest.
   * @param phase Dissolve progress in [0, 1]; the cleared fraction.
   * @details The band ahead of the frontier eases A/B toward rest over
   * DISSOLVE_FADE_FRACTION of the window. The swept set is re-cleared every
   * frame because B is autocatalytic and otherwise refills from its neighbors.
   */
  void convert_below(float phase) {
    for (int i = 0; i < RD_N; i++) {
      float h =
          math::hash01(static_cast<uint32_t>(i), transition.dissolve_seed);
      if (h < phase) {
        state.A[i] = 65535;
        state.B[i] = 0;
      } else if (h < phase + DISSOLVE_FADE_FRACTION) {
        float keep = (h - phase) * (1.0f / DISSOLVE_FADE_FRACTION);
        state.A[i] = static_cast<uint16_t>(
            65535.0f - (65535.0f - state.A[i]) * keep + 0.5f);
        state.B[i] =
            static_cast<uint16_t>(static_cast<float>(state.B[i]) * keep + 0.5f);
      }
    }
  }

  /**
   * @brief Ends a dissolve by seeding the next reaction at fresh cluster sites.
   * @details The field is already at rest (every node converted), so seeding
   * alone reproduces init()'s starting condition at new random sites;
   * feed/k are the user's and are left alone.
   */
  HS_COLD_MEMBER void start_reaction() {
    transition.dissolve_frames = -1;
    transition.grow_frames = 0;
    transition.stable_frames = 0;
    seed_reaction();
  }

  /**
   * @brief Starts a dissolve at a fresh node ordering.
   */
  void begin_dissolve() {
    transition.dissolve_frames = 0;
    transition.dissolve_seed = static_cast<uint32_t>(hs::random()());
  }

  /**
   * @brief Reports whether the reaction constants moved since the last frame.
   * @return True on the first frame after any of feed/k/dA/dB changes.
   * @details Latches the current values, so a slider drag reports once per
   * distinct value. Speed is excluded: it sets the integration rate, not the
   * reaction.
   */
  bool reaction_edited() {
    bool changed =
        params.feed != transition.last_feed || params.k != transition.last_k ||
        params.d_a != transition.last_d_a || params.d_b != transition.last_d_b;
    transition.last_feed = params.feed;
    transition.last_k = params.k;
    transition.last_d_a = params.d_a;
    transition.last_d_b = params.d_b;
    return changed;
  }

  /**
   * @brief Runs the reaction lifecycle: edit and stabilization detection, then
   *        the dissolve.
   * @param mean_db Mean per-node |dB| between the frame's start and end states.
   *        Substep motion that cancels over the frame reads as zero.
   * @details Dissolves when the user edits the reaction, or once the field has
   * stalled for STABLE_HOLD_FRAMES. Edits mid-dissolve are absorbed by the
   * in-flight one, which reseeds into whatever the constants read at its end.
   */
  void advance_transition(float mean_db) {
    if (transition.dissolve_frames >= 0) {
      transition.dissolve_frames++;
      convert_below(static_cast<float>(transition.dissolve_frames) /
                    DISSOLVE_FRAMES);
      // Latch edits made mid-dissolve; this dissolve already covers them.
      reaction_edited();
      if (transition.dissolve_frames >= DISSOLVE_FRAMES)
        start_reaction();
      return;
    }
    if (reaction_edited()) {
      begin_dissolve();
      return;
    }
    transition.grow_frames++;
    if (transition.grow_frames < MIN_GROW_FRAMES)
      return;
    // Scale the calibrated floor with the requested timestep.
    const float floor_db = MEAN_DB_STABLE * (params.dt * (1.0f / DEFAULT_DT)) *
                           (static_cast<float>(EVOLUTION_STEPS_PER_FRAME) /
                            BASELINE_STEPS_PER_FRAME);
    transition.stable_frames =
        mean_db < floor_db ? transition.stable_frames + 1 : 0;
    if (transition.stable_frames >= STABLE_HOLD_FRAMES)
      begin_dissolve();
  }

  /**
   * @brief Advances one Gray-Scott substep into the next buffers (Jacobi).
   * @param c_a Current A field (read-only), float in [0, 1] per node.
   * @param c_b Current B field (read-only), float in [0, 1] per node.
   * @param n_a Next A field (write target), float in [0, 1] per node.
   * @param n_b Next B field (write target), float in [0, 1] per node.
   * @details Gray-Scott: dA/dt = dA·∇²A - A·B² + feed·(1-A);
   * dB/dt = dB·∇²B + A·B² - (k+feed)·B. Double-buffered Jacobi: reads current
   * buffers, writes next. The
   * [0, 1] clamp saturates explicit-Euler overshoot past the stability bound
   * (see "Speed"); substeps stay in float so the Q16 state quantizes once per
   * frame, not once per substep.
   */
  HS_O3_FN void step_physics(const float *__restrict c_a,
                             const float *__restrict c_b, float *__restrict n_a,
                             float *__restrict n_b) {
    step_physics_nodes(c_a, c_b, [&](int i, float a, float b) {
      n_a[i] = a;
      n_b[i] = b;
    });
  }

  static constexpr int PHYSICS_NEIGHBOR_REACH = 144;

  HS_COLD_MEMBER static void validate_physics_neighbors(
      const ReactionGraph::NeighborRun *runs = ReactionGraph::neighbor_runs,
      unsigned count = ReactionGraph::NEIGHBOR_RUN_COUNT) {
    for (unsigned r = 0; r < count; ++r)
      for (int delta : runs[r].delta)
        HS_CHECK(delta >= -PHYSICS_NEIGHBOR_REACH,
                 "GS neighbor exceeds delayed-write history");
  }

  /** @brief Advances float A/B in place after their last stencil read. */
  HS_O3_FN void step_physics_inplace(float *a, float *b) {
    constexpr int HISTORY_SIZE = PHYSICS_NEIGHBOR_REACH + 1;
    std::array<float, HISTORY_SIZE> pending_a, pending_b;
    step_physics_nodes(a, b, [&](int i, float next_a, float next_b) {
      pending_a[i % HISTORY_SIZE] = next_a;
      pending_b[i % HISTORY_SIZE] = next_b;
      // Node i is the last possible reader of i - PHYSICS_NEIGHBOR_REACH.
      if (i >= PHYSICS_NEIGHBOR_REACH) {
        int done = i - PHYSICS_NEIGHBOR_REACH;
        a[done] = pending_a[done % HISTORY_SIZE];
        b[done] = pending_b[done % HISTORY_SIZE];
      }
    });
    for (int i = RD_N - PHYSICS_NEIGHBOR_REACH; i < RD_N; ++i) {
      a[i] = pending_a[i % HISTORY_SIZE];
      b[i] = pending_b[i % HISTORY_SIZE];
    }
  }

  template <typename StoreFn>
  HS_O3_FN void step_physics_nodes(const float *c_a, const float *c_b,
                                   StoreFn &&store) {
    const float feed = params.feed;
    const float KILL_RATE = params.k;
    const float d_a = params.d_a;
    const float d_b = params.d_b;
    const float dt = params.dt * STEP_DT_SCALE;
    int i = 0;
    for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT; ++r) {
      const auto &run = ReactionGraph::neighbor_runs[r];
      for (; i < run.end; ++i) {
        float a = c_a[i];
        float b = c_b[i];

        float l_a, l_b;
#if defined(__arm__)
        l_a = __builtin_fmaf(-static_cast<float>(RD_K), a,
                             c_a[i + run.delta[0]] + c_a[i + run.delta[1]]);
        l_b = __builtin_fmaf(-static_cast<float>(RD_K), b,
                             c_b[i + run.delta[0]] + c_b[i + run.delta[1]]);
        for (int k = 2; k < RD_K; ++k) {
          int ni = i + run.delta[k];
          l_a += c_a[ni];
          l_b += c_b[ni];
          // Preserve the Cortex-M7 Laplacian's floating-point order.
          if (k < RD_K - 1)
            asm("" : "+t"(l_a), "+t"(l_b));
        }
#elif defined(__FAST_MATH__) || (defined(_M_FP_FAST) && _M_FP_FAST)
        l_a = -RD_K * a;
        l_b = -RD_K * b;
        for (int k = 0; k < RD_K; ++k) {
          int ni = i + run.delta[k];
          l_a += c_a[ni];
          l_b += c_b[ni];
        }
#else
        float sum_a = 0.0f, sum_b = 0.0f;
        for (int k = 0; k < RD_K; ++k) {
          int ni = i + run.delta[k];
          sum_a += c_a[ni];
          sum_b += c_b[ni];
        }
        l_a = sum_a - RD_K * a;
        l_b = sum_b - RD_K * b;
#endif

        float abb = a * b * b;
        float next_a = hs::clamp(a + (d_a * l_a - abb + feed * (1.0f - a)) * dt,
                                 0.0f, 1.0f);
        float next_b = hs::clamp(
            b + (d_b * l_b + abb - (KILL_RATE + feed) * b) * dt, 0.0f, 1.0f);
        store(i, next_a, next_b);
      }
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
   * @details Off the render path: shade_pixel gathers the stencil once per
   * pixel and re-weights it inline. This one-sample form is the oracle
   * tests/test_effects.h bounds that shared stencil against.
   */
  float interpolate_b(const math::Vector &p, int seed,
                      const math::Vector *nodes) const {
    float tw = 0, wb = 0;
    refine_and_accumulate(p, nodes, seed, [&](int i, float w) {
      wb += from_q16(state.B[i]) * w;
      tw += w;
    });
    // Zero remains cullable by the test oracle's b threshold.
    if (tw <= Base::KERNEL_MIN_TOTAL_WEIGHT)
      return 0.0f;
    return wb / tw;
  }

  /**
   * @brief Shades one pixel's four sub-samples through an inlinable typed path.
   * @tparam Grid Scan::Shader::SsaaGrid type supplying the sub-pixel offsets.
   * @param seed Cubemap-LUT seed node id, or -1 for a culled pixel.
   * @param center_rv World-space direction at the pixel center.
   * @param world_nodes Oriented lattice node positions.
   * @param grid Row's SSAA sub-pixel grid.
   * @param x Pixel column.
   * @return The finished alpha-premultiplied pixel.
   * @details Accepts seeds inside a proven nearest-node radius immediately;
   * boundary pixels check all six neighbors. The center stencil is shared
   * across the four sub-pixel samples. The row offset is
   * 0.25 * RADIANS_PER_ROW<H>; stencil reuse can exceed one node spacing
   * at low vertical resolutions.
   */
  template <typename Grid>
  HS_O3_FN Pixel shade_pixel(int seed, const math::Vector &center_rv,
                             const math::Vector *world_nodes, const Grid &grid,
                             int x) const {
    if (seed < 0)
      return Pixel(0, 0, 0);

    float noise_value = 0.0f;
    if (params.hue_shift != 0.0f || params.shimmer != 0.0f)
      noise_value =
          sample_color_noise(Base::inverse_orientation.apply(center_rv));
    const float HUE_SHIFT = noise_value * params.hue_shift;
    const float LIGHTNESS = std::max(0.0f, noise_value) * params.shimmer;
    int center =
        Base::template refine_render_center<true>(center_rv, world_nodes, seed);
    constexpr uint32_t SAMPLES = Grid::SAMPLES;
    static_assert(SAMPLES == 4);
    const float st = math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::sin_theta[x];
    const float ct = math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::cos_theta(x);
    // Each horizontal pair is m +/- d, with m dot d = 0.
    math::Vector midpoints[2];
    float offset_squared[2], cross_scale[2];
    for (int row = 0; row < 2; ++row) {
      float sp = grid.sin_phi[row];
      midpoints[row] =
          math::Vector(sp * ct * grid.cos_dtheta, grid.cos_phi[row],
                       sp * st * grid.cos_dtheta);
      float offset = sp * grid.sin_dtheta;
      float ox = -st * offset, oz = ct * offset;
      offset_squared[row] = ox * ox + oz * oz;
      cross_scale[row] = 2.0f * offset * Base::INV_R2;
    }
    float weights[SAMPLES] = {}, weighted_b[SAMPLES] = {};
    float pigment_weights[RD_K + 1][SAMPLES] = {};
    const auto &stencil_run =
        ReactionGraph::neighbor_runs[ReactionGraph::neighbor_run_index[center]];
    for (int j = 0; j < RD_K + 1; ++j) {
      int ni = j == 0 ? center : center + stencil_run.delta[j - 1];
      const math::Vector &p = world_nodes[ni];
      float b = state.B[ni];
      float tangent = st * p.x - ct * p.z;
      for (int row = 0; row < 2; ++row) {
        float dx = midpoints[row].x - p.x;
        float dy = midpoints[row].y - p.y;
        float dz = midpoints[row].z - p.z;
        float base = dx * dx + dy * dy + dz * dz + offset_squared[row];
        float base_u = 1.0f - base * Base::INV_R2;
        float cross = cross_scale[row] * tangent;
        for (int col = 0; col < 2; ++col) {
          int i = 2 * row + col;
          float u = col == 0 ? base_u - cross : base_u + cross;
          if (u > 0.0f) {
            float w = u * u;
            pigment_weights[j][i] = b * w;
            weighted_b[i] += b * w;
            weights[i] += w;
          }
        }
      }
    }
    float accum_r = 0.0f, accum_g = 0.0f, accum_b = 0.0f;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC unroll 1
#endif
    for (int i = 0; i < Grid::SAMPLES; ++i) {
      float tw = weights[i], wb = weighted_b[i];
      if (tw <= Base::KERNEL_MIN_TOTAL_WEIGHT)
        continue;
      float b = wb * (Q16_INV / tw);
      if (b < B_CULL_THRESHOLD)
        continue;

      float t = hs::clamp((b - B_COLOR_FLOOR) * B_COLOR_SCALE, 0.0f, 1.0f);
      float scale = 1.0f / (wb * SAMPLES);
      Pixel color_cache[NUM_SEED_CLUSTERS];
      uint32_t cached = 0;
      auto sample_palette = [&](int palette_id) {
        uint32_t bit = 1u << palette_id;
        if (!(cached & bit)) {
          color_cache[palette_id] =
              modified_palette_color(palette_id, t, HUE_SHIFT, LIGHTNESS);
          cached |= bit;
        }
        return color_cache[palette_id];
      };
      for (int j = 0; j < RD_K + 1; ++j) {
        float weight = pigment_weights[j][i] * scale;
        if (weight <= 0.0f)
          continue;
        int ni = j == 0 ? center : center + stencil_run.delta[j - 1];
        Pixel rgb = mix_pigment(state.pigment[ni], sample_palette);
        accum_r += rgb.r * weight;
        accum_g += rgb.g * weight;
        accum_b += rgb.b * weight;
      }
    }
    return Pixel(
        static_cast<uint16_t>(hs::clamp(accum_r + 0.5f, 0.0f, 65535.0f)),
        static_cast<uint16_t>(hs::clamp(accum_g + 0.5f, 0.0f, 65535.0f)),
        static_cast<uint16_t>(hs::clamp(accum_b + 0.5f, 0.0f, 65535.0f)));
  }

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
   * culling on !hot2[seed] is exact, not approximate.
   */
  HS_O3_FN static void fill_hot_flags(const uint16_t *b, uint8_t *hot1,
                                      uint8_t *hot2, int count,
                                      uint16_t threshold) {
    int i = 0;
    for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT && i < count;
         ++r) {
      const auto &run = ReactionGraph::neighbor_runs[r];
      const int end = std::min<int>(run.end, count);
      for (; i < end; ++i) {
        bool hot = b[i] >= threshold;
        for (int k = 0; k < RD_K && !hot; ++k)
          hot = b[i + run.delta[k]] >= threshold;
        hot1[i] = hot;
      }
    }
    i = 0;
    for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT && i < count;
         ++r) {
      const auto &run = ReactionGraph::neighbor_runs[r];
      const int end = std::min<int>(run.end, count);
      for (; i < end; ++i) {
        bool hot = hot1[i];
        for (int k = 0; k < RD_K && !hot; ++k)
          hot = hot1[i + run.delta[k]];
        hot2[i] = hot;
      }
    }
  }

  /**
   * @brief Advances the sim STEPS_PER_FRAME substeps and rasterizes the B field.
   * @param canvas Destination canvas to draw the sphere into.
   * @details Rasterizes the B field onto the sphere via the orientation-aware
   * SSAA shader pipeline after advancing the simulation.
   */
  void render(Canvas &canvas) {
    HS_PROFILE(grd_render);
    advance_color_noise();
    ScratchScope frame_guard(scratch_arena_a);
    float mean_db = 0.0f;
    {
      // Q16 quantization occurs once per frame.
      HS_PROFILE(grd_simulate);
      ScratchScope physics_guard(scratch_arena_a);
      float *cur_a = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *cur_b = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      uint16_t *next_pigment = static_cast<uint16_t *>(
          scratch_arena_a.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));

      for (int i = 0; i < RD_N; i++) {
        cur_a[i] = from_q16(state.A[i]);
        cur_b[i] = from_q16(state.B[i]);
      }
      for (int step = 0; step < STEPS_PER_FRAME; ++step) {
        step_pigment(cur_a, cur_b, next_pigment);
        step_physics_inplace(cur_a, cur_b);
      }
      uint32_t db_sum_q16 = 0;
      for (int i = 0; i < RD_N; i++) {
        state.A[i] = to_q16(cur_a[i]);
        uint16_t next_b = to_q16(cur_b[i]);
        int db = static_cast<int>(next_b) - state.B[i];
        db_sum_q16 += static_cast<uint32_t>(db < 0 ? -db : db);
        state.B[i] = next_b;
      }
      mean_db = static_cast<float>(db_sum_q16) * (Q16_INV / RD_N);
    }
    advance_transition(mean_db);

    // Physics scratch is popped; the raster phase reuses the arena for the
    // oriented lattice so the kernel walks stay in world space, plus the
    // two-ring cull flags.
    HS_PROFILE(grd_rasterize);
    auto lattice = [this] {
      HS_PROFILE(grd_orient);
      return orient_lattice();
    }();
    math::Vector *world_nodes = lattice.get();
    uint8_t *hot1 = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *hot2 = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    {
      HS_PROFILE(grd_cull_flags);
      fill_hot_flags(state.B, hot1, hot2, RD_N, to_q16(B_CULL_THRESHOLD));
    }

    // Seed the cubemap lookup once per pixel center; a seed whose two-ring
    // sits below the render floor is culled for the whole pixel (v0 = -1).
    auto vertex_shader = [&](Fragment &frag) {
      if (!hot2[static_cast<int>(frag.v0)])
        frag.v0 = -1.0f;
    };

    auto pixel_shader = [&](Fragment &frag, const auto &grid, int x) -> Pixel {
      return shade_pixel(static_cast<int>(frag.v0), frag.pos, world_nodes, grid,
                         x);
    };

    {
      HS_PROFILE(grd_shader_draw);
      rasterize_lattice(canvas, vertex_shader, pixel_shader);
    }
  }

  /**
   * @brief Persistent Q16 state buffers for the two species.
   */
  struct {
    uint16_t *A = nullptr,
             *B = nullptr; /**< Per-node A/B concentrations, Q16. */
    uint16_t *pigment = nullptr;
  } state;

  /**
   * @brief Auto-transition state: current reaction's lifetime and dissolve
   *        progress.
   */
  struct {
    int grow_frames = 0;        /**< Frames since this reaction was seeded. */
    int stable_frames = 0;      /**< Consecutive sub-floor frames. */
    int dissolve_frames = -1;   /**< Dissolve progress; -1 when inactive. */
    uint32_t dissolve_seed = 0; /**< Per-transition node-order hash seed. */
    /** Reaction constants as of the last frame; reaction_edited() latches them
     *  to spot a user edit, seeded from params in init() so frame 1 is
     *  clean. */
    float last_feed = 0.0f, last_k = 0.0f, last_d_a = 0.0f, last_d_b = 0.0f;
  } transition;

  /** @brief Per-seed linear RGB ramps sampled by B concentration. */
  Pixel *palettes = nullptr;
  int8_t *color_noise_lut = nullptr;
  FastNoiseLite color_noise;
  HueNoiseBakeCache color_noise_cache;
  float color_noise_phase = 0.0f;

  /**
   * @brief GUI-tunable Gray-Scott parameters.
   */
  struct Params {
    float feed = 0.04f; /**< Feed rate of A. */
    float k = 0.06f;    /**< Kill rate of B. */
    float d_a = 0.02f;  /**< Diffusion coefficient of A. */
    float d_b = 0.01f;  /**< Diffusion coefficient of B. */
    float dt = 2.5f;    /**< Integration timestep (Speed slider). */
    float noise_speed = 0.00040800002f;
    float noise_scale = 4.941984f;
    float hue_shift = 0.2f;
    float shimmer = 0.4f;
  } params;
  static_assert(Params{}.dt == DEFAULT_DT,
                "the stabilization floor is calibrated at the Speed default");
};
