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

#include "core/platform/effect_identity.h"
#include <algorithm>
#include <bit>
#include <cmath>
#include <utility>
#include "effects/common/palette_recipes.h"
#include "core/color/noise_shimmer_palette.h"
#include "core/engine/engine.h"
#include "effects/common/reaction_diffusion.h"

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
 * feed/kill), producing spots/stripes/mazes. Persistent state is Q16;
 * substeps integrate in float and quantize back once per frame.
 *
 * A reaction runs until its field has all but stopped moving, then dissolves
 * and reseeds at new sites. Editing the constants dissolves the current field
 * too. Each node carries its two strongest seed palettes and their mixing
 * weight.
 */
template <int W, int H>
class GSReactionDiffusion
    : public ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H>;
  friend Base;

  using Base::for_each_neighbor;
  using Base::from_q16;
  using Base::orient_lattice;
  using Base::prepare_lattice_basis;
  using Base::Q16_INV;
  using Base::RD_K;
  using Base::RD_N;
  using Base::register_param;
  using Base::to_q16;
  using typename Base::LatticeBasis;
  using typename Base::FloatRgb;

public:
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::GSReactionDiffusion::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  GSReactionDiffusion() = default;

  /**
   * @brief One-time setup: arenas, GUI params, A/B state, cubemap LUT, lattice.
   */
  void init() override {
    constexpr size_t PALETTE_BYTES =
        NUM_SEED_CLUSTERS *
            (PALETTE_SIZE * sizeof(Pixel) +
             COLOR_VALUE_STEPS * COLOR_NOISE_STEPS * sizeof(FloatRgb)) +
        2 * alignof(Pixel) + 2 * HueNoiseLutView::SIZE * sizeof(int8_t);
    constexpr size_t PERSISTENT_BYTES = 192 * 1024 + 32;
    constexpr size_t PHYSICS_SCRATCH_BYTES =
        2u * RD_N * sizeof(float) + RD_N * sizeof(uint16_t) +
        2u * PHYSICS_HISTORY_SIZE * sizeof(float);
    constexpr size_t RASTER_SCRATCH_BYTES =
        RD_N * sizeof(math::Vector) + 2u * RD_N * sizeof(uint8_t);
    // Slab tenants: basis, marked, hot1, hot2, cur_a, cur_b, next_pigment,
    // pending_a, pending_b and start_b.
    constexpr size_t SLAB_SCRATCH_BYTES =
        sizeof(LatticeBasis) + 3u * HOT_WORDS * sizeof(uint32_t) +
        2u * RD_N * sizeof(float) + 2u * RD_N * sizeof(uint16_t) +
        2u * PHYSICS_HISTORY_SIZE * sizeof(float);
    constexpr size_t SCRATCH_BYTES = std::max(
        {PHYSICS_SCRATCH_BYTES, RASTER_SCRATCH_BYTES, SLAB_SCRATCH_BYTES});
    constexpr size_t SCRATCH_TENANTS = 10;
    Base::template configure_rd_arenas<uint16_t, 3, PERSISTENT_BYTES,
                                       PALETTE_BYTES, SCRATCH_BYTES,
                                       SCRATCH_TENANTS, true>();

    register_param("Feed", &params.feed, 0.0f, 0.1f);
    register_param("Kill", &params.k, 0.0f, 0.1f);
    // Each substep clamps to [0,1], so the joint Speed/diffusion maximum may
    // exceed the linear diffusion bound.
    register_param("dA", &params.d_a, 0.0f, 0.05f);
    register_param("dB", &params.d_b, 0.0f, 0.05f);
    register_param("Speed", &params.dt, 0.1f, 3.0f);
    register_param("Noise Speed", &params.noise_speed, -0.001f, 0.001f);
    register_param("Noise Scale", &params.noise_scale, 1.0f / 64.0f, 8.0f);
    register_param("Hue Shift", &params.hue_shift, -4.0f, 4.0f);
    register_param("Shimmer", &params.shimmer, 0.0f, 1.0f);

    state.A = persistent_arena.allocate_n<uint16_t>(RD_N);
    state.B = persistent_arena.allocate_n<uint16_t>(RD_N);
    state.pigment = persistent_arena.allocate_n<uint16_t>(RD_N);
    palettes =
        persistent_arena.allocate_n<Pixel>(NUM_SEED_CLUSTERS * PALETTE_SIZE);
    shimmer_cache.attach(persistent_arena.allocate_n<FloatRgb>(
        NUM_SEED_CLUSTERS * COLOR_VALUE_STEPS * COLOR_NOISE_STEPS));
    color_noise_lut =
        persistent_arena.allocate_n<int8_t>(HueNoiseLutView::SIZE);
    previous_noise_lut =
        persistent_arena.allocate_n<int8_t>(HueNoiseLutView::SIZE);
    color_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    color_noise.SetSeed(6047);
    color_noise.SetFrequency(1.0f);
    refresh_color_noise();

    validate_physics_neighbors();
    Base::template init_lattice<true>();
    seed_reaction();
    refresh_color_palettes(true);
    reaction_edited(); // latch the defaults; frame 1 is not an edit
  }

private:
  friend struct ::hs_test::effects_tests::GSWhiteBox;

  /**
   * @brief Fully-saturated B blobs seeded per reaction.
   * @details Nucleation sites for the GS instability to grow from; without them
   * the uniform A=1/B=0 field never moves.
   */
  static constexpr int NUM_SEED_CLUSTERS = 30;
  static constexpr int PALETTE_SIZE = 32;
  static constexpr int COLOR_NOISE_STEPS = 15;
  static constexpr int COLOR_VALUE_STEPS = 16;
  static constexpr float CACHED_HUE_LIMIT = 0.25f;
  static constexpr float CACHED_SHIMMER_LIMIT = 0.4f;

  /**
   * @brief Per-node pigment codec: two 5-bit seed-palette ids and a 6-bit
   * quantized mass of the first palette.
   */
  struct Pigment {
    static constexpr int ID_BITS = 5;
    static constexpr int ID_MASK = (1 << ID_BITS) - 1;
    static constexpr int MASS_SHIFT = 2 * ID_BITS;
    /** @brief Mass of a pigment that is entirely its first palette. */
    static constexpr int FULL_MASS = 63;
    /** @brief Keeps the first id and the mass; clears the second id. */
    static constexpr uint16_t SOLID_KEY_MASK =
        static_cast<uint16_t>(~(ID_MASK << ID_BITS));
    static_assert(NUM_SEED_CLUSTERS <= ID_MASK + 1);
    static_assert(FULL_MASS << MASS_SHIFT <= 0xffff);

    /** @brief Decoded pigment fields. */
    struct Unpacked {
      int first;  ///< First seed-palette id.
      int second; ///< Second seed-palette id.
      int mass;   ///< First palette's share, in [0, FULL_MASS].
    };

    __attribute__((always_inline)) static constexpr uint16_t
    pack(int first, int second, int mass) {
      return static_cast<uint16_t>(first | (second << ID_BITS) |
                                   (mass << MASS_SHIFT));
    }
    __attribute__((always_inline)) static constexpr Unpacked
    unpack(uint16_t pigment) {
      return {pigment & ID_MASK, (pigment >> ID_BITS) & ID_MASK,
              pigment >> MASS_SHIFT};
    }
    /** @brief The first palette's share of a pigment with `mass`, in [0, 1]. */
    __attribute__((always_inline)) static constexpr float first_mass(int mass) {
      return static_cast<float>(mass) * (1.0f / FULL_MASS);
    }
    /** @brief Quantizes `first` as a share of `total`. */
    __attribute__((always_inline)) static int quantize_mass(float first,
                                                            float total) {
      return static_cast<int>(static_cast<float>(FULL_MASS) * first / total +
                              0.5f);
    }
  };
  /** @brief Substep budget used to calibrate the stabilization threshold. */
  static constexpr int BASELINE_STEPS_PER_FRAME = 16;
  /** @brief Rendered frames the dissolve takes to convert every node back to
   * rest. */
  static constexpr int DISSOLVE_FRAMES = 103;
  /**
   * @brief Rendered frames a fresh reaction runs before the stabilization
   * detector arms.
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
   * BASELINE_STEPS_PER_FRAME.
   * @details A calibration heuristic: Q16 quantization makes the low-Speed
   * response nonlinear.
   */
  static constexpr float MEAN_DB_STABLE = 2.0e-4f;
  /** @brief Speed the stabilization floor is calibrated at. */
  static constexpr float DEFAULT_DT = 2.5f;
  /** @brief Base-dt substep equivalents advanced per rendered frame. */
  static constexpr int EVOLUTION_STEPS_PER_FRAME = 10;
  /** @brief Euler integrations performed per rendered frame. */
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
   * @brief Cull threshold; coincides with the color floor, so a pixel is
   * either fully transparent or on the gradient.
   */
  static constexpr float B_CULL_THRESHOLD = B_COLOR_FLOOR;

  HS_COLD_MEMBER static GenerativePalette make_palette() {
    return GenerativePalette{EffectPaletteRecipes::gs_reaction_diffusion(
        EffectPaletteRecipes::random_base_turns())};
  }

  static void add_pigment(float *weights, uint16_t pigment, float mass) {
    const auto PIGMENT = Pigment::unpack(pigment);
    float first = Pigment::first_mass(PIGMENT.mass);
    weights[PIGMENT.first] += mass * first;
    weights[PIGMENT.second] += mass * (1.0f - first);
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
    int mix = mass > 0.0f ? Pigment::quantize_mass(weights[first], mass)
                          : Pigment::FULL_MASS;
    return Pigment::pack(first, second, mix);
  }

  HS_COLD_MEMBER void seed_reaction(int first = 0,
                                    int end = NUM_SEED_CLUSTERS) {
    HS_PROFILE(grd_seed_reaction);
    if (first == 0) {
      for (int i = 0; i < RD_N; ++i) {
        state.A[i] = 65535;
        state.B[i] = 0;
        state.pigment[i] = Pigment::pack(0, 0, Pigment::FULL_MASS);
      }
      shimmer_cache.invalidate();
    }
    for (int seed = first; seed < end; ++seed) {
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

  template <int PIGMENT_STEPS = 1>
  HS_HOT_FLASH_MEMBER void step_pigment(const float *a, const float *b,
                                        uint16_t *next) {
    HS_PROFILE(grd_pigment);
    const float DT = params.dt * STEP_DT_SCALE;
    static_assert(RD_K == 6, "GS stability bound assumes a 6-NN lattice");
    const float DIFFUSION = params.d_b * DT;
    float weights[NUM_SEED_CLUSTERS] = {};
    int i = 0;
    for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT; ++r) {
      const auto &run = ReactionGraph::neighbor_runs[r];
      for (; i < run.end; ++i) {
        float retained = fmaxf(0.0f, b[i] * (1.0f - RD_K * DIFFUSION -
                                             (params.k + params.feed) * DT) +
                                         a[i] * b[i] * b[i] * DT);
        float diffusion = DIFFUSION;
        if constexpr (PIGMENT_STEPS > 1) {
          float incoming = 0.0f;
          for (int k = 0; k < RD_K; ++k)
            incoming += b[i + run.delta[k]] * DIFFUSION;
          float total = retained + incoming;
          float r = total > 0.0f ? retained / total : 0.0f;
          float power = 1.0f, scale = 1.0f;
          for (int step = 1; step < PIGMENT_STEPS; ++step) {
            power *= r;
            scale += power;
          }
          retained *= power;
          diffusion *= scale;
        }
        uint16_t uniform = state.pigment[i] & Pigment::SOLID_KEY_MASK;
        bool same = (uniform >> Pigment::MASS_SHIFT) == Pigment::FULL_MASS;
        for (int k = 0; k < RD_K && same; ++k)
          same = (state.pigment[i + run.delta[k]] & Pigment::SOLID_KEY_MASK) ==
                 uniform;
        if (same) {
          float mass = retained;
          for (int k = 0; k < RD_K; ++k)
            mass += b[i + run.delta[k]] * diffusion;
          int first = mass > 0.0f ? Pigment::unpack(uniform).first : 0;
          int second = first == 0 ? 1 : 0;
          next[i] = Pigment::pack(first, second, Pigment::FULL_MASS);
          continue;
        }
        uint32_t touched = 0;
        auto add = [&](int id, float mass) __attribute__((always_inline)) {
          weights[id] += mass;
          touched |= 1u << id;
        };
        auto gather = [&](int node, float mass) __attribute__((always_inline)) {
          const auto PIGMENT = Pigment::unpack(state.pigment[node]);
          float first = Pigment::first_mass(PIGMENT.mass);
          add(PIGMENT.first, mass * first);
          add(PIGMENT.second, mass * (1.0f - first));
        };
        gather(i, retained);
        for (int k = 0; k < RD_K; ++k) {
          int nb = i + run.delta[k];
          gather(nb, b[nb] * diffusion);
        }
        int first = 0, second = 1;
        // Nonnegative IEEE-754 weights sort by their unsigned representations.
        auto weight_bits = [&](int id) __attribute__((always_inline)) {
          return std::bit_cast<uint32_t>(weights[id]) & 0x7fffffffu;
        };
        uint32_t first_bits = weight_bits(0);
        uint32_t second_bits = weight_bits(1);
        weights[0] = weights[1] = 0.0f;
        if (second_bits > first_bits) {
          std::swap(first, second);
          std::swap(first_bits, second_bits);
        }
        uint32_t remaining = touched & ~3u;
        while (remaining) {
          int id = __builtin_ctz(remaining);
          remaining &= remaining - 1;
          uint32_t bits = weight_bits(id);
          weights[id] = 0.0f;
          if (bits > first_bits) {
            second = first;
            second_bits = first_bits;
            first = id;
            first_bits = bits;
          } else if (bits > second_bits) {
            second = id;
            second_bits = bits;
          }
        }
        float first_mass = std::bit_cast<float>(first_bits);
        float second_mass = std::bit_cast<float>(second_bits);
        float mass = first_mass + second_mass;
        int mix = (first_bits | second_bits) != 0
                      ? Pigment::quantize_mass(first_mass, mass)
                      : Pigment::FULL_MASS;
        next[i] = Pigment::pack(first, second, mix);
      }
    }
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

  HS_HOT_FLASH_MEMBER Pixel modified_palette_color(int seed, float t,
                                                   float shift,
                                                   float lightness) const {
    if (shift == 0.0f && lightness == 0.0f)
      return palette_color(seed, t);
    const Color4 SOURCE(palette_color(seed, t), 1.0f);
    const Color4 SHIFTED =
        shift == 0.0f ? SOURCE : hue_rotate_lut_gamut(SOURCE, shift);
    return NoiseShimmerPalette<SeedPalette>::lift_color(SHIFTED, lightness)
        .color;
  }

  /**
   * @brief Hue-shifted, shimmer-lifted seed palettes baked over a grid of noise
   * rows, keyed on the latched hue and shimmer.
   * @details Rows bake centre-out across frames. A row pair that is not baked
   * yet, or a hue or shimmer past the cached limits, reads as exact.
   */
  class ShimmerCache {
  public:
    /** @brief Noise row pair to blend and whether to bypass the cache. */
    struct Selection {
      int row;    ///< Lower noise row; the pair is row and row + 1.
      bool exact; ///< True when the pair is unbaked or past the cached limits.
    };

    /** @brief Binds the bake storage, one row per seed and noise step. */
    void attach(FloatRgb *storage) { palettes = storage; }

    /** @brief Marks every row stale after the source palettes change. */
    void invalidate() { valid = false; }

    /** @brief True when the latched hue or shimmer is nonzero. */
    bool noise_enabled() const { return enabled; }

    /**
     * @brief Latches `hue` and `shimmer` and bakes the next rows.
     * @param hue Hue shift to latch; zero with zero `shimmer` disables.
     * @param shimmer Shimmer amount to latch.
     * @param complete Bakes every remaining row.
     * @param bake Source color for (seed, t, hue shift, lightness).
     */
    template <typename Bake>
    __attribute__((always_inline)) void refresh(float hue, float shimmer,
                                                bool complete, Bake &&bake) {
      enabled = hue != 0.0f || shimmer != 0.0f;
      exact = !enabled || fabsf(hue) > CACHED_HUE_LIMIT ||
              shimmer > CACHED_SHIMMER_LIMIT;
      const bool RESTART =
          !valid || latched_hue != hue || latched_shimmer != shimmer;
      if (RESTART) {
        latched_hue = hue;
        latched_shimmer = shimmer;
        valid = true;
        rows = 0;
        next_row = 0;
      }
      if (next_row == COLOR_NOISE_STEPS || exact)
        return;
      HS_PROFILE(grd_color_palette);
      int count = complete ? COLOR_NOISE_STEPS : RESTART ? 5 : 2;
      for (; count > 0 && next_row < COLOR_NOISE_STEPS; --count) {
        int index = next_row++;
        int n = COLOR_NOISE_STEPS / 2 +
                ((index & 1) ? -(index + 1) / 2 : index / 2);
        float noise = -1.0f + 2.0f * n / (COLOR_NOISE_STEPS - 1);
        for (int seed = 0; seed < NUM_SEED_CLUSTERS; ++seed)
          for (int t = 0; t < COLOR_VALUE_STEPS; ++t)
            palettes[(seed * COLOR_NOISE_STEPS + n) * COLOR_VALUE_STEPS + t] =
                FloatRgb(bake(seed,
                              static_cast<float>(t) / (COLOR_VALUE_STEPS - 1),
                              noise * hue, fmaxf(noise, 0.0f) * shimmer));
        rows |= 1u << n;
      }
    }

    /** @brief The noise row pair below `position`, and whether to shade it
     * exactly. */
    __attribute__((always_inline)) Selection selection(float position) const {
      const int ROW =
          std::min(static_cast<int>(position), COLOR_NOISE_STEPS - 2);
      return {ROW, !valid || (rows & (3u << ROW)) != (3u << ROW) || exact};
    }

    /** @brief The baked colors of `seed` at noise row `noise_row`. */
    __attribute__((always_inline)) const FloatRgb *row(int seed,
                                                       int noise_row) const {
      return palettes +
             (seed * COLOR_NOISE_STEPS + noise_row) * COLOR_VALUE_STEPS;
    }

  private:
    friend struct ::hs_test::effects_tests::GSWhiteBox;

    FloatRgb *palettes = nullptr;
    bool valid = false;
    bool exact = true;
    bool enabled = false;
    uint16_t rows = 0;
    uint8_t next_row = 0;
    float latched_hue = 0.0f;
    float latched_shimmer = 0.0f;
  };

  HS_HOT_FLASH_MEMBER void refresh_color_palettes(bool complete = false) {
    shimmer_cache.refresh(
        params.hue_shift, params.shimmer, complete,
        [this](int seed, float t, float shift, float lightness)
            __attribute__((always_inline)) {
              return modified_palette_color(seed, t, shift, lightness);
            });
  }

  HS_FLASH_INLINE void refresh_color_noise() {
    color_noise_cache.refresh<true>(std::span<int8_t, HueNoiseLutView::SIZE>(
                                        color_noise_lut, HueNoiseLutView::SIZE),
                                    color_noise, params.noise_scale,
                                    color_noise_phase);
  }

  void advance_color_noise() {
    const float PREVIOUS_PHASE = color_noise_phase;
    color_noise_phase = math::wrap_t(color_noise_phase + params.noise_speed);
    if (params.hue_shift != 0.0f || params.shimmer != 0.0f) {
      if (!Spacetime::flat_surface())
        retain_color_noise(PREVIOUS_PHASE);
      refresh_color_noise();
    }
  }

  /**
   * @brief Keeps the previous frame's noise table before the next bake.
   * @param previous_phase Noise phase of the previous frame.
   * @details The table is history only when it was baked at that phase and
   * the current scale; otherwise slabs shade the current table alone.
   */
  void retain_color_noise(float previous_phase) {
    noise_history = color_noise_cache.phase == previous_phase &&
                    color_noise_cache.scale == params.noise_scale;
    if (noise_history)
      std::copy_n(color_noise_lut, HueNoiseLutView::SIZE, previous_noise_lut);
  }

  float sample_color_noise(const math::Vector &direction) const {
    return sample_hue_noise_lut({color_noise_lut, true}, direction);
  }

  HS_O3_FN __attribute__((always_inline)) float sample_color_noise(
      const ReactionGraph::CubemapLUT::Projection &projection) const {
    const int FACE = projection.face;
    const float U = FACE < 2 ? -projection.u : projection.u;
    const float V = FACE >= 2 && FACE < 4 ? -projection.v : projection.v;
    const HueNoiseSample SAMPLE = prepare_hue_noise_face(FACE, U, V);
    const float NOW = SAMPLE.sample({color_noise_lut, true});
    if (!slab_noise.previous)
      return NOW;
    const float BEFORE = SAMPLE.sample({slab_noise.previous, true});
    return BEFORE + slab_noise.weight * (NOW - BEFORE);
  }

  /**
   * @brief Remaining fraction ahead of the dissolve frontier.
   * @param i Node id.
   * @param phase Dissolve progress in [0, 1].
   * @return Fraction of the node's departure from rest to retain.
   */
  __attribute__((always_inline)) float dissolve_keep(int i, float phase) const {
    const float h =
        math::hash01(static_cast<uint32_t>(i), transition.dissolve_seed);
    if (h < phase)
      return 0.0f;
    if (h < phase + DISSOLVE_FADE_FRACTION)
      return (h - phase) * (1.0f / DISSOLVE_FADE_FRACTION);
    return 1.0f;
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
      const float keep = dissolve_keep(i, phase);
      if (keep < 1.0f) {
        state.A[i] = static_cast<uint16_t>(
            65535.0f - (65535.0f - state.A[i]) * keep + 0.5f);
        state.B[i] =
            static_cast<uint16_t>(static_cast<float>(state.B[i]) * keep + 0.5f);
      }
    }
  }

  /**
   * @brief A node's Q16 B as convert_below() leaves it.
   * @param i Node id.
   * @param phase Dissolve progress in [0, 1].
   * @param b The node's Q16 B before the conversion.
   * @return The converted B.
   */
  __attribute__((always_inline)) uint16_t dissolved_b(int i, float phase,
                                                      uint16_t b) const {
    const float keep = dissolve_keep(i, phase);
    return keep < 1.0f
               ? static_cast<uint16_t>(static_cast<float>(b) * keep + 0.5f)
               : b;
  }

  /**
   * @brief Ends a dissolve by seeding the next reaction at fresh cluster sites.
   * @param staged Splits seeding across a black frame and a complete seed frame.
   * @details Chemistry resumes after all seed clusters have been planted.
   */
  HS_COLD_MEMBER void start_reaction(bool staged = false) {
    transition.dissolve_frames = -1;
    transition.grow_frames = 0;
    transition.stable_frames = 0;
    transition.next_seed = staged ? NUM_SEED_CLUSTERS / 2 : NUM_SEED_CLUSTERS;
    seed_reaction(0, transition.next_seed);
  }

  /**
   * @brief Starts a dissolve at a fresh node ordering.
   */
  void begin_dissolve() {
    transition.dissolve_frames = 0;
    transition.dissolve_seed = static_cast<uint32_t>(hs::random()());
  }

  /**
   * @brief Plants the clusters a staged reseed deferred and latches edits made
   * while it was pending.
   */
  __attribute__((always_inline)) void finish_reseed() {
    seed_reaction(transition.next_seed, NUM_SEED_CLUSTERS);
    transition.next_seed = NUM_SEED_CLUSTERS;
    reaction_edited();
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
   * stalled for STABLE_HOLD_FRAMES.
   * Edits mid-dissolve are absorbed by the
   * in-flight one, which reseeds into whatever the constants read at its end.
   */
  void advance_transition(float mean_db) {
    const TransitionPlan PLAN = plan_transition();
    if (PLAN.dissolve >= 0.0f)
      convert_below(PLAN.dissolve);
    finish_transition(PLAN, mean_db);
  }

  /** @brief The parts of a frame's lifecycle step settled before its
   *         substeps. */
  struct TransitionPlan {
    bool edited;    ///< The reaction constants moved since the last frame.
    float dissolve; ///< Dissolve phase converted this frame, or -1.
    bool reseed; ///< The dissolve ends: the field reseeds, the frame is black.
  };

  /**
   * @brief Latches this frame's edits and advances an in-flight dissolve.
   * @return The frame's lifecycle outcome, but for stabilization.
   */
  TransitionPlan plan_transition() {
    const bool EDITED = reaction_edited();
    if (transition.dissolve_frames < 0)
      return {EDITED, -1.0f, false};
    transition.dissolve_frames++;
    return {EDITED,
            static_cast<float>(transition.dissolve_frames) / DISSOLVE_FRAMES,
            transition.dissolve_frames >= DISSOLVE_FRAMES};
  }

  /**
   * @brief Completes a frame's lifecycle step after its substeps.
   * @param plan The frame's plan_transition() outcome; a dissolve's
   *        conversion is already applied.
   * @param mean_db Mean per-node |dB| between the frame's start and end states.
   */
  void finish_transition(const TransitionPlan &plan, float mean_db) {
    if (plan.dissolve >= 0.0f) {
      if (plan.reseed)
        start_reaction(true);
      return;
    }
    if (plan.edited) {
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

  static constexpr int PHYSICS_NEIGHBOR_REACH = 144;
  static constexpr int PHYSICS_HISTORY_SIZE = 256;
  static_assert(PHYSICS_HISTORY_SIZE > PHYSICS_NEIGHBOR_REACH &&
                (PHYSICS_HISTORY_SIZE & (PHYSICS_HISTORY_SIZE - 1)) == 0);

  HS_COLD_MEMBER static void validate_physics_neighbors(
      const ReactionGraph::NeighborRun *runs = ReactionGraph::neighbor_runs,
      unsigned count = ReactionGraph::NEIGHBOR_RUN_COUNT) {
    for (unsigned r = 0; r < count; ++r)
      for (int delta : runs[r].delta)
        HS_CHECK(delta >= -PHYSICS_NEIGHBOR_REACH,
                 "GS neighbor exceeds delayed-write history");
  }

  /** @brief Advances float A/B in place after their last stencil read. */
  HS_O3_FN void step_physics_inplace(float *a, float *b, float *pending_a,
                                     float *pending_b) {
    HS_PROFILE(grd_physics);
    constexpr int HISTORY_MASK = PHYSICS_HISTORY_SIZE - 1;
    step_physics_nodes(a, b, [&](int i, float next_a, float next_b) {
      pending_a[i & HISTORY_MASK] = next_a;
      pending_b[i & HISTORY_MASK] = next_b;
      // Node i is the last possible reader of i - PHYSICS_NEIGHBOR_REACH.
      if (i >= PHYSICS_NEIGHBOR_REACH) {
        int done = i - PHYSICS_NEIGHBOR_REACH;
        a[done] = pending_a[done & HISTORY_MASK];
        b[done] = pending_b[done & HISTORY_MASK];
      }
    });
    for (int i = RD_N - PHYSICS_NEIGHBOR_REACH; i < RD_N; ++i) {
      a[i] = pending_a[i & HISTORY_MASK];
      b[i] = pending_b[i & HISTORY_MASK];
    }
  }

  template <typename StoreFn>
  HS_HOT_FLASH_MEMBER void
  step_physics_nodes(const float *c_a, const float *c_b, StoreFn &&store) {
    const float feed = params.feed;
    const float KILL_RATE = params.k;
    const float d_a = params.d_a;
    const float d_b = params.d_b;
    const float dt = params.dt * STEP_DT_SCALE;
    int i = 0;
    for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT; ++r) {
      const auto &run = ReactionGraph::neighbor_runs[r];
      auto calculate = [&](int i) __attribute__((always_inline)) {
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
        return std::pair{next_a, next_b};
      };
      for (; i + 3 < run.end; i += 4) {
        const auto first = calculate(i);
        const auto second = calculate(i + 1);
        const auto third = calculate(i + 2);
        const auto fourth = calculate(i + 3);
        store(i, first.first, first.second);
        store(i + 1, second.first, second.second);
        store(i + 2, third.first, third.second);
        store(i + 3, fourth.first, fourth.second);
      }
      for (; i + 1 < run.end; i += 2) {
        const auto first = calculate(i);
        const auto second = calculate(i + 1);
        store(i, first.first, first.second);
        store(i + 1, second.first, second.second);
      }
      if (i < run.end) {
        const auto next = calculate(i);
        store(i, next.first, next.second);
        ++i;
      }
    }
  }

  /// Zero words either side of a node bit set; they cover the widest
  /// neighbour offset, 144 nodes, plus the word a shifted read spills into.
  static constexpr int HOT_GUARD_WORDS = 6;
  /// Words in a node bit set, guards included.
  static constexpr int HOT_WORDS = (RD_N + 31) / 32 + 2 * HOT_GUARD_WORDS;
  static_assert(RD_N % 32 == 0, "node bit sets occupy whole words");

  /** @brief Per-node flags packed 32 to a word, read like a flag array. */
  struct HotBits {
    const uint32_t *words; ///< Word holding nodes 0 to 31.

    /** @brief Whether node @p i is flagged. */
    HS_HOT_INLINE bool operator[](int i) const {
      return (words[i >> 5] >> (i & 31)) & 1u;
    }
  };

  /** @brief The flags of a guarded bit set. */
  static HotBits hot_bits(const uint32_t *bits) {
    return {bits + HOT_GUARD_WORDS};
  }

  /** @brief Stencil terms of a world-space node seen from a fixed-meridian
   *         column. */
  struct MeridianColumn {
    float cos_theta; ///< Cosine of the column's azimuth.
    float sin_theta; ///< Sine of the column's azimuth.

    /** @brief The node's component along the column's azimuth. */
    HS_HOT_INLINE float radial(const math::Vector &p) const {
      return cos_theta * p.x + sin_theta * p.z;
    }
    /** @brief The node's component against increasing azimuth. */
    HS_HOT_INLINE float tangent(const math::Vector &p) const {
      return sin_theta * p.x - cos_theta * p.z;
    }
    /** @brief The node's component along the north axis. */
    HS_HOT_INLINE float up(const math::Vector &p) const { return p.y; }
  };

  /** @brief Stencil terms of a lattice-space node seen from a column at its
   *         sweep time. */
  struct LatticeColumn {
    const Scan::ObjectColumn &frame; ///< The column's lattice-space frame.

    /** @brief The node's component along the column's azimuth. */
    HS_HOT_INLINE float radial(const math::Vector &p) const {
      return math::dot(frame.radial, p);
    }
    /** @brief The node's component against increasing azimuth. */
    HS_HOT_INLINE float tangent(const math::Vector &p) const {
      return -math::dot(frame.tangent, p);
    }
    /** @brief The node's component along the column's north. */
    HS_HOT_INLINE float up(const math::Vector &p) const {
      return math::dot(frame.up, p);
    }
  };

  /** @brief The fixed-meridian stencil terms of column `x`. */
  template <typename Grid>
  static HS_HOT_INLINE MeridianColumn meridian_column(int x) {
    return {math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::cos_theta(x),
            math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::sin_theta[x]};
  }

  template <typename Grid, typename Column, typename OnNode>
  static HS_HOT_INLINE void
  accumulate_column_stencil(const Grid &grid, const Column &column, int center,
                            const math::Vector *nodes, OnNode &&on_node) {
    const auto &run =
        ReactionGraph::neighbor_runs[ReactionGraph::neighbor_run_index[center]];
    float radial_scale[2], cross_scale[2];
    for (int row = 0; row < 2; ++row) {
      radial_scale[row] = grid.sin_phi[row] * grid.cos_dtheta;
      cross_scale[row] =
          2.0f * grid.sin_phi[row] * grid.sin_dtheta * Base::INV_R2;
    }
    for (int slot = 0; slot < RD_K + 1; ++slot) {
      const int node = slot == 0 ? center : center + run.delta[slot - 1];
      const math::Vector &p = nodes[node];
      auto on_weight = on_node(slot, node);
      const float radial = column.radial(p);
      const float tangent = column.tangent(p);
      const float up = column.up(p);
      for (int row = 0; row < 2; ++row) {
        const float dot = radial_scale[row] * radial + grid.cos_phi[row] * up;
        const float base_u = 1.0f - (2.0f - 2.0f * dot) * Base::INV_R2;
        const float cross = cross_scale[row] * tangent;
        const float LEFT = fmaxf(0.0f, base_u - cross);
        const float RIGHT = fmaxf(0.0f, base_u + cross);
        on_weight(2 * row, LEFT * LEFT);
        on_weight(2 * row + 1, RIGHT * RIGHT);
      }
    }
  }

  template <typename Grid, typename OnNode>
  static __attribute__((always_inline)) void
  accumulate_render_stencil(const Grid &grid, int x, int center,
                            const math::Vector *nodes, OnNode &&on_node) {
    accumulate_column_stencil(grid, meridian_column<Grid>(x), center, nodes,
                              on_node);
  }

  HS_O3_FN __attribute__((always_inline)) Pixel shade_pigment(
      uint16_t pigment, float t, float scale, float noise_value) const {
    HS_PROFILE_DEEP(grd_shader_palette);
    const auto PIGMENT = Pigment::unpack(pigment);
    const float FIRST_MASS = Pigment::first_mass(PIGMENT.mass);
    const int PALETTE_COUNT = PIGMENT.mass == Pigment::FULL_MASS ? 1 : 2;
    const float NOISE_POSITION =
        hs::clamp((noise_value + 1.0f) * (0.5f * (COLOR_NOISE_STEPS - 1)), 0.0f,
                  static_cast<float>(COLOR_NOISE_STEPS - 1));
    const auto NOISE_SELECTION = shimmer_cache.selection(NOISE_POSITION);
    const int NOISE_ROW = NOISE_SELECTION.row;
    const bool EXACT = NOISE_SELECTION.exact;
    const float VALUE_POSITION = t * (COLOR_VALUE_STEPS - 1);
    const int VALUE_COLUMN =
        std::min(static_cast<int>(VALUE_POSITION), COLOR_VALUE_STEPS - 2);
    const float NOISE_WEIGHT = NOISE_POSITION - NOISE_ROW;
    const float VALUE_WEIGHT = VALUE_POSITION - VALUE_COLUMN;
    const float SCALE00 = (1.0f - VALUE_WEIGHT) * (1.0f - NOISE_WEIGHT) * scale;
    const float SCALE01 = VALUE_WEIGHT * (1.0f - NOISE_WEIGHT) * scale;
    const float SCALE10 = (1.0f - VALUE_WEIGHT) * NOISE_WEIGHT * scale;
    const float SCALE11 = VALUE_WEIGHT * NOISE_WEIGHT * scale;
    float accum_r = 0, accum_g = 0, accum_b = 0;
    for (int component = 0; component < PALETTE_COUNT; ++component) {
      int palette = component ? PIGMENT.second : PIGMENT.first;
      float mass = component ? 1.0f - FIRST_MASS : FIRST_MASS;
      if (EXACT) {
        Pixel rgb =
            modified_palette_color(palette, t, noise_value * params.hue_shift,
                                   fmaxf(noise_value, 0.0f) * params.shimmer);
        float weight = mass * scale;
        accum_r += rgb.r * weight;
        accum_g += rgb.g * weight;
        accum_b += rgb.b * weight;
        continue;
      }
      const FloatRgb *row =
          shimmer_cache.row(palette, NOISE_ROW) + VALUE_COLUMN;
      float w00 = mass * SCALE00, w01 = mass * SCALE01;
      float w10 = mass * SCALE10, w11 = mass * SCALE11;
      accum_r += row[0].r * w00 + row[1].r * w01 +
                 row[COLOR_VALUE_STEPS].r * w10 +
                 row[COLOR_VALUE_STEPS + 1].r * w11;
      accum_g += row[0].g * w00 + row[1].g * w01 +
                 row[COLOR_VALUE_STEPS].g * w10 +
                 row[COLOR_VALUE_STEPS + 1].g * w11;
      accum_b += row[0].b * w00 + row[1].b * w01 +
                 row[COLOR_VALUE_STEPS].b * w10 +
                 row[COLOR_VALUE_STEPS + 1].b * w11;
    }
    return Pixel(round_linear_channel(accum_r), round_linear_channel(accum_g),
                 round_linear_channel(accum_b));
  }

  HS_HOT_FLASH_MEMBER Pixel shade_render(
      uint16_t pigment, float t, float scale, const math::Vector &center_rv,
      const ReactionGraph::CubemapLUT::Projection *projection) const {
    float noise_value = 0;
    {
      HS_PROFILE_DEEP(grd_shader_noise);
      if (shimmer_cache.noise_enabled())
        noise_value = projection != nullptr
                          ? sample_color_noise(*projection)
                          : sample_color_noise(
                                Base::inverse_orientation.apply(center_rv));
    }
    return shade_pigment(pigment, t, scale, noise_value);
  }

  template <typename Grid, typename Column, typename Flags>
  HS_HOT_FLASH_MEMBER Pixel shade_column_full(
      int seed, const math::Vector &center_rv, const math::Vector *nodes,
      const Grid &grid, const Column &column, const Flags &hot_flags,
      const ReactionGraph::CubemapLUT::Projection *projection) const {
    if (seed < 0)
      return Pixel(0, 0, 0);

    constexpr uint32_t SAMPLES = Grid::SAMPLES;
    float weights[SAMPLES] = {}, weighted_b[SAMPLES] = {};
    uint16_t pigment;
    {
      HS_PROFILE_DEEP(grd_shader_stencil);
      int center =
          Base::template refine_render_center<true>(center_rv, nodes, seed);
      if (culled(hot_flags, center))
        return Pixel(0, 0, 0);
      pigment = state.pigment[center];
      accumulate_column_stencil(
          grid, column, center, nodes,
          [&](int, int node) __attribute__((always_inline)) {
            float b = state.B[node];
            return [&, b](int sample, float weight)
                       __attribute__((always_inline)) {
                         weights[sample] += weight;
                         weighted_b[sample] += b * weight;
                       };
          });
    }
    int covered = 0;
    float t_sum = 0.0f;
    for (int i = 0; i < Grid::SAMPLES; ++i) {
      if (weights[i] <= Base::KERNEL_MIN_TOTAL_WEIGHT)
        continue;
      constexpr float CULL_MASS_SCALE =
          (B_CULL_THRESHOLD / Q16_INV) * (1.0f - 1e-6f);
      constexpr float SATURATED_MASS_SCALE =
          ((B_COLOR_FLOOR + 1.0f / B_COLOR_SCALE) / Q16_INV) * (1.0f + 1e-6f);
      if (weighted_b[i] < CULL_MASS_SCALE * weights[i])
        continue;
      if (weighted_b[i] >= SATURATED_MASS_SCALE * weights[i]) {
        ++covered;
        t_sum += 1.0f;
        continue;
      }
      float b = weighted_b[i] * (Q16_INV / weights[i]);
      if (b < B_CULL_THRESHOLD)
        continue;
      ++covered;
      t_sum += hs::clamp((b - B_COLOR_FLOOR) * B_COLOR_SCALE, 0.0f, 1.0f);
    }
    if (!covered)
      return Pixel(0, 0, 0);
    return shade_render(pigment, t_sum / covered, covered * (1.0f / SAMPLES),
                        center_rv, projection);
  }

  template <typename Grid>
  Pixel shade_pixel_full(
      int seed, const math::Vector &center_rv, const math::Vector *world_nodes,
      const Grid &grid, int x, const uint8_t *hot_flags = nullptr,
      const ReactionGraph::CubemapLUT::Projection *projection = nullptr) const {
    return shade_column_full(seed, center_rv, world_nodes, grid,
                             meridian_column<Grid>(x), hot_flags, projection);
  }

  template <typename Grid>
  static __attribute__((always_inline)) float render_support_limit() {
    return Base::KERNEL_R * 0.99f -
           0.25f * (math::RADIANS_PER_COLUMN<Grid::WIDTH> +
                    math::RADIANS_PER_ROW<Grid::HEIGHT>);
  }

  /**
   * @brief Shades one pixel's four sub-samples through an inlinable typed path.
   * @tparam Grid Scan::Shader::SsaaGrid type supplying the sub-pixel offsets.
   * @tparam Column MeridianColumn or LatticeColumn.
   * @tparam Flags `const uint8_t *` or HotBits.
   * @param seed Cubemap-LUT seed node id, or -1 for a culled pixel.
   * @param center_rv Direction at the pixel center, in the frame of @p nodes.
   * @param nodes Lattice node positions in the frame of @p column.
   * @param grid Row's SSAA sub-pixel grid.
   * @param column The pixel column's stencil terms.
   * @param hot_flags Per-node flags, a byte array or HotBits; a pixel whose
   *        center node is unflagged shades black. A null array shades every
   *        pixel.
   * @param projection Pixel's precomputed cubemap projection for color noise,
   *        or null to project `center_rv`.
   * @return The finished alpha-premultiplied pixel.
   * @details Accepts seeds inside a proven nearest-node radius immediately;
   * boundary pixels check the neighbors. Sub-samples share the center
   * stencil, which can reach past one node spacing at low vertical resolution.
   */
  template <typename Grid, typename Column, typename Flags>
  HS_HOT_INLINE Pixel shade_column_pixel(
      int seed, const math::Vector &center_rv, const math::Vector *nodes,
      const Grid &grid, const Column &column, const Flags &hot_flags,
      const ReactionGraph::CubemapLUT::Projection *projection) const {
    if (seed < 0)
      return Pixel(0, 0, 0);

    static_assert(B_CULL_THRESHOLD == B_COLOR_FLOOR);
    constexpr uint32_t SAMPLES = Grid::SAMPLES;
    float w00 = 0, w11 = 0;
    float b00 = 0, b01 = 0, b10 = 0, b11 = 0;
    uint16_t pigment;
    {
      HS_PROFILE_DEEP(grd_shader_stencil);
      int center =
          Base::template refine_render_center<true>(center_rv, nodes, seed);
      if (culled(hot_flags, center))
        return Pixel(0, 0, 0);
      const float LIMIT = render_support_limit<Grid>();
      if (LIMIT <= 0.0f ||
          Base::dist2(center_rv, nodes[center]) > LIMIT * LIMIT)
        return shade_column_full(seed, center_rv, nodes, grid, column,
                                 hot_flags, projection);
      pigment = state.pigment[center];
      const auto &run = ReactionGraph::neighbor_runs
          [ReactionGraph::neighbor_run_index[center]];
      const float RADIAL0 = grid.sin_phi[0] * grid.cos_dtheta;
      const float RADIAL1 = grid.sin_phi[1] * grid.cos_dtheta;
      const float CROSS0 =
          2.0f * grid.sin_phi[0] * grid.sin_dtheta * Base::INV_R2;
      const float CROSS1 =
          2.0f * grid.sin_phi[1] * grid.sin_dtheta * Base::INV_R2;
      // An object-frame column folds its frame into one direction per term,
      // so each node costs a dot product per term.
      constexpr bool OBJECT_FRAME = std::is_same_v<Column, LatticeColumn>;
      [[maybe_unused]] math::Vector base_dir0, base_dir1, side0, side1;
      if constexpr (OBJECT_FRAME) {
        const Scan::ObjectColumn &frame = column.frame;
        const float SCALE = 2.0f * Base::INV_R2;
        base_dir0 =
            (frame.radial * RADIAL0 + frame.up * grid.cos_phi[0]) * SCALE;
        base_dir1 =
            (frame.radial * RADIAL1 + frame.up * grid.cos_phi[1]) * SCALE;
        side0 = frame.tangent * -CROSS0;
        side1 = frame.tangent * -CROSS1;
      }
      constexpr float BASE_OFFSET = 1.0f - 2.0f * Base::INV_R2;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC unroll 2
#endif
      for (int slot = 0; slot < RD_K + 1; ++slot) {
        const int NODE = slot == 0 ? center : center + run.delta[slot - 1];
        const math::Vector &p = nodes[NODE];
        const float B = state.B[NODE] - B_CULL_THRESHOLD / Q16_INV;
        float BASE0, BASE1, D0, D1;
        if constexpr (OBJECT_FRAME) {
          BASE0 = BASE_OFFSET + math::dot(base_dir0, p);
          BASE1 = BASE_OFFSET + math::dot(base_dir1, p);
          D0 = math::dot(side0, p);
          D1 = math::dot(side1, p);
        } else {
          const float RADIAL = column.radial(p);
          const float TANGENT = column.tangent(p);
          const float UP = column.up(p);
          const float DOT0 = RADIAL0 * RADIAL + grid.cos_phi[0] * UP;
          const float DOT1 = RADIAL1 * RADIAL + grid.cos_phi[1] * UP;
          BASE0 = 1.0f - (2.0f - 2.0f * DOT0) * Base::INV_R2;
          BASE1 = 1.0f - (2.0f - 2.0f * DOT1) * Base::INV_R2;
          D0 = CROSS0 * TANGENT;
          D1 = CROSS1 * TANGENT;
        }
        float left = fmaxf(0.0f, BASE0 - D0), right = fmaxf(0.0f, BASE0 + D0);
        left *= left;
        right *= right;
        w00 += left;
        b00 += B * left;
        b01 += B * right;
        left = fmaxf(0.0f, BASE1 - D1);
        right = fmaxf(0.0f, BASE1 + D1);
        left *= left;
        right *= right;
        w11 += right;
        b10 += B * left;
        b11 += B * right;
      }
    }
    // Within the guard, Q16 rounding can flip a corner's sign.
    constexpr float COVERAGE_ROUNDING_GUARD = 0.125f;
    if (fabsf(b00) < COVERAGE_ROUNDING_GUARD ||
        fabsf(b01) < COVERAGE_ROUNDING_GUARD ||
        fabsf(b10) < COVERAGE_ROUNDING_GUARD ||
        fabsf(b11) < COVERAGE_ROUNDING_GUARD)
      return shade_column_full(seed, center_rv, nodes, grid, column, hot_flags,
                               projection);
    int covered = (b00 >= 0.0f) + (b01 >= 0.0f) + (b10 >= 0.0f) + (b11 >= 0.0f);
    if (!covered)
      return Pixel(0, 0, 0);
    const float INVERSE0 = Q16_INV * B_COLOR_SCALE / w00;
    const float INVERSE1 = Q16_INV * B_COLOR_SCALE / w11;
    const float INVERSE_CROSS = 0.5f * (INVERSE0 + INVERSE1);
    float t_sum = hs::clamp(b00 * INVERSE0, 0.0f, 1.0f) +
                  hs::clamp(b11 * INVERSE1, 0.0f, 1.0f) +
                  hs::clamp(b01 * INVERSE_CROSS, 0.0f, 1.0f) +
                  hs::clamp(b10 * INVERSE_CROSS, 0.0f, 1.0f);
    return shade_render(pigment, t_sum / covered, covered * (1.0f / SAMPLES),
                        center_rv, projection);
  }

  /**
   * @brief Shades one pixel against the oriented lattice at the fixed meridian.
   * @tparam Grid Scan::Shader::SsaaGrid type supplying the sub-pixel offsets.
   * @param seed Cubemap-LUT seed node id, or -1 for a culled pixel.
   * @param center_rv World-space direction at the pixel center.
   * @param world_nodes Oriented lattice node positions.
   * @param grid Row's SSAA sub-pixel grid.
   * @param x Pixel column.
   * @param hot_flags Per-node flags; a pixel whose center node is unflagged
   *        shades black. Null shades every pixel.
   * @param projection Pixel's precomputed cubemap projection for color noise,
   *        or null to project `center_rv`.
   * @return The finished alpha-premultiplied pixel.
   */
  template <typename Grid>
  HS_O3_FN Pixel shade_pixel(
      int seed, const math::Vector &center_rv, const math::Vector *world_nodes,
      const Grid &grid, int x, const uint8_t *hot_flags = nullptr,
      const ReactionGraph::CubemapLUT::Projection *projection = nullptr) const {
    return shade_column_pixel(seed, center_rv, world_nodes, grid,
                              meridian_column<Grid>(x), hot_flags, projection);
  }

  /**
   * @brief Shades one pixel against the unrotated lattice from its column's
   *        frame at its sweep time.
   * @tparam Grid Scan::Shader::SsaaGrid type supplying the sub-pixel offsets.
   * @param seed Cubemap-LUT seed node id, or -1 for a culled pixel.
   * @param center Lattice-space direction at the pixel center.
   * @param grid Row's SSAA sub-pixel grid.
   * @param frame The pixel column's lattice-space frame.
   * @param hot_flags One-ring flags; a pixel whose center node is unflagged
   *        shades black.
   * @param projection Pixel's cubemap projection for color noise.
   * @return The finished alpha-premultiplied pixel.
   */
  template <typename Grid>
  HS_HOT_FLASH_MEMBER Pixel shade_lattice_pixel(
      int seed, const math::Vector &center, const Grid &grid,
      const Scan::ObjectColumn &frame, const HotBits &hot_flags,
      const ReactionGraph::CubemapLUT::Projection *projection) const {
    return shade_column_pixel(seed, center, Base::nodes, grid,
                              LatticeColumn{frame}, hot_flags, projection);
  }

  /**
   * @brief Builds per-node two-ring "renderable" flags for the B field.
   * @param b Per-node B concentrations, Q16.
   * @param hot1 Output: per-node flag, set when any of {node, neighbors}
   *        reaches the threshold; read by the per-center shader cull.
   * @param hot2 Output: per-node flag, set when any node within two hops
   *        reaches the threshold.
   * @param count Node count.
   * @param threshold Q16 render floor.
   * @details The stencil reaches two hops from the seed via the refined
   * center, so culling on !hot2[seed] is exact.
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

  /** @brief Whether a flag array culls a node; a null byte array culls none.
   */
  static HS_HOT_INLINE bool culled(const uint8_t *flags, int node) {
    return flags && !flags[node];
  }
  /** @brief Whether a bit set culls a node. */
  static HS_HOT_INLINE bool culled(const HotBits &flags, int node) {
    return !flags[node];
  }

  /**
   * @brief Sets the bit of each node whose B reaches a threshold.
   * @param b Per-node B concentrations, Q16.
   * @param bits HOT_WORDS words receiving the set, guards zeroed.
   * @param threshold Q16 render floor.
   */
  HS_O3_FN HS_NOINLINE_NOCLONE static void
  mark_hot_bits(const uint16_t *b, uint32_t *bits, uint16_t threshold) {
    std::fill_n(bits, HOT_WORDS, 0u);
    uint32_t *words = bits + HOT_GUARD_WORDS;
    for (int w = 0; w < RD_N / 32; ++w) {
      uint32_t word = 0;
      for (int j = 0; j < 32; ++j)
        word |= static_cast<uint32_t>(b[32 * w + j] >= threshold) << j;
      words[w] = word;
    }
  }

  /**
   * @brief The 32 bits of a guarded bit set starting at a node index.
   * @param words Word of node 0.
   * @param bit First node, at least -32 * (HOT_GUARD_WORDS - 1).
   * @return Bits [bit, bit + 32).
   */
  static HS_HOT_INLINE uint32_t bits_from(const uint32_t *words, int bit) {
    const int word = bit >> 5;
    const int shift = bit & 31;
    const uint32_t low = words[word] >> shift;
    return shift ? low | (words[word + 1] << (32 - shift)) : low;
  }

  /**
   * @brief Sets the bit of each node that is set, or has a set neighbour, in
   *        another bit set.
   * @param in Input set, HOT_WORDS words, guards zeroed.
   * @param out HOT_WORDS words receiving the output set, guards zeroed.
   * @details Each neighbour run shares its offsets, so a run's nodes take
   * their neighbours 32 at a time as shifted words of @p in.
   */
  HS_O3_FN HS_NOINLINE_NOCLONE static void dilate_hot_bits(const uint32_t *in,
                                                           uint32_t *out) {
    std::fill_n(out, HOT_WORDS, 0u);
    const uint32_t *source = in + HOT_GUARD_WORDS;
    uint32_t *target = out + HOT_GUARD_WORDS;
    int start = 0;
    for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT; ++r) {
      const auto &run = ReactionGraph::neighbor_runs[r];
      const int end = run.end;
      for (int w = start >> 5; w <= (end - 1) >> 5; ++w) {
        const int base = w << 5;
        uint32_t hot = source[w];
        for (int k = 0; k < RD_K; ++k)
          hot |= bits_from(source, base + run.delta[k]);
        const int lo = std::max(start - base, 0);
        const int hi = std::min(end - base, 32);
        const uint32_t below_hi = hi == 32 ? ~0u : (1u << hi) - 1u;
        target[w] |= hot & below_hi & ~((1u << lo) - 1u);
      }
      start = end;
    }
  }

  HS_O3_FN void draw_lattice(Canvas &canvas, const math::Vector *world_nodes,
                             const uint8_t *hot1, const uint8_t *hot2) const {
    Scan::Shader::walk_grid<W, H>(
        canvas, [&](const math::Vector &center, const auto &grid, int x) {
          const math::Vector object_direction =
              Base::inverse_orientation.apply(center);
          const auto PROJECTION =
              ReactionGraph::CubemapLUT::project(object_direction);
          const int seed = Base::cube_lut.lookup(PROJECTION);
          return hot2[seed] ? shade_pixel(seed, center, world_nodes, grid, x,
                                          hot1, &PROJECTION)
                            : Pixel(0, 0, 0);
        });
  }

  /**
   * @brief Shades the B field, each column against the unrotated lattice from
   *        its frame at its sweep time.
   * @param canvas Destination canvas; an open TimeWindow limits the columns.
   * @param basis The frame's per-column lattice frames.
   * @param hot1 One-ring renderable flags.
   * @param hot2 Two-ring renderable flags.
   */
  HS_HOT_FLASH_MEMBER void draw_lattice_columns(Canvas &canvas,
                                                const LatticeBasis &basis,
                                                const HotBits &hot1,
                                                const HotBits &hot2) const {
    Scan::Shader::walk_object_grid<W, H>(
        canvas, basis,
        [&](const math::Vector &center, const Scan::ObjectColumn &frame,
            const auto &grid, int) HS_HOT_INLINE {
          const auto PROJECTION = ReactionGraph::CubemapLUT::project(center);
          const int seed = Base::cube_lut.lookup(PROJECTION);
          return hot2[seed] ? shade_lattice_pixel(seed, center, grid, frame,
                                                  hot1, &PROJECTION)
                            : Pixel(0, 0, 0);
        });
  }

  /**
   * @brief Quantizes the float B field into the Q16 field the shader reads,
   *        marking the nodes that reach the render floor.
   * @param b Float B per node.
   * @param dissolve The frame's dissolve phase, or negative for none.
   * @param marked HOT_WORDS words whose node bits receive the marks, as
   *        mark_hot_bits() sets them; the guard words are left alone.
   * @details Quantizes as to_q16(), with the clamp branch-free.
   */
  HS_O3_FN HS_NOINLINE_NOCLONE void publish_b(const float *b, float dissolve,
                                              uint32_t *marked) {
    HS_PROFILE(grd_publish);
    const uint16_t THRESHOLD = to_q16(B_CULL_THRESHOLD);
    uint32_t *words = marked + HOT_GUARD_WORDS;
    const auto publish = [&](auto &&adjust) HS_HOT_INLINE {
      for (int w = 0; w < RD_N / 32; ++w) {
        uint32_t word = 0;
        for (int j = 0; j < 32; ++j) {
          const int i = 32 * w + j;
          const uint16_t q =
              adjust(i, static_cast<uint16_t>(
                            fminf(fmaxf(b[i], 0.0f), 1.0f) * 65535.0f + 0.5f));
          state.B[i] = q;
          word |= static_cast<uint32_t>(q >= THRESHOLD) << j;
        }
        words[w] = word;
      }
    };
    if (dissolve < 0.0f)
      publish([](int, uint16_t q) HS_HOT_INLINE { return q; });
    else
      publish([&](int i, uint16_t q)
                  HS_HOT_INLINE { return dissolved_b(i, dissolve, q); });
  }

  /**
   * @brief Advances the sim and shades it in sweep slabs, one per substep.
   * @param canvas Destination canvas to draw the sphere into.
   * @details The frame's pigment and lifecycle outcome are settled before the
   * first substep. Each slab shades the field after its substep through the
   * dissolve mask at the slab's end, with its own cull flags, and colour noise
   * blended between the frame's two noise tables at the slab's end; slabs
   * ending in the first half of the frame shade the frame-start pigment and
   * the rest the frame-end pigment. A frame that ends a dissolve stays black;
   * a reseed frame shades the seeded field whole. The Q16 state ends the
   * frame as render() leaves it with temporal sampling disabled.
   */
  HS_HOT_FLASH_MEMBER void render_slabs(Canvas &canvas) {
    const LatticeBasis &basis = *prepare_lattice_basis();
    uint32_t *marked = scratch_arena_a.allocate_n<uint32_t>(HOT_WORDS);
    uint32_t *hot1 = scratch_arena_a.allocate_n<uint32_t>(HOT_WORDS);
    uint32_t *hot2 = scratch_arena_a.allocate_n<uint32_t>(HOT_WORDS);
    std::fill_n(marked, HOT_WORDS, 0u);
    auto shade = [&]() HS_HOT_INLINE {
      {
        HS_PROFILE(grd_cull_flags);
        dilate_hot_bits(marked, hot1);
        dilate_hot_bits(hot1, hot2);
      }
      HS_PROFILE(grd_shader_draw);
      draw_lattice_columns(canvas, basis, hot_bits(hot1), hot_bits(hot2));
    };
    if (transition.reseed_pending()) {
      finish_reseed();
      refresh_color_palettes();
      mark_hot_bits(state.B, marked, to_q16(B_CULL_THRESHOLD));
      shade();
      return;
    }
    const TransitionPlan PLAN = plan_transition();
    if (!PLAN.reseed)
      refresh_color_palettes();
    float *cur_a = scratch_arena_a.allocate_n<float>(RD_N);
    float *cur_b = scratch_arena_a.allocate_n<float>(RD_N);
    for (int i = 0; i < RD_N; i++) {
      cur_a[i] = from_q16(state.A[i]);
      cur_b[i] = from_q16(state.B[i]);
    }
    uint16_t *next_pigment = scratch_arena_a.allocate_n<uint16_t>(RD_N);
    step_pigment<STEPS_PER_FRAME>(cur_a, cur_b, next_pigment);
    float *pending_a = scratch_arena_a.allocate_n<float>(PHYSICS_HISTORY_SIZE);
    float *pending_b = scratch_arena_a.allocate_n<float>(PHYSICS_HISTORY_SIZE);
    uint16_t *start_b = scratch_arena_a.allocate_n<uint16_t>(RD_N);
    std::copy_n(state.B, RD_N, start_b);
    int slab = 0;
    Spacetime::simulate_in_slabs<STEPS_PER_FRAME, STEPS_PER_FRAME>(
        canvas,
        [&](int substeps) HS_HOT_INLINE {
          for (int step = 0; step < substeps; ++step)
            step_physics_inplace(cur_a, cur_b, pending_a, pending_b);
        },
        [&]() HS_HOT_INLINE {
          const int ENDS = ++slab;
          if (ENDS == STEPS_PER_FRAME / 2 + 1)
            std::copy_n(next_pigment, RD_N, state.pigment);
          if (PLAN.reseed)
            return;
          const float dissolve =
              PLAN.dissolve < 0.0f
                  ? PLAN.dissolve
                  : (static_cast<float>(transition.dissolve_frames - 1) +
                     static_cast<float>(ENDS) / STEPS_PER_FRAME) /
                        DISSOLVE_FRAMES;
          publish_b(cur_b, dissolve, marked);
          slab_noise = {noise_history ? previous_noise_lut : nullptr,
                        static_cast<float>(ENDS) / STEPS_PER_FRAME};
          shade();
        });
    slab_noise = {};
    uint32_t db_sum_q16 = 0;
    for (int i = 0; i < RD_N; i++) {
      state.A[i] = to_q16(cur_a[i]);
      uint16_t next_b = to_q16(cur_b[i]);
      int db = static_cast<int>(next_b) - start_b[i];
      db_sum_q16 += static_cast<uint32_t>(db < 0 ? -db : db);
      state.B[i] = next_b;
    }
    if (PLAN.dissolve >= 0.0f)
      convert_below(PLAN.dissolve);
    finish_transition(PLAN, static_cast<float>(db_sum_q16) * (Q16_INV / RD_N));
  }

  /**
   * @brief Advances the sim STEPS_PER_FRAME substeps and rasterizes the B field.
   * @param canvas Destination canvas to draw the sphere into.
   * @details With temporal sampling enabled, render_slabs() advances and
   * shades the frame. Reseeding uses one black frame followed by a complete
   * seed frame, both without partial-state chemistry or rendering.
   */
  void render(Canvas &canvas) {
    HS_PROFILE(grd_render);
    {
      HS_PROFILE(grd_color_noise);
      advance_color_noise();
    }
    ScratchScope frame_guard(scratch_arena_a);
    if (!Spacetime::flat_surface()) {
      render_slabs(canvas);
      return;
    }
    float mean_db = 0.0f;
    if (transition.reseed_pending()) {
      finish_reseed();
    } else {
      {
        // Q16 quantization occurs once per frame.
        HS_PROFILE(grd_simulate);
        ScratchScope physics_guard(scratch_arena_a);
        float *cur_a = scratch_arena_a.allocate_n<float>(RD_N);
        float *cur_b = scratch_arena_a.allocate_n<float>(RD_N);
        uint16_t *next_pigment = scratch_arena_a.allocate_n<uint16_t>(RD_N);
        float *pending_a =
            scratch_arena_a.allocate_n<float>(PHYSICS_HISTORY_SIZE);
        float *pending_b =
            scratch_arena_a.allocate_n<float>(PHYSICS_HISTORY_SIZE);

        for (int i = 0; i < RD_N; i++) {
          cur_a[i] = from_q16(state.A[i]);
          cur_b[i] = from_q16(state.B[i]);
        }
        step_pigment<STEPS_PER_FRAME>(cur_a, cur_b, next_pigment);
        std::copy_n(next_pigment, RD_N, state.pigment);
        for (int step = 0; step < STEPS_PER_FRAME; ++step) {
          step_physics_inplace(cur_a, cur_b, pending_a, pending_b);
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
      if (transition.reseed_pending())
        return;
    }
    refresh_color_palettes();

    // Physics scratch is popped; the oriented lattice and cull flags reuse it.
    HS_PROFILE(grd_rasterize);
    auto lattice = [this] {
      HS_PROFILE(grd_orient);
      return orient_lattice();
    }();
    math::Vector *world_nodes = lattice.get();
    uint8_t *hot1 = scratch_arena_a.allocate_n<uint8_t>(RD_N);
    uint8_t *hot2 = scratch_arena_a.allocate_n<uint8_t>(RD_N);
    {
      HS_PROFILE(grd_cull_flags);
      fill_hot_flags(state.B, hot1, hot2, RD_N, to_q16(B_CULL_THRESHOLD));
    }

    {
      HS_PROFILE(grd_shader_draw);
      draw_lattice(canvas, world_nodes, hot1, hot2);
    }
  }

  /**
   * @brief Persistent Q16 A/B concentrations and packed two-palette pigment.
   */
  struct {
    uint16_t *A = nullptr,
             *B = nullptr;       /**< Per-node A/B concentrations, Q16. */
    uint16_t *pigment = nullptr; /**< Per-node packed `Pigment`. */
  } state;

  /**
   * @brief Auto-transition state: current reaction's lifetime and dissolve
   *        progress.
   */
  struct Transition {
    int next_seed = NUM_SEED_CLUSTERS; /**< First cluster awaiting placement. */
    int grow_frames = 0;        /**< Frames since this reaction was seeded. */
    int stable_frames = 0;      /**< Consecutive sub-floor frames. */
    int dissolve_frames = -1;   /**< Dissolve progress; -1 when inactive. */
    uint32_t dissolve_seed = 0; /**< Per-transition node-order hash seed. */
    /** Reaction constants as of the last frame, latched by
     *  reaction_edited(). */
    float last_feed = 0.0f, last_k = 0.0f, last_d_a = 0.0f, last_d_b = 0.0f;

    /** @brief True while a staged reseed has clusters left to plant. */
    bool reseed_pending() const { return next_seed < NUM_SEED_CLUSTERS; }
  } transition;

  /** @brief Per-seed linear RGB ramps sampled by B concentration. */
  Pixel *palettes = nullptr;
  ShimmerCache shimmer_cache;
  int8_t *color_noise_lut = nullptr;
  /** @brief The previous frame's noise table, under temporal sampling. */
  int8_t *previous_noise_lut = nullptr;
  /** @brief `previous_noise_lut` holds the previous frame's bake. */
  bool noise_history = false;
  /** @brief Noise blend of the slab being shaded. */
  struct SlabNoise {
    /// Table the blend starts from, or null for the current table alone.
    const int8_t *previous = nullptr;
    float weight = 1.0f; ///< Weight of the current table.
  } slab_noise;
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
