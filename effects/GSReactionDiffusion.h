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

#include <algorithm>
#include <bit>
#include <cmath>
#include <utility>
#include "core/color/effect_palette_recipes.h"
#include "core/color/noise_shimmer_palette.h"
#include "core/engine/engine.h"
#include "effects/ReactionDiffusionBase.h"

// Unit-test accessor for private state the smoke harness cannot pin.
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
  using Base::for_each_neighbor;
  using Base::from_q16;
  using Base::init_lattice;
  using Base::orient_lattice;
  using Base::Q16_INV;
  using Base::RD_K;
  using Base::RD_N;
  using Base::refine_and_accumulate;
  using Base::refine_render_center;
  using Base::register_param;
  using Base::to_q16;

public:
  static constexpr const char *EFFECT_ID = "GSReactionDiffusion";

  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  GSReactionDiffusion() = default;

  /**
   * @brief One-time setup: arenas, GUI params, A/B state, cubemap LUT, lattice.
   * @details Carves the persistent arena, registers the GUI params, allocates
   * the A/B/pigment state, seed palettes and colour-noise LUT, seeds the first
   * reaction, binds the flash lattice and builds the cubemap LUT once.
   */
  void init() override {
    constexpr size_t PALETTE_BYTES =
        NUM_SEED_CLUSTERS *
            (PALETTE_SIZE * sizeof(Pixel) +
             COLOR_VALUE_STEPS * COLOR_NOISE_STEPS * sizeof(FloatColor)) +
        2 * alignof(Pixel) + HueNoiseLutView::SIZE * sizeof(int8_t);
    constexpr size_t PERSISTENT_BYTES = 192 * 1024 + 32;
    constexpr size_t PHYSICS_SCRATCH_BYTES =
        2u * RD_N * sizeof(float) + RD_N * sizeof(uint16_t) +
        2u * PHYSICS_HISTORY_SIZE * sizeof(float);
    constexpr size_t RASTER_SCRATCH_BYTES =
        RD_N * sizeof(math::Vector) + 2u * RD_N * sizeof(uint8_t);
    constexpr size_t SCRATCH_BYTES =
        std::max(PHYSICS_SCRATCH_BYTES, RASTER_SCRATCH_BYTES);
    Base::template configure_rd_arenas<uint16_t, 3, PERSISTENT_BYTES,
                                       PALETTE_BYTES, SCRATCH_BYTES, 3, true>();

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

    state.A = persistent_arena.allocate_n<uint16_t>(RD_N);
    state.B = persistent_arena.allocate_n<uint16_t>(RD_N);
    state.pigment = persistent_arena.allocate_n<uint16_t>(RD_N);
    palettes =
        persistent_arena.allocate_n<Pixel>(NUM_SEED_CLUSTERS * PALETTE_SIZE);
    modified_palettes = persistent_arena.allocate_n<FloatColor>(
        NUM_SEED_CLUSTERS * COLOR_VALUE_STEPS * COLOR_NOISE_STEPS);
    color_noise_lut =
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
  // Test seam for private state the smoke harness cannot pin.
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
  struct FloatColor {
    float r, g, b;
    FloatColor &operator=(const Pixel &p) {
      r = p.r;
      g = p.g;
      b = p.b;
      return *this;
    }
    Pixel pixel() const { return Pixel(r, g, b); }
  };
  static constexpr float CACHED_HUE_LIMIT = 0.25f;
  static constexpr float CACHED_SHIMMER_LIMIT = 0.4f;
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

  HS_COLD_MEMBER void seed_reaction(int first = 0,
                                    int end = NUM_SEED_CLUSTERS) {
    HS_PROFILE(grd_seed_reaction);
    if (first == 0) {
      for (int i = 0; i < RD_N; ++i) {
        state.A[i] = 65535;
        state.B[i] = 0;
        state.pigment[i] = 63u << 10;
      }
      color_palette_valid = false;
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
        uint16_t uniform = state.pigment[i] & 0xfc1fu;
        bool same = (uniform >> 10) == 63;
        for (int k = 0; k < RD_K && same; ++k)
          same = (state.pigment[i + run.delta[k]] & 0xfc1fu) == uniform;
        if (same) {
          float mass = retained;
          for (int k = 0; k < RD_K; ++k)
            mass += b[i + run.delta[k]] * diffusion;
          int first = mass > 0.0f ? uniform & 31 : 0;
          int second = first == 0 ? 1 : 0;
          next[i] = static_cast<uint16_t>(first | (second << 5) | (63u << 10));
          continue;
        }
        uint32_t touched = 0;
        auto add = [&](int id, float mass) __attribute__((always_inline)) {
          weights[id] += mass;
          touched |= 1u << id;
        };
        auto gather = [&](int node, float mass) __attribute__((always_inline)) {
          uint16_t pigment = state.pigment[node];
          float first = static_cast<float>(pigment >> 10) * (1.0f / 63.0f);
          add(pigment & 31, mass * first);
          add((pigment >> 5) & 31, mass * (1.0f - first));
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
                      ? static_cast<int>(63.0f * first_mass / mass + 0.5f)
                      : 63;
        next[i] = static_cast<uint16_t>(first | (second << 5) | (mix << 10));
      }
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

  HS_HOT_FLASH_MEMBER void refresh_color_palettes(bool complete = false) {
    color_noise_enabled = params.hue_shift != 0.0f || params.shimmer != 0.0f;
    color_palette_exact = !color_noise_enabled ||
                          fabsf(params.hue_shift) > CACHED_HUE_LIMIT ||
                          params.shimmer > CACHED_SHIMMER_LIMIT;
    const bool RESTART = !color_palette_valid ||
                         color_palette_hue != params.hue_shift ||
                         color_palette_shimmer != params.shimmer;
    if (RESTART) {
      color_palette_hue = params.hue_shift;
      color_palette_shimmer = params.shimmer;
      color_palette_valid = true;
      color_palette_rows = 0;
      color_palette_next_row = 0;
    }
    if (color_palette_next_row == COLOR_NOISE_STEPS ||
        fabsf(params.hue_shift) > CACHED_HUE_LIMIT ||
        params.shimmer > CACHED_SHIMMER_LIMIT ||
        (params.hue_shift == 0.0f && params.shimmer == 0.0f))
      return;
    HS_PROFILE(grd_color_palette);
    int count = complete ? COLOR_NOISE_STEPS : RESTART ? 5 : 2;
    for (; count > 0 && color_palette_next_row < COLOR_NOISE_STEPS; --count) {
      int index = color_palette_next_row++;
      int n =
          COLOR_NOISE_STEPS / 2 + ((index & 1) ? -(index + 1) / 2 : index / 2);
      float noise = -1.0f + 2.0f * n / (COLOR_NOISE_STEPS - 1);
      for (int seed = 0; seed < NUM_SEED_CLUSTERS; ++seed)
        for (int t = 0; t < COLOR_VALUE_STEPS; ++t)
          modified_palettes[(seed * COLOR_NOISE_STEPS + n) * COLOR_VALUE_STEPS +
                            t] =
              modified_palette_color(
                  seed, static_cast<float>(t) / (COLOR_VALUE_STEPS - 1),
                  noise * params.hue_shift,
                  fmaxf(noise, 0.0f) * params.shimmer);
      color_palette_rows |= 1u << n;
    }
  }

  struct ColorNoiseSample {
    int row;
    uint16_t weight;
    float value;
    bool exact;
  };

  struct ColorValueSample {
    int column;
    uint16_t weight;
    float value;
  };

  __attribute__((always_inline)) ColorNoiseSample
  color_noise_sample(float noise) const {
    float position =
        hs::clamp((noise + 1.0f) * (0.5f * (COLOR_NOISE_STEPS - 1)), 0.0f,
                  static_cast<float>(COLOR_NOISE_STEPS - 1));
    int row = std::min(static_cast<int>(position), COLOR_NOISE_STEPS - 2);
    return {row, static_cast<uint16_t>((position - row) * 65535.0f), noise,
            !color_palette_valid ||
                (color_palette_rows & (3u << row)) != (3u << row) ||
                color_palette_exact};
  }

  __attribute__((always_inline)) static ColorValueSample
  color_value_sample(float t) {
    float position = t * (COLOR_VALUE_STEPS - 1);
    int column = std::min(static_cast<int>(position), COLOR_VALUE_STEPS - 2);
    return {column, lut_index_weight(position, column), t};
  }

  __attribute__((always_inline)) Pixel
  cached_palette_color(int seed, const ColorValueSample &value,
                       const ColorNoiseSample &noise) const {
    if (noise.exact)
      return modified_palette_color(seed, value.value,
                                    noise.value * params.hue_shift,
                                    fmaxf(noise.value, 0.0f) * params.shimmer);
    const FloatColor *row =
        modified_palettes +
        (seed * COLOR_NOISE_STEPS + noise.row) * COLOR_VALUE_STEPS +
        value.column;
    Pixel first = row[0].pixel().lerp16(row[1].pixel(), value.weight);
    Pixel second = row[COLOR_VALUE_STEPS].pixel().lerp16(
        row[COLOR_VALUE_STEPS + 1].pixel(), value.weight);
    return first.lerp16(second, noise.weight);
  }

  Pixel cached_palette_color(int seed, float t, float noise) const {
    return cached_palette_color(seed, color_value_sample(t),
                                color_noise_sample(noise));
  }

  HS_FLASH_INLINE void refresh_color_noise() {
    color_noise_cache.refresh<true>(std::span<int8_t, HueNoiseLutView::SIZE>(
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

  HS_O3_FN __attribute__((always_inline)) float sample_color_noise(
      const ReactionGraph::CubemapLUT::Projection &projection) const {
    const int FACE = projection.face;
    const float U = FACE < 2 ? -projection.u : projection.u;
    const float V = FACE >= 2 && FACE < 4 ? -projection.v : projection.v;
    return sample_hue_noise_face({color_noise_lut, true}, FACE, U, V);
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
        start_reaction(true);
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
   * @brief Jacobi reference substep for the in-place physics oracle tests.
   * @details The render path uses step_physics_inplace.
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

  /**
   * @brief Kernel-weighted sample of the B concentration at a point.
   * @param p Query point on the sphere.
   * @param seed Seed node id from the cubemap LUT; the fused stencil walk
   * selects the nearest node among the seed and its direct neighbors.
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

  template <typename Grid, typename OnNode>
  static __attribute__((always_inline)) void
  accumulate_render_stencil(const Grid &grid, int x, int center,
                            const math::Vector *nodes, OnNode &&on_node) {
    const float st = math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::sin_theta[x];
    const float ct = math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::cos_theta(x);
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
      const float radial = ct * p.x + st * p.z;
      const float tangent = st * p.x - ct * p.z;
      for (int row = 0; row < 2; ++row) {
        const float dot = radial_scale[row] * radial + grid.cos_phi[row] * p.y;
        const float base_u = 1.0f - (2.0f - 2.0f * dot) * Base::INV_R2;
        const float cross = cross_scale[row] * tangent;
        const float LEFT = fmaxf(0.0f, base_u - cross);
        const float RIGHT = fmaxf(0.0f, base_u + cross);
        on_weight(2 * row, LEFT * LEFT);
        on_weight(2 * row + 1, RIGHT * RIGHT);
      }
    }
  }

  HS_O3_FN __attribute__((always_inline)) Pixel shade_pigment(
      uint16_t pigment, float t, float scale, float noise_value) const {
    HS_PROFILE_DEEP(grd_shader_palette);
    const float FIRST_MASS = (pigment >> 10) * (1.0f / 63.0f);
    const int PALETTE_COUNT = (pigment >> 10) == 63 ? 1 : 2;
    const float NOISE_POSITION =
        hs::clamp((noise_value + 1.0f) * (0.5f * (COLOR_NOISE_STEPS - 1)), 0.0f,
                  static_cast<float>(COLOR_NOISE_STEPS - 1));
    const int NOISE_ROW =
        std::min(static_cast<int>(NOISE_POSITION), COLOR_NOISE_STEPS - 2);
    const bool EXACT =
        !color_palette_valid ||
        (color_palette_rows & (3u << NOISE_ROW)) != (3u << NOISE_ROW) ||
        color_palette_exact;
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
      int palette = component ? (pigment >> 5) & 31 : pigment & 31;
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
      const FloatColor *row =
          modified_palettes +
          (palette * COLOR_NOISE_STEPS + NOISE_ROW) * COLOR_VALUE_STEPS +
          VALUE_COLUMN;
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
    return Pixel(
        static_cast<uint16_t>(hs::clamp(accum_r + 0.5f, 0.0f, 65535.0f)),
        static_cast<uint16_t>(hs::clamp(accum_g + 0.5f, 0.0f, 65535.0f)),
        static_cast<uint16_t>(hs::clamp(accum_b + 0.5f, 0.0f, 65535.0f)));
  }

  HS_HOT_FLASH_MEMBER Pixel shade_render(
      uint16_t pigment, float t, float scale, const math::Vector &center_rv,
      const ReactionGraph::CubemapLUT::Projection *projection) const {
    float noise_value = 0;
    {
      HS_PROFILE_DEEP(grd_shader_noise);
      if (color_noise_enabled)
        noise_value = projection != nullptr
                          ? sample_color_noise(*projection)
                          : sample_color_noise(
                                Base::inverse_orientation.apply(center_rv));
    }
    return shade_pigment(pigment, t, scale, noise_value);
  }

  template <typename Grid>
  HS_HOT_FLASH_MEMBER Pixel shade_pixel_full(
      int seed, const math::Vector &center_rv, const math::Vector *world_nodes,
      const Grid &grid, int x, const uint8_t *hot_flags = nullptr,
      const ReactionGraph::CubemapLUT::Projection *projection = nullptr) const {
    if (seed < 0)
      return Pixel(0, 0, 0);

    constexpr uint32_t SAMPLES = Grid::SAMPLES;
    float weights[SAMPLES] = {}, weighted_b[SAMPLES] = {};
    uint16_t pigment;
    {
      HS_PROFILE_DEEP(grd_shader_stencil);
      int center = Base::template refine_render_center<true>(center_rv,
                                                             world_nodes, seed);
      if (hot_flags && !hot_flags[center])
        return Pixel(0, 0, 0);
      pigment = state.pigment[center];
      accumulate_render_stencil(
          grid, x, center, world_nodes,
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
  static __attribute__((always_inline)) float render_support_limit() {
    return Base::KERNEL_R * 0.99f -
           0.25f * (math::RADIANS_PER_COLUMN<Grid::WIDTH> +
                    math::RADIANS_PER_ROW<Grid::HEIGHT>);
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
  HS_O3_FN Pixel shade_pixel(
      int seed, const math::Vector &center_rv, const math::Vector *world_nodes,
      const Grid &grid, int x, const uint8_t *hot_flags = nullptr,
      const ReactionGraph::CubemapLUT::Projection *projection = nullptr) const {
    if (seed < 0)
      return Pixel(0, 0, 0);

    static_assert(B_CULL_THRESHOLD == B_COLOR_FLOOR);
    constexpr uint32_t SAMPLES = Grid::SAMPLES;
    float w00 = 0, w11 = 0;
    float b00 = 0, b01 = 0, b10 = 0, b11 = 0;
    uint16_t pigment;
    {
      HS_PROFILE_DEEP(grd_shader_stencil);
      int center = Base::template refine_render_center<true>(center_rv,
                                                             world_nodes, seed);
      if (hot_flags && !hot_flags[center])
        return Pixel(0, 0, 0);
      const float LIMIT = render_support_limit<Grid>();
      if (LIMIT <= 0.0f ||
          Base::dist2(center_rv, world_nodes[center]) > LIMIT * LIMIT)
        return shade_pixel_full(seed, center_rv, world_nodes, grid, x,
                                hot_flags, projection);
      pigment = state.pigment[center];
      const float ST = math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::sin_theta[x];
      const float CT = math::TrigLUT<Grid::WIDTH, Grid::HEIGHT>::cos_theta(x);
      const auto &run = ReactionGraph::neighbor_runs
          [ReactionGraph::neighbor_run_index[center]];
      const float RADIAL0 = grid.sin_phi[0] * grid.cos_dtheta;
      const float RADIAL1 = grid.sin_phi[1] * grid.cos_dtheta;
      const float CROSS0 =
          2.0f * grid.sin_phi[0] * grid.sin_dtheta * Base::INV_R2;
      const float CROSS1 =
          2.0f * grid.sin_phi[1] * grid.sin_dtheta * Base::INV_R2;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC unroll 2
#endif
      for (int slot = 0; slot < RD_K + 1; ++slot) {
        const int NODE = slot == 0 ? center : center + run.delta[slot - 1];
        const math::Vector &p = world_nodes[NODE];
        const float B = state.B[NODE] - B_CULL_THRESHOLD / Q16_INV;
        const float RADIAL = CT * p.x + ST * p.z;
        const float TANGENT = ST * p.x - CT * p.z;
        const float DOT0 = RADIAL0 * RADIAL + grid.cos_phi[0] * p.y;
        const float DOT1 = RADIAL1 * RADIAL + grid.cos_phi[1] * p.y;
        const float BASE0 = 1.0f - (2.0f - 2.0f * DOT0) * Base::INV_R2;
        const float BASE1 = 1.0f - (2.0f - 2.0f * DOT1) * Base::INV_R2;
        const float D0 = CROSS0 * TANGENT, D1 = CROSS1 * TANGENT;
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
    // Seven weighted Q16 terms have less than 0.08 signed-mass rounding error.
    constexpr float COVERAGE_ROUNDING_GUARD = 0.125f;
    if (fabsf(b00) < COVERAGE_ROUNDING_GUARD ||
        fabsf(b01) < COVERAGE_ROUNDING_GUARD ||
        fabsf(b10) < COVERAGE_ROUNDING_GUARD ||
        fabsf(b11) < COVERAGE_ROUNDING_GUARD)
      return shade_pixel_full(seed, center_rv, world_nodes, grid, x, hot_flags,
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
   * @brief Advances the sim STEPS_PER_FRAME substeps and rasterizes the B field.
   * @param canvas Destination canvas to draw the sphere into.
   * @details Rasterizes the B field onto the sphere via the orientation-aware
   * SSAA shader pipeline after advancing the simulation. Reseeding uses one
   * black frame followed by a complete seed frame, both without partial-state
   * chemistry or rendering.
   */
  void render(Canvas &canvas) {
    HS_PROFILE(grd_render);
    {
      HS_PROFILE(grd_color_noise);
      advance_color_noise();
    }
    ScratchScope frame_guard(scratch_arena_a);
    float mean_db = 0.0f;
    if (transition.next_seed < NUM_SEED_CLUSTERS) {
      seed_reaction(transition.next_seed, NUM_SEED_CLUSTERS);
      transition.next_seed = NUM_SEED_CLUSTERS;
      reaction_edited();
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
      if (transition.next_seed < NUM_SEED_CLUSTERS)
        return;
    }
    refresh_color_palettes();

    // Physics scratch is popped; the raster phase reuses the arena for the
    // oriented lattice so the kernel walks stay in world space, plus the
    // two-ring cull flags.
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
    int next_seed = NUM_SEED_CLUSTERS; /**< First cluster awaiting placement. */
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
  FloatColor *modified_palettes = nullptr;
  bool color_palette_valid = false;
  bool color_palette_exact = true;
  bool color_noise_enabled = false;
  uint16_t color_palette_rows = 0;
  uint8_t color_palette_next_row = 0;
  float color_palette_hue = 0.0f;
  float color_palette_shimmer = 0.0f;
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
