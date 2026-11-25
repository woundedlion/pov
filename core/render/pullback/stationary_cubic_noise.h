/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/memory.h"
#include "core/render/pullback/modulation.h"

/** @file stationary_cubic_noise.h
 * @brief Cubic patches over the byte-cube grid for stationary direct noise. */

namespace Pullback::Modulation {

/** @brief Immutable cubic coefficients for each cell on a cube face. */
struct StationaryCubicNoise : ImmutableTable {
  /// Sample nodes per cube face edge.
  static constexpr int N = Pattern::SphereLutView::FACE_STEPS;
  /// Cubic intervals per cube face edge.
  static constexpr int CELLS_PER_EDGE = N - 1;
  /// Signed normalized coefficients for one cubic cell.
  using Cell = std::array<int16_t, 10>;
  /// Edge samples at the nodes and at thirds within each interval.
  using Edges = std::array<std::array<float, 3 * N - 2>, 4>;
  /// Resident cubic coefficients indexed by face and cell.
  std::array<Cell, Pattern::SphereLutView::FACE_COUNT * CELLS_PER_EDGE *
                       CELLS_PER_EDGE>
      cells;
  /// Highest spatial frequency resolved by CELLS_PER_EDGE cells per face edge.
  static constexpr float MAX_SCALE = 1.25f;
  /// Noise recipe and loop offset identifying a bake.
  struct Key {
    /// Generator seed.
    int32_t seed;
    /// Spatial frequency.
    float scale;
    /** @var Key::y
     * @brief Y component of the loop offset. */
    /** @var Key::z
     * @brief Z component of the loop offset. */
    /// X component of the loop offset.
    float x, y, z;
    bool operator==(const Key &) const = default;
  };
  /** @var StationaryCubicNoise::pending
   * @brief Recipe queued for the next bake. */
  /// Identity of the currently baked coefficient table.
  std::optional<Key> resident, pending;

  StationaryCubicNoise() : ImmutableTable{sample, sample_byte_tap} {}

  /** @brief Looks up a direction in the fitted cube.
   * @param p Prepared source borrowing the coefficient table.
   * @param direction Normalized source direction.
   * @return Interpolated noise value. */
  static float sample(const Prepared &p,
                      const Pattern::SphereDirection &direction) {
    return sample_byte_tap(p, Pattern::prepare_sphere_lut_sample(direction));
  }

  /** @brief Evaluates the cell's bilinear part and two cubic edge corrections.
   * @param p Prepared source borrowing the coefficient table.
   * @param tap Face and interpolation coordinates.
   * @return Interpolated noise value. */
  HS_O3_FN HS_NOINLINE_NOCLONE static float
  sample_byte_tap(const Prepared &p, const Pattern::SphereLutSample &tap) {
    const auto &cache =
        *static_cast<const StationaryCubicNoise *>(p.immutable_table);
    const auto &c = cache.cells[tap.cell_offset];
    const float x = tap.x_fraction, y = tap.y_fraction;
    const float base = c[0] + x * c[1] + y * (c[2] + x * c[3]);
    const float horizontal = c[4] + x * c[5] + y * c[6];
    const float vertical = c[7] + y * c[8] + x * c[9];
    return (base + x * (1.0f - x) * horizontal + y * (1.0f - y) * vertical) *
           (1.0f / 32767.0f);
  }

  /** @brief Fits the curvature between two ends from samples at 1/3 and 2/3.
   * @param start Value at the first endpoint.
   * @param first Value one third of the way across.
   * @param second Value two thirds of the way across.
   * @param end Value at the last endpoint.
   * @return Coefficients for the curved part of the edge. */
  static std::array<float, 2> curvature(float start, float first, float second,
                                        float end) {
    const float low = 4.5f * (first - hs::lerp(start, end, 1.f / 3.f));
    const float high = 4.5f * (second - hs::lerp(start, end, 2.f / 3.f));
    return {2.0f * low - high, 3.0f * (high - low)};
  }

  /** @brief Stores the total-degree-three fit in signed normalized coefficients.
   * @param cell Destination coefficients.
   * @param edges Lower and upper row samples.
   * @param x First sample index of the cell. */
  HS_SIZE_FLASH_MEMBER static void fit(Cell &cell, const Edges &edges, int x) {
    const float a = edges[0][x], b = edges[0][x + 3];
    const float c = edges[3][x], d = edges[3][x + 3];
    const auto bottom = curvature(a, edges[0][x + 1], edges[0][x + 2], b);
    const auto top = curvature(c, edges[3][x + 1], edges[3][x + 2], d);
    const auto left = curvature(a, edges[1][x], edges[2][x], c);
    const auto right = curvature(b, edges[1][x + 3], edges[2][x + 3], d);
    const std::array coefficients = {a,
                                     b - a,
                                     c - a,
                                     d - c - b + a,
                                     bottom[0],
                                     0.5f * (bottom[1] + top[1]),
                                     top[0] - bottom[0],
                                     left[0],
                                     0.5f * (left[1] + right[1]),
                                     right[0] - left[0]};
    for (size_t i = 0; i < cell.size(); ++i)
      cell[i] = static_cast<int16_t>(
          roundf(hs::clamp(coefficients[i], -1.f, 1.f) * 32767.f));
  }

  /** @brief Admits stationary single-octave simplex at the cubic cell resolution.
   * @param config Candidate source recipe.
   * @return Whether the recipe is supported by this cubic table. */
  static constexpr bool supports(const Config &config) {
    return config.enabled && config.source == Pattern::SourceKind::NOISE &&
           config.evaluation == Pattern::Evaluation::DIRECT &&
           config.noise.variant == Pattern::NoiseVariant::CONFIGURED_SINGLE &&
           config.noise.speed == 0.0f && config.noise.scale <= MAX_SCALE;
  }

  /** @brief Fills the immutable coefficients from the canonical direct generator.
   * @param p Prepared canonical noise generator and spatial controls. */
  HS_SIZE_FLASH_MEMBER void bake(const Prepared &p) {
    constexpr int intervals = 3 * (N - 1);
    ScratchScope scratch(scratch_arena_a);
    auto &edges = *scratch_arena_a.make_default_initialized<Edges>();
    for (int face = 0; face < 6; ++face) {
      for (int y = 0; y < N - 1; ++y) {
        for (int row = y == 0 ? 0 : 1; row < 4; ++row) {
          const float v = float(2 * (3 * y + row) - intervals) / intervals;
          const int step = row == 1 || row == 2 ? 3 : 1;
          for (int x = 0; x <= intervals; x += step) {
            const float u = float(2 * x - intervals) / intervals;
            edges[row][x] = Pattern::SphereNoise::sample(
                Pattern::sphere_face_coordinate(face, u, v),
                std::get<0>(p.source));
          }
        }
        for (int x = 0; x < N - 1; ++x)
          fit(cells[(face * CELLS_PER_EDGE + y) * CELLS_PER_EDGE + x], edges,
              3 * x);
        edges[0] = edges[3];
      }
    }
  }

  /** @brief Binds a stable generator once its previous bake is no longer borrowed.
   * @param p Freshly prepared canonical endpoint.
   * @param retained Endpoint still in use by the retained frame, if any. */
  HS_HOT_FLASH_MEMBER void bind(Prepared &p, const Prepared *retained) {
    if (!supports(p.config) || p.source.index() != 0 ||
        (p.evaluate != Detail::DIRECT_SAMPLERS[0] &&
         p.evaluate != Detail::sample_unit_noise_origin)) {
      pending.reset();
      return;
    }
    const auto &noise = std::get<0>(p.source);
    FastNoiseLite generator;
    init_effect_noise(generator, p.config.seed);
    if (!noise.noise.same_noise_single_3d(generator)) {
      pending.reset();
      return;
    }
    const Key key{p.config.seed, noise.scale, noise.loop_offset.x,
                  noise.loop_offset.y, noise.loop_offset.z};
    if (resident != key) {
      if (resident &&
          (pending != key || (retained && retained->immutable_table == this))) {
        pending = key;
        return;
      }
      bake(p);
      resident = key;
    }
    pending.reset();
    p.immutable_table = this;
    p.evaluate = sample;
  }
};
} // namespace Pullback::Modulation
