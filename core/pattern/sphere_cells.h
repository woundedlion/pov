/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_cells.h
 * @brief Scalar spherical Voronoi patterns over canonical point layouts. */

#include <array>
#include <span>

#include "core/math/spherical.h"
#include "core/mesh/solid_tables.h"
#include "core/pattern/field.h"

namespace Pattern {

/** @brief Seed layout of a spherical cell pattern. */
enum class SphereCellLayout : uint8_t {
  TETRAHEDRON,
  CUBE,
  OCTAHEDRON,
  DODECAHEDRON,
  ICOSAHEDRON,
  FIBONACCI,
  SCATTERED
};

/** @brief Scalar drawn inside each spherical Voronoi cell. */
enum class SphereCellOutput : uint8_t { BORDERS, CENTERS, CELL_VALUE };

/** @brief Passive spherical cell controls. */
struct SphereCellsParams {
  float count = 42.0f; ///< Rounded seed count for generated layouts.
  float seed = 0.0f; ///< Rounded seed for scattered positions and cell values.
  float line_thickness = 0.06f; ///< Border half-width in mean cell spacings.
  float line_softness = 0.025f; ///< Border fade width in mean cell spacings.
  uint8_t layout =
      static_cast<uint8_t>(SphereCellLayout::ICOSAHEDRON); ///< Seed layout.
  uint8_t output =
      static_cast<uint8_t>(SphereCellOutput::BORDERS); ///< Scalar output.
};

/** @brief Exact spherical Voronoi cells with frame-prepared seed positions. */
struct SphereCells {
  using Coordinate = SphereDirection;      ///< Admitted unit direction.
  using Params = SphereCellsParams;        ///< Seed layout and scalar controls.
  static constexpr size_t MAX_CELLS = 128; ///< Maximum prepared seed count.

  /** @brief Immutable seed positions and output controls. */
  struct Prepared {
    std::array<math::Vector, MAX_CELLS> centers{}; ///< Unit seed positions.
    uint32_t count = 0;                            ///< Active seed positions.
    uint32_t seed = 0;                             ///< Cell-value hash stream.
    float distance_scale = 1.0f; ///< Inverse square root of mean cell area.
    float line_thickness = 0.0f; ///< Border half-width in mean cell spacings.
    float line_softness = 0.0f;  ///< Border fade width in mean cell spacings.
    SphereCellOutput output = SphereCellOutput::BORDERS; ///< Scalar output.
  };

  /** @brief Checks the finite bounded cell recipe.
   * @param params Seed and output controls.
   * @return Whether preparation admits the recipe. */
  static bool valid(const Params &params) {
    return parameter_in_range(params.count, 4.0f,
                              static_cast<float>(MAX_CELLS)) &&
           parameter_in_range(params.seed, 0.0f, 65535.0f) &&
           parameter_in_range(params.line_thickness, 0.0f, 1.0f) &&
           parameter_in_range(params.line_softness, 1.0f / 1024.0f, 1.0f) &&
           params.layout <= static_cast<uint8_t>(SphereCellLayout::SCATTERED) &&
           params.output <= static_cast<uint8_t>(SphereCellOutput::CELL_VALUE);
  }

  /** @brief Resolves the canonical or generated unit seed positions.
   * @param params Seed and output controls.
   * @return Immutable cell pattern. */
  static Prepared prepare(const Params &params) {
    HS_CHECK(valid(params), "sphere cells: invalid parameters");
    Prepared prepared;
    prepared.seed = static_cast<uint32_t>(lroundf(params.seed));
    prepared.output = static_cast<SphereCellOutput>(params.output);
    prepared.line_thickness = params.line_thickness;
    prepared.line_softness = params.line_softness;
    std::span<const math::Vector> solid;
    switch (static_cast<SphereCellLayout>(params.layout)) {
    case SphereCellLayout::TETRAHEDRON:
      solid = Solids::Tetrahedron::vertices;
      break;
    case SphereCellLayout::CUBE:
      solid = Solids::Cube::vertices;
      break;
    case SphereCellLayout::OCTAHEDRON:
      solid = Solids::Octahedron::vertices;
      break;
    case SphereCellLayout::DODECAHEDRON:
      solid = Solids::Dodecahedron::vertices;
      break;
    case SphereCellLayout::ICOSAHEDRON:
      solid = Solids::Icosahedron::vertices;
      break;
    case SphereCellLayout::FIBONACCI:
    case SphereCellLayout::SCATTERED:
      break;
    }
    if (!solid.empty()) {
      prepared.count = static_cast<uint32_t>(solid.size());
      std::copy(solid.begin(), solid.end(), prepared.centers.begin());
    } else {
      prepared.count = static_cast<uint32_t>(lroundf(params.count));
      hs::Pcg32 random(prepared.seed);
      for (size_t index = 0; index < prepared.count; ++index) {
        prepared.centers[index] =
            params.layout == static_cast<uint8_t>(SphereCellLayout::FIBONACCI)
                ? math::fib_spiral(static_cast<int>(prepared.count), 0.5f,
                                   static_cast<int>(index))
                : math::random_vector(random);
      }
    }
    prepared.distance_scale =
        sqrtf(static_cast<float>(prepared.count) / (4.0f * math::PI_F));
    return prepared;
  }

  /** @brief Samples the cell borders, center falloff or stable cell value.
   * @param coordinate Unit sphere direction.
   * @param prepared Immutable cell pattern.
   * @return Scalar in [-1, 1]. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    assert(prepared.count >= 4 && prepared.count <= MAX_CELLS);
    const auto &point = coordinate.value();
    size_t nearest = 0;
    float nearest_distance = math::distance_squared(point, prepared.centers[0]);
    for (size_t index = 1; index < prepared.count; ++index) {
      const float DISTANCE =
          math::distance_squared(point, prepared.centers[index]);
      if (DISTANCE < nearest_distance) {
        nearest = index;
        nearest_distance = DISTANCE;
      }
    }
    if (prepared.output == SphereCellOutput::CELL_VALUE)
      return 2.0f *
                 math::hash01(static_cast<uint32_t>(nearest), prepared.seed) -
             1.0f;
    if (prepared.output == SphereCellOutput::CENTERS) {
      const float DISTANCE =
          math::angle_between(point, prepared.centers[nearest]) *
          prepared.distance_scale;
      return 1.0f - 2.0f * hs::clamp(DISTANCE, 0.0f, 1.0f);
    }
    float edge_sine = 1.0f;
    for (size_t index = 0; index < prepared.count; ++index) {
      if (index == nearest)
        continue;
      const auto normal = prepared.centers[nearest] - prepared.centers[index];
      const float LENGTH_SQUARED = math::dot(normal, normal);
      if (LENGTH_SQUARED < math::EPS_NORMALIZE_SQ)
        continue;
      const float SINE = math::dot(point, normal) / sqrtf(LENGTH_SQUARED);
      edge_sine = std::min(edge_sine, SINE);
    }
    const float DISTANCE =
        asinf(hs::clamp(edge_sine, 0.0f, 1.0f)) * prepared.distance_scale;
    return 1.0f - 2.0f * math::smooth_ramp(prepared.line_thickness,
                                           prepared.line_thickness +
                                               prepared.line_softness,
                                           DISTANCE);
  }
};

} // namespace Pattern
