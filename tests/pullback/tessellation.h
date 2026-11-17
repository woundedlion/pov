/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Tessellation fields against independent geometric references. */

#include "core/render/pullback.h"
#include "tests/test_harness.h"
#include "tests/test_fixture.h"

namespace hs_test::pullback_tests {

/** @brief Square cells evaluate the circle-to-box blend formula. */
inline void test_pullback_square_cells_formula() {
  const auto reference = [](const Pattern::PlanePosition &coordinate,
                            const Pattern::TessellationCellParams &params) {
    const float x =
        math::wrap_t(params.cell_scale * coordinate.value.re + 0.5f) - 0.5f;
    const float y =
        math::wrap_t(params.cell_scale * coordinate.value.im + 0.5f) - 0.5f;
    const float circle = sqrtf(x * x + y * y) - params.cell_radius;
    const float bx = fabsf(x) - params.cell_radius;
    const float by = fabsf(y) - params.cell_radius;
    const float square = sqrtf(fmaxf(bx, 0.0f) * fmaxf(bx, 0.0f) +
                               fmaxf(by, 0.0f) * fmaxf(by, 0.0f)) +
                         fminf(fmaxf(bx, by), 0.0f);
    return 1.0f - 2.0f * math::smooth_ramp(
                             -params.cell_softness, params.cell_softness,
                             hs::lerp(circle, square, params.cell_shape));
  };
  hs::Pcg32 random(0x12345678u);
  for (float scale : {1.0f / 64.0f, 0.98971874f, 1.0f, 8.0f})
    for (float radius : {1.0f / 64.0f, 0.25f, 0.37757313f, 0.49f})
      for (float softness : {1.0f / 1024.0f, 0.1f, 0.34563965f, 1.0f})
        for (float blend : {0.0f, 0.00001f, 0.25f, 0.99999f, 1.0f}) {
          Pattern::TessellationCellParams params;
          params.cell_scale = scale;
          params.cell_radius = radius;
          params.cell_softness = softness;
          params.cell_shape = blend;
          for (int sample = 0; sample < 256; ++sample) {
            const float x = rand_uniform(random, -1024.0f, 1024.0f);
            const float y = rand_uniform(random, -1024.0f, 1024.0f);
            const Pattern::PlanePosition coordinate{
                math::Complex(sample < 8 ? float(sample - 4) * 0.5f : x,
                              sample < 8 ? 0.0f : y)};
            const float expected = reference(coordinate, params);
            const float actual = Pattern::square_cells(coordinate, params);
            HS_EXPECT_NEAR(actual, expected,
                           8.0f * std::numeric_limits<float>::epsilon());
          }
        }
}

/** @brief Independent double-precision tessellation reference. */
struct TessellationReference {
  /** Half-turn angle, radians. */
  static constexpr double HALF_TURN = 3.14159265358979323846;
  Pattern::TessellationKind kind;           /**< Cell geometry. */
  Pattern::TessellationFill fill;           /**< Cell or edge fill. */
  Pattern::TessellationSourceParams params; /**< Field parameters. */

  /** @brief Signed distance to a regular polygon.
   * @param x Plane x coordinate.
   * @param y Plane y coordinate.
   * @param sides Polygon side count.
   * @param angle First outward normal angle, radians.
   * @param radius Polygon inradius.
   * @return Negative inside, positive outside. */
  static double polygon_distance(double x, double y, int sides, double angle,
                                 double radius) {
    const double step = 2.0 * HALF_TURN / sides;
    double inside = -1e30;
    for (int edge = 0; edge < sides; ++edge) {
      const double theta = angle + edge * step;
      inside =
          std::max(inside, x * std::cos(theta) + y * std::sin(theta) - radius);
    }
    if (inside <= 0.0)
      return inside;
    const double circumradius = radius / std::cos(HALF_TURN / sides);
    double outside = 1e30;
    for (int edge = 0; edge < sides; ++edge) {
      const double theta = angle + edge * step;
      const double ax = circumradius * std::cos(theta - 0.5 * step);
      const double ay = circumradius * std::sin(theta - 0.5 * step);
      const double bx = circumradius * std::cos(theta + 0.5 * step);
      const double by = circumradius * std::sin(theta + 0.5 * step);
      const double t =
          std::clamp(((x - ax) * (bx - ax) + (y - ay) * (by - ay)) /
                         ((bx - ax) * (bx - ax) + (by - ay) * (by - ay)),
                     0.0, 1.0);
      outside = std::min(
          outside, std::hypot(x - ax - t * (bx - ax), y - ay - t * (by - ay)));
    }
    return outside;
  }

  /** @brief Cubic smooth ramp.
   * @param edge0 Lower boundary.
   * @param edge1 Upper boundary.
   * @param value Sample coordinate.
   * @return Ramp value in [0, 1]. */
  static double ramp(double edge0, double edge1, double value) {
    const double t = std::clamp((value - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
  }

  /** @brief Evaluates the reference field.
   * @param point Plane sample position.
   * @return Field value in [-1, 1]. */
  double operator()(const math::Complex &point) const {
    using K = Pattern::TessellationKind;
    const double sqrt3 = std::sqrt(3.0);
    const double u = double(params.cell_scale) * point.re;
    const double v = double(params.cell_scale) * point.im;
    double x = 0.0, y = 0.0, angle = 0.0;
    int sides = 4;
    if (kind == K::SQUARE) {
      x = u - std::round(u);
      y = v - std::round(v);
    } else if (kind == K::HEXAGONAL) {
      sides = 6;
      const double row = std::round(v / (0.5 * sqrt3));
      const double column = std::round(u - 0.5 * row);
      double best = 1e30;
      for (double j = row - 2; j <= row + 2; ++j)
        for (double i = column - 2; i <= column + 2; ++i) {
          const double dx = u - (i + 0.5 * j);
          const double dy = v - 0.5 * sqrt3 * j;
          if (dx * dx + dy * dy < best) {
            best = dx * dx + dy * dy;
            x = dx;
            y = dy;
          }
        }
    } else {
      sides = 3;
      const double b = v / (0.5 * sqrt3);
      const double a = u - 0.5 * b;
      double best = -1e30;
      for (double j = std::floor(b) - 1; j <= std::floor(b) + 1; ++j)
        for (double i = std::floor(a) - 1; i <= std::floor(a) + 1; ++i)
          for (const bool up : {true, false}) {
            const double third = up ? 1.0 / 3.0 : 2.0 / 3.0;
            const double cx = (i + third) + 0.5 * (j + third);
            const double cy = 0.5 * sqrt3 * (j + third);
            const double dx = sqrt3 * (u - cx);
            const double dy = sqrt3 * (v - cy);
            const double facing = up ? -0.5 * HALF_TURN : 0.5 * HALF_TURN;
            const double depth = -polygon_distance(dx, dy, 3, facing, 0.5);
            if (depth > best) {
              best = depth;
              x = dx;
              y = up ? dy : -dy;
            }
          }
      angle = -0.5 * HALF_TURN;
    }
    if (fill == Pattern::TessellationFill::EDGES) {
      const double distance = -polygon_distance(x, y, sides, angle, 0.5);
      return 1.0 -
             2.0 * ramp(params.line_thickness,
                        double(params.line_thickness) + params.line_softness,
                        distance);
    }
    const double circle = std::hypot(x, y) - params.cell_radius;
    const double polygon =
        polygon_distance(x, y, sides, angle, params.cell_radius);
    const double distance = circle + (polygon - circle) * params.cell_shape;
    return 1.0 - 2.0 * ramp(-double(params.cell_softness),
                            double(params.cell_softness), distance);
  }
};

/** @brief Every tessellation kind and fill matches an independent polygon
    reference. */
inline void test_pullback_tessellation_reference() {
  using K = Pattern::TessellationKind;
  using F = Pattern::TessellationFill;
  hs::Pcg32 random(0x9e3779b9u);
  const auto next = [&] { return rand_uniform(random, -3.0f, 3.0f); };
  for (const K kind : {K::TRIANGULAR, K::SQUARE, K::HEXAGONAL})
    for (const F fill : {F::CELLS, F::EDGES})
      for (const float shape : {0.0f, 0.35f, 1.0f}) {
        HS_CONTEXT("kind / fill", static_cast<int>(kind),
                   static_cast<int>(fill));
        HS_CONTEXT("shape percent", static_cast<int>(shape * 100.0f));
        TessellationReference reference{kind, fill, {}};
        reference.params.cell_scale = 1.7f;
        reference.params.cell_shape = shape;
        reference.params.cell_radius = 0.3f;
        reference.params.cell_softness = 0.08f;
        reference.params.line_thickness = 0.05f;
        reference.params.line_softness = 0.04f;
        int inside = 0;
        int outside = 0;
        for (int sample = 0; sample < 2000; ++sample) {
          const float x = next();
          const float y = next();
          const math::Complex point(x, y);
          const float actual = Pattern::tessellation(
              Pattern::PlanePosition{point}, reference.params, kind, fill);
          const double expected = reference(point);
          HS_EXPECT_NEAR(actual, expected, 2e-4);
          inside += expected > 0.99;
          outside += expected < -0.99;
        }
        HS_EXPECT_GT(inside, 100);
        HS_EXPECT_GT(outside, 100);
      }
}

/** @brief Hexagonal cells pack circles and hexagons on the hexagonal lattice,
    with each hexagon's edges facing its six neighbours. */
inline void test_pullback_hexagonal_cells_pack_hex_lattice() {
  using K = Pattern::TessellationKind;
  using F = Pattern::TessellationFill;
  Pattern::TessellationCellParams params;
  params.cell_scale = 2.5f;
  params.cell_radius = 0.4f;
  params.cell_softness = 0.05f;
  const float PERIOD = 0.4f;
  const float HALF_SQRT3 = 0.8660254038f;
  const auto at = [&](float u, float v) {
    return Pattern::tessellation<K::HEXAGONAL, F::CELLS>(
        Pattern::PlanePosition{math::Complex(PERIOD * u, PERIOD * v)}, params);
  };
  for (const float shape : {0.0f, 1.0f}) {
    HS_CONTEXT("shape", static_cast<int>(shape));
    params.cell_shape = shape;
    for (int j = -3; j <= 3; ++j)
      for (int i = -3; i <= 3; ++i) {
        const float u = i + 0.5f * j;
        const float v = HALF_SQRT3 * j;
        HS_EXPECT_EQ(at(u, v), 1.0f);
        HS_EXPECT_EQ(at(u + 0.5f, v), -1.0f);
        HS_EXPECT_EQ(at(u + 0.25f, v + 0.5f * HALF_SQRT3), -1.0f);
        HS_EXPECT_EQ(at(u - 0.25f, v + 0.5f * HALF_SQRT3), -1.0f);
        HS_EXPECT_LT(at(u + 0.42f, v), 0.0f);
        const float corner = at(u + 0.44f * HALF_SQRT3, v + 0.44f * 0.5f);
        if (shape == 0.0f)
          HS_EXPECT_LT(corner, 0.0f);
        else
          HS_EXPECT_GT(corner, 0.0f);
      }
  }
}

/** @brief Field of one fixed tessellation at a plane point.
 * @tparam Kind Cell geometry.
 * @tparam Fill Cell or edge fill.
 * @param params Field parameters.
 * @param point Plane sample position.
 * @return Field value in [-1, 1]. */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill>
float fixed_tessellation(
    const Pullback::Source::FixedTessellationParams<Kind, Fill> &params,
    const math::Complex &point) {
  return Pattern::tessellation<Kind, Fill>(Pattern::PlanePosition{point},
                                           params);
}

/** @brief One fixed tessellation is invariant under its reported basis and
    not under half of it.
 * @tparam Kind Cell geometry.
 * @tparam Fill Cell or edge fill. */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill>
void check_tessellation_translation_lattice() {
  HS_CONTEXT("kind / fill", static_cast<int>(Kind), static_cast<int>(Fill));
  Pullback::Source::FixedTessellationParams<Kind, Fill> params;
  params.cell_scale = 1.7f;
  if constexpr (Fill == Pattern::TessellationFill::CELLS) {
    params.cell_shape = 0.8f;
    params.cell_radius = 0.35f;
    params.cell_softness = 0.1f;
  } else {
    params.line_thickness = 0.06f;
    params.line_softness = 0.05f;
  }
  const Pullback::TranslationLattice lattice =
      Pullback::Source::translation_lattice(params);
  const float PERIOD = 1.0f / 1.7f;
  const bool square = Kind == Pattern::TessellationKind::SQUARE;
  HS_EXPECT_EQ(lattice.a.re, PERIOD);
  HS_EXPECT_EQ(lattice.a.im, 0.0f);
  HS_EXPECT_NEAR(lattice.b.re, square ? 0.0f : 0.5f * PERIOD, 1e-7f);
  HS_EXPECT_NEAR(lattice.b.im, square ? PERIOD : 0.8660254038f * PERIOD, 1e-7f);
  bool half_differs = false;
  hs::Pcg32 random(0x2545f491u);
  for (int sample = 0; sample < 300; ++sample) {
    const float x = rand_uniform(random, -2.0f, 2.0f);
    const float y = rand_uniform(random, -2.0f, 2.0f);
    const math::Complex point(x, y);
    const float base = fixed_tessellation(params, point);
    for (const math::Complex shift :
         {lattice.a, lattice.b, lattice.a - lattice.b,
          math::Complex() - lattice.a - lattice.b})
      HS_EXPECT_NEAR(fixed_tessellation(params, point + shift), base, 1e-4f);
    const math::Complex half(0.5f * lattice.a.re, 0.5f * lattice.a.im);
    half_differs |=
        fabsf(fixed_tessellation(params, point + half) - base) > 0.5f;
  }
  HS_EXPECT_TRUE(half_differs);
}

/** @brief Each tessellation kind reports a translation lattice that leaves
    its field unchanged. */
inline void test_pullback_tessellation_translation_lattice() {
  using K = Pattern::TessellationKind;
  using F = Pattern::TessellationFill;
  check_tessellation_translation_lattice<K::TRIANGULAR, F::CELLS>();
  check_tessellation_translation_lattice<K::SQUARE, F::CELLS>();
  check_tessellation_translation_lattice<K::HEXAGONAL, F::CELLS>();
  check_tessellation_translation_lattice<K::TRIANGULAR, F::EDGES>();
  check_tessellation_translation_lattice<K::SQUARE, F::EDGES>();
  check_tessellation_translation_lattice<K::HEXAGONAL, F::EDGES>();
}

} // namespace hs_test::pullback_tests
