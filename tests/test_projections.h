/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Projection coordinates, cuts, edge identities, and Peirce layout contracts.
 */
#pragma once

#include <array>
#include <bit>
#include <limits>

#include "core/math/projections.h"
#include "core/render/pullback/projection.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace projections_tests {

using namespace projections;

inline constexpr float INV_SQRT_TWO = 0.7071067811865475f;

constexpr float PEIRCE_QUARTER_PERIOD = 1.8540746773013719f;
static_assert(std::size(AIROCEAN_CENTERS) == AIROCEAN_FACE_COUNT);
static_assert(std::size(AIROCEAN_NORMALS) == AIROCEAN_FACE_COUNT);
static_assert(std::size(AIROCEAN_PLANAR_FACES) == AIROCEAN_FACE_COUNT);
static_assert(std::size(AIROCEAN_TRANSFORMS) == AIROCEAN_FACE_COUNT);
static_assert(std::size(AIROCEAN_CUT_MASKS) == AIROCEAN_FACE_COUNT);

inline math::Vector direction(float latitude, float longitude) {
  const float radius = cosf(latitude);
  return math::Vector(radius * cosf(longitude), sinf(latitude),
                      radius * sinf(longitude));
}

/** @brief The kernel's y-up to z-up permutation, as airocean_projection does it
 *         at a zero central meridian. */
inline AiroceanVector airocean_axes(const math::Vector &v) {
  return AiroceanVector{v.x, v.z, v.y};
}

inline bool planar_points_equal(const AiroceanPoint &a,
                                const AiroceanPoint &b) {
  return a.x == b.x && a.y == b.y;
}

inline bool planar_segments_equal(const AiroceanPoint &a,
                                  const AiroceanPoint &b,
                                  const AiroceanPoint &c,
                                  const AiroceanPoint &d) {
  return (planar_points_equal(a, c) && planar_points_equal(b, d)) ||
         (planar_points_equal(a, d) && planar_points_equal(b, c));
}

inline bool spherical_points_equal(const AiroceanVector &a,
                                   const AiroceanVector &b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool spherical_segments_equal(const AiroceanVector &a,
                                     const AiroceanVector &b,
                                     const AiroceanVector &c,
                                     const AiroceanVector &d) {
  return (spherical_points_equal(a, c) && spherical_points_equal(b, d)) ||
         (spherical_points_equal(a, d) && spherical_points_equal(b, c));
}

inline const AiroceanPoint &planar_vertex(size_t face, size_t vertex) {
  return AIROCEAN_PLANAR_FACES[face][vertex % 3];
}

/** @brief Twice the signed area of the planar triangle a-b-c. */
inline float planar_cross(const AiroceanPoint &a, const AiroceanPoint &b,
                          const AiroceanPoint &c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

/** @brief Gnomonic projection onto a face's plane followed by that face's
 *         affine unfold, spelled out independently of the kernel. */
inline AiroceanPoint airocean_unfold(size_t face, const AiroceanVector &p) {
  const AiroceanVector &center = AIROCEAN_CENTERS[face];
  const AiroceanVector &normal = AIROCEAN_NORMALS[face];
  const float plane =
      center.x * normal.x + center.y * normal.y + center.z * normal.z;
  const float ray = p.x * normal.x + p.y * normal.y + p.z * normal.z;
  const float scale = plane / ray;
  const float (&transform)[2][4] = AIROCEAN_TRANSFORMS[face];
  return AiroceanPoint{
      transform[0][0] * p.x * scale + transform[0][1] * p.y * scale +
          transform[0][2] * p.z * scale + transform[0][3],
      transform[1][0] * p.x * scale + transform[1][1] * p.y * scale +
          transform[1][2] * p.z * scale + transform[1][3]};
}

/** @brief Longitude wrapping is idempotent, preserves direction and maps pi to negative pi. */
inline void test_wrap_longitude_range() {
  for (int step = -40; step <= 40; ++step) {
    const float raw = step * 0.4f;
    const float wrapped = wrap_longitude(raw);
    HS_EXPECT_GE(wrapped, -math::PI_F);
    HS_EXPECT_LT(wrapped, math::PI_F);
    HS_EXPECT_NEAR(cosf(wrapped), cosf(raw), 2e-6f);
    HS_EXPECT_NEAR(sinf(wrapped), sinf(raw), 2e-6f);
    HS_EXPECT_NEAR(wrap_longitude(wrapped), wrapped, 1e-6f);
  }
  HS_EXPECT_EQ(wrap_longitude(math::PI_F), -math::PI_F);
  HS_EXPECT_EQ(wrap_longitude(-math::PI_F), -math::PI_F);
  HS_EXPECT_EQ(wrap_longitude(3.0f * math::PI_F), -math::PI_F);
  const float below_cut = wrap_longitude(
      std::nextafter(-math::PI_F, -std::numeric_limits<float>::infinity()));
  HS_EXPECT_LT(below_cut, math::PI_F);
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(below_cut),
               std::bit_cast<uint32_t>(math::PI_F) - 2u);
}

/** @brief Longitude reduction retains fmod rounding at cuts and large inputs. */
inline void test_wrap_longitude_fmod_contract() {
  const auto reference = [](float raw) {
    float reduced = std::fmod(raw + math::PI_F, math::TWO_PI_F);
    if (reduced < 0.0f)
      reduced += math::TWO_PI_F;
    return reduced - math::PI_F;
  };
  const auto check = [&](float raw) {
    const float expected = reference(raw);
    const float actual = wrap_longitude(raw);
    if (std::isnan(expected))
      HS_EXPECT_TRUE(std::isnan(actual));
    else
      HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual),
                   std::bit_cast<uint32_t>(expected));
  };
  const float INFINITY_VALUE = std::numeric_limits<float>::infinity();
  for (float value : {-3.0f * math::PI_F, -math::PI_F, -0.0f, 0.0f, math::PI_F,
                      3.0f * math::PI_F}) {
    check(std::nextafter(value, -INFINITY_VALUE));
    check(value);
    check(std::nextafter(value, INFINITY_VALUE));
  }
  for (float value : {std::numeric_limits<float>::lowest(),
                      std::numeric_limits<float>::max(), -INFINITY_VALUE,
                      INFINITY_VALUE, std::numeric_limits<float>::quiet_NaN()})
    check(value);
  uint32_t state = 0x8e2bd471u;
  for (int index = 0; index < 4096; ++index) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    check(std::bit_cast<float>(state));
  }
}

/** @brief Bonne converges to sinusoidal coordinates as the standard parallel approaches zero. */
inline void test_bonne_sinusoidal_limit() {
  for (const float pole : {-1.0f, 1.0f}) {
    const auto result =
        bonne_projection(math::Vector(0.0f, pole, 0.0f), 0.7f, 0.0f);
    HS_EXPECT_EQ(result.coords.re, 0.0f);
    HS_EXPECT_EQ(result.fade_edge_distance, 0.0f);
  }
  for (int latitude_step = -12; latitude_step <= 12; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * math::PI_F / 12.0f);
    for (int longitude_step = 0; longitude_step < 16; ++longitude_step) {
      const float longitude =
          longitude_step * (math::TWO_PI_F / 16.0f) - math::PI_F + 0.1f;
      const math::Vector v = direction(latitude, longitude);
      const ProjectionKernelResult flat = bonne_projection(v, 0.0f, 0.0f);
      const float measured_latitude = asinf(hs::clamp(v.y, -1.0f, 1.0f));
      const float measured_longitude = wrap_longitude(atan2f(v.z, v.x));
      HS_EXPECT_NEAR(flat.coords.re,
                     measured_longitude * cosf(measured_latitude), 2e-6f);
      HS_EXPECT_NEAR(flat.coords.im, measured_latitude, 2e-6f);
      const ProjectionKernelResult near_flat = bonne_projection(v, 0.0f, 1e-3f);
      HS_EXPECT_NEAR(near_flat.coords.re, flat.coords.re, 0.02f);
      HS_EXPECT_NEAR(near_flat.coords.im, flat.coords.im, 0.02f);
    }
  }
}

/** @brief General and polar Bonne coordinates match Snyder’s equations. */
inline void test_bonne_general_parallel_matches_snyder() {
  for (float parallel : {0.4f, -0.4f, 0.7f, -0.7f, 1.2f, -1.2f,
                         math::PI_F * .5f, -math::PI_F * .5f}) {
    for (float latitude : {-1.3f, -.6f, .2f, 1.1f})
      for (float longitude : {-2.4f, -.7f, .4f, 2.2f}) {
        const auto v = direction(latitude, longitude);
        const auto actual = bonne_projection(v, 0.0f, parallel);
        const double cot = fabsf(parallel) == math::PI_F * .5f
                               ? 0.0
                               : 1.0 / std::tan(double(parallel));
        const double rho = cot + parallel - latitude;
        const double angle = longitude * std::cos(double(latitude)) / rho;
        HS_EXPECT_NEAR(actual.coords.re, rho * std::sin(angle), 1e-5);
        HS_EXPECT_NEAR(actual.coords.im, cot - rho * std::cos(angle), 1e-5);
      }
    const auto origin =
        bonne_projection(direction(parallel, 0.0f), 0.0f, parallel);
    HS_EXPECT_NEAR(origin.coords.re, 0.0f, 1e-5f);
    HS_EXPECT_NEAR(origin.coords.im, 0.0f, 1e-5f);
  }
}

/** @brief Bonne edge metadata measures the nearest longitude cut. */
inline void test_bonne_general_parallel_edge_metadata() {
  for (float parallel :
       {.4f, -.4f, .7f, -.7f, 1.2f, -1.2f, math::PI_F * .5f, -math::PI_F * .5f})
    for (float latitude : {-1.3f, -.6f, .2f, 1.1f})
      for (float longitude : {-2.4f, -.7f, .4f, 2.2f}) {
        const auto actual =
            bonne_projection(direction(latitude, longitude), 0, parallel);
        HS_EXPECT_NEAR(actual.fade_edge_distance,
                       std::cos(double(latitude)) *
                           (math::PI_F - std::abs(double(longitude))),
                       1e-5);
        HS_EXPECT_EQ(actual.region_id, static_cast<uint32_t>(longitude < 0));
      }
}

/** @brief Bonne coordinates are invariant under a shared meridian rotation. */
inline void test_bonne_general_parallel_meridian_invariance() {
  for (float parallel :
       {.4f, -.4f, .7f, -.7f, 1.2f, -1.2f, math::PI_F * .5f, -math::PI_F * .5f})
    for (float latitude : {-1.3f, -.6f, .2f, 1.1f})
      for (float longitude : {-2.4f, -.7f, .4f, 2.2f}) {
        const auto actual =
            bonne_projection(direction(latitude, longitude), 0, parallel);
        constexpr float MERIDIAN = .37f;
        const auto shifted = bonne_projection(
            direction(latitude, longitude + MERIDIAN), MERIDIAN, parallel);
        HS_EXPECT_NEAR(shifted.coords.re, actual.coords.re, 1e-5);
        HS_EXPECT_NEAR(shifted.coords.im, actual.coords.im, 1e-5);
      }
}

/** @brief Peirce meridians and strip scroll preserve direction and magnitude. */
inline void test_peirce_meridian_and_scroll_coordinates() {
  constexpr float PERIOD = 4.0f * PEIRCE_QUARTER_PERIOD;
  for (const auto layout : {PeirceLayout::DIAMOND, PeirceLayout::SQUARE,
                            PeirceLayout::HORIZONTAL, PeirceLayout::VERTICAL})
    for (float latitude : {-1.0f, -.3f, .4f, 1.1f})
      for (float longitude : {-2.0f, -.4f, .3f, 1.7f}) {
        const auto v = direction(latitude, longitude);
        const auto base = peirce_projection(v, 0.0f, layout, 0.0f);
        constexpr float MERIDIAN = .23f;
        const auto shifted = peirce_projection(
            direction(latitude, longitude + MERIDIAN), MERIDIAN, layout, 0.0f);
        HS_EXPECT_NEAR(shifted.coords.re, base.coords.re, 1e-5f);
        HS_EXPECT_NEAR(shifted.coords.im, base.coords.im, 1e-5f);
        if (layout == PeirceLayout::HORIZONTAL ||
            layout == PeirceLayout::VERTICAL)
          for (float scroll : {.125f, .375f}) {
            const auto moved = peirce_projection(v, 0.0f, layout, scroll);
            const bool horizontal = layout == PeirceLayout::HORIZONTAL;
            const float delta = horizontal ? moved.coords.re - base.coords.re
                                           : moved.coords.im - base.coords.im;
            HS_EXPECT_NEAR(std::remainder(delta, PERIOD),
                           std::remainder(scroll * PERIOD, PERIOD), 1e-5f);
            HS_EXPECT_NEAR(horizontal ? moved.coords.im : moved.coords.re,
                           horizontal ? base.coords.im : base.coords.re, 1e-5f);
          }
      }
}

/** @brief The elliptic integral is odd, monotone and matches reference endpoint values. */
inline void test_peirce_elliptic_integral_shape() {
  HS_EXPECT_NEAR(peirce_elliptic_integral(0.25f * math::PI_F),
                 0.8260178762492452f, 1e-6f);
  HS_EXPECT_NEAR(peirce_elliptic_integral(0.0f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(peirce_elliptic_integral(0.5f * math::PI_F),
                 PEIRCE_QUARTER_PERIOD, 1e-6f);
  float previous = peirce_elliptic_integral(-0.5f * math::PI_F);
  HS_EXPECT_NEAR(previous, -PEIRCE_QUARTER_PERIOD, 1e-6f);
  for (int step = -63; step <= 64; ++step) {
    const float phi = step * (0.5f * math::PI_F / 64.0f);
    const float value = peirce_elliptic_integral(phi);
    HS_EXPECT_GT(value, previous);
    HS_EXPECT_NEAR(peirce_elliptic_integral(-phi), -value, 1e-6f);
    previous = value;
  }
}

/** @brief Sector boundaries snap within their tie band while interior longitudes remain unchanged. */
inline void test_peirce_sector_longitude_snapping() {
  HS_EXPECT_EQ(peirce_sector_longitude(math::Vector(0.0f, 1.0f, 0.0f), 0.0f),
               0.5f * math::PI_F);
  HS_EXPECT_EQ(peirce_sector_longitude(math::Vector(0.0f, -1.0f, 0.0f), 0.0f),
               0.5f * math::PI_F);
  constexpr std::array<float, 4> BOUNDARIES = {
      -0.75f * math::PI_F, -0.25f * math::PI_F, 0.25f * math::PI_F,
      0.75f * math::PI_F};
  for (float boundary : BOUNDARIES)
    for (float offset : {-1e-6f, 0.0f, 1e-6f}) {
      const float angle = boundary + offset;
      const math::Vector v(cosf(angle), 0.0f, sinf(angle));
      HS_EXPECT_EQ(peirce_sector_longitude(v, 0.0f), boundary);
    }
  const math::Vector interior(cosf(0.5f), 0.0f, sinf(0.5f));
  HS_EXPECT_NEAR(peirce_sector_longitude(interior, 0.0f), 0.5f, 1e-6f);
}

inline void check_peirce_fast_square_matches_exact(const math::Vector &v) {
  const ProjectionKernelResult exact =
      peirce_projection(v, 0.0f, projections::PeirceLayout::SQUARE, 0.0f);
  const ProjectionKernelResult fast = peirce_projection_fast_square(v);
  HS_EXPECT_EQ(fast.region_id, exact.region_id);
  HS_EXPECT_EQ(fast.edge_class, exact.edge_class);
  HS_EXPECT_EQ(fast.flags, exact.flags);
  HS_EXPECT_EQ(fast.traits, exact.traits);
  HS_EXPECT_EQ(fast.boundary_flags, exact.boundary_flags);
  HS_EXPECT_EQ(fast.component_id, exact.component_id);
  constexpr auto &METRICS = Pullback::Projection::PEIRCE_FAST_SQUARE_METRICS;
  static_assert(METRICS[0].domain ==
                    Pullback::ApproximationDomain::PROJECTED_COORDINATE &&
                METRICS[1].domain ==
                    Pullback::ApproximationDomain::PROJECTED_EDGE_DISTANCE);
  HS_EXPECT_NEAR(fast.coords.re, exact.coords.re, METRICS[0].limit);
  HS_EXPECT_NEAR(fast.coords.im, exact.coords.im, METRICS[0].limit);
  HS_EXPECT_NEAR(fast.fade_edge_distance, exact.fade_edge_distance,
                 METRICS[1].limit);
}

/** @brief The fast Peirce kernel matches exact coordinates, edge distance and metadata across the sphere. */
inline void test_peirce_fast_square_matches_exact() {
  // Half-step azimuths: the sector seams sit at even multiples of pi/96, so
  // every off-pole sample lands in a sector interior.
  for (int latitude_step = 0; latitude_step <= 96; ++latitude_step) {
    const float y = -1.0f + 2.0f * latitude_step / 96.0f;
    const float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
    for (int longitude_step = 0; longitude_step < 96; ++longitude_step) {
      const float angle =
          math::PI_F * (2 * longitude_step + 1) / 96.0f - math::PI_F;
      check_peirce_fast_square_matches_exact(
          math::Vector(radius * cosf(angle), y, radius * sinf(angle)));
    }
  }
}

/** @brief Rounded pole heights must preserve the exact kernel's edge identity. */
inline void test_peirce_fast_square_rounded_pole_cap() {
  for (float sign : {-1.0f, 1.0f})
    for (float radius : {1e-6f, 1e-5f, 1e-4f, 2e-4f, 3e-4f, 1e-3f})
      for (int longitude = 0; longitude < 64; ++longitude) {
        const float angle = (longitude + 0.5f) * math::TWO_PI_F / 64.0f;
        const math::Vector v(radius * cosf(angle),
                             sign * sqrtf(1.0f - radius * radius),
                             radius * sinf(angle));
        const auto exact =
            peirce_projection(v, 0.0f, projections::PeirceLayout::SQUARE, 0.0f);
        const auto fast = peirce_projection_fast_square(v);
        HS_EXPECT_EQ(fast.edge_class, exact.edge_class);
      }
}

/** @brief The fast Peirce kernel matches exact values and metadata at seams and poles. */
inline void test_peirce_fast_square_on_seams_and_poles() {
  check_peirce_fast_square_matches_exact(math::Vector(0.0f, 1.0f, 0.0f));
  check_peirce_fast_square_matches_exact(math::Vector(0.0f, -1.0f, 0.0f));
  for (int latitude_step = -7; latitude_step <= 7; ++latitude_step) {
    const float y = latitude_step / 8.0f;
    const float radius = sqrtf(1.0f - y * y);
    const float diagonal = radius * INV_SQRT_TWO;
    check_peirce_fast_square_matches_exact(math::Vector(diagonal, y, diagonal));
    check_peirce_fast_square_matches_exact(
        math::Vector(-diagonal, y, diagonal));
    check_peirce_fast_square_matches_exact(
        math::Vector(-diagonal, y, -diagonal));
    check_peirce_fast_square_matches_exact(
        math::Vector(diagonal, y, -diagonal));
    check_peirce_fast_square_matches_exact(math::Vector(radius, y, 0.0f));
    check_peirce_fast_square_matches_exact(math::Vector(0.0f, y, radius));
    check_peirce_fast_square_matches_exact(math::Vector(-radius, y, 0.0f));
    check_peirce_fast_square_matches_exact(math::Vector(0.0f, y, -radius));
  }
}

/** @brief Diagonal tie-band samples retain their seam region and edge identity. */
inline void test_peirce_fast_square_ties_the_diagonal_band_to_its_seam() {
  // Inside peirce_sector_longitude's snap band the exact kernel folds the
  // azimuth onto the boundary; diagonal_tie labels those samples as the seam.
  constexpr float OFFSETS[] = {-1.0e-6f, -1.0e-7f, 1.0e-7f, 1.0e-6f};
  constexpr int QUADRANT_X[] = {1, -1, -1, 1};
  constexpr int QUADRANT_Z[] = {1, 1, -1, -1};
  for (int latitude_step = -6; latitude_step <= 6; ++latitude_step) {
    const float y = latitude_step / 8.0f;
    const float diagonal = sqrtf(1.0f - y * y) * INV_SQRT_TWO;
    for (int quadrant = 0; quadrant < 4; ++quadrant) {
      const float sx = static_cast<float>(QUADRANT_X[quadrant]);
      const float sz = static_cast<float>(QUADRANT_Z[quadrant]);
      const ProjectionKernelResult seam = peirce_projection_fast_square(
          math::Vector(sx * diagonal, y, sz * diagonal));
      for (float offset : OFFSETS) {
        const ProjectionKernelResult banded = peirce_projection_fast_square(
            math::Vector(sx * diagonal * (1.0f + offset), y,
                         sz * diagonal * (1.0f - offset)));
        HS_EXPECT_EQ(banded.region_id, seam.region_id);
        HS_EXPECT_EQ(banded.edge_class, seam.edge_class);
      }
    }
  }
}

/** @brief Peirce edge distance vanishes at the map's four singularities. */
inline void test_peirce_edge_distance_locates_the_singularities() {
  // acos saturates at the singularity, so a float direction there lands one
  // sqrt-of-epsilon step out rather than exactly on zero.
  constexpr float SATURATION = 1e-3f;
  for (float x : {INV_SQRT_TWO, -INV_SQRT_TWO})
    for (float z : {INV_SQRT_TWO, -INV_SQRT_TWO}) {
      const math::Vector singular(x, 0.0f, z);
      for (PeirceLayout layout :
           {PeirceLayout::DIAMOND, PeirceLayout::SQUARE,
            PeirceLayout::HORIZONTAL, PeirceLayout::VERTICAL})
        HS_EXPECT_NEAR(
            peirce_projection(singular, 0.0f, layout, 0.0f).fade_edge_distance,
            0.0f, SATURATION);
      HS_EXPECT_NEAR(peirce_projection_fast_square(singular).fade_edge_distance,
                     0.0f, SATURATION);
    }
  // The four regular equatorial points sit a half-quadrant from the nearest
  // singularity, and read the same distance on both sides of the equator.
  for (int quadrant = 0; quadrant < 4; ++quadrant) {
    const float longitude = quadrant * (0.5f * math::PI_F);
    for (float latitude : {-1e-4f, 1e-4f})
      HS_EXPECT_NEAR(peirce_projection(direction(latitude, longitude), 0.0f,
                                       projections::PeirceLayout::SQUARE, 0.0f)
                         .fade_edge_distance,
                     0.25f * math::PI_F, 1e-4f);
  }
  // Distance rises monotonically as the equator walks away from the singularity.
  float previous = -1.0f;
  for (int step = 0; step <= 16; ++step) {
    const float longitude = 0.25f * math::PI_F * (1.0f - step / 16.0f);
    const float distance =
        peirce_projection(direction(0.0f, longitude), 0.0f,
                          projections::PeirceLayout::SQUARE, 0.0f)
            .fade_edge_distance;
    HS_EXPECT_GT(distance, previous);
    previous = distance;
  }
}

/** @brief Square coordinates are the diamond rotated by pi/4 with the same region. */
inline void test_peirce_square_is_the_rotated_diamond() {
  for (int latitude_step = -8; latitude_step <= 8; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * math::PI_F / 8.0f);
    for (int longitude_step = 0; longitude_step < 12; ++longitude_step) {
      const math::Vector v =
          direction(latitude, longitude_step * (math::TWO_PI_F / 12.0f) -
                                  math::PI_F + 0.2f);
      const ProjectionKernelResult diamond =
          peirce_projection(v, 0.0f, PeirceLayout::DIAMOND, 0.0f);
      const ProjectionKernelResult square =
          peirce_projection(v, 0.0f, projections::PeirceLayout::SQUARE, 0.0f);
      HS_EXPECT_NEAR(square.coords.re,
                     INV_SQRT_TWO * (diamond.coords.re - diamond.coords.im),
                     1e-6f);
      HS_EXPECT_NEAR(square.coords.im,
                     INV_SQRT_TWO * (diamond.coords.re + diamond.coords.im),
                     1e-6f);
      HS_EXPECT_EQ(square.region_id, diamond.region_id);
    }
  }
}

/** @brief A full scroll turn preserves strip coordinates and periodic edge
 *  metadata; each strip layout reports one edge class of its own, outside the
 *  square net's four edge classes. */
inline void test_peirce_strip_scroll_is_periodic() {
  constexpr uint8_t SQUARE_EDGE_CLASSES = 4;
  int strip_class[2] = {-1, -1};
  for (PeirceLayout layout : {PeirceLayout::HORIZONTAL, PeirceLayout::VERTICAL})
    for (int latitude_step = -6; latitude_step <= 6; ++latitude_step) {
      const float latitude = latitude_step * (0.5f * math::PI_F / 6.0f);
      for (int longitude_step = 0; longitude_step < 8; ++longitude_step) {
        const math::Vector v =
            direction(latitude, longitude_step * (math::TWO_PI_F / 8.0f) -
                                    math::PI_F + 0.15f);
        const ProjectionKernelResult base =
            peirce_projection(v, 0.0f, layout, 0.125f);
        const ProjectionKernelResult wrapped =
            peirce_projection(v, 0.0f, layout, 1.125f);
        HS_EXPECT_NEAR(wrapped.coords.re, base.coords.re, 1e-4f);
        HS_EXPECT_NEAR(wrapped.coords.im, base.coords.im, 1e-4f);
        HS_EXPECT_EQ(wrapped.region_id, base.region_id);
        HS_EXPECT_EQ(wrapped.component_id, base.component_id);
        HS_EXPECT_EQ(wrapped.boundary_flags, base.boundary_flags);
        HS_EXPECT_NEAR(wrapped.fade_edge_distance, base.fade_edge_distance,
                       1e-6f);
        HS_EXPECT_EQ(wrapped.flags, base.flags);
        HS_EXPECT_EQ(wrapped.traits, base.traits);
        HS_EXPECT_EQ(wrapped.edge_class, base.edge_class);
        HS_EXPECT_GE(base.edge_class, SQUARE_EDGE_CLASSES);
        int &expected = strip_class[layout == PeirceLayout::VERTICAL];
        if (expected < 0)
          expected = base.edge_class;
        HS_EXPECT_EQ(int(base.edge_class), expected);
      }
    }
  HS_EXPECT_NE(strip_class[0], strip_class[1]);
}

/** @brief A Peirce strip layout reports the equator tear its reflection leaves
 *         open, and keeps the opposite quarters glued modulo the strip period. */
inline void test_peirce_strip_tears_the_unglued_equator() {
  constexpr float LATITUDE = 1e-4f;
  for (PeirceLayout layout :
       {PeirceLayout::HORIZONTAL, PeirceLayout::VERTICAL}) {
    const float torn =
        layout == PeirceLayout::HORIZONTAL ? 0.0f : 0.5f * math::PI_F;
    const float glued =
        layout == PeirceLayout::HORIZONTAL ? 0.5f * math::PI_F : 0.0f;
    const ProjectionKernelResult north =
        peirce_projection(direction(LATITUDE, torn), 0.0f, layout, 0.0f);
    const ProjectionKernelResult south =
        peirce_projection(direction(-LATITUDE, torn), 0.0f, layout, 0.0f);
    HS_EXPECT_EQ(north.traits & projection_traits(ProjectionTrait::CUT),
                 projection_traits(ProjectionTrait::CUT));
    const float jump = std::max(std::fabs(north.coords.re - south.coords.re),
                                std::fabs(north.coords.im - south.coords.im));
    HS_EXPECT_GT(jump, 2.0f * PEIRCE_QUARTER_PERIOD - 1e-3f);
    HS_EXPECT_EQ(north.boundary_flags,
                 projection_boundary(ProjectionBoundary::CUT));
    HS_EXPECT_NEAR(north.fade_edge_distance, LATITUDE, 1e-6f);
    HS_EXPECT_NEAR(south.fade_edge_distance, LATITUDE, 1e-6f);
    const ProjectionKernelResult seam_north =
        peirce_projection(direction(LATITUDE, glued), 0.0f, layout, 0.0f);
    const ProjectionKernelResult seam_south =
        peirce_projection(direction(-LATITUDE, glued), 0.0f, layout, 0.0f);
    float delta_x = seam_north.coords.re - seam_south.coords.re;
    float delta_y = seam_north.coords.im - seam_south.coords.im;
    if (layout == PeirceLayout::HORIZONTAL)
      delta_x = std::remainder(delta_x, 4.0f * PEIRCE_QUARTER_PERIOD);
    else
      delta_y = std::remainder(delta_y, 4.0f * PEIRCE_QUARTER_PERIOD);
    HS_EXPECT_NEAR(delta_x, 0.0f, 1e-3f);
    HS_EXPECT_NEAR(delta_y, 0.0f, 1e-3f);
    HS_EXPECT_EQ(seam_north.boundary_flags,
                 projection_boundary(ProjectionBoundary::SINGULAR));
  }
  for (PeirceLayout layout : {PeirceLayout::DIAMOND, PeirceLayout::SQUARE})
    HS_EXPECT_EQ(
        peirce_projection(direction(LATITUDE, 0.0f), 0.0f, layout, 0.0f)
                .traits &
            projection_traits(ProjectionTrait::CUT),
        0);
}

/** @brief Paired half-edges share the minimum canonical edge identity. */
inline void test_airocean_edge_identity_is_the_canonical_half_edge() {
  for (size_t face = 0; face < AIROCEAN_FACE_COUNT; ++face)
    for (size_t edge = 0; edge < 3; ++edge) {
      size_t canonical = face * 3 + edge;
      for (size_t other = 0; other < AIROCEAN_FACE_COUNT; ++other)
        for (size_t other_edge = 0; other_edge < 3; ++other_edge) {
          if (other * 3 + other_edge >= canonical)
            continue;
          if (planar_segments_equal(planar_vertex(face, edge),
                                    planar_vertex(face, edge + 1),
                                    planar_vertex(other, other_edge),
                                    planar_vertex(other, other_edge + 1)))
            canonical = other * 3 + other_edge;
        }
      HS_EXPECT_EQ(size_t(airocean_edge_identity(uint8_t(face), uint8_t(edge))),
                   canonical);
    }
}

/** @brief Glued planar edge endpoints are bit-identical and share an edge identity. */
inline void test_airocean_glued_edges_are_bit_identical() {
  const auto point_bits_equal = [](const AiroceanPoint &a,
                                   const AiroceanPoint &b) {
    return std::bit_cast<uint32_t>(a.x) == std::bit_cast<uint32_t>(b.x) &&
           std::bit_cast<uint32_t>(a.y) == std::bit_cast<uint32_t>(b.y);
  };
  const auto edge_bits_equal =
      [&](const AiroceanPoint &a, const AiroceanPoint &b,
          const AiroceanPoint &c, const AiroceanPoint &d) {
        return (point_bits_equal(a, c) && point_bits_equal(b, d)) ||
               (point_bits_equal(a, d) && point_bits_equal(b, c));
      };

  int glued = 0;
  int shared = 0;
  for (size_t face = 0; face < AIROCEAN_FACE_COUNT; ++face)
    for (size_t edge = 0; edge < 3; ++edge)
      for (size_t other = face + 1; other < AIROCEAN_FACE_COUNT; ++other)
        for (size_t other_edge = 0; other_edge < 3; ++other_edge) {
          if (!spherical_segments_equal(
                  AIROCEAN_FACES[face][edge],
                  AIROCEAN_FACES[face][(edge + 1) % 3],
                  AIROCEAN_FACES[other][other_edge],
                  AIROCEAN_FACES[other][(other_edge + 1) % 3]))
            continue;
          ++shared;
          const bool cut = airocean_edge_is_cut(uint8_t(face), uint8_t(edge));
          HS_EXPECT_EQ(
              cut, airocean_edge_is_cut(uint8_t(other), uint8_t(other_edge)));
          if (cut)
            continue;
          ++glued;
          HS_EXPECT_TRUE(edge_bits_equal(planar_vertex(face, edge),
                                         planar_vertex(face, edge + 1),
                                         planar_vertex(other, other_edge),
                                         planar_vertex(other, other_edge + 1)));
          HS_EXPECT_EQ(
              airocean_edge_identity(uint8_t(face), uint8_t(edge)),
              airocean_edge_identity(uint8_t(other), uint8_t(other_edge)));
        }
  HS_EXPECT_EQ(shared, 33);
  HS_EXPECT_EQ(glued, 21);
}

/** @brief Face vertices share a positive plane offset and their centroid equals the face center. */
inline void test_airocean_face_planes() {
  for (size_t face = 0; face < AIROCEAN_FACE_COUNT; ++face) {
    const AiroceanVector &normal = AIROCEAN_NORMALS[face];
    HS_EXPECT_NEAR(
        sqrtf(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z),
        1.0f, 1e-6f);
    const AiroceanVector &center = AIROCEAN_CENTERS[face];
    const float offset =
        center.x * normal.x + center.y * normal.y + center.z * normal.z;
    HS_EXPECT_GT(offset, 0.0f);
    AiroceanVector centroid{0.0f, 0.0f, 0.0f};
    for (size_t vertex = 0; vertex < 3; ++vertex) {
      const AiroceanVector &v = AIROCEAN_FACES[face][vertex];
      HS_EXPECT_NEAR(v.x * normal.x + v.y * normal.y + v.z * normal.z, offset,
                     1e-6f);
      HS_EXPECT_TRUE(airocean_contains(v, face));
      centroid.x += v.x / 3.0f;
      centroid.y += v.y / 3.0f;
      centroid.z += v.z / 3.0f;
    }
    HS_EXPECT_NEAR(centroid.x, center.x, 1e-6f);
    HS_EXPECT_NEAR(centroid.y, center.y, 1e-6f);
    HS_EXPECT_NEAR(centroid.z, center.z, 1e-6f);
  }
}

/** @brief Unfolding maps each spherical vertex to its clockwise planar face vertex. */
inline void test_airocean_unfold_sends_vertices_to_planar_vertices() {
  for (size_t face = 0; face < AIROCEAN_FACE_COUNT; ++face) {
    HS_EXPECT_LT(planar_cross(AIROCEAN_PLANAR_FACES[face][0],
                              AIROCEAN_PLANAR_FACES[face][1],
                              AIROCEAN_PLANAR_FACES[face][2]),
                 0.0f);
    for (size_t vertex = 0; vertex < 3; ++vertex) {
      const AiroceanPoint mapped =
          airocean_unfold(face, AIROCEAN_FACES[face][vertex]);
      HS_EXPECT_NEAR(mapped.x, AIROCEAN_PLANAR_FACES[face][vertex].x, 1e-5f);
      HS_EXPECT_NEAR(mapped.y, AIROCEAN_PLANAR_FACES[face][vertex].y, 1e-5f);
    }
  }
}

/** @brief Samples spherical directions and their selected Airocean faces. */
template <typename Check> inline void for_each_airocean_sample(Check check) {
  for (int latitude_step = 0; latitude_step <= 48; ++latitude_step) {
    const float y = -1.0f + 2.0f * latitude_step / 48.0f;
    const float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
    for (int longitude_step = 0; longitude_step < 64; ++longitude_step) {
      const float angle =
          math::TWO_PI_F * (longitude_step + 0.31f) / 64.0f - math::PI_F;
      const math::Vector v(radius * cosf(angle), y, radius * sinf(angle));
      const ProjectionKernelResult net =
          airocean_projection_meridian(v, 0.0f, false);
      check(v, net);
    }
  }
}

/** @brief Projected points stay in their selected faces with consistent edge metadata. */
inline void test_airocean_projection_stays_inside_its_face() {
  for_each_airocean_sample(
      [](const math::Vector &v, const ProjectionKernelResult &net) {
        const size_t face = net.region_id;
        HS_EXPECT_LT(face, AIROCEAN_FACE_COUNT);
        if (face >= AIROCEAN_FACE_COUNT)
          return;
        HS_EXPECT_TRUE(airocean_contains(airocean_axes(v), face));
        HS_EXPECT_EQ(airocean_outside_score(airocean_axes(v), face), 0.0f);
        const auto expected = airocean_unfold(face, airocean_axes(v));
        HS_EXPECT_NEAR(net.coords.re, expected.x, 1e-5f);
        HS_EXPECT_NEAR(net.coords.im, expected.y, 1e-5f);
        if (AIROCEAN_CUT_MASKS[face]) {
          float cut_distance = INFINITY;
          for (uint8_t edge = 0; edge < 3; ++edge)
            if (AIROCEAN_CUT_MASKS[face] & (1u << edge)) {
              const auto a = planar_vertex(face, edge);
              const auto b = face == 14 && edge == 0
                                 ? AIROCEAN_PLANAR_FACES[18][1]
                                 : planar_vertex(face, edge + 1);
              cut_distance = std::min(
                  cut_distance,
                  point_segment_distance({net.coords.re, net.coords.im}, a, b));
            }
          HS_EXPECT_NEAR(net.fade_edge_distance, cut_distance, 1e-5f);
        }
        const AiroceanPoint point{net.coords.re, net.coords.im};
        for (size_t edge = 0; edge < 3; ++edge)
          HS_EXPECT_LE(planar_cross(planar_vertex(face, edge),
                                    planar_vertex(face, edge + 1), point),
                       1e-5f);
        HS_EXPECT_TRUE(
            net.edge_class == airocean_edge_identity(uint8_t(face), 0) ||
            net.edge_class == airocean_edge_identity(uint8_t(face), 1) ||
            net.edge_class == airocean_edge_identity(uint8_t(face), 2) ||
            (face == 14 && net.edge_class == airocean_edge_identity(18, 0)));
        HS_EXPECT_EQ(net.boundary_flags,
                     AIROCEAN_CUT_MASKS[face] != 0
                         ? projection_boundary(ProjectionBoundary::CUT)
                         : projection_boundary(ProjectionBoundary::NONE));
      });
}

/** @brief Horizontal Airocean coordinates are a quarter-turn of the vertical net. */
inline void test_airocean_horizontal_is_a_quarter_turn() {
  for_each_airocean_sample(
      [](const math::Vector &v, const ProjectionKernelResult &net) {
        const auto horizontal = airocean_projection_meridian(v, 0.0f, true);
        HS_EXPECT_EQ(horizontal.region_id, net.region_id);
        HS_EXPECT_NEAR(horizontal.coords.re,
                       AIROCEAN_NET_HEIGHT - net.coords.im, 1e-5f);
        HS_EXPECT_NEAR(horizontal.coords.im, net.coords.re, 1e-5f);
      });
}

/** @brief Rotating the direction and central meridian preserves the Airocean net. */
inline void test_airocean_meridian_invariance() {
  for_each_airocean_sample([](const math::Vector &v,
                              const ProjectionKernelResult &net) {
    constexpr float MERIDIAN = .37f;
    const math::Vector rotated(v.x * cosf(MERIDIAN) - v.z * sinf(MERIDIAN), v.y,
                               v.x * sinf(MERIDIAN) + v.z * cosf(MERIDIAN));
    const auto shifted = airocean_projection_meridian(rotated, MERIDIAN, false);
    HS_EXPECT_EQ(shifted.region_id, net.region_id);
    HS_EXPECT_NEAR(shifted.coords.re, net.coords.re, 1e-5f);
    HS_EXPECT_NEAR(shifted.coords.im, net.coords.im, 1e-5f);
  });
}

/** @brief Face 14 selects cut or glued traits and edge identity on opposite edge halves. */
inline void test_airocean_face14_half_cut_edge() {
  const auto &a = AIROCEAN_FACES[14][0];
  const auto &b = AIROCEAN_FACES[14][1];
  const auto &c = AIROCEAN_FACES[14][2];
  for (const float t : {0.25f, 0.75f}) {
    constexpr float INSET = 1e-4f;
    const math::Vector point(
        (1.0f - INSET) * ((1.0f - t) * a.x + t * b.x) + INSET * c.x,
        (1.0f - INSET) * ((1.0f - t) * a.z + t * b.z) + INSET * c.z,
        (1.0f - INSET) * ((1.0f - t) * a.y + t * b.y) + INSET * c.y);
    const auto net =
        airocean_projection_meridian(point.normalized(), 0.0f, false);
    HS_EXPECT_EQ(net.region_id, 14);
    HS_EXPECT_EQ(net.traits,
                 projection_traits(t < 0.5f ? ProjectionTrait::CUT
                                            : ProjectionTrait::GLUED));
    HS_EXPECT_EQ(net.edge_class, airocean_edge_identity(t < 0.5f ? 14 : 18, 0));
  }
}

/** @brief Face selection remains bounded for huge, infinite and NaN directions. */
inline void test_airocean_projection_face_index_stays_in_range() {
  // Far off the unit sphere, every face scores above the fallback's sentinel.
  const float MAGNITUDES[] = {1.0f, 1e18f, 1e30f, 3.0e38f};
  for (float magnitude : MAGNITUDES)
    for (size_t face = 0; face < AIROCEAN_FACE_COUNT; ++face)
      for (size_t edge = 0; edge < 3; ++edge)
        for (int step = 0; step <= 8; ++step) {
          const AiroceanVector &a = AIROCEAN_FACES[face][edge];
          const AiroceanVector &b = AIROCEAN_FACES[face][(edge + 1) % 3];
          const float t = step / 8.0f;
          const AiroceanVector seam{a.x + t * (b.x - a.x),
                                    a.y + t * (b.y - a.y),
                                    a.z + t * (b.z - a.z)};
          const float scale =
              magnitude /
              sqrtf(seam.x * seam.x + seam.y * seam.y + seam.z * seam.z);
          // Inverse of airocean_axes: the kernel reads y-up, the faces z-up.
          const math::Vector v(seam.x * scale, seam.z * scale, seam.y * scale);
          HS_EXPECT_LT(
              size_t(airocean_projection_meridian(v, 0.0f, false).region_id),
              AIROCEAN_FACE_COUNT);
        }

  const float infinity = std::numeric_limits<float>::infinity();
  const float not_a_number = std::numeric_limits<float>::quiet_NaN();
  const float exotic[] = {infinity, -infinity, not_a_number, 0.0f, 1.0f};
  for (float x : exotic)
    for (float y : exotic)
      for (float z : exotic)
        HS_EXPECT_LT(size_t(airocean_projection_meridian(math::Vector(x, y, z),
                                                         0.4f, false)
                                .region_id),
                     AIROCEAN_FACE_COUNT);
}

/** @brief Segment distance clamps to endpoints and handles degenerate segments. */
inline void test_point_segment_distance() {
  const AiroceanPoint A{1.0f, 2.0f};
  const AiroceanPoint B{5.0f, 2.0f};
  HS_EXPECT_NEAR(point_segment_distance(A, A, B), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(B, A, B), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{3.0f, 5.0f}, A, B), 3.0f,
                 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{-3.0f, 2.0f}, A, B), 4.0f,
                 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{9.0f, 5.0f}, A, B), 5.0f,
                 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{1.0f, 6.0f}, A, A), 4.0f,
                 1e-6f);
  const float squared =
      point_segment_distance_squared(AiroceanPoint{3.0f, 5.0f}, A, B);
  HS_EXPECT_NEAR(squared, 9.0f, 1e-6f);
}

/** @brief Projection trait bits compose and match Bonne and Peirce kernel metadata. */
inline void test_projection_trait_packing() {
  HS_EXPECT_EQ(projection_traits(ProjectionTrait::NONE), 0);
  constexpr auto BIT = [](ProjectionTrait trait) {
    return static_cast<uint8_t>(trait);
  };
  HS_EXPECT_EQ(projection_traits(ProjectionTrait::CUT, ProjectionTrait::GLUED),
               BIT(ProjectionTrait::CUT) | BIT(ProjectionTrait::GLUED));
  HS_EXPECT_EQ(projection_traits(ProjectionTrait::CUT, ProjectionTrait::GLUED,
                                 ProjectionTrait::SINGULAR),
               BIT(ProjectionTrait::CUT) | BIT(ProjectionTrait::GLUED) |
                   BIT(ProjectionTrait::SINGULAR));
  const ProjectionKernelResult bonne =
      bonne_projection(math::Vector(1.0f, 0.0f, 0.0f), 0.0f, 0.4f);
  HS_EXPECT_EQ(bonne.traits, projection_traits(ProjectionTrait::CUT));
  const math::Vector EQUATOR(1.0f, 0.0f, 0.0f);
  const auto square =
      peirce_projection(EQUATOR, 0.0f, projections::PeirceLayout::SQUARE, 0.0f);
  const auto strip = peirce_projection(
      EQUATOR, 0.0f, projections::PeirceLayout::HORIZONTAL, 0.0f);
  const auto TRAITS =
      projection_traits(ProjectionTrait::GLUED, ProjectionTrait::FOLDED,
                        ProjectionTrait::PERIODIC, ProjectionTrait::SINGULAR);
  HS_EXPECT_EQ(square.traits, TRAITS);
  HS_EXPECT_EQ(strip.traits, TRAITS | projection_traits(ProjectionTrait::CUT));
}

/** @brief Equirectangular and folded sinusoidal coordinates match longitude and latitude formulas. */
inline void test_cylindrical_coordinates() {
  for (float latitude : {-0.7f, 0.0f, 0.7f}) {
    for (float longitude : {-2.0f, -0.5f, 0.5f, 2.0f}) {
      const math::Vector v = direction(latitude, longitude);
      const auto rectangular = equirectangular(v, 0.2f);
      const auto folded = folded_sinusoidal(v, 0.2f);
      HS_EXPECT_NEAR(rectangular.re, longitude - 0.2f, 0.005f);
      HS_EXPECT_NEAR(rectangular.im, latitude, 0.005f);
      HS_EXPECT_NEAR(folded.re, fabsf(longitude - 0.2f) * cosf(latitude),
                     0.005f);
      HS_EXPECT_NEAR(folded.im, latitude, 0.005f);
    }
  }
}

/** @brief Signed distance from an atlas point to a tile's nearest side.
 *  @return Positive inside the tile. */
template <size_t F, size_t K, size_t P>
inline float atlas_tile_distance(const PolyhedralAtlas<F, K, P> &atlas,
                                 size_t face, const math::Complex &point) {
  const math::Complex local = point - atlas.centers[face];
  float distance = NO_EDGE_DISTANCE;
  for (size_t k = 0; k < atlas.side_counts[face]; ++k) {
    const AtlasSide &line = atlas.sides[face][k];
    distance =
        std::min(distance, line.x * local.re + line.y * local.im + line.offset);
  }
  return distance;
}

template <typename Solid> inline void expect_polyhedral_atlas_sampling() {
  auto *volatile build_atlas = &make_polyhedral_atlas<Solid>;
  const auto atlas = build_atlas();
  HS_EXPECT_NEAR(atlas.radius, POLYHEDRAL_ATLAS<Solid>.radius, 2e-6f);
  constexpr size_t F = Solid::face_counts.size();
  constexpr int SAMPLES = 2048;
  for (int index = 0; index < SAMPLES; ++index) {
    const math::Vector v = math::fib_spiral(SAMPLES, 0.5f, index);
    const auto actual = polyhedral_atlas_projection(atlas, v, 1.0f, 0.0f, 0.0f);
    size_t best = 0;
    for (size_t f = 1; f < F; ++f)
      if (math::dot(v, atlas.selectors[f]) >
          math::dot(v, atlas.selectors[best]))
        best = f;
    HS_EXPECT_EQ(actual.region_id, best);
    const float inside =
        atlas_tile_distance(atlas, actual.region_id, actual.coords);
    HS_EXPECT_GE(inside, -1e-5f);
    HS_EXPECT_NEAR(actual.fade_edge_distance, std::max(0.0f, inside), 1e-6f);
    for (size_t f = 0; f < F; ++f)
      if (f != actual.region_id)
        HS_EXPECT_LT(atlas_tile_distance(atlas, f, actual.coords), 1e-5f);
    HS_EXPECT_LE(std::hypot(actual.coords.re, actual.coords.im),
                 atlas.radius * (1.0f + 1e-5f));
    HS_EXPECT_EQ(actual.edge_class / atlas.sides[actual.region_id].size(),
                 actual.region_id);
    HS_EXPECT_EQ(actual.traits, projection_traits(ProjectionTrait::CUT));
    HS_EXPECT_EQ(actual.boundary_flags,
                 projection_boundary(ProjectionBoundary::CUT));
    const float c = cosf(0.53f), sn = sinf(0.53f);
    const math::Vector rotated(v.x * c - v.z * sn, v.y, v.z * c + v.x * sn);
    const auto turned =
        polyhedral_atlas_projection(atlas, rotated, c, sn, 0.0f);
    HS_EXPECT_EQ(turned.region_id, actual.region_id);
    HS_EXPECT_NEAR(turned.coords.re, actual.coords.re, 2e-5f);
    HS_EXPECT_NEAR(turned.coords.im, actual.coords.im, 2e-5f);
  }
}

template <typename Solid> inline void expect_polyhedral_atlas_correction() {
  const auto &atlas = POLYHEDRAL_ATLAS<Solid>;
  constexpr size_t F = Solid::face_counts.size();
  constexpr int SAMPLES = 1024;
  constexpr float STEP = 2e-3f;
  for (int index = 0; index < SAMPLES; ++index) {
    const math::Vector v = math::fib_spiral(SAMPLES, 0.5f, index);
    const auto base = polyhedral_atlas_projection(atlas, v, 1.0f, 0.0f, 0.0f);
    const size_t face = base.region_id;
    const math::Vector normal = atlas.selectors[face].normalized();
    const float theta =
        atan2f(math::cross(v, normal).length(), math::dot(v, normal));
    const float lambert =
        theta < 1e-4f ? 1.0f : 2.0f * sinf(0.5f * theta) / tanf(theta);
    const math::Complex center = atlas.centers[face];
    for (float correction : {0.35f, 1.0f}) {
      const auto actual =
          polyhedral_atlas_projection(atlas, v, 1.0f, 0.0f, correction);
      const float shrink = 1.0f + correction * (lambert - 1.0f);
      HS_EXPECT_EQ(actual.region_id, face);
      HS_EXPECT_NEAR(actual.coords.re - center.re,
                     (base.coords.re - center.re) * shrink, 2e-5f);
      HS_EXPECT_NEAR(actual.coords.im - center.im,
                     (base.coords.im - center.im) * shrink, 2e-5f);
      HS_EXPECT_NEAR(actual.fade_edge_distance,
                     base.fade_edge_distance * shrink, 2e-5f);
      HS_EXPECT_GE(atlas_tile_distance(atlas, face, actual.coords), -1e-5f);
      for (size_t f = 0; f < F; ++f)
        if (f != face)
          HS_EXPECT_LT(atlas_tile_distance(atlas, f, actual.coords), 1e-5f);
    }
    const math::Vector east =
        math::cross(math::Vector(0.0f, 1.0f, 0.0f), v).length() > 0.1f
            ? math::cross(math::Vector(0.0f, 1.0f, 0.0f), v).normalized()
            : math::cross(math::Vector(1.0f, 0.0f, 0.0f), v).normalized();
    const math::Vector north = math::cross(v, east);
    const math::Vector a = v, b = (v + east * STEP).normalized(),
                       c = (v + north * STEP).normalized();
    const auto pa = polyhedral_atlas_projection(atlas, a, 1.0f, 0.0f, 1.0f);
    const auto pb = polyhedral_atlas_projection(atlas, b, 1.0f, 0.0f, 1.0f);
    const auto pc = polyhedral_atlas_projection(atlas, c, 1.0f, 0.0f, 1.0f);
    if (pb.region_id != face || pc.region_id != face)
      continue;
    const float planar =
        fabsf((pb.coords.re - pa.coords.re) * (pc.coords.im - pa.coords.im) -
              (pb.coords.im - pa.coords.im) * (pc.coords.re - pa.coords.re));
    const float spherical = math::cross(b - a, c - a).length();
    const float height = atlas.heights[face] * atlas.x_axes[face].length();
    HS_EXPECT_NEAR(planar / spherical / (height * height), 1.0f, 1e-2f);
  }
}

/** @brief Area correction pulls each point toward its tile center by the
 *  equal-area azimuthal ratio and makes full correction area-preserving. */
inline void test_polyhedral_atlas_area_correction() {
  expect_polyhedral_atlas_correction<Solids::Tetrahedron>();
  expect_polyhedral_atlas_correction<Solids::Cube>();
  expect_polyhedral_atlas_correction<Solids::Octahedron>();
  expect_polyhedral_atlas_correction<Solids::Dodecahedron>();
  expect_polyhedral_atlas_correction<Solids::Icosahedron>();
}

/** @brief Every direction lands inside the tile of the face its ray crosses,
 *  outside all other tiles, with its distance to that tile's nearest side. */
inline void test_polyhedral_atlas_sampling() {
  expect_polyhedral_atlas_sampling<Solids::Tetrahedron>();
  expect_polyhedral_atlas_sampling<Solids::Cube>();
  expect_polyhedral_atlas_sampling<Solids::Octahedron>();
  expect_polyhedral_atlas_sampling<Solids::Dodecahedron>();
  expect_polyhedral_atlas_sampling<Solids::Icosahedron>();
}

/** @brief A truncated tetrahedron has unequal-height antipodal face planes. */
struct TruncatedTetrahedronAtlasFixture {
  static constexpr auto vertices = [] {
    std::array<math::Vector, 12> result{};
    size_t index = 0;
    for (size_t from = 0; from < 4; ++from)
      for (size_t to = 0; to < 4; ++to)
        if (from != to)
          result[index++] = polyhedral_atlas_detail::unit(
              Solids::Tetrahedron::vertices[from] * 2.0f +
              Solids::Tetrahedron::vertices[to]);
    return result;
  }();
  static constexpr std::array<uint8_t, 8> face_counts{3, 3, 3, 3, 6, 6, 6, 6};
  static constexpr std::array<int, 36> faces{
      1, 2, 0, 5,  4, 3, 7, 8, 6, 11, 10, 9, 2, 9,  10, 5, 3, 0,
      1, 6, 8, 11, 9, 2, 0, 3, 4, 7,  6,  1, 5, 10, 11, 8, 7, 4};
};

/** @brief Unequal antipodal face heights select the full face scan. */
inline void test_polyhedral_atlas_unequal_antipodal_heights() {
  static_assert(polyhedral_atlas_detail::antipodal_pairs<
                    TruncatedTetrahedronAtlasFixture>() == 0);
  HS_EXPECT_EQ(
      POLYHEDRAL_ATLAS<TruncatedTetrahedronAtlasFixture>.antipodes.size(), 0u);
  expect_polyhedral_atlas_sampling<TruncatedTetrahedronAtlasFixture>();
}

template <typename Solid> inline void expect_polyhedral_atlas_tiles() {
  const auto &atlas = POLYHEDRAL_ATLAS<Solid>;
  constexpr size_t F = Solid::face_counts.size();
  constexpr size_t COLUMNS = (F + 1) / 2;
  const size_t sides = Solid::face_counts[0];
  const float inradius = 0.5f / tanf(math::PI_F / static_cast<float>(sides));
  size_t north = F;
  float upper_south = 1.0f, lower_north = -1.0f;
  for (size_t f = 0; f < F; ++f) {
    const math::Vector normal = atlas.selectors[f].normalized();
    const auto center =
        polyhedral_atlas_projection(atlas, normal, 1.0f, 0.0f, 0.0f);
    HS_EXPECT_EQ(center.region_id, f);
    HS_EXPECT_NEAR(center.coords.re, atlas.centers[f].re, 1e-5f);
    HS_EXPECT_NEAR(center.coords.im, atlas.centers[f].im, 1e-5f);
    HS_EXPECT_NEAR(center.fade_edge_distance, inradius, 1e-5f);
    float &row_extreme = atlas.centers[f].im > 0.0f ? upper_south : lower_north;
    row_extreme = atlas.centers[f].im > 0.0f ? std::min(row_extreme, normal.y)
                                             : std::max(row_extreme, normal.y);
    if (normal.y > 0.999f)
      north = f;
    if (fabsf(normal.y) > 0.99f)
      continue;
    const math::Vector up =
        (normal + math::Vector(0.0f, 1.0f, 0.0f) * 0.05f).normalized();
    const auto above = polyhedral_atlas_projection(atlas, up, 1.0f, 0.0f, 0.0f);
    HS_EXPECT_EQ(above.region_id, f);
    HS_EXPECT_NEAR(above.coords.re, center.coords.re, 1e-5f);
    HS_EXPECT_GT(above.coords.im, center.coords.im);
  }
  HS_EXPECT_GE(upper_south, lower_north - 1e-3f);
  HS_EXPECT_LT(north, F);
  if (north < F) {
    const float cell = atlas.centers[north].im * 2.0f;
    HS_EXPECT_NEAR(atlas.centers[north].re,
                   -0.5f * cell * static_cast<float>(COLUMNS - 1), 1e-5f);
  }
  for (size_t f = 0; f < F; ++f) {
    const math::Vector normal = atlas.selectors[f].normalized();
    HS_CONTEXT("face", f);
    const math::Vector up = atlas.y_axes[f].normalized();
    const math::Vector east = atlas.x_axes[f].normalized();
    const math::Vector p0 = (normal + up * 0.15f - east * 0.2f).normalized();
    const math::Vector p1 = (normal + up * 0.15f + east * 0.2f).normalized();
    const math::Vector mid = (p0 + p1).normalized();
    const auto q0 = polyhedral_atlas_projection(atlas, p0, 1.0f, 0.0f, 0.0f);
    const auto q1 = polyhedral_atlas_projection(atlas, p1, 1.0f, 0.0f, 0.0f);
    const auto qm = polyhedral_atlas_projection(atlas, mid, 1.0f, 0.0f, 0.0f);
    HS_EXPECT_EQ(q0.region_id, f);
    HS_EXPECT_EQ(q1.region_id, f);
    HS_EXPECT_EQ(qm.region_id, f);
    const float cross =
        (q1.coords.re - q0.coords.re) * (qm.coords.im - q0.coords.im) -
        (q1.coords.im - q0.coords.im) * (qm.coords.re - q0.coords.re);
    HS_EXPECT_NEAR(cross, 0.0f, 1e-5f);
  }
}

/** @brief Atlas tiles have unit sides, sit north up in latitude-ordered rows,
 *  and keep great circles straight. */
inline void test_polyhedral_atlas_tiles() {
  expect_polyhedral_atlas_tiles<Solids::Tetrahedron>();
  expect_polyhedral_atlas_tiles<Solids::Cube>();
  expect_polyhedral_atlas_tiles<Solids::Octahedron>();
  expect_polyhedral_atlas_tiles<Solids::Dodecahedron>();
  expect_polyhedral_atlas_tiles<Solids::Icosahedron>();
}

/** @brief Skipping edge measurements preserves projection coordinates and identity. */
inline void test_projection_without_edge_distance() {
  auto compare = [](const ProjectionKernelResult &measured,
                    const ProjectionKernelResult &skipped) {
    HS_EXPECT_EQ(skipped.coords.re, measured.coords.re);
    HS_EXPECT_EQ(skipped.coords.im, measured.coords.im);
    HS_EXPECT_EQ(skipped.region_id, measured.region_id);
    HS_EXPECT_EQ(skipped.flags, measured.flags);
    HS_EXPECT_EQ(skipped.traits, measured.traits);
    HS_EXPECT_EQ(skipped.edge_class, measured.edge_class);
    HS_EXPECT_EQ(skipped.fade_edge_distance, NO_EDGE_DISTANCE);
    HS_EXPECT_EQ(skipped.boundary_flags,
                 projection_boundary(ProjectionBoundary::NONE));
  };
  for (float latitude : {-1.4f, -.7f, -.1f, .1f, .7f, 1.4f})
    for (float longitude : {-2.7f, -1.3f, -.4f, .2f, 1.1f, 2.6f}) {
      const auto v = direction(latitude, longitude);
      for (auto layout : {PeirceLayout::DIAMOND, PeirceLayout::SQUARE,
                          PeirceLayout::HORIZONTAL, PeirceLayout::VERTICAL})
        compare(peirce_projection(v, .17f, layout, .1f, true),
                peirce_projection(v, .17f, layout, .1f, false));
      for (bool horizontal : {false, true})
        compare(airocean_projection_meridian(v, .17f, horizontal, true),
                airocean_projection_meridian(v, .17f, horizontal, false));
    }
}

/** @brief Longitudes immediately outside the tie band retain their angle. */
inline void test_peirce_sector_outside_tie_band() {
  for (float boundary : {-.75f * math::PI_F, -.25f * math::PI_F,
                         .25f * math::PI_F, .75f * math::PI_F})
    for (float offset : {-5e-6f, 5e-6f})
      for (int index = -14; index <= 14; ++index) {
        const float angle = boundary + offset;
        const auto v = direction(index * .1f, angle);
        const float actual = peirce_sector_longitude(v, 0.0f);
        HS_EXPECT_NE(actual, boundary);
        HS_EXPECT_NEAR(actual, angle, 1e-6f);
        check_peirce_fast_square_matches_exact(v);
      }
}

/** @brief Folded Peirce layouts glue the equator and center the north pole. */
inline void test_peirce_folded_equator_gluing() {
  for (auto layout : {PeirceLayout::DIAMOND, PeirceLayout::SQUARE}) {
    const auto pole =
        peirce_projection(direction(.5f * math::PI_F, 0), 0, layout, 0);
    HS_EXPECT_NEAR(pole.coords.re, 0.0f, 1e-6f);
    HS_EXPECT_NEAR(pole.coords.im, 0.0f, 1e-6f);
    for (int index = 0; index < 64; ++index) {
      const float longitude =
          (index / 16) * (.5f * math::PI_F) + .25f * math::PI_F + .1f +
          (index % 16 + .5f) / 16 * (.5f * math::PI_F - .2f);
      const auto north =
          peirce_projection(direction(1e-4f, longitude), 0, layout, 0);
      const auto south =
          peirce_projection(direction(-1e-4f, longitude), 0, layout, 0);
      HS_EXPECT_NEAR(north.coords.re, south.coords.re, 2e-3f);
      HS_EXPECT_NEAR(north.coords.im, south.coords.im, 2e-3f);
    }
  }
}

inline int run_projections_tests() {
  hs_test::ModuleFixture fixture("projections");
  HS_RUN(test_projection_without_edge_distance);
  HS_RUN(test_peirce_sector_outside_tie_band);
  HS_RUN(test_peirce_folded_equator_gluing);
  HS_RUN(test_wrap_longitude_range);
  HS_RUN(test_wrap_longitude_fmod_contract);
  HS_RUN(test_cylindrical_coordinates);
  HS_RUN(test_polyhedral_atlas_sampling);
  HS_RUN(test_polyhedral_atlas_unequal_antipodal_heights);
  HS_RUN(test_polyhedral_atlas_tiles);
  HS_RUN(test_polyhedral_atlas_area_correction);
  HS_RUN(test_bonne_sinusoidal_limit);
  HS_RUN(test_bonne_general_parallel_matches_snyder);
  HS_RUN(test_bonne_general_parallel_edge_metadata);
  HS_RUN(test_bonne_general_parallel_meridian_invariance);
  HS_RUN(test_peirce_elliptic_integral_shape);
  HS_RUN(test_peirce_sector_longitude_snapping);
  HS_RUN(test_peirce_fast_square_matches_exact);
  HS_RUN(test_peirce_fast_square_on_seams_and_poles);
  HS_RUN(test_peirce_fast_square_rounded_pole_cap);
  HS_RUN(test_peirce_fast_square_ties_the_diagonal_band_to_its_seam);
  HS_RUN(test_peirce_edge_distance_locates_the_singularities);
  HS_RUN(test_peirce_square_is_the_rotated_diamond);
  HS_RUN(test_peirce_strip_scroll_is_periodic);
  HS_RUN(test_peirce_meridian_and_scroll_coordinates);
  HS_RUN(test_peirce_strip_tears_the_unglued_equator);
  HS_RUN(test_airocean_edge_identity_is_the_canonical_half_edge);
  HS_RUN(test_airocean_glued_edges_are_bit_identical);
  HS_RUN(test_airocean_face_planes);
  HS_RUN(test_airocean_unfold_sends_vertices_to_planar_vertices);
  HS_RUN(test_airocean_projection_stays_inside_its_face);
  HS_RUN(test_airocean_horizontal_is_a_quarter_turn);
  HS_RUN(test_airocean_meridian_invariance);
  HS_RUN(test_airocean_face14_half_cut_edge);
  HS_RUN(test_airocean_projection_face_index_stays_in_range);
  HS_RUN(test_point_segment_distance);
  HS_RUN(test_projection_trait_packing);
  return fixture.result();
}

} // namespace projections_tests
} // namespace hs_test
