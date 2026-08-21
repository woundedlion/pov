/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <limits>

#include "core/math/projections.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace projections_tests {

using namespace projections;

constexpr float PEIRCE_QUARTER_PERIOD = 1.8540746773013719f;
constexpr size_t AIROCEAN_FACE_COUNT = 23;

inline Vector direction(float latitude, float longitude) {
  const float radius = cosf(latitude);
  return Vector(radius * cosf(longitude), sinf(latitude),
                radius * sinf(longitude));
}

/** @brief The kernel's y-up to z-up permutation, as airocean_projection does it
 *         at a zero central meridian. */
inline AiroceanVector airocean_axes(const Vector &v) {
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

inline void test_wrap_longitude_range() {
  for (int step = -40; step <= 40; ++step) {
    const float raw = step * 0.4f;
    const float wrapped = wrap_longitude(raw);
    HS_EXPECT_GT(wrapped, -PI_F - 1e-5f);
    HS_EXPECT_LE(wrapped, PI_F + 1e-5f);
    HS_EXPECT_NEAR(cosf(wrapped), cosf(raw), 2e-6f);
    HS_EXPECT_NEAR(sinf(wrapped), sinf(raw), 2e-6f);
    HS_EXPECT_NEAR(wrap_longitude(wrapped), wrapped, 1e-6f);
  }
}

inline void test_bonne_sinusoidal_limit() {
  for (int latitude_step = -12; latitude_step <= 12; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * PI_F / 12.0f);
    for (int longitude_step = 0; longitude_step < 16; ++longitude_step) {
      const float longitude = longitude_step * (TWO_PI_F / 16.0f) - PI_F + 0.1f;
      const Vector v = direction(latitude, longitude);
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

inline void test_bonne_polar_limit_is_finite() {
  for (float standard_parallel : {0.5f * PI_F, -0.5f * PI_F}) {
    for (int step = -8; step <= 8; ++step) {
      const Vector v = direction(step * (0.5f * PI_F / 8.0f), 0.7f);
      const ProjectionKernelResult werner =
          bonne_projection(v, 0.0f, standard_parallel);
      HS_EXPECT_TRUE(std::isfinite(werner.coords.re));
      HS_EXPECT_TRUE(std::isfinite(werner.coords.im));
      HS_EXPECT_GE(werner.fade_edge_distance, 0.0f);
    }
  }
}

inline void test_peirce_elliptic_integral_shape() {
  HS_EXPECT_NEAR(peirce_elliptic_integral(0.0f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(peirce_elliptic_integral(0.5f * PI_F), PEIRCE_QUARTER_PERIOD,
                 1e-6f);
  float previous = peirce_elliptic_integral(-0.5f * PI_F);
  HS_EXPECT_NEAR(previous, -PEIRCE_QUARTER_PERIOD, 1e-6f);
  for (int step = -63; step <= 64; ++step) {
    const float phi = step * (0.5f * PI_F / 64.0f);
    const float value = peirce_elliptic_integral(phi);
    HS_EXPECT_GT(value, previous);
    HS_EXPECT_NEAR(peirce_elliptic_integral(-phi), -value, 1e-6f);
    previous = value;
  }
}

inline void test_peirce_sector_longitude_snapping() {
  HS_EXPECT_EQ(peirce_sector_longitude(Vector(0.0f, 1.0f, 0.0f), 0.0f),
               0.5f * PI_F);
  HS_EXPECT_EQ(peirce_sector_longitude(Vector(0.0f, -1.0f, 0.0f), 0.0f),
               0.5f * PI_F);
  constexpr std::array<float, 4> BOUNDARIES = {-0.75f * PI_F, -0.25f * PI_F,
                                               0.25f * PI_F, 0.75f * PI_F};
  for (float boundary : BOUNDARIES)
    for (float offset : {-1e-6f, 0.0f, 1e-6f}) {
      const float angle = boundary + offset;
      const Vector v(cosf(angle), 0.0f, sinf(angle));
      HS_EXPECT_EQ(peirce_sector_longitude(v, 0.0f), boundary);
    }
  const Vector interior(cosf(0.5f), 0.0f, sinf(0.5f));
  HS_EXPECT_NEAR(peirce_sector_longitude(interior, 0.0f), 0.5f, 1e-6f);
}

inline void check_peirce_fast_square_matches_exact(const Vector &v) {
  const ProjectionKernelResult exact = peirce_projection(v, 0.0f, 1, 0.0f);
  const ProjectionKernelResult fast = peirce_projection_fast_square(v);
  HS_EXPECT_NEAR(fast.coords.re, exact.coords.re, 1e-3f);
  HS_EXPECT_NEAR(fast.coords.im, exact.coords.im, 1e-3f);
  HS_EXPECT_NEAR(fast.fade_edge_distance, exact.fade_edge_distance, 1e-3f);
  HS_EXPECT_EQ(fast.region_id, exact.region_id);
  HS_EXPECT_EQ(fast.edge_class, exact.edge_class);
  HS_EXPECT_EQ(fast.flags, exact.flags);
  HS_EXPECT_EQ(fast.traits, exact.traits);
}

inline void test_peirce_fast_square_matches_exact() {
  // Half-step azimuths: the sector seams sit at even multiples of pi/96, so
  // every sample lands in a sector interior. The snap band itself is covered by
  // test_peirce_fast_square_ties_the_diagonal_band_to_its_seam.
  for (int latitude_step = 0; latitude_step <= 96; ++latitude_step) {
    const float y = -1.0f + 2.0f * latitude_step / 96.0f;
    const float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
    for (int longitude_step = 0; longitude_step < 96; ++longitude_step) {
      const float angle = PI_F * (2 * longitude_step + 1) / 96.0f - PI_F;
      check_peirce_fast_square_matches_exact(
          Vector(radius * cosf(angle), y, radius * sinf(angle)));
    }
  }
}

inline void test_peirce_fast_square_on_seams_and_poles() {
  constexpr float INV_SQRT_TWO = 0.7071067811865475f;
  check_peirce_fast_square_matches_exact(Vector(0.0f, 1.0f, 0.0f));
  check_peirce_fast_square_matches_exact(Vector(0.0f, -1.0f, 0.0f));
  for (int latitude_step = -7; latitude_step <= 7; ++latitude_step) {
    const float y = latitude_step / 8.0f;
    const float radius = sqrtf(1.0f - y * y);
    const float diagonal = radius * INV_SQRT_TWO;
    check_peirce_fast_square_matches_exact(Vector(diagonal, y, diagonal));
    check_peirce_fast_square_matches_exact(Vector(-diagonal, y, diagonal));
    check_peirce_fast_square_matches_exact(Vector(-diagonal, y, -diagonal));
    check_peirce_fast_square_matches_exact(Vector(diagonal, y, -diagonal));
    check_peirce_fast_square_matches_exact(Vector(radius, y, 0.0f));
    check_peirce_fast_square_matches_exact(Vector(0.0f, y, radius));
    check_peirce_fast_square_matches_exact(Vector(-radius, y, 0.0f));
    check_peirce_fast_square_matches_exact(Vector(0.0f, y, -radius));
  }
}

inline void test_peirce_fast_square_ties_the_diagonal_band_to_its_seam() {
  // Inside peirce_sector_longitude's snap band the exact kernel folds the
  // azimuth onto the boundary; the fast kernel's diagonal_tie must label
  // those samples as the seam itself, so a glued edge reads the same from
  // either side.
  constexpr float INV_SQRT_TWO = 0.7071067811865475f;
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
          Vector(sx * diagonal, y, sz * diagonal));
      for (float offset : OFFSETS) {
        const ProjectionKernelResult banded = peirce_projection_fast_square(
            Vector(sx * diagonal * (1.0f + offset), y,
                   sz * diagonal * (1.0f - offset)));
        HS_EXPECT_EQ(banded.region_id, seam.region_id);
        HS_EXPECT_EQ(banded.edge_class, seam.edge_class);
      }
    }
  }
}

inline void test_peirce_square_is_the_rotated_diamond() {
  constexpr float INV_SQRT_TWO = 0.7071067811865475f;
  for (int latitude_step = -8; latitude_step <= 8; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * PI_F / 8.0f);
    for (int longitude_step = 0; longitude_step < 12; ++longitude_step) {
      const Vector v = direction(latitude, longitude_step * (TWO_PI_F / 12.0f) -
                                               PI_F + 0.2f);
      const ProjectionKernelResult diamond =
          peirce_projection(v, 0.0f, 0, 0.0f);
      const ProjectionKernelResult square = peirce_projection(v, 0.0f, 1, 0.0f);
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

inline void test_peirce_strip_scroll_is_periodic() {
  for (uint8_t layout : {uint8_t(2), uint8_t(3)})
    for (int latitude_step = -6; latitude_step <= 6; ++latitude_step) {
      const float latitude = latitude_step * (0.5f * PI_F / 6.0f);
      for (int longitude_step = 0; longitude_step < 8; ++longitude_step) {
        const Vector v = direction(
            latitude, longitude_step * (TWO_PI_F / 8.0f) - PI_F + 0.15f);
        const ProjectionKernelResult base =
            peirce_projection(v, 0.0f, layout, 0.125f);
        const ProjectionKernelResult wrapped =
            peirce_projection(v, 0.0f, layout, 1.125f);
        HS_EXPECT_NEAR(wrapped.coords.re, base.coords.re, 1e-4f);
        HS_EXPECT_NEAR(wrapped.coords.im, base.coords.im, 1e-4f);
        HS_EXPECT_EQ(wrapped.edge_class, layout == 2 ? 4 : 5);
      }
    }
}

inline void test_airocean_cut_masks_match_the_edge_lists() {
  constexpr size_t CUT_HALF_EDGES =
      sizeof(AIROCEAN_CUT_FACES) / sizeof(AIROCEAN_CUT_FACES[0]);
  HS_EXPECT_EQ(sizeof(AIROCEAN_CUT_EDGES) / sizeof(AIROCEAN_CUT_EDGES[0]),
               CUT_HALF_EDGES);
  std::array<uint8_t, AIROCEAN_FACE_COUNT> rebuilt{};
  for (size_t i = 0; i < CUT_HALF_EDGES; ++i) {
    HS_EXPECT_LT(AIROCEAN_CUT_FACES[i], AIROCEAN_FACE_COUNT);
    HS_EXPECT_LT(AIROCEAN_CUT_EDGES[i], 3);
    rebuilt[AIROCEAN_CUT_FACES[i]] |= uint8_t(1U << AIROCEAN_CUT_EDGES[i]);
  }
  for (size_t face = 0; face < AIROCEAN_FACE_COUNT; ++face) {
    HS_EXPECT_EQ(rebuilt[face], AIROCEAN_CUT_MASKS[face]);
    for (uint8_t edge = 0; edge < 3; ++edge)
      HS_EXPECT_EQ(airocean_edge_is_cut(uint8_t(face), edge),
                   (rebuilt[face] & (1U << edge)) != 0);
  }
}

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

inline void test_airocean_glued_edges_are_bit_identical() {
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
          HS_EXPECT_TRUE(planar_segments_equal(
              planar_vertex(face, edge), planar_vertex(face, edge + 1),
              planar_vertex(other, other_edge),
              planar_vertex(other, other_edge + 1)));
          HS_EXPECT_EQ(
              airocean_edge_identity(uint8_t(face), uint8_t(edge)),
              airocean_edge_identity(uint8_t(other), uint8_t(other_edge)));
        }
  HS_EXPECT_EQ(shared, 33);
  HS_EXPECT_EQ(glued, 21);
}

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
      HS_EXPECT_TRUE(airocean_contains(v, AIROCEAN_FACES[face]));
      centroid.x += v.x / 3.0f;
      centroid.y += v.y / 3.0f;
      centroid.z += v.z / 3.0f;
    }
    HS_EXPECT_NEAR(centroid.x, center.x, 1e-6f);
    HS_EXPECT_NEAR(centroid.y, center.y, 1e-6f);
    HS_EXPECT_NEAR(centroid.z, center.z, 1e-6f);
  }
}

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

inline void test_airocean_projection_stays_inside_its_face() {
  for (int latitude_step = 0; latitude_step <= 48; ++latitude_step) {
    const float y = -1.0f + 2.0f * latitude_step / 48.0f;
    const float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
    for (int longitude_step = 0; longitude_step < 64; ++longitude_step) {
      const float angle = TWO_PI_F * (longitude_step + 0.31f) / 64.0f - PI_F;
      const Vector v(radius * cosf(angle), y, radius * sinf(angle));
      const ProjectionKernelResult net = airocean_projection(v, 0.0f, false);
      const size_t face = net.region_id;
      HS_EXPECT_LT(face, AIROCEAN_FACE_COUNT);
      HS_EXPECT_TRUE(airocean_contains(airocean_axes(v), AIROCEAN_FACES[face]));
      HS_EXPECT_EQ(
          airocean_outside_score(airocean_axes(v), AIROCEAN_FACES[face]), 0.0f);
      const AiroceanPoint point{net.coords.re, net.coords.im};
      for (size_t edge = 0; edge < 3; ++edge)
        HS_EXPECT_LE(planar_cross(planar_vertex(face, edge),
                                  planar_vertex(face, edge + 1), point),
                     1e-5f);
      HS_EXPECT_TRUE(
          net.edge_class == airocean_edge_identity(uint8_t(face), 0) ||
          net.edge_class == airocean_edge_identity(uint8_t(face), 1) ||
          net.edge_class == airocean_edge_identity(uint8_t(face), 2) ||
          net.edge_class == airocean_edge_identity(18, 0));
      HS_EXPECT_EQ(net.boundary_flags != 0,
                   (net.traits & projection_traits(ProjectionTrait::CUT)) != 0);
    }
  }
}

inline void test_airocean_projection_face_index_stays_in_range() {
  // A direction this far off the unit sphere scores above the fallback's
  // 65536 sentinel on every face, so no candidate is ever taken.
  const float magnitudes[] = {1.0f, 1e18f, 1e30f, 3.0e38f};
  for (float magnitude : magnitudes)
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
          const Vector v(seam.x * scale, seam.z * scale, seam.y * scale);
          HS_EXPECT_LT(size_t(airocean_projection(v, 0.0f, false).region_id),
                       AIROCEAN_FACE_COUNT);
        }

  const float infinity = std::numeric_limits<float>::infinity();
  const float not_a_number = std::numeric_limits<float>::quiet_NaN();
  const float exotic[] = {infinity, -infinity, not_a_number, 0.0f, 1.0f};
  for (float x : exotic)
    for (float y : exotic)
      for (float z : exotic)
        HS_EXPECT_LT(
            size_t(airocean_projection(Vector(x, y, z), 0.4f, false).region_id),
            AIROCEAN_FACE_COUNT);
}

inline void test_point_segment_distance() {
  const AiroceanPoint a{1.0f, 2.0f};
  const AiroceanPoint b{5.0f, 2.0f};
  HS_EXPECT_NEAR(point_segment_distance(a, a, b), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(b, a, b), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{3.0f, 5.0f}, a, b), 3.0f,
                 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{-3.0f, 2.0f}, a, b), 4.0f,
                 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{9.0f, 5.0f}, a, b), 5.0f,
                 1e-6f);
  HS_EXPECT_NEAR(point_segment_distance(AiroceanPoint{1.0f, 6.0f}, a, a), 4.0f,
                 1e-6f);
  const float squared =
      point_segment_distance_squared(AiroceanPoint{3.0f, 5.0f}, a, b);
  HS_EXPECT_NEAR(squared, 9.0f, 1e-6f);
}

inline void test_projection_trait_packing() {
  HS_EXPECT_EQ(projection_traits(ProjectionTrait::NONE), 0);
  HS_EXPECT_EQ(projection_traits(ProjectionTrait::CUT, ProjectionTrait::GLUED),
               3);
  HS_EXPECT_EQ(projection_traits(ProjectionTrait::CUT, ProjectionTrait::GLUED,
                                 ProjectionTrait::SINGULAR),
               19);
  const ProjectionKernelResult bonne =
      bonne_projection(Vector(1.0f, 0.0f, 0.0f), 0.0f, 0.4f);
  HS_EXPECT_EQ(bonne.traits, projection_traits(ProjectionTrait::CUT));
  const ProjectionKernelResult peirce =
      peirce_projection(Vector(1.0f, 0.0f, 0.0f), 0.0f, 1, 0.0f);
  HS_EXPECT_EQ(peirce.traits & projection_traits(ProjectionTrait::GLUED,
                                                 ProjectionTrait::FOLDED,
                                                 ProjectionTrait::PERIODIC),
               projection_traits(ProjectionTrait::GLUED,
                                 ProjectionTrait::FOLDED,
                                 ProjectionTrait::PERIODIC));
}

inline int run_projections_tests() {
  hs_test::ModuleFixture fixture("projections");
  test_wrap_longitude_range();
  test_bonne_sinusoidal_limit();
  test_bonne_polar_limit_is_finite();
  test_peirce_elliptic_integral_shape();
  test_peirce_sector_longitude_snapping();
  test_peirce_fast_square_matches_exact();
  test_peirce_fast_square_on_seams_and_poles();
  test_peirce_fast_square_ties_the_diagonal_band_to_its_seam();
  test_peirce_square_is_the_rotated_diamond();
  test_peirce_strip_scroll_is_periodic();
  test_airocean_cut_masks_match_the_edge_lists();
  test_airocean_edge_identity_is_the_canonical_half_edge();
  test_airocean_glued_edges_are_bit_identical();
  test_airocean_face_planes();
  test_airocean_unfold_sends_vertices_to_planar_vertices();
  test_airocean_projection_stays_inside_its_face();
  test_airocean_projection_face_index_stays_in_range();
  test_point_segment_distance();
  test_projection_trait_packing();
  return fixture.result();
}

} // namespace projections_tests
} // namespace hs_test
