/*
 * Airocean constants/transforms and Peirce coefficients derived from PROJ
 * commit 71eb06bfa597ed1e6d673723fc509e7604fe72af.
 *
 * Copyright (c) 2000, Frank Warmerdam
 * Copyright (c) 2019, 2020, 2021, 2022, 2023, 2024 PROJ contributors
 * Copyright (c) 2005, 2006, 2009 Gerald I. Evenden
 * Copyright (c) 2020 Kristian Evers
 * Copyright (c) 2021 Toby C Wilkinson
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */
#pragma once

/**
 * @file projections.h
 * @brief Sphere-to-plane projection kernels: folded sinusoidal,
 *        equirectangular, Bonne, Peirce quincuncial, and Fuller Airocean.
 * @details Each kernel maps a unit direction to plane coordinates; the
 * interrupted ones add the seam metadata a shader needs to fade a cut and to
 * keep a glued edge continuous. Kernels are pure and frame-independent; the
 * caller applies the coordinate scale and the pole attenuation. Constants
 * derive from PROJ at the commit named in the header above.
 */

#include "math/3dmath.h"

namespace projections {

/** @brief Topological properties a projection's image carries at its edges. */
enum class ProjectionTrait : uint8_t {
  NONE = 0,
  /** Image is torn along a seam; the two sides do not join. */
  CUT = 1U << 0,
  /** Seam sides are geometrically continuous. */
  GLUED = 1U << 1,
  /** Image tiles along an axis. */
  PERIODIC = 1U << 2,
  /** A hemisphere is reflected into the other's image. */
  FOLDED = 1U << 3,
  /** Image has at least one point of infinite scale. */
  SINGULAR = 1U << 4
};

/**
 * @brief Widens one trait to the packed trait mask.
 * @param a Trait to encode.
 * @return The mask holding just that trait.
 */
constexpr uint8_t projection_traits(ProjectionTrait a) {
  return static_cast<uint8_t>(a);
}

/**
 * @brief Packs two or three traits into one mask.
 * @param a First trait.
 * @param b Second trait.
 * @param c Optional third trait.
 * @return The union of the given traits.
 */
constexpr uint8_t projection_traits(ProjectionTrait a, ProjectionTrait b,
                                    ProjectionTrait c = ProjectionTrait::NONE) {
  return projection_traits(a) | projection_traits(b) | projection_traits(c);
}

/**
 * @brief Boundary kind a kernel reports alongside `fade_edge_distance`.
 * @details Classifies the boundary that distance is measured toward, not
 * whether the sample has reached it. A kernel whose image carries one boundary
 * kind everywhere reports it on every sample (Bonne CUT, Peirce SINGULAR);
 * only Airocean varies the mask per point, by whether the nearest face edge is
 * cut or glued.
 */
enum class ProjectionBoundary : uint8_t {
  NONE = 0,
  /** The measured edge is a cut. */
  CUT = 1U << 0,
  /** The measured edge carries a locus of infinite scale. */
  SINGULAR = 1U << 1
};

/**
 * @brief Widens one boundary kind to the packed boundary mask.
 * @param a Boundary kind to encode.
 * @return The mask holding just that kind.
 */
constexpr uint8_t projection_boundary(ProjectionBoundary a) {
  return static_cast<uint8_t>(a);
}

/**
 * @brief `fade_edge_distance` when no cut distance was measured. Larger than
 *        any real distance a kernel produces, so it also seeds their nearest-
 *        edge searches.
 */
inline constexpr float NO_EDGE_DISTANCE = 65536.0f;
inline constexpr float NO_EDGE_DISTANCE_SQUARED =
    NO_EDGE_DISTANCE * NO_EDGE_DISTANCE;

/** @brief One projection kernel's plane coordinates plus its seam metadata. */
struct ProjectionKernelResult {
  /** Plane position in the kernel's native units. */
  Complex coords{};
  /** Which sheet of an interrupted image the point fell in. */
  uint8_t region_id = 0;
  /** Disconnected component within the region. */
  uint8_t component_id = 0;
  /** ProjectionBoundary mask for the edge `fade_edge_distance` measures to. */
  uint8_t boundary_flags = 0;
  /** Distance to the nearest cut in the kernel's own units; NO_EDGE_DISTANCE
   *  when the kernel was asked to skip it or found no cut. */
  float fade_edge_distance = NO_EDGE_DISTANCE;
  /** Kernel-specific per-point flags. */
  uint8_t flags = 0;
  /** ProjectionTrait mask for the image as a whole. */
  uint8_t traits = 0;
  /** Identity of the nearest edge, shared by both sides of a glued seam. */
  uint8_t edge_class = 0;
};

/**
 * @brief Wraps a longitude into [-pi, pi].
 * @param longitude Angle in radians, unbounded.
 * @return The same direction expressed in [-pi, pi]. Both endpoints are
 *         reachable — the antimeridian can land on either sign — so a caller
 *         partitioning the range must accept it on both sides.
 */
inline float wrap_longitude(float longitude) {
  float wrapped = fmodf(longitude + PI_F, TWO_PI_F);
  if (wrapped < 0.0f)
    wrapped += TWO_PI_F;
  return wrapped - PI_F;
}

/**
 * @brief Folded sinusoidal (Sanson-Flamsteed) pseudocylindrical projection.
 * @param v Unit direction on the sphere.
 * @param central_meridian Longitude placed at the image's axis, in radians.
 * @return Plane coordinates in radians with approximate angular terms:
 *         absolute azimuth tapered by cos(latitude), against latitude.
 * @details Folding the azimuth about the central meridian maps both
 * hemispheres onto one image, so the antimeridian carries no seam; the
 * cos(latitude) taper collapses each pole to a point.
 */
HS_FLASH_INLINE inline Complex
folded_sinusoidal(const Vector &v, float central_meridian = 0.0f) {
  const float radius = sqrtf(v.x * v.x + v.z * v.z);
  return {std::fabs(wrap_longitude(fast_atan2(v.z, v.x) - central_meridian)) *
              radius,
          0.5f * PI_F - fast_acos(v.y)};
}

/**
 * @brief Equirectangular (plate carree) cylindrical projection.
 * @param v Unit direction on the sphere.
 * @param central_meridian Longitude placed at the image's axis, in radians.
 * @return Plane coordinates in radians with approximate angular terms:
 *         wrapped longitude against latitude.
 * @details Longitude is periodic, so the image is cut at the antimeridian and
 * each pole spreads across a full image row.
 */
HS_FLASH_INLINE inline Complex equirectangular(const Vector &v,
                                               float central_meridian = 0.0f) {
  return {wrap_longitude(fast_atan2(v.z, v.x) - central_meridian),
          0.5f * PI_F - fast_acos(v.y)};
}

/**
 * @brief Bonne pseudoconical equal-area projection.
 * @param v Unit direction on the sphere.
 * @param central_meridian Longitude placed at the image's axis, in radians.
 * @param standard_parallel Latitude of the parallel held true to scale, in
 *        radians; the sign selects the hemisphere the cone opens toward, 0
 *        degenerates to the sinusoidal limit and +/-pi/2 to the polar (Werner)
 *        limit.
 * @return Plane coordinates in radians, with `fade_edge_distance` set to the
 *         angular distance to the antimeridian cut.
 */
HS_FLASH_INLINE inline ProjectionKernelResult
bonne_projection(const Vector &v, float central_meridian,
                 float standard_parallel) {
  const float longitude = wrap_longitude(atan2f(v.z, v.x) - central_meridian);
  const float latitude = asinf(hs::clamp(v.y, -1.0f, 1.0f));
  const float cos_latitude = cosf(latitude);
  const float cut_distance = cos_latitude * (PI_F - fabsf(longitude));
  Complex coords;
  if (standard_parallel == 0.0f) {
    // cot(phi1) diverges here; the limit is the sinusoidal projection.
    coords = Complex(longitude * cos_latitude, latitude);
  } else {
    const float c = fabsf(standard_parallel) == 0.5f * PI_F
                        ? 0.0f
                        : cosf(standard_parallel) / sinf(standard_parallel);
    const float rho = c + standard_parallel - latitude;
    float e = 0.0f;
    if (fabsf(rho) > 1e-7f)
      e = longitude * cos_latitude / rho;
    const float half_sine = sinf(0.5f * e);
    coords = Complex(rho * sinf(e), latitude - standard_parallel +
                                        2.0f * rho * half_sine * half_sine);
  }
  return {.coords = coords,
          .region_id = static_cast<uint8_t>(longitude < 0.0f),
          .component_id = 0,
          .boundary_flags = projection_boundary(ProjectionBoundary::CUT),
          .fade_edge_distance = std::max(0.0f, cut_distance),
          .flags = 0,
          .traits = projection_traits(ProjectionTrait::CUT),
          .edge_class = 0};
}

/**
 * @brief Incomplete elliptic integral of the first kind at modulus 1/sqrt(2).
 * @param phi Amplitude in radians, |phi| <= pi/2.
 * @return F(phi, 1/sqrt(2)); F(pi/2) is the quarter period
 *         K = 1.8540746773013719 that peirce_projection tiles with.
 * @details Clenshaw recurrence over an eight-term Chebyshev expansion in
 * y = 2*(2*phi/pi)^2 - 1. `C` holds the coefficients in descending order
 * (C[0] is the highest) and the order-zero coefficient C0 is applied halved,
 * per the Clenshaw convention. The coefficients come from PROJ at the commit
 * named in the file header; they are a fit, not a closed form, and cannot be
 * rederived from anything else here.
 */
inline float peirce_elliptic_integral(float phi) {
  constexpr float C0 = 2.19174570831038f;
  constexpr float C[] = {-8.58691003636495e-7f, 2.02692115653689e-7f,
                         3.12960480765314e-5f,  5.30394739921063e-5f,
                         -0.0012804644680613f,  -0.00575574836830288f,
                         0.0914203033408211f};
  float y = phi * (2.0f / PI_F);
  y = 2.0f * y * y - 1.0f;
  const float y2 = 2.0f * y;
  float d1 = 0.0f;
  float d2 = 0.0f;
  for (float coefficient : C) {
    const float previous = d1;
    d1 = y2 * d1 - d2 + coefficient;
    d2 = previous;
  }
  return phi * (y * d1 - d2 + 0.5f * C0);
}

/**
 * @brief Longitude of a direction, snapped onto the quincuncial sector seams.
 * @param v Unit direction on the sphere.
 * @param central_meridian Longitude placed at the image's axis, in radians.
 * @return The wrapped longitude, or the exact sector boundary when within
 *         2e-6 rad of one so both sides of a seam agree; pi/2 on the poles.
 */
inline float peirce_sector_longitude(const Vector &v, float central_meridian) {
  if (v.x == 0.0f && v.z == 0.0f)
    return 0.5f * PI_F;
  float longitude = wrap_longitude(atan2f(v.z, v.x) - central_meridian);
  constexpr float BOUNDARIES[] = {-0.75f * PI_F, -0.25f * PI_F, 0.25f * PI_F,
                                  0.75f * PI_F};
  constexpr float TIE_EPSILON = 2e-6f;
  for (float boundary : BOUNDARIES)
    if (fabsf(longitude - boundary) <= TIE_EPSILON)
      return boundary;
  return longitude;
}

/**
 * @brief Peirce quincuncial projection, conformal except at four singularities.
 * @param v Unit direction on the sphere.
 * @param central_meridian Longitude placed at the image's axis, in radians.
 * @param layout 0 diamond, 1 square (the diamond turned 45 degrees), 2
 *        horizontal strip, 3 vertical strip. Layouts 0 and 1 fold the southern
 *        hemisphere into four triangles around the northern square; 2 and 3
 *        instead lay the hemispheres side by side and tile.
 * @param scroll Fraction of a full period to translate a strip layout by;
 *        ignored for layouts 0 and 1.
 * @param calculate_edge_distance When false, `fade_edge_distance` is left at
 *        NO_EDGE_DISTANCE and the inverse-trig calls that compute it are
 *        skipped.
 * @return Plane coordinates in units of the quarter period
 *         K = 1.8540746773013719; the southern fold reflects about 2K and the
 *         strip layouts repeat every 4K.
 */
HS_FLASH_INLINE inline ProjectionKernelResult
peirce_projection(const Vector &v, float central_meridian, uint8_t layout,
                  float scroll, bool calculate_edge_distance = true) {
  constexpr float INV_SQRT_TWO = 0.7071067811865475f;
  constexpr float K = 1.8540746773013719f;
  constexpr float SHIFT = 2.0f * K;
  const float longitude = peirce_sector_longitude(v, central_meridian);
  const float sl = sinf(longitude);
  const float cl = cosf(longitude);
  const float y = hs::clamp(v.y, -1.0f, 1.0f);
  const float cp = sqrtf(std::max(0.0f, 1.0f - y * y));
  const float cos_a = hs::clamp(cp * (sl + cl) * INV_SQRT_TWO, -1.0f, 1.0f);
  const float cos_b = hs::clamp(cp * (sl - cl) * INV_SQRT_TWO, -1.0f, 1.0f);
  const float sin_product =
      sqrtf(std::max(0.0f, (1.0f - cos_a * cos_a) * (1.0f - cos_b * cos_b)));
  const float cos_sum = hs::clamp(cos_a * cos_b - sin_product, -1.0f, 1.0f);
  const float cos_difference =
      hs::clamp(cos_a * cos_b + sin_product, -1.0f, 1.0f);
  float m = asinf(sqrtf(std::max(0.0f, 1.0f + std::min(0.0f, cos_sum))));
  float n = asinf(sqrtf(fabsf(1.0f - std::max(0.0f, cos_difference))));
  if (sl < 0.0f)
    m = -m;
  if (cl > 0.0f)
    n = -n;
  float x = peirce_elliptic_integral(m);
  float projected_y = peirce_elliptic_integral(n);
  uint8_t region = 0;
  uint8_t flags = 0;
  if (v.y < 0.0f) {
    flags = 1;
    if (longitude < -0.75f * PI_F || longitude >= 0.75f * PI_F)
      region = 1;
    else if (longitude < -0.25f * PI_F)
      region = 2;
    else if (longitude < 0.25f * PI_F)
      region = 3;
    else
      region = 4;
    if (layout <= 1) {
      if (longitude < -0.75f * PI_F) {
        projected_y = SHIFT - projected_y;
      } else if (longitude < -0.25f * PI_F) {
        x = -SHIFT - x;
      } else if (longitude < 0.25f * PI_F) {
        projected_y = -SHIFT - projected_y;
      } else if (longitude < 0.75f * PI_F) {
        x = SHIFT - x;
      } else {
        projected_y = SHIFT - projected_y;
      }
    }
  }
  if (layout == 1) {
    const float old_x = x;
    x = INV_SQRT_TWO * (x - projected_y);
    projected_y = INV_SQRT_TWO * (old_x + projected_y);
  } else if (layout == 2) {
    if (v.y < 0.0f)
      x = SHIFT - x;
    x -= K;
    x += scroll * 2.0f * SHIFT;
    x = fmodf(x + SHIFT, 2.0f * SHIFT);
    if (x < 0.0f)
      x += 2.0f * SHIFT;
    x -= SHIFT;
  } else if (layout == 3) {
    if (v.y < 0.0f)
      projected_y = SHIFT - projected_y;
    projected_y -= K;
    projected_y += scroll * 2.0f * SHIFT;
    projected_y = fmodf(projected_y + SHIFT, 2.0f * SHIFT);
    if (projected_y < 0.0f)
      projected_y += 2.0f * SHIFT;
    projected_y -= SHIFT;
  }
  const float rotated_x = cp * cl;
  const float rotated_z = cp * sl;
  float edge = NO_EDGE_DISTANCE;
  if (calculate_edge_distance) {
    edge = acosf(
        hs::clamp(std::max(fabsf(rotated_x), fabsf(rotated_z)), 0.0f, 1.0f));
    if (v.y < 0.0f && layout <= 1) {
      const float fold_sine = cp * fabsf(fabsf(sl) - fabsf(cl)) * INV_SQRT_TWO;
      edge = std::min(edge, asinf(hs::clamp(fold_sine, 0.0f, 1.0f)));
    }
  }
  uint8_t edge_class = fabsf(rotated_x) >= fabsf(rotated_z)
                           ? static_cast<uint8_t>(rotated_x < 0.0f)
                           : static_cast<uint8_t>(2 + (rotated_z < 0.0f));
  if (layout == 2)
    edge_class = 4;
  else if (layout == 3)
    edge_class = 5;
  return {.coords = Complex(x, projected_y),
          .region_id = region,
          .component_id = 0,
          .boundary_flags = projection_boundary(ProjectionBoundary::SINGULAR),
          .fade_edge_distance = edge,
          .flags = flags,
          .traits =
              projection_traits(ProjectionTrait::GLUED, ProjectionTrait::FOLDED,
                                ProjectionTrait::PERIODIC) |
              projection_traits(ProjectionTrait::SINGULAR),
          .edge_class = edge_class};
}

/**
 * @brief Fast square-layout Peirce projection for the zero-meridian renderer.
 * @param v Unit direction on the sphere.
 * @return Square-layout coordinates and seam metadata with approximate angular
 *         terms.
 */
HS_FLASH_INLINE inline ProjectionKernelResult
peirce_projection_fast_square(const Vector &v) {
  constexpr float INV_SQRT_TWO = 0.7071067811865475f;
  constexpr float K = 1.8540746773013719f;
  constexpr float SHIFT = 2.0f * K;
  const float cos_a = hs::clamp((v.z + v.x) * INV_SQRT_TWO, -1.0f, 1.0f);
  const float cos_b = hs::clamp((v.z - v.x) * INV_SQRT_TWO, -1.0f, 1.0f);
  const float sin_product =
      sqrtf(std::max(0.0f, (1.0f - cos_a * cos_a) * (1.0f - cos_b * cos_b)));
  const float cos_sum = hs::clamp(cos_a * cos_b - sin_product, -1.0f, 1.0f);
  const float cos_difference =
      hs::clamp(cos_a * cos_b + sin_product, -1.0f, 1.0f);
  float m = 0.5f * PI_F -
            fast_acos(sqrtf(std::max(0.0f, 1.0f + std::min(0.0f, cos_sum))));
  float n = 0.5f * PI_F -
            fast_acos(sqrtf(fabsf(1.0f - std::max(0.0f, cos_difference))));
  uint8_t sector = 0;
  uint8_t edge_class = 0;
  const float horizontal_sq = v.x * v.x + v.z * v.z;
  if (horizontal_sq < 1e-12f) {
    const float longitude = peirce_sector_longitude(v, 0.0f);
    if (longitude < 0.0f)
      m = -m;
    if (longitude > -0.5f * PI_F && longitude < 0.5f * PI_F)
      n = -n;
    if (longitude < -0.75f * PI_F || longitude >= 0.75f * PI_F)
      sector = 1;
    else if (longitude < -0.25f * PI_F)
      sector = 2;
    else if (longitude < 0.25f * PI_F)
      sector = 3;
    else
      sector = 4;
  } else {
    if (v.z < 0.0f)
      m = -m;
    if (v.x > 0.0f)
      n = -n;
    const float diagonal_delta = fabsf(v.x) - fabsf(v.z);
    const bool diagonal_tie =
        diagonal_delta * diagonal_delta <= horizontal_sq * 9.0e-12f;
    if (v.x < 0.0f) {
      if (v.z < 0.0f) {
        sector = v.z < v.x || diagonal_tie ? 2 : 1;
        edge_class = v.z > v.x || diagonal_tie ? 1 : 3;
      } else {
        sector = v.z + v.x <= 0.0f || diagonal_tie ? 1 : 4;
        edge_class = v.z + v.x > 0.0f && !diagonal_tie ? 2 : 1;
      }
    } else if (v.z < 0.0f) {
      sector = v.z + v.x >= 0.0f || diagonal_tie ? 3 : 2;
      edge_class = v.z + v.x < 0.0f && !diagonal_tie ? 3 : 0;
    } else {
      sector = v.z >= v.x || diagonal_tie ? 4 : 3;
      edge_class = v.z > v.x && !diagonal_tie ? 2 : 0;
    }
  }
  float x = peirce_elliptic_integral(m);
  float projected_y = peirce_elliptic_integral(n);
  uint8_t region = 0;
  uint8_t flags = 0;
  if (v.y < 0.0f) {
    flags = 1;
    region = sector;
    if (sector == 1) {
      projected_y = SHIFT - projected_y;
    } else if (sector == 2) {
      x = -SHIFT - x;
    } else if (sector == 3) {
      projected_y = -SHIFT - projected_y;
    } else {
      x = SHIFT - x;
    }
  }
  const float old_x = x;
  x = INV_SQRT_TWO * (x - projected_y);
  projected_y = INV_SQRT_TWO * (old_x + projected_y);
  float edge =
      fast_acos(hs::clamp(std::max(fabsf(v.x), fabsf(v.z)), 0.0f, 1.0f));
  if (v.y < 0.0f) {
    const float fold_sine = fabsf(fabsf(v.z) - fabsf(v.x)) * INV_SQRT_TWO;
    edge = std::min(edge,
                    0.5f * PI_F - fast_acos(hs::clamp(fold_sine, 0.0f, 1.0f)));
  }
  return {.coords = Complex(x, projected_y),
          .region_id = region,
          .component_id = 0,
          .boundary_flags = projection_boundary(ProjectionBoundary::SINGULAR),
          .fade_edge_distance = edge,
          .flags = flags,
          .traits =
              projection_traits(ProjectionTrait::GLUED, ProjectionTrait::FOLDED,
                                ProjectionTrait::PERIODIC) |
              projection_traits(ProjectionTrait::SINGULAR),
          .edge_class = edge_class};
}

/**
 * @brief A direction in the Airocean kernel's own axis convention.
 * @details PROJ's icosahedron tables are authored with z up, so the kernel
 * permutes the engine's y-up Vector into this type on entry.
 */
struct AiroceanVector {
  float x;
  float y;
  float z;
};

/** @brief A point in the Airocean unfolded plane. */
struct AiroceanPoint {
  float x;
  float y;
};

/**
 * @brief Vertices of the 23 spherical triangles the Airocean net unfolds.
 * @details The icosahedron has 20 faces; two of them are subdivided so the net
 * can cut through ocean rather than land, giving 23 entries. Faces 18-19 share
 * one plane and faces 20-22 share another, which is why AIROCEAN_NORMALS
 * repeats those rows. Vertex order carries the sign convention
 * airocean_contains tests against; reordering a face inverts it.
 */
inline constexpr AiroceanVector AIROCEAN_FACES[23][3] = {
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {-0.5154559599f, -0.3817168983f, 0.7672009925f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f}},
    {{0.9950094394f, -0.09134779528f, 0.04014717588f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f}},
    {{0.5154559599f, 0.3817168983f, -0.7672009925f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {-0.3557814025f, 0.8435800025f, -0.4022342266f}},
    {{-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}},
    {{-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}},
    {{-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f}},
    {{-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f}},
    {{-0.5188367303f, -0.8354203804f, -0.1813318376f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f}},
    {{0.4146822253f, -0.6559624054f, -0.6306758079f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f}},
    {{0.5154559599f, 0.3817168983f, -0.7672009925f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f}},
    {{-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.3879669146f, 0.3827173765f, -0.6531583886f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f},
     {-0.3879669146f, 0.3827173765f, -0.6531583886f}},
    {{-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.5884910224f, 0.5302967344f, 0.0627648018f}},
    {{-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f},
     {-0.5884910224f, 0.5302967344f, 0.0627648018f}},
    {{-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.5884910224f, 0.5302967344f, 0.0627648018f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}}};

constexpr AiroceanVector airocean_cross(const AiroceanVector &a,
                                        const AiroceanVector &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

struct AiroceanEdgeNormals {
  AiroceanVector values[23][3]{};

  constexpr AiroceanEdgeNormals() {
    for (int face = 0; face < 23; ++face) {
      values[face][0] =
          airocean_cross(AIROCEAN_FACES[face][1], AIROCEAN_FACES[face][2]);
      values[face][1] =
          airocean_cross(AIROCEAN_FACES[face][2], AIROCEAN_FACES[face][0]);
      values[face][2] =
          airocean_cross(AIROCEAN_FACES[face][0], AIROCEAN_FACES[face][1]);
    }
  }
};

inline constexpr AiroceanEdgeNormals AIROCEAN_EDGE_NORMALS{};
inline constexpr float AIROCEAN_CONTAINS_EPS = 1e-7f;
/** Floor on the gnomonic ray's plane component; see airocean_projection. */
inline constexpr float AIROCEAN_RAY_EPS = 1e-6f;
/**
 * @brief Each face's vertex centroid, unnormalized.
 * @details Supplies the plane offset the gnomonic ray is scaled onto.
 */
inline constexpr AiroceanVector AIROCEAN_CENTERS[23] = {
    {0.6446661988f, 0.2740726115f, 0.375187188f},
    {0.1747689772f, 0.5231760117f, 0.5720300654f},
    {-0.1699952529f, 0.1174635855f, 0.7673197837f},
    {0.08682595643f, -0.3823838838f, 0.6911725899f},
    {0.5903144229f, -0.2855941828f, 0.4488213177f},
    {0.6764340432f, 0.3752631611f, -0.1819073264f},
    {0.2261704292f, 0.6869057604f, -0.3293677939f},
    {-0.08387563251f, 0.7783209294f, 0.1365911396f},
    {-0.6417158749f, 0.1218644341f, 0.4525765415f},
    {-0.6764340432f, -0.3752631611f, 0.1819073264f},
    {-0.2261704292f, -0.6869057604f, 0.3293677939f},
    {0.08387563251f, -0.7783209294f, -0.1365911396f},
    {0.5884910224f, -0.5302967344f, -0.0627648018f},
    {0.6417158749f, -0.1218644341f, -0.4525765415f},
    {-0.5903144229f, 0.2855941828f, -0.4488213177f},
    {-0.6446661988f, -0.2740726115f, -0.375187188f},
    {-0.1747689772f, -0.5231760117f, -0.5720300654f},
    {0.1699952529f, -0.1174635855f, -0.7673197837f},
    {-0.0760974524f, 0.5360047591f, -0.6075312026f},
    {-0.09755446046f, 0.2287630085f, -0.7748139772f},
    {-0.6464272881f, 0.4884081774f, -0.1265388669f},
    {-0.4529848834f, 0.6766130474f, 0.09706879436f},
    {-0.6660608957f, 0.4258689784f, 0.2177644779f}};
/**
 * @brief Each face's unit plane normal; coplanar sub-triangles repeat a row.
 */
inline constexpr AiroceanVector AIROCEAN_NORMALS[23] = {
    {0.8112534709f, 0.3448953238f, 0.4721387736f},
    {0.2199307791f, 0.658369178f, 0.7198475379f},
    {-0.2139234835f, 0.147817183f, 0.9656017935f},
    {0.1092625279f, -0.4811951573f, 0.8697775121f},
    {0.7428567302f, -0.3593941678f, 0.5648005937f},
    {0.8512303986f, 0.4722343789f, -0.2289137389f},
    {0.284614807f, 0.8644080973f, -0.4144792552f},
    {-0.105549815f, 0.9794457296f, 0.171887461f},
    {-0.807540758f, 0.1533552486f, 0.5695261995f},
    {-0.8512303986f, -0.4722343789f, 0.2289137389f},
    {-0.284614807f, -0.8644080973f, 0.4144792552f},
    {0.105549815f, -0.9794457296f, -0.171887461f},
    {0.7405621474f, -0.6673299565f, -0.07898376463f},
    {0.807540758f, -0.1533552486f, -0.5695261995f},
    {-0.7428567302f, 0.3593941678f, -0.5648005937f},
    {-0.8112534709f, -0.3448953238f, -0.4721387736f},
    {-0.2199307791f, -0.658369178f, -0.7198475379f},
    {0.2139234835f, -0.147817183f, -0.9656017935f},
    {-0.1092625279f, 0.4811951573f, -0.8697775121f},
    {-0.1092625279f, 0.4811951573f, -0.8697775121f},
    {-0.7405621474f, 0.6673299565f, 0.07898376463f},
    {-0.7405621474f, 0.6673299565f, 0.07898376463f},
    {-0.7405621474f, 0.6673299565f, 0.07898376463f}};
/**
 * @brief Each face's vertices once unfolded into the plane, in AIROCEAN_FACES
 *        vertex order.
 * @details Two faces that share an edge carry bit-identical endpoints for it,
 * so edges match by coordinate equality rather than by tolerance.
 */
inline constexpr AiroceanPoint AIROCEAN_PLANAR_FACES[23][3] = {
    {{1.821185995f, 3.154386673f},
     {1.821185995f, 4.205848897f},
     {2.731778992f, 3.680117785f}},
    {{1.821185995f, 3.154386673f},
     {0.9105929973f, 3.680117785f},
     {1.821185995f, 4.205848897f}},
    {{1.821185995f, 3.154386673f},
     {0.9105929973f, 2.628655561f},
     {0.9105929973f, 3.680117785f}},
    {{1.821185995f, 3.154386673f},
     {1.821185995f, 2.102924448f},
     {0.9105929973f, 2.628655561f}},
    {{1.821185995f, 3.154386673f},
     {2.731778992f, 3.680117785f},
     {2.731778992f, 2.628655561f}},
    {{2.731778992f, 3.680117785f},
     {1.821185995f, 4.205848897f},
     {2.731778992f, 4.731580009f}},
    {{1.821185995f, 5.257311121f},
     {1.821185995f, 4.205848897f},
     {0.9105929973f, 4.731580009f}},
    {{0.9105929973f, 4.731580009f},
     {1.821185995f, 4.205848897f},
     {0.9105929973f, 3.680117785f}},
    {{0.9105929973f, 2.628655561f},
     {0.0f, 3.154386673f},
     {0.9105929973f, 3.680117785f}},
    {{0.9105929973f, 2.628655561f},
     {0.9105929973f, 1.577193336f},
     {0.0f, 2.102924448f}},
    {{0.9105929973f, 2.628655561f},
     {1.821185995f, 2.102924448f},
     {0.9105929973f, 1.577193336f}},
    {{0.9105929973f, 1.577193336f},
     {1.821185995f, 2.102924448f},
     {1.821185995f, 1.051462224f}},
    {{1.821185995f, 1.051462224f},
     {1.821185995f, 2.102924448f},
     {2.731778992f, 1.577193336f}},
    {{1.821185995f, 0.0f},
     {1.821185995f, 1.051462224f},
     {2.731778992f, 0.5257311121f}},
    {{0.0f, 5.257311121f}, {0.9105929973f, 4.731580009f}, {0.0f, 4.205848897f}},
    {{0.0f, 1.051462224f}, {0.0f, 2.102924448f}, {0.9105929973f, 1.577193336f}},
    {{0.9105929973f, 0.5257311121f},
     {0.9105929973f, 1.577193336f},
     {1.821185995f, 1.051462224f}},
    {{0.9105929973f, 0.5257311121f},
     {1.821185995f, 1.051462224f},
     {1.821185995f, 0.0f}},
    {{0.9105929973f, 4.731580009f},
     {0.4552964987f, 4.994445565f},
     {0.9105929973f, 5.783042233f}},
    {{0.9105929973f, 0.5257311121f},
     {1.821185995f, 0.0f},
     {0.9105929973f, 0.0f}},
    {{0.0f, 4.205848897f},
     {0.9105929973f, 4.731580009f},
     {0.6070619982f, 4.205848897f}},
    {{0.9105929973f, 4.731580009f},
     {0.9105929973f, 3.680117785f},
     {0.6070619982f, 4.205848897f}},
    {{0.0f, 3.154386673f},
     {0.3035309991f, 3.680117785f},
     {0.9105929973f, 3.680117785f}}};
/**
 * @brief Per-face affine map from the gnomonic point on the face plane to the
 *        unfolded plane.
 * @details Row `r` holds `{ax, ay, az, b}`, evaluated as
 * `out[r] = ax*q.x + ay*q.y + az*q.z + b`.
 */
inline constexpr float AIROCEAN_TRANSFORMS[23][2][4] = {
    {{0.5771127853f, -0.6019490725f, -0.5519041105f, 2.124716994f},
     {0.09385435001f, 0.7202114479f, -0.6873767753f, 3.680117785f}},
    {{0.9709901201f, -0.2187361325f, -0.09660585362f, 1.517654996f},
     {0.09385435001f, 0.7202114479f, -0.6873767753f, 3.680117785f}},
    {{0.9721374115f, -0.06476823822f, 0.2252863255f, 1.214123996f},
     {0.09584151699f, 0.9868916636f, -0.1298431665f, 3.154386673f}},
    {{0.9921258754f, -0.001098710673f, -0.1252399307f, 1.517654996f},
     {0.061220482f, 0.876612807f, 0.4772861187f, 2.628655561f}},
    {{0.280304148f, -0.5991800397f, -0.7499419075f, 2.428247993f},
     {0.6079419899f, 0.7154153424f, -0.3443652491f, 3.154386673f}},
    {{0.2596066191f, -0.7580069591f, -0.5983559587f, 2.428247993f},
     {-0.4560824616f, 0.4499112594f, -0.7678337365f, 4.205848897f}},
    {{0.9586365701f, -0.2580860646f, 0.1200312867f, 1.517654996f},
     {-0.003215303703f, -0.4314976531f, -0.902108329f, 4.731580009f}},
    {{0.9928349404f, 0.09405868119f, 0.07370037669f, 1.214123996f},
     {0.05601801133f, 0.1784349382f, -0.9823558191f, 4.205848897f}},
    {{0.5819727896f, 0.05026939416f, 0.8116530418f, 0.6070619982f},
     {0.09584151699f, 0.9868916636f, -0.1298431665f, 3.154386673f}},
    {{0.5247823075f, -0.7686380597f, 0.3657855423f, 0.6070619982f},
     {0.003215303703f, 0.4314976531f, 0.902108329f, 2.102924448f}},
    {{0.9586365701f, -0.2580860646f, 0.1200312867f, 1.214123996f},
     {0.003215303703f, 0.4314976531f, 0.902108329f, 2.102924448f}},
    {{0.9928349404f, 0.09405868119f, 0.07370037669f, 1.517654996f},
     {-0.05601801133f, -0.1784349382f, 0.9823558191f, 1.577193336f}},
    {{0.6696489291f, 0.7230710214f, 0.1695246581f, 2.124716994f},
     {-0.05601801133f, -0.1784349382f, 0.9823558191f, 1.577193336f}},
    {{0.5819727896f, 0.05026939416f, 0.8116530418f, 2.124716994f},
     {-0.09584151699f, -0.9868916636f, 0.1298431665f, 0.5257311121f}},
    {{0.3863411333f, 0.9191578806f, 0.07674189989f, 0.3035309991f},
     {0.5467215079f, -0.1611974646f, -0.8216513678f, 4.731580009f}},
    {{0.2072761413f, -0.9246959463f, 0.3193336941f, 0.3035309991f},
     {-0.5467215079f, 0.1611974646f, 0.8216513678f, 1.577193336f}},
    {{0.9709901201f, -0.2187361325f, -0.09660585362f, 1.214123996f},
     {-0.09385435001f, -0.7202114479f, 0.6873767753f, 1.051462224f}},
    {{0.9721374115f, -0.06476823822f, 0.2252863255f, 1.517654996f},
     {-0.09584151699f, -0.9868916636f, 0.1298431665f, 0.5257311121f}},
    {{0.5490814303f, 0.7586196048f, 0.3507219383f, 0.6070619982f},
     {0.8285959708f, -0.4392579149f, -0.347104021f, 5.257311121f}},
    {{0.9921258754f, -0.001098710673f, -0.1252399307f, 1.214123996f},
     {-0.061220482f, -0.876612807f, -0.4772861187f, 0.0f}},
    {{0.6696489291f, 0.7230710214f, 0.1695246581f, 0.6070619982f},
     {0.05601801133f, 0.1784349382f, -0.9823558191f, 4.205848897f}},
    {{0.6696489291f, 0.7230710214f, 0.1695246581f, 0.6070619982f},
     {0.05601801133f, 0.1784349382f, -0.9823558191f, 4.205848897f}},
    {{0.2863114437f, 0.2070063213f, 0.9355074239f, 0.3035309991f},
     {0.6079419899f, 0.7154153424f, -0.3443652491f, 3.680117785f}}};

/** @brief Extent of the unfolded net along its second axis, the axis the
 *         quarter-turn layout reflects about. */
inline constexpr float AIROCEAN_NET_HEIGHT = 5.78304223331047f;

/**
 * @brief The net's torn edges as parallel (face, edge) arrays.
 * @details Authoring form of AIROCEAN_CUT_MASKS: 26 half-edges, one entry per
 * index in both arrays. AIROCEAN_CUT_MASKS is the per-face bitset the kernel
 * actually reads.
 */
inline constexpr uint8_t AIROCEAN_CUT_FACES[] = {
    3,  4,  4,  5,  5,  6,  6,  8,  9,  12, 12, 13, 13,
    14, 14, 15, 15, 16, 18, 18, 19, 19, 20, 21, 22, 22};
/** @brief Edge index within the face at the same position in
 *         AIROCEAN_CUT_FACES. */
inline constexpr uint8_t AIROCEAN_CUT_EDGES[] = {0, 1, 2, 1, 2, 0, 2, 0, 2,
                                                 1, 2, 1, 2, 0, 2, 0, 2, 0,
                                                 1, 2, 1, 2, 2, 1, 0, 1};

/** @brief Per face, a bit per edge index set when that edge is torn. */
inline constexpr uint8_t AIROCEAN_CUT_MASKS[23] = {
    0, 0, 0, 1, 6, 6, 5, 0, 1, 4, 0, 0, 6, 6, 5, 5, 1, 0, 6, 6, 4, 2, 3};

/**
 * @brief Per (face, edge), the shared identity of that geometric edge.
 * @details The identity is `canonical_face * 3 + canonical_edge`, where the
 * canonical half-edge is the lowest-numbered face carrying the same planar
 * segment. Both halves of a seam therefore report the same value, which lets
 * the shader treat the two sides asymmetrically without knowing which face it
 * landed on.
 */
inline constexpr uint8_t AIROCEAN_EDGE_IDENTITIES[23][3] = {
    {0, 1, 2},    {3, 4, 0},    {6, 7, 3},    {9, 10, 6},   {2, 13, 14},
    {1, 16, 17},  {18, 19, 20}, {19, 4, 23},  {24, 25, 7},  {27, 28, 29},
    {10, 31, 27}, {31, 34, 35}, {34, 37, 38}, {39, 40, 41}, {42, 43, 44},
    {45, 28, 47}, {48, 35, 50}, {50, 39, 53}, {54, 55, 56}, {53, 58, 59},
    {43, 61, 62}, {23, 64, 61}, {66, 67, 25}};

/**
 * @brief Reports whether a face's edge is torn rather than glued.
 * @param face Face index in [0, 23).
 * @param edge Edge index in [0, 3).
 * @return True when the net tears along that edge.
 */
inline bool airocean_edge_is_cut(uint8_t face, uint8_t edge) {
  return (AIROCEAN_CUT_MASKS[face] & (1U << edge)) != 0;
}

/**
 * @brief Shared identity of a face's edge.
 * @param face Face index in [0, 23).
 * @param edge Edge index in [0, 3).
 * @return The identity both halves of that seam report.
 */
inline uint8_t airocean_edge_identity(uint8_t face, uint8_t edge) {
  return AIROCEAN_EDGE_IDENTITIES[face][edge];
}

inline float airocean_edge_halfspace(const AiroceanVector &p, uint8_t face,
                                     uint8_t edge) {
  const AiroceanVector &normal = AIROCEAN_EDGE_NORMALS.values[face][edge];
  return p.x * normal.x + p.y * normal.y + p.z * normal.z;
}

/**
 * @brief Tests a direction against a spherical triangle.
 * @param p Direction to test; need not be normalized.
 * @param face Face index in [0, 23).
 * @return True when `p` lies in the triangle's cone, boundary included.
 */
inline bool airocean_contains(const AiroceanVector &p, uint8_t face) {
  return airocean_edge_halfspace(p, face, 0) <= AIROCEAN_CONTAINS_EPS &&
         airocean_edge_halfspace(p, face, 1) <= AIROCEAN_CONTAINS_EPS &&
         airocean_edge_halfspace(p, face, 2) <= AIROCEAN_CONTAINS_EPS;
}

/**
 * @brief How far outside a spherical triangle a direction falls.
 * @param p Direction to test; need not be normalized.
 * @param face Face index in [0, 23).
 * @return 0 when contained, otherwise the largest violated half-space
 *         determinant. Used to pick the least-wrong face when rounding leaves
 *         a direction in no triangle at all.
 */
inline float airocean_outside_score(const AiroceanVector &p, uint8_t face) {
  return std::max(0.0f,
                  std::max(airocean_edge_halfspace(p, face, 0),
                           std::max(airocean_edge_halfspace(p, face, 1),
                                    airocean_edge_halfspace(p, face, 2))));
}

/**
 * @brief Squared distance from a point to a segment.
 * @param p Query point.
 * @param a Segment start.
 * @param b Segment end; equal to `a` degenerates to the point distance.
 * @return The squared distance, in the unfolded plane's units.
 */
inline float point_segment_distance_squared(const AiroceanPoint &p,
                                            const AiroceanPoint &a,
                                            const AiroceanPoint &b) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float denom = dx * dx + dy * dy;
  const float t = denom == 0.0f
                      ? 0.0f
                      : hs::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy) / denom,
                                  0.0f, 1.0f);
  const float ex = p.x - (a.x + t * dx);
  const float ey = p.y - (a.y + t * dy);
  return ex * ex + ey * ey;
}

/**
 * @brief Distance from a point to a segment.
 * @param p Query point.
 * @param a Segment start.
 * @param b Segment end.
 * @return The distance, in the unfolded plane's units.
 */
inline float point_segment_distance(const AiroceanPoint &p,
                                    const AiroceanPoint &a,
                                    const AiroceanPoint &b) {
  return sqrtf(point_segment_distance_squared(p, a, b));
}

/**
 * @brief Fuller Airocean (Dymaxion) projection onto the unfolded icosahedral
 *        net.
 * @param v Unit direction on the sphere.
 * @param central_meridian Longitude rotated onto the net's axis, in radians.
 * @param horizontal Rotates the finished net a quarter turn.
 * @param calculate_edge_distance When false, `fade_edge_distance` is left at
 *        NO_EDGE_DISTANCE and the per-edge cut distances are skipped.
 * @return Plane coordinates on the net, with `region_id` set to the face and
 *         `edge_class` to the nearest edge's shared identity.
 * @details Each direction is assigned to the spherical triangle that contains
 * it, gnomonically projected onto that face's plane, then mapped by the face's
 * affine transform. Face 14's first edge is half torn and half glued, so its
 * cut distance runs to face 18's vertex and its identity switches partway
 * along.
 */
HS_FLASH_INLINE inline ProjectionKernelResult
airocean_projection(const Vector &v, float central_meridian, bool horizontal,
                    bool calculate_edge_distance = true) {
  const float c = cosf(central_meridian);
  const float s = sinf(central_meridian);
  const AiroceanVector p{v.x * c + v.z * s, v.z * c - v.x * s, v.y};
  uint8_t face = 0;
  for (; face < 23; ++face)
    if (airocean_contains(p, face))
      break;
  if (face == 23) {
    // A far-from-unit direction can score above the sentinel on every face.
    face = 0;
    float best_score = 65536.0f;
    for (uint8_t candidate = 0; candidate < 23; ++candidate) {
      const float score = airocean_outside_score(p, candidate);
      if (score < best_score) {
        best_score = score;
        face = candidate;
      }
    }
  }
  const AiroceanVector &center = AIROCEAN_CENTERS[face];
  const AiroceanVector &normal = AIROCEAN_NORMALS[face];
  const float plane =
      center.x * normal.x + center.y * normal.y + center.z * normal.z;
  const float ray = p.x * normal.x + p.y * normal.y + p.z * normal.z;
  // The outside-score fallback can settle on a face the ray does not cross, so
  // floor the divisor; copysignf keys on the sign bit, so -0.0f floors negative.
  const float divisor =
      fabsf(ray) < AIROCEAN_RAY_EPS ? copysignf(AIROCEAN_RAY_EPS, ray) : ray;
  const float scale = plane / divisor;
  const AiroceanVector q{p.x * scale, p.y * scale, p.z * scale};
  const float (&transform)[2][4] = AIROCEAN_TRANSFORMS[face];
  AiroceanPoint output{transform[0][0] * q.x + transform[0][1] * q.y +
                           transform[0][2] * q.z + transform[0][3],
                       transform[1][0] * q.x + transform[1][1] * q.y +
                           transform[1][2] * q.z + transform[1][3]};
  uint8_t edge_class = 0;
  float cut_edge_distance_squared = NO_EDGE_DISTANCE_SQUARED;
  float face_edge_distance_squared = NO_EDGE_DISTANCE_SQUARED;
  for (uint8_t candidate = 0; candidate < 3; ++candidate) {
    const AiroceanPoint &a = AIROCEAN_PLANAR_FACES[face][candidate];
    const AiroceanPoint &b = AIROCEAN_PLANAR_FACES[face][(candidate + 1) % 3];
    const float distance_squared = point_segment_distance_squared(output, a, b);
    if (distance_squared < face_edge_distance_squared) {
      face_edge_distance_squared = distance_squared;
      edge_class = candidate;
    }
    if (calculate_edge_distance && airocean_edge_is_cut(face, candidate)) {
      const float cut_distance_squared =
          face == 14 && candidate == 0
              ? point_segment_distance_squared(output, a,
                                               AIROCEAN_PLANAR_FACES[18][1])
              : distance_squared;
      cut_edge_distance_squared =
          std::min(cut_edge_distance_squared, cut_distance_squared);
    }
  }
  const float edge = calculate_edge_distance && cut_edge_distance_squared <
                                                    NO_EDGE_DISTANCE_SQUARED
                         ? sqrtf(cut_edge_distance_squared)
                         : NO_EDGE_DISTANCE;
  bool cut_edge = airocean_edge_is_cut(face, edge_class);
  uint8_t edge_identity = airocean_edge_identity(face, edge_class);
  if (face == 14 && edge_class == 0) {
    const AiroceanPoint &a = AIROCEAN_PLANAR_FACES[14][0];
    const AiroceanPoint &b = AIROCEAN_PLANAR_FACES[14][1];
    const float dx = b.x - a.x;
    const float dy = b.y - a.y;
    const float t =
        ((output.x - a.x) * dx + (output.y - a.y) * dy) / (dx * dx + dy * dy);
    cut_edge = t <= 0.5f;
    if (!cut_edge)
      edge_identity = airocean_edge_identity(18, 0);
  }
  if (horizontal)
    output = {AIROCEAN_NET_HEIGHT - output.y, output.x};
  return {.coords = Complex(output.x, output.y),
          .region_id = face,
          .component_id = 0,
          .boundary_flags = static_cast<uint8_t>(
              cut_edge ? projection_boundary(ProjectionBoundary::CUT)
                       : projection_boundary(ProjectionBoundary::NONE)),
          .fade_edge_distance = edge,
          .flags = 0,
          .traits = projection_traits(cut_edge ? ProjectionTrait::CUT
                                               : ProjectionTrait::GLUED),
          .edge_class = edge_identity};
}

} // namespace projections
