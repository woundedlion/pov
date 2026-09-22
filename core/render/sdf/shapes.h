/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include "math/geometry.h"
#include "platform/constants.h"
#include "render/sdf/common.h"

/**
 * @file shapes.h
 * @brief The polygon, star, flower and line leaves.
 */

namespace SDF {

/**
 * @brief Calculates signed distance to a planar polygon.
 * @details Register semantics: the DistanceResult table (row: PlanarPolygon).
 */
struct PlanarPolygon {
  const math::Basis
      &basis;         /**< Orientation frame (v = polygon axis); retained by
                         reference, so it must outlive the shape. */
  float circumradius; /**< Angular radius from center to vertex (radians). */
  int sides;          /**< Number of polygon sides. */
  float phase;        /**< Azimuth phase offset (radians). */
  float sector;       /**< Angular width of one polygon sector. */
  float reciprocal_sector; /**< Reciprocal angular sector width. */
  float apothem;           /**< Precomputed inradius (radians). */
  float ny, r_val,
      alpha_angle; /**< Axis y-component, XZ projection length and azimuth. */
  float cos_cap, sin_cap; /**< Circumradius trig for the scanline cap pad. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float sign;             /**< +1 fills the polygon, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a planar polygon from its basis, radius, side count, phase.
   * @param b Orientation frame (v = polygon axis); retained by reference and
   *          read by every distance() call, so it must outlive the shape.
   * @param radius Circumradius as a fraction of the hemisphere.
   * @param s Number of polygon sides (must be >= 3).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  PlanarPolygon(const math::Basis &b, float radius, int s, float ph,
                bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides >= 3, "SDF PlanarPolygon: sides must be at least 3");
    HS_CHECK(radius > 0.0f, "SDF PlanarPolygon: radius must be positive");
    // arc_stretch<PlanarPolygon> = 2 holds only within a hemisphere; a wider
    // shape must be built inverted, about its antipode.
    HS_CHECK(radius <= 1.0f, "SDF PlanarPolygon: radius exceeds unit sphere");
    circumradius = radius * (math::PI_F / 2.0f);
    sector = math::TWO_PI_F / sides;
    reciprocal_sector = static_cast<float>(sides) / math::TWO_PI_F;
    apothem = circumradius * cosf(math::PI_F / sides);

    CapBounds cb = cap_bounds(basis.v, circumradius, invert);
    ny = cb.ny;
    r_val = cb.r_val;
    alpha_angle = cb.alpha_angle;
    phi_min = cb.phi_min;
    phi_max = cb.phi_max;
    cos_cap = cb.cos_radius;
    sin_cap = cb.sin_radius;
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The polygon retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  PlanarPolygon(const math::Basis &&, float, int, float, bool = false) = delete;

  /**
   * @brief Maps the polygon's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the polygon plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the bounding-cap interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   * @details Covers the polygon body and its whole AA fringe: distance()
   *   clamps against the circumscribed disc, so nothing past circumradius +
   *   pixel_width shades.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    return emit_padded_cap_row<W, H>(sign, cos_cap, sin_cap, ny, r_val,
                                     alpha_angle, y, out);
  }

  /**
   * @brief Signed distance to the planar polygon edge, writing into res.
   * @tparam ComputeUVs When true, also computes the normalized radial t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = polar*cos(local) - apothem clamped below
   *        by polar - circumradius, raw_dist = polar angle from center, t =
   *        polar/circumradius when ComputeUVs.
   * @note `polar*cos(local) - apothem` is the azimuthal-equidistant chart's
   *       distance to the edge line: exact along the apothem inside the
   *       circumscribed disc, under-estimating near the sector corners
   *       (gradient cos(PI/sides) there), and over-reporting beyond the disc,
   *       where the chart stretches tangential distances by polar/sin(polar)
   *       and the nearest boundary point is a vertex (2.02 vs 1.57 rad near
   *       the antipode of a square of circumradius PI/2). The disc term
   *       `polar - circumradius` is a lower bound of the true distance; the
   *       max of the two is not, so for scanline shading, not a march-safe
   *       metric.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float polar =
        math::fast_acos(hs::clamp(math::dot(p, basis.v), -1.0f, 1.0f));

    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local = centered_sector_angle(azimuth, sector, reciprocal_sector);

    float dist_edge = std::max(polar * math::fast_cosf(local) - apothem,
                               polar - circumradius);
    float t_val = 0.0f;
    if constexpr (ComputeUVs)
      t_val = polar / circumradius;

    res = DistanceResult(sign * dist_edge, t_val, polar, 0.0f, apothem);
  }
};

/**
 * @brief Calculates signed distance to a spherical polygon (great-circle
 * edges).
 * @details Uses sector folding plus a precomputed great-circle plane normal for
 * O(1) per-pixel distance, with exact angular distances for smooth AA.
 * Register semantics: the DistanceResult table (row: SphericalPolygon).
 */
struct SphericalPolygon {
  const math::Basis &basis; /**< Orientation frame (v = polygon axis); retained
                              by reference, so it must outlive the shape. */
  int sides;                /**< Number of polygon sides. */
  float phase;              /**< Azimuth phase offset (radians). */
  float sector;             /**< Angular width of one polygon sector. */
  float reciprocal_sector;  /**< Reciprocal angular sector width. */
  float circumradius; /**< Angular distance from center to vertex (radians). */
  float edge_nv;      /**< Edge normal dotted with the center axis. */
  float edge_nu;      /**< Edge normal dotted with the u-axis. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float ny, r_val,
      alpha_angle; /**< Axis y-component, XZ projection length and azimuth. */
  float cos_cap, sin_cap; /**< Circumradius trig for the scanline cap pad. */
  float sign;             /**< +1 fills the polygon, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a spherical polygon from its basis, radius, side count,
   * phase.
   * @param b Orientation frame (v = polygon axis); retained by reference and
   *          read by every distance() call, so it must outlive the shape.
   * @param radius Polygon radius as a fraction of the hemisphere.
   * @param s Number of polygon sides (must be >= 3).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  SphericalPolygon(const math::Basis &b, float radius, int s, float ph,
                   bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides >= 3, "SDF SphericalPolygon: sides must be at least 3");
    HS_CHECK(radius > 0.0f, "SDF SphericalPolygon: radius must be positive");
    // A shape wider than a hemisphere must be built inverted, about its
    // antipode.
    HS_CHECK(radius <= 1.0f,
             "SDF SphericalPolygon: radius exceeds unit sphere");
    sector = math::TWO_PI_F / sides;
    reciprocal_sector = static_cast<float>(sides) / math::TWO_PI_F;
    circumradius = radius * (math::PI_F / 2.0f);

    // Build canonical edge: between vertices at azimuth ±π/n from
    // the sector bisector (u-axis), at angular distance circumradius
    float half_step = math::PI_F / sides;
    float sin_r = sinf(circumradius);
    float cos_r = cosf(circumradius);
    float cos_hs = cosf(half_step);
    float sin_hs = sinf(half_step);

    math::Vector v1 =
        basis.v * cos_r + (basis.u * cos_hs + basis.w * sin_hs) * sin_r;
    math::Vector v2 =
        basis.v * cos_r + (basis.u * cos_hs - basis.w * sin_hs) * sin_r;

    // Normal pointing outward (away from polygon interior)
    math::Vector en = math::cross(v2, v1);
    float len = en.magnitude();
    if (len > 1e-9f) {
      en = en * (1.0f / len);
    } else {
      // Degenerate canonical edge (circumradius near 0 or PI): both limits
      // drive en toward ±u, so substitute the sector bisector as a defined unit
      // normal.
      en = basis.u;
    }
    // Ensure outward: dot(center, n) should be negative
    if (math::dot(en, basis.v) > 0)
      en = -en;

    edge_nv = math::dot(en, basis.v);
    edge_nu = math::dot(en, basis.u);

    CapBounds cb = cap_bounds(basis.v, circumradius, invert);
    ny = cb.ny;
    r_val = cb.r_val;
    alpha_angle = cb.alpha_angle;
    phi_min = cb.phi_min;
    phi_max = cb.phi_max;
    cos_cap = cb.cos_radius;
    sin_cap = cb.sin_radius;
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The polygon retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  SphericalPolygon(const math::Basis &&, float, int, float,
                   bool = false) = delete;

  /**
   * @brief Maps the polygon's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the polygon plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the bounding-cap interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   * @details Covers the polygon body and its whole AA fringe: distance()
   *   clamps against the circumscribed disc, so nothing past circumradius +
   *   pixel_width shades.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    return emit_padded_cap_row<W, H>(sign, cos_cap, sin_cap, ny, r_val,
                                     alpha_angle, y, out);
  }

  /**
   * @brief Signed distance to the spherical polygon, writing into res.
   * @tparam ComputeUVs When true, also computes the normalized radial t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed angular distance to the nearest
   *        great-circle edge clamped below by polar - circumradius, raw_dist =
   *        polar angle, t = polar/circumradius when ComputeUVs.
   * @note The folded edge dot measures distance to the whole great circle, not
   *       the edge arc, so past a vertex it under-estimates with a radial
   *       gradient of cos(PI/sides). The circumscribed-disc distance
   *       `polar - circumradius` is the tighter bound out there; both are lower
   *       bounds of the true distance, so the max of the two is too.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float cos_p = hs::clamp(math::dot(p, basis.v), -1.0f, 1.0f);
    float polar = math::fast_acos(cos_p);

    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local = centered_sector_angle(azimuth, sector, reciprocal_sector);

    // Angular distance to the nearest great circle edge via precomputed normal
    // cos(local) is even, so sector folding works automatically
    float sin_p = sqrtf(std::max(0.0f, 1.0f - cos_p * cos_p));
    float dp = edge_nv * cos_p + edge_nu * math::fast_cosf(local) * sin_p;
    float dist_edge =
        std::max(asinf(hs::clamp(dp, -1.0f, 1.0f)), polar - circumradius);

    float t_val = 0.0f;
    if constexpr (ComputeUVs)
      t_val = polar / circumradius;
    res = DistanceResult(sign * dist_edge, t_val, polar, 0.0f, circumradius);
  }

  /**
   * @brief Signed edge-plane dot for sine-domain solid antialiasing.
   * @param p Point on sphere (normalized).
   * @return Signed sine-domain distance, with a hemisphere bound away from
   *         the edge.
   * @note Carries distance()'s circumscribed-disc clamp in the sine domain:
   *       sin(polar - circumradius) expands to sin_p*cos_cap - cos_p*sin_cap,
   *       so the tighter bound past a vertex costs no transcendental here.
   *       The hemisphere bound can dominate only beyond angular distance PI/2.
   */
  float sine_distance(const math::Vector &p) const {
    float cos_p = hs::clamp(math::dot(p, basis.v), -1.0f, 1.0f);
    float sin_p = sqrtf(std::max(0.0f, 1.0f - cos_p * cos_p));

    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local = centered_sector_angle(azimuth, sector, reciprocal_sector);
    float dp = edge_nv * cos_p + edge_nu * math::fast_cosf(local) * sin_p;
    float disc = sin_p * cos_cap - cos_p * sin_cap;
    return sign * std::max(hs::clamp(dp, -cos_p, 1.0f), disc);
  }
};

/**
 * @brief Calculates signed distance to a star shape.
 * @details Register semantics: the DistanceResult table (row: Star).
 */
struct Star {
  const math::Basis &basis; /**< Orientation frame (v = star axis); retained by
                              reference, so it must outlive the shape. */
  int sides;                /**< Number of star points. */
  float phase;              /**< Azimuth phase offset (radians). */
  float sector;             /**< Angular width of one star sector. */
  float reciprocal_sector;  /**< Reciprocal angular sector width. */
  static constexpr bool is_solid =
      true; /**< Star renders as a filled region. */

  float edge_nx, edge_ny,
      plane_d;        /**< 2D edge plane (normal and offset) for one point. */
  float circumradius; /**< Angular radius from center to point tip (radians). */

  float ny, r_val,
      alpha_angle; /**< Axis y-component, XZ projection length and azimuth. */
  float cos_cap, sin_cap; /**< Circumradius trig for the scanline cap pad. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float sign;             /**< +1 fills the star, -1 its complement. */

  /**
   * @brief Builds a star from its basis, radius, point count, and phase.
   * @param b Orientation frame (v = star axis); retained by reference and read
   *          by every distance() call, so it must outlive the shape.
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of star points (must be >= 3).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  Star(const math::Basis &b, float radius, int s, float ph, bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides >= 3, "SDF Star: sides must be at least 3");
    HS_CHECK(radius > 0.0f, "SDF Star: radius must be positive");
    // arc_stretch<Star> = 2 holds only within a hemisphere; a wider shape must
    // be built inverted, about its antipode.
    HS_CHECK(radius <= 1.0f, "SDF Star: radius exceeds unit sphere");
    sector = math::TWO_PI_F / sides;
    reciprocal_sector = static_cast<float>(sides) / math::TWO_PI_F;
    float outer_radius = radius * (math::PI_F / 2.0f);
    float inner_radius = outer_radius * STAR_INNER_RATIO;
    float angle_step = math::PI_F / sides;

    float v_t = outer_radius;
    float v_vx = inner_radius * cosf(angle_step);
    float v_vy = inner_radius * sinf(angle_step);

    float dx = v_vx - v_t;
    float dy = v_vy;
    float len = sqrtf(dx * dx + dy * dy);
    edge_nx = -dy / len;
    edge_ny = dx / len;
    plane_d = -(edge_nx * v_t);
    circumradius = outer_radius;

    CapBounds cb = cap_bounds(basis.v, outer_radius, invert);
    ny = cb.ny;
    r_val = cb.r_val;
    alpha_angle = cb.alpha_angle;
    phi_min = cb.phi_min;
    phi_max = cb.phi_max;
    cos_cap = cb.cos_radius;
    sin_cap = cb.sin_radius;
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The star retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  Star(const math::Basis &&, float, int, float, bool = false) = delete;

  /**
   * @brief Maps the star's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the star plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the bounding-circle interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   * @details Covers the star body and its whole AA fringe: distance() clamps
   *   against the circumscribed disc, so nothing past circumradius +
   *   pixel_width shades.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    return emit_padded_cap_row<W, H>(sign, cos_cap, sin_cap, ny, r_val,
                                     alpha_angle, y, out);
  }

  /**
   * @brief Signed distance to the star, writing into res.
   * @tparam ComputeUVs When true, also stores the normalized azimuth t. The
   *        azimuth itself is always computed (the petal-sector geometry needs
   *        it); only the final t store is gated, matching Ring/Flower.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance to the nearest point edge,
   *        raw_dist = polar distance from center, t = normalized azimuth when
   *        ComputeUVs (0 otherwise).
   * @note Folding a sector onto one edge half-plane gives a radial gradient of
   *       |edge_nx| at a tip (0.309 at 5 points, 0.220 at 8), which reads a
   *       fringe several pixels past the tip; the circumscribed-disc distance
   *       `scan_dist - circumradius` is the tighter bound there. The edge term
   *       is a chart-plane distance, so beyond the disc it over-reports the
   *       true distance where the chart stretches tangential distances by
   *       scan_dist/sin(scan_dist) and the nearest boundary point is a tip
   *       (0.89 vs 0.78 rad at 1.47 rad along a notch of a 5-point star of
   *       radius 0.6). The disc term is a lower bound of the true distance;
   *       the max of the two is not, so for scanline shading, not a
   *       march-safe metric.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float scan_dist =
        math::fast_acos(hs::clamp(math::dot(p, basis.v), -1.0f, 1.0f));
    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local_azimuth =
        centered_sector_angle(azimuth, sector, reciprocal_sector);
    local_azimuth = std::abs(local_azimuth);

    float px = scan_dist * math::fast_cosf(local_azimuth);
    float py = scan_dist * math::fast_sinf(local_azimuth);

    float dist_edge = std::max(-(px * edge_nx + py * edge_ny + plane_d),
                               scan_dist - circumradius);

    float t = 0.0f;
    if constexpr (ComputeUVs)
      t = math::wrap_t(azimuth / math::TWO_PI_F);
    res = DistanceResult(sign * dist_edge, t, scan_dist, 0.0f, circumradius);
  }
};

/**
 * @brief Calculates signed distance to a flower shape.
 * @details Register semantics: the DistanceResult table (row: Flower).
 */
struct Flower {
  const math::Basis
      &basis;   /**< Orientation frame (-v = flower center); retained by
                              reference, so it must outlive the shape. */
  int sides;    /**< Number of petals. */
  float phase;  /**< Azimuth phase offset (radians). */
  float sector; /**< Angular width of one flower sector. */
  float reciprocal_sector; /**< Reciprocal angular sector width. */
  float circumradius;      /**< Angular radius from the antipode to petal tip
                         (radians). */
  float apothem;           /**< Petal inradius offset (PI - outer radius). */
  math::Vector antipode;   /**< Antipode of the flower axis (scan origin). */
  float ny, r_val, alpha_angle; /**< Antipode y-component, XZ projection length
                                     and azimuth. */
  float cos_cap, sin_cap; /**< Circumradius trig for the scanline cap pad. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float sign;             /**< +1 fills the flower, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Flower renders as a filled region. */

  /**
   * @brief Builds a flower from its basis, radius, petal count, and phase.
   * @param b Orientation frame (-v = flower center); retained by reference and
   *          read by every distance() call, so it must outlive the shape.
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of petals (must be >= 3).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  Flower(const math::Basis &b, float radius, int s, float ph,
         bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides >= 3, "SDF Flower: sides must be at least 3");
    HS_CHECK(radius > 0.0f, "SDF Flower: radius must be positive");
    // A shape wider than a hemisphere must be built inverted, about its
    // antipode.
    HS_CHECK(radius <= 1.0f, "SDF Flower: radius exceeds unit sphere");
    sector = math::TWO_PI_F / sides;
    reciprocal_sector = static_cast<float>(sides) / math::TWO_PI_F;
    float outer = radius * (math::PI_F / 2.0f);
    apothem = math::PI_F - outer;
    circumradius = outer;
    antipode = -basis.v;

    CapBounds cb = cap_bounds(antipode, circumradius, invert);
    ny = cb.ny;
    r_val = cb.r_val;
    alpha_angle = cb.alpha_angle;
    phi_min = cb.phi_min;
    phi_max = cb.phi_max;
    cos_cap = cb.cos_radius;
    sin_cap = cb.sin_radius;
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The flower retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  Flower(const math::Basis &&, float, int, float, bool = false) = delete;

  /**
   * @brief Maps the flower's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the flower plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the bounding-circle interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   * @details A petal tip is the shape's radial extreme and its gradient there
   *   is 1, so the fringe reaches one pixel_width past the tip and the cap pad
   *   covers it with no disc clamp. The cap is centred on the antipode, where
   *   the fill sits.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    return emit_padded_cap_row<W, H>(sign, cos_cap, sin_cap, ny, r_val,
                                     alpha_angle, y, out);
  }

  /**
   * @brief Signed distance to the flower, writing into res.
   * @tparam ComputeUVs When true, also computes the normalized scan-distance t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance from the petal edge,
   *        raw_dist = scan distance from the antipode, t =
   *        scan_dist/circumradius when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float scan_dist =
        math::fast_acos(hs::clamp(math::dot(p, antipode), -1.0f, 1.0f));
    float polar = math::PI_F - scan_dist;

    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local = centered_sector_angle(azimuth, sector, reciprocal_sector);

    float dist_edge = polar * math::fast_cosf(local) - apothem;
    float t_val = 0.0f;
    if constexpr (ComputeUVs)
      t_val = scan_dist / circumradius;

    res =
        DistanceResult(sign * -dist_edge, t_val, scan_dist, 0.0f, circumradius);
  }
};

/**
 * @brief Signed distance to a great-circle arc segment of given thickness.
 * @details Register semantics: the DistanceResult table (stroke row: Line).
 *
 * Antipodal endpoints select no arc; the shape degrades to the full great
 * circle through them, bounded and culled as such.
 */
struct Line {
  math::Vector a, b; /**< Arc endpoints (unit vectors). */
  float thickness;   /**< Half-width of the stroke (radians). */

  math::Vector n;         /**< Great-circle plane normal of the arc. */
  float len;              /**< Arc length (radians). */
  float phi_min, phi_max; /**< Precomputed vertical bounds (radians). */

  // Bounding-cap geometry, loop-invariant across scanlines (precomputed in
  // ctor).
  math::Vector mid; /**< Arc midpoint axis (bounding-cap center). */
  float mid_ny = 0.0f, mid_r = 0.0f,
        mid_alpha = 0.0f; /**< Midpoint y, XZ projection length, azimuth. */
  float cap_D_min = 0.0f; /**< Cosine of the bounding-cap radius. */
  bool cap_horiz_valid =
      false; /**< False when the cap yields no horizontal cull: a ~vertical cap
                axis, or antipodal endpoints (whose rendered great circle no cap
                bounds). */
  static constexpr bool is_solid = false; /**< Line renders as a stroke. */

  /**
   * @brief Builds a great-circle arc segment from two endpoints.
   * @param start First endpoint (unit vector).
   * @param end Second endpoint (unit vector).
   * @param th Half-width of the stroke (radians).
   */
  Line(const math::Vector &start, const math::Vector &end, float th)
      : a(start), b(end), thickness(th) {
    // A negative half-width inverts the band, culling every probe.
    HS_CHECK(thickness >= 0.0f, "Line: negative stroke half-width");
    len = math::angle_between(a, b);
    bool antipodal = false;
    math::Vector cr = math::cross(a, b);
    // EPS_CROSS_SQ, not EPS_NORMALIZE_SQ: the bound is on the direction the
    // cross carries, not on whether it normalizes at all. |cross| = sin of
    // the endpoint separation, so 1e-8 names the same band an angular 1e-4
    // does; above it the components' ~1e-7 of rounding leaves the arc plane
    // and the bounding-cap axis under ~2e-3 rad of direction error, a tenth
    // of a pixel at W = 288.
    if (cr.x * cr.x + cr.y * cr.y + cr.z * cr.z < math::EPS_CROSS_SQ) {
      // sin collapses at both ends of the range, so the dot's sign separates
      // them; angle_between cannot, its acos having no resolution there.
      if (math::dot(a, b) < 0.0f) {
        // Antipodal endpoints (len ~ π) leave the arc plane undefined: any great
        // circle through them serves, and distance() then measures the whole
        // circle rather than an arc (every projected point sits at
        // ang_a + ang_b == len). Pick a plane the endpoints actually lie in.
        antipodal = true;
        math::Vector ref =
            (fabsf(a.y) < 0.9f) ? math::Vector(0, 1, 0) : math::Vector(1, 0, 0);
        n = math::cross(a, ref).normalized();
      } else {
        // Coincident endpoints: the arc is the point `a`. len is zeroed so
        // distance() takes the same branch.
        len = 0.0f;
        n = math::Vector(0, 0, 0);
      }
    } else {
      n = cr.normalized();
    }

    // A great-circle arc bulges toward a pole between its endpoints, so its
    // peak latitude can lie outside [phi_a, phi_b]; extend by the interior
    // extremum.
    float phi_a = acosf(std::max(-1.0f, std::min(1.0f, a.y)));
    float phi_b = acosf(std::max(-1.0f, std::min(1.0f, b.y)));
    float phi_lo = std::min(phi_a, phi_b);
    float phi_hi = std::max(phi_a, phi_b);
    // The latitude turns inside the arc iff the forward tangent's y-component
    // (cross(n, p).y) flips sign between endpoints; the extremum is
    // ±sqrt(1-n.y²).
    float t0 = math::cross(n, a).y;
    float t1 = math::cross(n, b).y;
    if ((t0 > 0.0f) != (t1 > 0.0f)) {
      float peak = sqrtf(std::max(0.0f, 1.0f - n.y * n.y));
      float phi_ext = acosf(t0 > 0.0f ? peak : -peak);
      phi_lo = std::min(phi_lo, phi_ext);
      phi_hi = std::max(phi_hi, phi_ext);
    }
    float margin = thickness + BOUNDS_MARGIN;
    phi_min = std::max(0.0f, phi_lo - margin);
    phi_max = std::min(math::PI_F, phi_hi + margin);

    // Bounding cap centered on the segment midpoint. Antipodal endpoints sum to
    // ~0 (no defined midpoint); guard the normalize.
    mid = math::normalized_or(a + b, math::Vector(1, 0, 0));
    mid_ny = mid.y;
    mid_r = sqrtf(mid.x * mid.x + mid.z * mid.z);
    mid_alpha = atan2f(mid.z, mid.x);
    // A cap wider than pi covers the sphere; cos turns back up past pi, so an
    // unclamped radius would report a tighter cap than the arc occupies.
    float cap_radius = std::min(len * 0.5f + thickness, math::PI_F);
    cap_D_min = cosf(cap_radius);
    cap_horiz_valid = mid_r >= MIN_HORIZONTAL_PROJ;

    if (antipodal) {
      // The rendered geometry is the full great circle of `n`, which spans the
      // latitudes between its two extrema and every azimuth: no half-circle cap
      // bounds it, so drop the horizontal cull.
      float peak = sqrtf(std::max(0.0f, 1.0f - n.y * n.y));
      phi_min = std::max(0.0f, acosf(peak) - margin);
      phi_max = std::min(math::PI_F, acosf(-peak) + margin);
      cap_horiz_valid = false;
    }
  }

  /**
   * @brief Maps the arc's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the arc plus thickness margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Signed distance to the arc segment, writing into res.
   * @tparam ComputeUVs Accepted for interface parity; the arc stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = unsigned segment distance minus thickness,
   *        raw_dist = unsigned angular distance to the segment.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    if (len < 1e-6f) {
      float dist = math::angle_between(p, a);
      res = DistanceResult(dist - thickness, 0.0f, dist, 0.0f, thickness);
      return;
    }

    float d_plane = math::dot(p, n);
    math::Vector p_proj_plane = p - n * d_plane;
    float proj_mag = p_proj_plane.magnitude();

    if (proj_mag < 1e-6f) {
      float d_a = math::angle_between(p, a);
      float d_b = math::angle_between(p, b);
      float dist = std::min(d_a, d_b);
      res = DistanceResult(dist - thickness, 0.0f, dist, 0.0f, thickness);
      return;
    }

    math::Vector p_proj = p_proj_plane / proj_mag;
    float ang_a = math::angle_between(a, p_proj);
    float ang_b = math::angle_between(b, p_proj);

    float dist_seg = 0.0f;
    // Two geodesic metrics, C0-continuous at the join: in-segment the foot lies
    // on the arc so asinf(|d_plane|) is exact; past an endpoint the closest
    // point is that cap (asinf(|d_plane|) alone would measure the whole great
    // circle).
    if (ang_a + ang_b <= len + 1e-4f) {
      dist_seg = asinf(std::abs(d_plane));
    } else {
      float d_a = math::angle_between(p, a);
      float d_b = math::angle_between(p, b);
      dist_seg = std::min(d_a, d_b);
    }

    res = DistanceResult(dist_seg - thickness, 0.0f, dist_seg, 0.0f, thickness);
  }

  /**
   * @brief Emits the bounding-cap interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan (also when
   *         the cap has no horizontal projection).
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!cap_horiz_valid)
      return false;

    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    float cos_phi = math::TrigLUT<W, H>::cos_phi[y];
    float sin_phi = math::TrigLUT<W, H>::sin_phi[y];

    return emit_cap_interval<W>(cap_D_min, mid_ny, mid_r, mid_alpha, cos_phi,
                                sin_phi, out);
  }
};

// Leaf roster for the CSG composition contract: a leaf that stops satisfying
// SDFShape fails here rather than at whichever composition happens to use it.
static_assert(SDFShape<PlanarPolygon>);
static_assert(SDFShape<SphericalPolygon>);
static_assert(SDFShape<Star>);
static_assert(SDFShape<Flower>);
static_assert(SDFShape<Line>);

} // namespace SDF
