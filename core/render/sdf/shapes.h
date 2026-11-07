/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include "core/math/geometry.h"
#include "core/platform/constants.h"
#include "core/render/sdf/common.h"

/**
 * @file shapes.h
 * @brief The polygon, star, flower and line leaves.
 */

namespace SDF {

/** @brief Cap bounds shared by solid leaves.
 * @tparam Derived Leaf carrying cap and sign. */
template <typename Derived> struct CapBoundedLeaf {
  /** @brief Maps the cap latitude band to inclusive rows.
   * @tparam H Canvas height.
   * @return Padded row bounds. */
  template <int H> Bounds get_vertical_bounds() const {
    const auto &cap = static_cast<const Derived &>(*this).cap;
    return phi_bounds_to_rows<H>(cap.phi_min, cap.phi_max);
  }

  /** @brief Copies the leaf with bounds covering an angular sweep.
   * @param arc Largest point displacement in radians.
   * @return Leaf with expanded cap bounds. */
  Derived swept(float arc) const {
    const auto &leaf = static_cast<const Derived &>(*this);
    Derived out = leaf;
    out.cap = swept_cap(leaf.cap, arc);
    return out;
  }

  /** @brief Emits one padded bounding-cap row interval.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @tparam OutputIt Interval sink.
   * @param y Row index.
   * @param out Sink accepting start and end.
   * @return True for a handled row; false requests a full scan. */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    const auto &leaf = static_cast<const Derived &>(*this);
    const auto &cap = leaf.cap;
    return emit_padded_cap_row<W, H>(leaf.sign, cap.cos_radius, cap.sin_radius,
                                     cap.ny, cap.r_val, cap.alpha_angle, y,
                                     out);
  }
};

/**
 * @brief Calculates signed distance to a planar polygon.
 * @details Register semantics: the DistanceResult table (row: PlanarPolygon).
 */
struct PlanarPolygon : CapBoundedLeaf<PlanarPolygon> {
  /** Orientation frame (v = polygon axis); retained by reference, so
   * it must outlive the shape. */
  const math::Basis &basis;
  float circumradius; /**< Angular radius from center to vertex (radians). */
  int sides;          /**< Number of polygon sides. */
  float phase;        /**< Azimuth phase offset (radians). */
  float sector;       /**< Angular width of one polygon sector. */
  float reciprocal_sector; /**< Reciprocal angular sector width. */
  float apothem;           /**< Precomputed inradius (radians). */
  CapBounds cap; /**< Circumscribed cap: axis projection, row band, radius. */
  float sign;    /**< +1 fills the polygon, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a planar polygon from its basis, radius, side count, phase.
   * @param b Orientation frame (v = polygon axis); retained by reference and
   *          read by every distance() call, so it must outlive the shape.
   * @param radius Circumradius as a fraction of the hemisphere.
   * @param s Number of polygon sides (must be >= 3).
   * @param ph Azimuth phase offset (radians), added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
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
    HS_CHECK(radius <= 1.0f,
             "SDF PlanarPolygon: radius exceeds a hemisphere (build inverted)");
    circumradius = radius * (math::PI_F / 2.0f);
    sector = math::TWO_PI_F / sides;
    reciprocal_sector = static_cast<float>(sides) / math::TWO_PI_F;
    apothem = circumradius * cosf(math::PI_F / sides);

    cap = cap_bounds(basis.v, circumradius, invert);
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The polygon retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  PlanarPolygon(const math::Basis &&, float, int, float, bool = false) = delete;

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
   *       (gradient cos(PI/sides) there), and over-reporting beyond the disc.
   *       The disc term `polar - circumradius` is a lower bound of the true
   *       distance; the max of the two is not, so it is for scanline shading,
   *       not a march-safe metric.
   * @note The circumscribed cap includes the disc-clamped AA fringe.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float polar =
        math::fast_acos(hs::clamp(math::dot(p, basis.v), -1.0f, 1.0f));

    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local = centered_sector_angle(azimuth, sector, reciprocal_sector);

    float dist_edge =
        fmaxf(polar * math::fast_cosf(local) - apothem, polar - circumradius);
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
 * O(1) per-pixel distance to the edge great circle. Past a vertex it is a
 * lower bound; see distance().
 * Register semantics: the DistanceResult table (row: SphericalPolygon).
 */
struct SphericalPolygon : CapBoundedLeaf<SphericalPolygon> {
  const math::Basis &basis; /**< Orientation frame (v = polygon axis); retained
                              by reference, so it must outlive the shape. */
  int sides;                /**< Number of polygon sides. */
  float phase;              /**< Azimuth phase offset (radians). */
  float sector;             /**< Angular width of one polygon sector. */
  float reciprocal_sector;  /**< Reciprocal angular sector width. */
  float circumradius; /**< Angular distance from center to vertex (radians). */
  float edge_nv;      /**< Edge normal dotted with the center axis. */
  float edge_nu;      /**< Edge normal dotted with the u-axis. */
  CapBounds cap; /**< Circumscribed cap: axis projection, row band, radius. */
  float sign;    /**< +1 fills the polygon, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a spherical polygon from its basis, radius, side count,
   * phase.
   * @param b Orientation frame (v = polygon axis); retained by reference and
   *          read by every distance() call, so it must outlive the shape.
   * @param radius Polygon radius as a fraction of the hemisphere.
   * @param s Number of polygon sides (must be >= 3).
   * @param ph Azimuth phase offset (radians), added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
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
    HS_CHECK(
        radius <= 1.0f,
        "SDF SphericalPolygon: radius exceeds a hemisphere (build inverted)");
    sector = math::TWO_PI_F / sides;
    reciprocal_sector = static_cast<float>(sides) / math::TWO_PI_F;
    circumradius = radius * (math::PI_F / 2.0f);

    // Canonical edge: vertices at azimuth ±π/n from the sector bisector
    // (u-axis), at angular distance circumradius.
    float half_step = math::PI_F / sides;
    float sin_r = sinf(circumradius);
    float cos_r = cosf(circumradius);
    float cos_hs = cosf(half_step);
    float sin_hs = sinf(half_step);

    math::Vector v1 =
        basis.v * cos_r + (basis.u * cos_hs + basis.w * sin_hs) * sin_r;
    math::Vector v2 =
        basis.v * cos_r + (basis.u * cos_hs - basis.w * sin_hs) * sin_r;

    math::Vector en = math::cross(v2, v1);
    float len = en.magnitude();
    if (len > 1e-9f) {
      en = en * (1.0f / len);
    } else {
      // Degenerate canonical edge near zero circumradius: substitute the
      // sector bisector as a defined unit normal.
      en = basis.u;
    }
    // Ensure outward: dot(center, n) should be negative
    if (math::dot(en, basis.v) > 0)
      en = -en;

    edge_nv = math::dot(en, basis.v);
    edge_nu = math::dot(en, basis.u);

    cap = cap_bounds(basis.v, circumradius, invert);
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The polygon retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  SphericalPolygon(const math::Basis &&, float, int, float,
                   bool = false) = delete;

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
   * @note The circumscribed cap includes the disc-clamped AA fringe.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float cos_p = hs::clamp(math::dot(p, basis.v), -1.0f, 1.0f);
    float polar = math::fast_acos(cos_p);

    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local = centered_sector_angle(azimuth, sector, reciprocal_sector);

    // Angular distance to the nearest great circle edge via precomputed normal
    // cos(local) is even, so sector folding works automatically
    float sin_p = sqrtf(fmaxf(0.0f, 1.0f - cos_p * cos_p));
    float dp = edge_nv * cos_p + edge_nu * math::fast_cosf(local) * sin_p;
    float dist_edge =
        fmaxf(asinf(hs::clamp(dp, -1.0f, 1.0f)), polar - circumradius);

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
    float sin_p = sqrtf(fmaxf(0.0f, 1.0f - cos_p * cos_p));

    float azimuth = basis_azimuth(p, basis.u, basis.w, phase);

    float local = centered_sector_angle(azimuth, sector, reciprocal_sector);
    float dp = edge_nv * cos_p + edge_nu * math::fast_cosf(local) * sin_p;
    float disc = sin_p * cap.cos_radius - cos_p * cap.sin_radius;
    return sign * fmaxf(hs::clamp(dp, -cos_p, 1.0f), disc);
  }
};

/**
 * @brief Calculates signed distance to a star shape.
 * @details Register semantics: the DistanceResult table (row: Star).
 */
struct Star : CapBoundedLeaf<Star> {
  const math::Basis &basis; /**< Orientation frame (v = star axis); retained by
                              reference, so it must outlive the shape. */
  int sides;                /**< Number of star points. */
  float phase;              /**< Azimuth phase offset (radians). */
  float sector;             /**< Angular width of one star sector. */
  float reciprocal_sector;  /**< Reciprocal angular sector width. */
  static constexpr bool is_solid =
      true; /**< Star renders as a filled region. */

  float edge_nx;      ///< Edge plane normal x-component.
  float edge_ny;      ///< Edge plane normal y-component.
  float plane_d;      ///< Edge plane offset.
  float circumradius; /**< Angular radius from center to point tip (radians). */
  CapBounds cap; /**< Circumscribed cap: axis projection, row band, radius. */
  float sign;    /**< +1 fills the star, -1 its complement. */

  /**
   * @brief Builds a star from its basis, radius, point count, and phase.
   * @param b Orientation frame (v = star axis); retained by reference and read
   *          by every distance() call, so it must outlive the shape.
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of star points (must be >= 3).
   * @param ph Azimuth phase offset (radians), added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  Star(const math::Basis &b, float radius, int s, float ph, bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides >= 3, "SDF Star: sides must be at least 3");
    HS_CHECK(radius > 0.0f, "SDF Star: radius must be positive");
    // arc_stretch<Star> = 2 holds only within a hemisphere; a wider shape must
    // be built inverted, about its antipode.
    HS_CHECK(radius <= 1.0f,
             "SDF Star: radius exceeds a hemisphere (build inverted)");
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

    cap = cap_bounds(basis.v, outer_radius, invert);
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The star retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  Star(const math::Basis &&, float, int, float, bool = false) = delete;

  /**
   * @brief Cap bounding the star, for bounds re-centred on another axis.
   * @param axis Receives the cap axis.
   * @param radius Receives the cap's angular radius (radians).
   * @return False for a complement star, which no cap bounds.
   */
  bool bounding_cap(math::Vector &axis, float &radius) const {
    axis = basis.v;
    radius = circumradius;
    return sign > 0.0f;
  }

  /**
   * @brief Bounds proxy covering the star turned onto another axis and swept
   *        by an arc.
   * @param axis Axis the star's centre is turned onto (unit).
   * @param arc Largest arc any point of the turned star moves (radians).
   * @return Copy whose bounds and intervals cover the turned, swept star.
   */
  Star moved(const math::Vector &axis, float arc) const {
    Star out = *this;
    out.cap =
        moved_cap(cap, circumradius, basis.v, axis, arc + MOVED_AXIS_SLACK);
    return out;
  }

  /**
   * @brief Signed distance to the star, writing into res.
   * @tparam ComputeUVs When true, also stores the normalized azimuth t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance to the nearest point edge,
   *        raw_dist = polar distance from center, t = normalized azimuth when
   *        ComputeUVs (0 otherwise).
   * @note The folded edge half-plane has radial gradient |edge_nx| at a tip;
   *       the circumscribed-disc distance `scan_dist - circumradius` is the
   *       tighter bound there. The edge term is a chart-plane distance and
   *       over-reports beyond the disc. The disc term is a lower bound of the
   *       true distance; the max of the two is not, so it is for scanline
   *       shading, not a march-safe metric.
   * @note The circumscribed cap includes the disc-clamped AA fringe.
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

    float dist_edge = fmaxf(-(px * edge_nx + py * edge_ny + plane_d),
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
struct Flower : CapBoundedLeaf<Flower> {
  /** Orientation frame (-v = flower center); retained by reference, so
   * it must outlive the shape. */
  const math::Basis &basis;
  int sides;               /**< Number of petals. */
  float phase;             /**< Azimuth phase offset (radians). */
  float sector;            /**< Angular width of one flower sector. */
  float reciprocal_sector; /**< Reciprocal angular sector width. */
  float circumradius;      /**< Angular radius from the antipode to petal tip
                         (radians). */
  float apothem;           /**< Petal inradius offset (PI - outer radius). */
  math::Vector antipode;   /**< Antipode of the flower axis (scan origin). */
  CapBounds cap; /**< Cap about the antipode: projection, row band, radius. */
  float sign;    /**< +1 fills the flower, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Flower renders as a filled region. */

  /**
   * @brief Builds a flower from its basis, radius, petal count, and phase.
   * @param b Orientation frame (-v = flower center); retained by reference and
   *          read by every distance() call, so it must outlive the shape.
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of petals (must be >= 3).
   * @param ph Azimuth phase offset (radians), added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
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
    HS_CHECK(radius <= 1.0f,
             "SDF Flower: radius exceeds a hemisphere (build inverted)");
    sector = math::TWO_PI_F / sides;
    reciprocal_sector = static_cast<float>(sides) / math::TWO_PI_F;
    float outer = radius * (math::PI_F / 2.0f);
    apothem = math::PI_F - outer;
    circumradius = outer;
    antipode = -basis.v;

    cap = cap_bounds(antipode, circumradius, invert);
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The flower retains its basis by reference and reads it in every
   * distance() call, so binding a temporary would leave those reads dangling.
   */
  Flower(const math::Basis &&, float, int, float, bool = false) = delete;

  /**
   * @brief Signed distance to the flower, writing into res.
   * @tparam ComputeUVs When true, also computes the normalized scan-distance t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance from the petal edge,
   *        raw_dist = scan distance from the antipode, t =
   *        scan_dist/circumradius when ComputeUVs.
   * @note The antipodal cap includes the unit-gradient petal-tip AA fringe.
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
  math::Vector a;  ///< First arc endpoint (unit vector).
  math::Vector b;  ///< Second arc endpoint (unit vector).
  float thickness; /**< Half-width of the stroke (radians). */

  math::Vector n;      /**< Great-circle plane normal of the arc. */
  math::Vector edge_a; ///< Inward arc-start half-plane normal.
  math::Vector edge_b; ///< Inward arc-end half-plane normal.
  float len;           /**< Arc length (radians). */
  float phi_min;       ///< Lower colatitude bound (radians).
  float phi_max;       ///< Upper colatitude bound (radians).

  // Bounding-cap geometry, precomputed in the ctor.
  math::Vector mid;       /**< Arc midpoint axis (bounding-cap center). */
  float mid_ny = 0.0f;    ///< Midpoint axis y-component.
  float mid_r = 0.0f;     ///< Length of the midpoint axis' XZ projection.
  float mid_alpha = 0.0f; ///< Azimuth of the midpoint axis' XZ projection.
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
    // |cross| = sin of the endpoint separation; EPS_CROSS_SQ bounds the arc
    // plane's direction error, not just whether the cross normalizes.
    if (cr.x * cr.x + cr.y * cr.y + cr.z * cr.z < math::EPS_CROSS_SQ) {
      // sin collapses at both ends of the range; the dot's sign separates them.
      if (math::dot(a, b) < 0.0f) {
        // Antipodal endpoints (len ~ π) leave the arc plane undefined: any great
        // circle through them serves, and distance() then measures the whole
        // circle rather than an arc. Pick a plane through the endpoints.
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
    edge_a = math::cross(n, a);
    edge_b = math::cross(b, n);

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
    float t0 = edge_a.y;
    float t1 = -edge_b.y;
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
      len = 2.0f * math::PI_F;
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
   * @brief Bounds proxy covering the arc swept by an arc.
   * @param arc Largest arc any point of the segment moves (radians).
   * @return Copy whose latitude band and bounding cap are grown by @p arc;
   *         its distance() is the unswept segment's.
   */
  Line swept(float arc) const {
    Line out = *this;
    out.phi_min = std::max(0.0f, phi_min - arc);
    out.phi_max = std::min(math::PI_F, phi_max + arc);
    out.cap_D_min = cosf(std::min(len * 0.5f + thickness + arc, math::PI_F));
    return out;
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

    float dist_seg = 0.0f;
    if (len > math::PI_F ||
        (math::dot(p, edge_a) >= 0.0f && math::dot(p, edge_b) >= 0.0f)) {
      const float D_PLANE = math::dot(p, n);
      dist_seg = asinf(fminf(fabsf(D_PLANE), 1.0f));
    } else {
      float d_a = math::angle_between(p, a);
      float d_b = math::angle_between(p, b);
      dist_seg = fminf(d_a, d_b);
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

// Leaf roster for the CSG composition contract.
static_assert(SDFShape<PlanarPolygon>);
static_assert(SDFShape<SphericalPolygon>);
static_assert(SDFShape<Star>);
static_assert(SDFShape<Flower>);
static_assert(SDFShape<Line>);

} // namespace SDF
