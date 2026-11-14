/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief White-box access to SDF::PreparedFace state.
 */
#pragma once

#include "core/render/sdf/face.h"

#include <concepts>
#include <type_traits>

namespace hs_test {

/** @brief Face type accepted by PreparedFaceTestAccess, const or mutable. */
template <typename F>
concept PreparedFaceRef =
    std::same_as<std::remove_const_t<F>, SDF::PreparedFace>;

/**
 * @brief White-box accessor for SDF::PreparedFace internals.
 * @details Each accessor returns a reference with the constness of its
 *          argument.
 */
struct PreparedFaceTestAccess {
  /** @brief Projected polygon with its wrap entry. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &poly_2d(F &f) { return f.poly_2d; }
  /** @brief Packed per-edge records. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &packed_edges(F &f) {
    return f.packed_edges;
  }
  /** @brief Convex edge half-planes. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &half_planes(F &f) {
    return f.half_planes;
  }
  /** @brief Cull radius in plane units. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &max_dist(F &f) {
    return f.max_dist;
  }
  /** @brief Whether distances report in plane units. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &linear_dist(F &f) {
    return f.linear_dist;
  }
  /** @brief Whether the sector walk is usable. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &sector_ok(F &f) {
    return f.sector_ok;
  }
  /** @brief Unwrapped sector angle keys. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &sector_keys(F &f) {
    return f.sector_keys;
  }
  /** @brief Sector winding sign. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &sector_sgn(F &f) {
    return f.sector_sgn;
  }
  /** @brief Squared origin-to-boundary bound. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &sector_min_radius_sq(F &f) {
    return f.sector_min_radius_sq;
  }
  /** @brief First row of the bounds. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &y_min(F &f) { return f.y_min; }
  /** @brief Last row of the bounds. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &y_max(F &f) { return f.y_max; }
  /** @brief Latitude mapping of the bounds. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &build_geometry(F &f) {
    return f.build_geometry;
  }
  /** @brief Optional row padding table. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &build_azimuth_pads(F &f) {
    return f.build_azimuth_pads;
  }
  /** @brief Azimuth coverage intervals. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &intervals(F &f) {
    return f.intervals;
  }
  /** @brief Whether the face spans all columns. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &full_width(F &f) {
    return f.full_width;
  }
  /** @brief LUT grid-x coefficient on px. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_ax(F &f) { return f.lut_ax; }
  /** @brief LUT grid-x coefficient on py. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_bx(F &f) { return f.lut_bx; }
  /** @brief LUT grid-x constant term. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_cx(F &f) { return f.lut_cx; }
  /** @brief LUT grid-y coefficient on px. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_ay(F &f) { return f.lut_ay; }
  /** @brief LUT grid-y coefficient on py. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_by(F &f) { return f.lut_by; }
  /** @brief LUT grid-y constant term. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_cy(F &f) { return f.lut_cy; }
  /** @brief LUT grid clamp bound. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_clamp(F &f) {
    return f.lut_clamp;
  }
  /** @brief LUT int16 to plane-unit scale. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_dequant(F &f) {
    return f.lut_dequant;
  }
  /** @brief LUT safe distance in quantized units. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &lut_q_safe(F &f) {
    return f.lut_q_safe;
  }
  /** @brief Squared report radius. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &report_reach_sq(F &f) {
    return f.report_reach_sq;
  }
  /** @brief Scratch buffer the spans view. @param f Face. @return Field. */
  template <PreparedFaceRef F> static auto &scratch_owner(F &f) {
    return f.scratch_owner;
  }

  /**
   * @brief Binds a congruence-class LUT to a bare prepared face.
   * @param f Face receiving the binding.
   * @param lut Canonical-frame LUT.
   * @param canon_xy Canonical centered polygon, x/y pairs.
   * @param vert_offset Cyclic vertex alignment offset.
   * @param reflected True for the mirror family.
   * @return Whether the LUT was bound.
   */
  static bool bind_class_lut(SDF::PreparedFace &f, const SDF::ClassLut *lut,
                             const float *canon_xy, int vert_offset,
                             bool reflected) {
    return f.bind_class_lut(lut, canon_xy, vert_offset, reflected);
  }
};

} // namespace hs_test
