/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cassert>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <array>
#include "math/geometry.h"
#include "render/shading.h"
#include "render/clip.h"
#include "render/filter/splat.h"
#include "memory.h"

/**
 * @file cull.h
 * @brief Edge samplers and the screen row/column span and clip-cull kernel the
 * curve rasterizer walks each segment through.
 */

namespace Plot {

/** @brief Arc length between unit sphere points, preserving short chords. */
__attribute__((always_inline)) inline float
unit_arc_length(const math::Vector &a, const math::Vector &b) {
  const math::Vector chord = a - b;
  const float length_sq = math::dot(chord, chord);
  // Below 1e-3 radians the chord differs from the arc by less than 4.2e-8 relative.
  return length_sq < 1e-6f ? sqrtf(length_sq) : math::angle_between(a, b);
}

/**
 * @brief Geodesic segment shorter than this (radians) collapses to a point.
 * @details A slerp-axis stability bound that picks the interpolation strategy,
 * independent of math::EPS_GEOMETRIC.
 */
inline constexpr float EPS_GEODESIC_SEGMENT = 0.001f;

/**
 * @brief Minimum |cross(a, b)|² for which the arc pole of a geodesic edge is
 *        taken from the cross product rather than a stable perpendicular.
 * @details Bounds the quantity the pole normalization consumes. |cross| =
 * sin(angle), so 1e-8 is the angular 1e-4 band without angle_between's ULP
 * amplification near a normalized dot of ±1.
 */
inline constexpr float EPS_ARC_POLE_SQ = 1e-8f;

/**
 * @brief Minimum |axis.y| for which a geodesic edge's endpoint columns bound
 *        its azimuth span.
 * @details Below this the great circle runs near the poles, where longitude is
 * ill-conditioned and the interior leaves the endpoint columns.
 */
inline constexpr float AXIS_Y_EPS = 1e-4f;

/**
 * @brief Minimum worst-case sin(φ) over a curve for which its plotted columns
 *        are trusted enough to cull by.
 * @details The azimuth Lipschitz bound scales as 1/sin(φ), and nearer the poles
 * the plotted column is float noise; below this the column cull bails.
 */
inline constexpr float MIN_SIN_PHI = 0.05f;

/**
 * @brief Floor on the adaptive sub-step length, as a fraction of base_step.
 * @details Lower clamp on the screen-velocity step, which tends to zero where
 * the azimuthal velocity diverges at the poles.
 */
inline constexpr float MIN_POLE_SCALE = 0.05f;

/**
 * @brief Target screen-space spacing (pixels) between adaptive sub-samples.
 * @details Slightly sub-pixel so neighbouring AntiAlias splats overlap without
 * holes.
 */
inline constexpr float SCREEN_STEP_PX = 0.9f;

/**
 * @brief Columns of slack added on each side of a culled column span.
 * @details Absorbs plot rounding and the AntiAlias tap spread.
 */
inline constexpr int COL_PAD = 2;

/**
 * @brief Columns a padded span reaches past its fractional end.
 * @details The pad plus the boundary column ceil() adds.
 */
inline constexpr int COL_FOOTPRINT = COL_PAD + 1;

/**
 * @brief Rows of slack added to the high end of a geodesic row span.
 * @details AntiAlias emits a sample into floor(row) and floor(row) + 1.
 */
inline constexpr float GEODESIC_ROW_AA_PAD = 1.0f;

/**
 * @brief Columns outside the render band a clip cut is placed at.
 * @details A piece ending on the band edge still overlaps it after the
 * COL_FOOTPRINT widening; the extra column absorbs fast-trig error in the cut.
 */
inline constexpr int CLIP_CUT_COL_PAD = COL_FOOTPRINT + 1;

/**
 * @brief Rows outside the render band a clip cut is placed at.
 * @details could_intersect_y takes the band edge itself as an intersection, and
 * the row map's fast-trig round trip moves a cut by a fraction of a row.
 */
inline constexpr int CLIP_CUT_ROW_PAD = 1;

#include "render/plot/cull/samplers.h"
#include "render/plot/cull/spans.h"
#include "render/plot/cull/pipeline.h"
#include "render/plot/cull/visibility.h"
} // namespace Plot
