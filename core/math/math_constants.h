/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file math_constants.h @brief Scalar mathematical constants and geometric tolerances. */

#include "core/platform/platform.h"

namespace math {

/** Exact 2^96 lift for divisors whose squared magnitude underflows the normal range. */
inline constexpr float COMPLEX_UNDERFLOW_LIFT = 0x1p96f;

/**
 * @brief The Golden Ratio constant (Phi).
 */
inline constexpr float PHI = 1.6180339887498948482045868343656f;
/**
 * @brief The inverse of the Golden Ratio (1/Phi).
 */
inline constexpr float INV_PHI = 1 / PHI;
/** @brief Generic floating-point comparison tolerance. */
inline constexpr float TOLERANCE = 1e-4f;
/// Geometric near-equality (positions, angles).
inline constexpr float EPS_GEOMETRIC = 1e-5f;
/// Squared edge length below which an edge is degenerate.
inline constexpr float EPS_LEN_SQ = 1e-6f;
/**
 * @brief Squared |cross(from, to)| below which make_rotation synthesizes a π
 *        turn instead of using the half-vector form.
 * @details Placed where the half-vector form's angular error
 * (~2*ULP(1)/|cross|) meets the synthesized turn's (|cross|).
 */
inline constexpr float EPS_ANTIPARALLEL_SQ = 4e-7f;
/// Squared cross-product magnitude below which it is degenerate.
inline constexpr float EPS_CROSS_SQ = 1e-8f;
/// Squared face-normal magnitude below which the face is degenerate.
inline constexpr float EPS_NORMAL_SQ = 1e-9f;
/// Squared length below which normalize() has no reliable direction.
inline constexpr float EPS_NORMALIZE_SQ = 1e-12f;
/** @brief Squared length below which a blend of two directions has cancelled. */
inline constexpr float EPS_BLEND_LEN_SQ = 1e-8f;
/// Slack on |q|^2 - 1 for is-unit quaternion assertions.
inline constexpr float EPS_UNIT_QUAT_SQ = 0.01f;
/// Slack on |v|^2 - 1 for is-unit vector assertions.
inline constexpr float EPS_UNIT_VEC_SQ = 0.02f;
/**
 * @brief Cosine above which a vector is treated as parallel to a reference axis.
 */
inline constexpr float COS_AXIS_PARALLEL = 1.0f - TOLERANCE;
/**
 * @brief Floating-point representation of PI.
 */
inline constexpr float PI_F = static_cast<float>(PI);
/**
 * @brief A full turn in radians (2π).
 * @details Suffixed _F because Arduino's wiring.h defines TWO_PI as a macro.
 */
inline constexpr float TWO_PI_F = 2.0f * PI_F;

} // namespace math
