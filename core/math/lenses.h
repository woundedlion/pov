/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file lenses.h
 * @brief Direction-domain sphere lenses: the glitch fold, the axial twist, the
 *        hexagonal kaleidoscope wedge, and the spherical reflection-group folds
 *        with their chamber mirror tables.
 */

#include <algorithm>
#include <array>

#include "platform/platform.h"
#include "math/3dmath.h"

/**
 * @brief Lenses that remap a sphere direction to another sphere direction.
 */
namespace lenses {

/** @brief Azimuthal rotation in radians applied per unit of height. */
inline constexpr float TWIST_RATE = 3.0f;
/** @brief Signed mirror distance treated as "on the chamber wall". */
inline constexpr float POLYHEDRAL_MIRROR_EPS = 1e-6f;
/** @brief Loop passes a chamber fold may take before it is declared stuck.
 *  @details The dodecahedral chamber is the worst case: H3's 15 positive roots
 *  bound the fold at 15 reflections, plus the pass that confirms the direction
 *  is inside. The excess is slack for a near-wall direction the mirror epsilon
 *  does not catch; a converging fold exits at its own pass count. */
inline constexpr int POLYHEDRAL_REFLECTION_LIMIT = 24;

/** @brief Inward mirror normals of the tetrahedral (*332) chamber. */
inline constexpr std::array<Vector, 3> TETRAHEDRAL_MIRRORS = {
    Vector(1.0f, 0.0f, 0.0f), Vector(-0.5f, 0.8660254038f, 0.0f),
    Vector(0.0f, -0.5773502692f, 0.8164965809f)};
/** @brief Inward mirror normals of the octahedral (*432) chamber. */
inline constexpr std::array<Vector, 3> OCTAHEDRAL_MIRRORS = {
    Vector(1.0f, 0.0f, 0.0f), Vector(-0.7071067812f, 0.7071067812f, 0.0f),
    Vector(0.0f, -0.7071067812f, 0.7071067812f)};
/** @brief Inward mirror normals of the dodecahedral (*532) chamber. */
inline constexpr std::array<Vector, 3> DODECAHEDRAL_MIRRORS = {
    Vector(1.0f, 0.0f, 0.0f), Vector(-0.8090169944f, 0.3090169944f, -0.5f),
    Vector(0.0f, 0.0f, 1.0f)};
/** @brief Inward mirror normals of the triangular-prism chamber. */
inline constexpr std::array<Vector, 3> TRIANGULAR_PRISM_MIRRORS = {
    Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
    Vector(0.8660254038f, 0.0f, -0.5f)};
/** @brief Inward mirror normals of the square-prism chamber. */
inline constexpr std::array<Vector, 3> SQUARE_PRISM_MIRRORS = {
    Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
    Vector(0.7071067812f, 0.0f, -0.7071067812f)};
/** @brief Inward mirror normals of the pentagonal-prism chamber. */
inline constexpr std::array<Vector, 3> PENTAGONAL_PRISM_MIRRORS = {
    Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
    Vector(0.5877852523f, 0.0f, -0.8090169944f)};
/** @brief Inward mirror normals of the hexagonal-prism chamber. */
inline constexpr std::array<Vector, 3> HEXAGONAL_PRISM_MIRRORS = {
    Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
    Vector(0.5f, 0.0f, -0.8660254038f)};
/** @brief Inward mirror normals of the octagonal-prism chamber. */
inline constexpr std::array<Vector, 3> OCTAGONAL_PRISM_MIRRORS = {
    Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
    Vector(0.3826834324f, 0.0f, -0.9238795325f)};

/**
 * @brief Applies a trig-free glitch lens to a sphere direction.
 * @param v Unit direction vector on the sphere.
 * @return Direction after latitude doubling and azimuth tripling; returns
 * the up vector near the lens axis.
 */
inline Vector glitch_lens(const Vector &v) {
  const float x2 = v.x * v.x;
  const float z2 = v.z * v.z;
  const float radius2 = x2 + z2;
  constexpr float MIN_AXIS_RADIUS2 = 1e-6f;
  if (radius2 < MIN_AXIS_RADIUS2)
    return Vector(0.0f, 1.0f, 0.0f);

  const float inverse_radius2 = 1.0f / radius2;
  const float double_y = 2.0f * v.y;
  return Vector(double_y * v.x * (4.0f * x2 * inverse_radius2 - 3.0f),
                2.0f * v.y * v.y - 1.0f,
                double_y * v.z * (3.0f - 4.0f * z2 * inverse_radius2));
}

/**
 * @brief Rotates a direction about the Y axis in proportion to its height.
 * @param v Unit direction on the sphere.
 * @return The direction rotated by TWIST_RATE * v.y radians.
 */
inline Vector twist_lens(const Vector &v) {
  const float angle = TWIST_RATE * v.y;
  const float c = fast_cosf(angle);
  const float s = fast_sinf(angle);
  return Vector(v.x * c - v.z * s, v.y, v.x * s + v.z * c);
}

/**
 * @brief Folds a direction's azimuth into one twelfth-turn (30 degree) wedge.
 * @param v Unit direction on the sphere.
 * @return A symmetry-equivalent direction inside the wedge.
 */
HS_FLASH_MEMBER inline Vector kaleidoscope_lens(const Vector &v) {
  constexpr float COS_60 = 0.5f;
  constexpr float SIN_60 = 0.8660254037844386f;
  constexpr float SQRT_3 = 1.7320508075688772f;
  float x = fabsf(v.x);
  float z = fabsf(v.z);
  if (z > SQRT_3 * x) {
    const float reflected_x = -COS_60 * x + SIN_60 * z;
    z = SIN_60 * x + COS_60 * z;
    x = reflected_x;
  }
  if (SQRT_3 * z > x) {
    const float reflected_x = COS_60 * x + SIN_60 * z;
    z = SIN_60 * x - COS_60 * z;
    x = reflected_x;
  }
  return Vector(x, v.y, z);
}

/**
 * @brief Folds a direction into a spherical reflection-group chamber.
 * @param v Unit direction on the sphere.
 * @param mirrors Inward unit normals of the chamber's simple mirrors.
 * @return A symmetry-equivalent direction inside the chamber.
 */
template <size_t N>
HS_O3_FN inline Vector
polyhedral_kaleidoscope_lens(Vector v, const std::array<Vector, N> &mirrors) {
  [[maybe_unused]] uint32_t reflections = 0;
  for (int reflection = 0; reflection < POLYHEDRAL_REFLECTION_LIMIT;
       ++reflection) {
    bool inside = true;
    for (const Vector &normal : mirrors) {
      const float distance = dot(v, normal);
      if (distance >= -POLYHEDRAL_MIRROR_EPS)
        continue;
      v.x -= 2.0f * distance * normal.x;
      v.y -= 2.0f * distance * normal.y;
      v.z -= 2.0f * distance * normal.z;
      HS_SB_STAGE_COUNT(++reflections);
      inside = false;
      break;
    }
    if (inside) {
      HS_SB_STAGE_COUNT(
          ++hs::g_shader_workbench_stage_cycles.polyhedral_pixels);
      HS_SB_STAGE_COUNT(
          hs::g_shader_workbench_stage_cycles.polyhedral_reflections +=
          reflections);
      HS_SB_STAGE_COUNT(
          hs::g_shader_workbench_stage_cycles.polyhedral_max_reflections =
              std::max(hs::g_shader_workbench_stage_cycles
                           .polyhedral_max_reflections,
                       reflections));
      return v;
    }
  }
  HS_CHECK(false, "polyhedral kaleidoscope fold did not converge");
  return v;
}

/**
 * @brief Folds a direction into the dodecahedral chamber.
 * @details Specializes polyhedral_kaleidoscope_lens for DODECAHEDRAL_MIRRORS:
 * the scan over the mirror array is unrolled and the two axis mirrors collapse
 * into sign flips, leaving only the fold iteration.
 * @param v Unit direction on the sphere.
 * @return A symmetry-equivalent direction inside the chamber.
 */
HS_O3_FN __attribute__((always_inline)) inline Vector
dodecahedral_kaleidoscope_lens(Vector v) {
  [[maybe_unused]] uint32_t reflections = 0;
  constexpr Vector OBLIQUE = DODECAHEDRAL_MIRRORS[1];
  for (int reflection = 0; reflection < POLYHEDRAL_REFLECTION_LIMIT;
       ++reflection) {
    if (v.x < -POLYHEDRAL_MIRROR_EPS) {
      v.x = -v.x;
    } else {
      const float distance = dot(v, OBLIQUE);
      if (distance < -POLYHEDRAL_MIRROR_EPS) {
        v.x -= 2.0f * distance * OBLIQUE.x;
        v.y -= 2.0f * distance * OBLIQUE.y;
        v.z -= 2.0f * distance * OBLIQUE.z;
      } else if (v.z < -POLYHEDRAL_MIRROR_EPS) {
        v.z = -v.z;
      } else {
        HS_SB_STAGE_COUNT(
            ++hs::g_shader_workbench_stage_cycles.polyhedral_pixels);
        HS_SB_STAGE_COUNT(
            hs::g_shader_workbench_stage_cycles.polyhedral_reflections +=
            reflections);
        HS_SB_STAGE_COUNT(
            hs::g_shader_workbench_stage_cycles.polyhedral_max_reflections =
                std::max(hs::g_shader_workbench_stage_cycles
                             .polyhedral_max_reflections,
                         reflections));
        return v;
      }
    }
    HS_SB_STAGE_COUNT(++reflections);
  }
  HS_CHECK(false, "dodecahedral kaleidoscope fold did not converge");
  return v;
}

} // namespace lenses
