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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "core/platform/platform.h"
#include "core/math/3dmath.h"

/**
 * @brief Lenses that remap a sphere direction to another sphere direction.
 */
namespace lenses {

/** @brief Azimuthal rotation in radians applied per unit of height. */
inline constexpr float TWIST_RATE = 3.0f;
/** @brief Signed mirror distance treated as "on the chamber wall". */
inline constexpr float POLYHEDRAL_MIRROR_EPS = 1e-6f;
/** @brief Loop passes a chamber fold may take before it is declared stuck. */
inline constexpr int POLYHEDRAL_REFLECTION_LIMIT = 24;

/** @brief Inward mirror normals of the tetrahedral (*332) chamber. */
inline constexpr std::array<math::Vector, 3> TETRAHEDRAL_MIRRORS = {
    math::Vector(1.0f, 0.0f, 0.0f), math::Vector(-0.5f, 0.8660254038f, 0.0f),
    math::Vector(0.0f, -0.5773502692f, 0.8164965809f)};
/** @brief Inward mirror normals of the octahedral (*432) chamber. */
inline constexpr std::array<math::Vector, 3> OCTAHEDRAL_MIRRORS = {
    math::Vector(1.0f, 0.0f, 0.0f),
    math::Vector(-0.7071067812f, 0.7071067812f, 0.0f),
    math::Vector(0.0f, -0.7071067812f, 0.7071067812f)};
/** @brief Inward mirror normals of the dodecahedral (*532) chamber. */
inline constexpr std::array<math::Vector, 3> DODECAHEDRAL_MIRRORS = {
    math::Vector(1.0f, 0.0f, 0.0f),
    math::Vector(-0.8090169944f, 0.3090169944f, -0.5f),
    math::Vector(0.0f, 0.0f, 1.0f)};
/** @brief Inward mirror normals of the triangular-prism chamber. */
inline constexpr std::array<math::Vector, 3> TRIANGULAR_PRISM_MIRRORS = {
    math::Vector(0.0f, 1.0f, 0.0f), math::Vector(0.0f, 0.0f, 1.0f),
    math::Vector(0.8660254038f, 0.0f, -0.5f)};
/** @brief Inward mirror normals of the square-prism chamber. */
inline constexpr std::array<math::Vector, 3> SQUARE_PRISM_MIRRORS = {
    math::Vector(0.0f, 1.0f, 0.0f), math::Vector(0.0f, 0.0f, 1.0f),
    math::Vector(0.7071067812f, 0.0f, -0.7071067812f)};
/** @brief Inward mirror normals of the pentagonal-prism chamber. */
inline constexpr std::array<math::Vector, 3> PENTAGONAL_PRISM_MIRRORS = {
    math::Vector(0.0f, 1.0f, 0.0f), math::Vector(0.0f, 0.0f, 1.0f),
    math::Vector(0.5877852523f, 0.0f, -0.8090169944f)};
/** @brief Inward mirror normals of the hexagonal-prism chamber. */
inline constexpr std::array<math::Vector, 3> HEXAGONAL_PRISM_MIRRORS = {
    math::Vector(0.0f, 1.0f, 0.0f), math::Vector(0.0f, 0.0f, 1.0f),
    math::Vector(0.5f, 0.0f, -0.8660254038f)};
/** @brief Inward mirror normals of the octagonal-prism chamber. */
inline constexpr std::array<math::Vector, 3> OCTAGONAL_PRISM_MIRRORS = {
    math::Vector(0.0f, 1.0f, 0.0f), math::Vector(0.0f, 0.0f, 1.0f),
    math::Vector(0.3826834324f, 0.0f, -0.9238795325f)};

/**
 * @brief Applies a trig-free glitch lens to a sphere direction.
 * @param v Unit direction vector on the sphere.
 * @return Direction after latitude doubling and azimuth tripling; returns
 * the up vector near the lens axis.
 */
inline math::Vector glitch_lens(const math::Vector &v) {
  const float x2 = v.x * v.x;
  const float z2 = v.z * v.z;
  const float radius2 = x2 + z2;
  constexpr float MIN_AXIS_RADIUS2 = 1e-6f;
  if (radius2 < MIN_AXIS_RADIUS2)
    return math::Vector(0.0f, 1.0f, 0.0f);

  const float inverse_radius2 = 1.0f / radius2;
  const float double_y = 2.0f * v.y;
  return math::Vector(double_y * v.x * (4.0f * x2 * inverse_radius2 - 3.0f),
                      2.0f * v.y * v.y - 1.0f,
                      double_y * v.z * (3.0f - 4.0f * z2 * inverse_radius2));
}

/**
 * @brief Rotates a direction about the Y axis in proportion to its height.
 * @param v Unit direction on the sphere.
 * @param rate Rotation in radians per unit height.
 * @return The direction rotated by rate * v.y radians.
 * @note The fast trigonometric approximation can move a unit input slightly
 *       off unit length; the result is not renormalized.
 */
inline math::Vector twist_lens(const math::Vector &v, float rate = TWIST_RATE) {
  const float angle = rate * v.y;
  const float c = math::fast_cosf(angle);
  const float s = math::fast_sinf(angle);
  return math::Vector(v.x * c - v.z * s, v.y, v.x * s + v.z * c);
}

/**
 * @brief Folds a direction's azimuth into one twelfth-turn (30 degree) wedge.
 * @param v Unit direction on the sphere.
 * @return A symmetry-equivalent direction inside the wedge.
 */
HS_FLASH_MEMBER inline math::Vector kaleidoscope_lens(const math::Vector &v) {
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
  return math::Vector(x, v.y, z);
}

/**
 * @brief Folds a direction into a spherical reflection-group chamber.
 * @param v Unit direction on the sphere.
 * @param mirrors Inward unit normals of the chamber's simple mirrors.
 * @return A symmetry-equivalent direction inside the chamber.
 */
template <size_t N>
HS_O3_FN inline math::Vector
polyhedral_kaleidoscope_lens(math::Vector v,
                             const std::array<math::Vector, N> &mirrors) {
  [[maybe_unused]] uint32_t reflections = 0;
  for (int reflection = 0; reflection < POLYHEDRAL_REFLECTION_LIMIT;
       ++reflection) {
    bool inside = true;
    for (const math::Vector &normal : mirrors) {
      const float distance = math::dot(v, normal);
      if (distance >= -POLYHEDRAL_MIRROR_EPS)
        continue;
      v.x -= 2.0f * distance * normal.x;
      v.y -= 2.0f * distance * normal.y;
      v.z -= 2.0f * distance * normal.z;
      HS_PULLBACK_PROJECTION_COUNT(++reflections);
      inside = false;
      break;
    }
    if (inside) {
      HS_PULLBACK_PROJECTION_COUNT(
          ++hs::g_pullback_projection_counts.polyhedral_pixels);
      HS_PULLBACK_PROJECTION_COUNT(
          hs::g_pullback_projection_counts.polyhedral_reflections +=
          reflections);
      HS_PULLBACK_PROJECTION_COUNT(
          hs::g_pullback_projection_counts.polyhedral_max_reflections =
              std::max(
                  hs::g_pullback_projection_counts.polyhedral_max_reflections,
                  reflections));
      return v;
    }
  }
  HS_CHECK(false, "polyhedral kaleidoscope fold did not converge");
  return v;
}

/**
 * @brief Applies the tetrahedral fold to a dodecahedral chamber direction.
 * @param v Unit direction inside the dodecahedral chamber.
 * @return The direction folded into the tetrahedral chamber.
 */
HS_HOT_INLINE inline math::Vector
tetrahedral_from_dodecahedral_chamber(math::Vector v) {
  constexpr math::Vector NORMAL = TETRAHEDRAL_MIRRORS[2];
  constexpr math::Vector SECOND = TETRAHEDRAL_MIRRORS[1];
  constexpr math::Vector OBLIQUE = DODECAHEDRAL_MIRRORS[1];
  constexpr float X_FROM_Y = -OBLIQUE.y / OBLIQUE.x;
  constexpr float X_FROM_Z = -OBLIQUE.z / OBLIQUE.x;
  constexpr float SECOND_Y =
      SECOND.y * (1.0f - 2.0f * NORMAL.y * NORMAL.y) + SECOND.x * X_FROM_Y;
  constexpr float SECOND_Z =
      -2.0f * SECOND.y * NORMAL.y * NORMAL.z + SECOND.x * X_FROM_Z;
  static_assert(SECOND_Y > 0.0f && SECOND_Z > 0.0f);
  // Height bound of the EPS-expanded dodecahedral chamber's other tetrahedral wall.
  constexpr float MIN_SAFE_Y = POLYHEDRAL_MIRROR_EPS *
                               (1.0f + SECOND.x / OBLIQUE.x + SECOND_Z) /
                               SECOND_Y;
  if (v.y <= MIN_SAFE_Y)
    return polyhedral_kaleidoscope_lens(v, TETRAHEDRAL_MIRRORS);
  const float DISTANCE = math::dot(v, NORMAL);
  if (DISTANCE < -POLYHEDRAL_MIRROR_EPS) {
    v.y -= 2.0f * DISTANCE * NORMAL.y;
    v.z -= 2.0f * DISTANCE * NORMAL.z;
  }
  HS_PULLBACK_PROJECTION_COUNT(
      ++hs::g_pullback_projection_counts.polyhedral_pixels);
  HS_PULLBACK_PROJECTION_COUNT(
      hs::g_pullback_projection_counts.polyhedral_reflections +=
      DISTANCE < -POLYHEDRAL_MIRROR_EPS ? 1 : 0);
  HS_PULLBACK_PROJECTION_COUNT(
      hs::g_pullback_projection_counts.polyhedral_max_reflections =
          std::max(hs::g_pullback_projection_counts.polyhedral_max_reflections,
                   uint32_t(DISTANCE < -POLYHEDRAL_MIRROR_EPS ? 1 : 0)));
  return v;
}

/**
 * @brief Folds a direction into the dodecahedral chamber.
 * @details Specialization of polyhedral_kaleidoscope_lens for
 * DODECAHEDRAL_MIRRORS.
 * @param v Unit direction on the sphere.
 * @return A symmetry-equivalent direction inside the chamber.
 */
HS_O3_FN __attribute__((always_inline)) inline math::Vector
dodecahedral_kaleidoscope_lens(math::Vector v) {
  [[maybe_unused]] uint32_t reflections = 0;
  constexpr math::Vector OBLIQUE = DODECAHEDRAL_MIRRORS[1];
  for (int remaining = POLYHEDRAL_REFLECTION_LIMIT; remaining > 0;
       --remaining) {
    if (!(v.x >= -POLYHEDRAL_MIRROR_EPS)) {
      v.x = -v.x;
      HS_PULLBACK_PROJECTION_COUNT(++reflections);
      if (--remaining == 0)
        break;
    }
    const float distance = math::dot(v, OBLIQUE);
    if (!(distance >= -POLYHEDRAL_MIRROR_EPS)) {
      v.x -= 2.0f * distance * OBLIQUE.x;
      v.y -= 2.0f * distance * OBLIQUE.y;
      v.z -= 2.0f * distance * OBLIQUE.z;
    } else if (!(v.z >= -POLYHEDRAL_MIRROR_EPS)) {
      v.z = -v.z;
    } else {
      HS_PULLBACK_PROJECTION_COUNT(
          ++hs::g_pullback_projection_counts.polyhedral_pixels);
      HS_PULLBACK_PROJECTION_COUNT(
          hs::g_pullback_projection_counts.polyhedral_reflections +=
          reflections);
      HS_PULLBACK_PROJECTION_COUNT(
          hs::g_pullback_projection_counts.polyhedral_max_reflections =
              std::max(
                  hs::g_pullback_projection_counts.polyhedral_max_reflections,
                  reflections));
      return v;
    }
    HS_PULLBACK_PROJECTION_COUNT(++reflections);
  }
  HS_CHECK(false, "dodecahedral kaleidoscope fold did not converge");
  return v;
}

/** @brief Canonical reflection-group assets and cube classifiers. */
namespace DodecahedralTable {
inline constexpr int FACE_STEPS = 64; ///< Classifier cells per cube-face axis.
inline constexpr float WALL_MARGIN = 1e-5f; ///< Strict chamber-interior margin.
static_assert(WALL_MARGIN > POLYHEDRAL_MIRROR_EPS);
using Matrix = std::array<std::array<float, 3>, 3>; ///< Linear group transform.
/** @brief Finite closure of the canonical chamber mirrors. */
struct Group {
  std::array<Matrix, 120> values{}; ///< Transforms in discovery order.
  uint8_t count = 0;                ///< Admitted transform count.
};
/** @brief Unit reflection-group transform.
 * @return Identity matrix. */
constexpr Matrix identity() {
  return {{{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}}};
}
/** @brief Reflects each transform column through one chamber mirror.
 * @param matrix Transform to reflect.
 * @param mirror Chamber-mirror index.
 * @return Reflected transform. */
constexpr Matrix reflect(Matrix matrix, int mirror) {
  const auto normal = DODECAHEDRAL_MIRRORS[mirror];
  const float n[] = {normal.x, normal.y, normal.z};
  for (int column = 0; column < 3; ++column) {
    const float distance = matrix[0][column] * n[0] + matrix[1][column] * n[1] +
                           matrix[2][column] * n[2];
    for (int row = 0; row < 3; ++row)
      matrix[row][column] -= 2.0f * distance * n[row];
  }
  return matrix;
}
/** @brief Compares transform coefficients within the group-generation tolerance.
 * @param a First transform.
 * @param b Second transform.
 * @return Whether all coefficients match. */
constexpr bool same(const Matrix &a, const Matrix &b) {
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column) {
      const float delta = a[row][column] - b[row][column];
      if (delta < -2e-5f || delta > 2e-5f)
        return false;
    }
  return true;
}
/** @brief Closes the chamber mirrors under reflection.
 * @return Group transforms in discovery order. */
constexpr Group make_group() {
  Group result;
  result.values[result.count++] = identity();
  for (int index = 0; index < result.count; ++index)
    for (int mirror = 0; mirror < 3; ++mirror) {
      const auto candidate = reflect(result.values[index], mirror);
      bool found = false;
      for (int other = 0; other < result.count; ++other)
        if (same(candidate, result.values[other])) {
          found = true;
          break;
        }
      if (found)
        continue;
      if (result.count == result.values.size()) {
        ++result.count;
        return result;
      }
      result.values[result.count++] = candidate;
    }
  return result;
}
/** @brief Canonical finite reflection group. */
inline constexpr auto
    GROUP HS_PROGMEM_UNIQUE(dodecahedral_fold_group) = make_group();
static_assert(GROUP.count == 120);
/** @brief Resolves each mirror transition from one group transform.
 * @param index Source group index.
 * @return Destination index for each chamber mirror. */
constexpr auto make_transition_row(size_t index) {
  std::array<uint8_t, 3> result{255, 255, 255};
  for (int mirror = 0; mirror < 3; ++mirror) {
    const auto candidate = reflect(GROUP.values[index], mirror);
    for (int next = 0; next < GROUP.count; ++next)
      if (same(candidate, GROUP.values[next])) {
        result[mirror] = static_cast<uint8_t>(next);
        break;
      }
  }
  return result;
}
/** @brief Cached mirror transitions for one group transform.
 * @tparam Index Source group index. */
template <size_t Index>
inline constexpr auto TRANSITION_ROW = make_transition_row(Index);
/** @brief Assembles the cached transition rows.
 * @tparam Indices Source group indices.
 * @return Group-index by mirror-index transitions. */
template <size_t... Indices>
constexpr auto make_transitions(std::index_sequence<Indices...>) {
  return std::array{TRANSITION_ROW<Indices>...};
}
/** @brief Mirror transitions between canonical group transforms. */
inline constexpr auto TRANSITIONS =
    make_transitions(std::make_index_sequence<GROUP.count>{});
static_assert([] {
  for (const auto &row : TRANSITIONS)
    for (uint8_t index : row)
      if (index >= GROUP.count)
        return false;
  return true;
}());
/** @brief Maps cube-face chart coordinates to a direction.
 * @param face Cube-face index.
 * @param u Horizontal chart coordinate.
 * @param v Vertical chart coordinate.
 * @return Unnormalized direction on the selected face. */
constexpr std::array<float, 3> face_direction(int face, float u, float v) {
  switch (face) {
  case 0:
    return {1.0f, v, u};
  case 1:
    return {-1.0f, v, -u};
  case 2:
    return {u, 1.0f, v};
  case 3:
    return {u, -1.0f, -v};
  case 4:
    return {u, v, 1.0f};
  default:
    return {-u, v, -1.0f};
  }
}
/** @brief Resolves the ordered canonical fold of one chart direction.
 * @param point Unnormalized chart direction.
 * @return Group index, or 255 when the fold exceeds its iteration limit. */
constexpr uint8_t classify(std::array<float, 3> point) {
  uint8_t matrix_index = 0;
  for (int step = 0; step < POLYHEDRAL_REFLECTION_LIMIT; ++step) {
    int mirror;
    const auto normal = DODECAHEDRAL_MIRRORS[1];
    const float distance =
        point[0] * normal.x + point[1] * normal.y + point[2] * normal.z;
    if (point[0] < -POLYHEDRAL_MIRROR_EPS)
      mirror = 0;
    else if (distance < -POLYHEDRAL_MIRROR_EPS)
      mirror = 1;
    else if (point[2] < -POLYHEDRAL_MIRROR_EPS)
      mirror = 2;
    else
      return matrix_index;
    matrix_index = TRANSITIONS[matrix_index][mirror];
    if (mirror == 0)
      point[0] = -point[0];
    else if (mirror == 2)
      point[2] = -point[2];
    else {
      point[0] -= 2.0f * distance * normal.x;
      point[1] -= 2.0f * distance * normal.y;
      point[2] -= 2.0f * distance * normal.z;
    }
  }
  return 255;
}
/** @brief Evaluates one cube-face classifier row.
 * @param face Cube-face index.
 * @param y Row index.
 * @return Classifier indices in horizontal order. */
constexpr auto make_row(int face, int y) {
  std::array<uint8_t, FACE_STEPS> result{};
  for (int x = 0; x < FACE_STEPS; ++x)
    result[x] =
        classify(face_direction(face, (x + 0.5f) * (2.0f / FACE_STEPS) - 1.0f,
                                (y + 0.5f) * (2.0f / FACE_STEPS) - 1.0f));
  return result;
}
/** @brief Independently evaluated classifier row.
 * @tparam Face Cube-face index.
 * @tparam Row Vertical row index. */
template <int Face, size_t Row> inline constexpr auto ROW = make_row(Face, Row);
/** @brief Flattens cached classifier rows into one face.
 * @tparam Face Cube-face index.
 * @tparam Rows Vertical row indices.
 * @return Row-major classifier map. */
template <int Face, size_t... Rows>
constexpr auto make_face(std::index_sequence<Rows...>) {
  const std::array rows{ROW<Face, Rows>...};
  std::array<uint8_t, FACE_STEPS * FACE_STEPS> result{};
  for (int y = 0; y < FACE_STEPS; ++y)
    for (int x = 0; x < FACE_STEPS; ++x)
      result[y * FACE_STEPS + x] = rows[y][x];
  return result;
}
/** @brief Evaluates one complete cube-face map.
 * @tparam Face Cube-face index.
 * @return Row-major classifier map. */
template <int Face> constexpr auto make_face() {
  return make_face<Face>(std::make_index_sequence<FACE_STEPS>{});
}
/** @brief Classifier map on the +X cube face. */
inline constexpr auto
    FACE_0 HS_PROGMEM_UNIQUE(dodecahedral_fold_face0) = make_face<0>();
/** @brief Classifier map on the -X cube face. */
inline constexpr auto
    FACE_1 HS_PROGMEM_UNIQUE(dodecahedral_fold_face1) = make_face<1>();
/** @brief Classifier map on the +Y cube face. */
inline constexpr auto
    FACE_2 HS_PROGMEM_UNIQUE(dodecahedral_fold_face2) = make_face<2>();
/** @brief Classifier map on the -Y cube face. */
inline constexpr auto
    FACE_3 HS_PROGMEM_UNIQUE(dodecahedral_fold_face3) = make_face<3>();
/** @brief Classifier map on the +Z cube face. */
inline constexpr auto
    FACE_4 HS_PROGMEM_UNIQUE(dodecahedral_fold_face4) = make_face<4>();
/** @brief Classifier map on the -Z cube face. */
inline constexpr auto
    FACE_5 HS_PROGMEM_UNIQUE(dodecahedral_fold_face5) = make_face<5>();
/** @brief Classifier maps in cube-face orientation order. */
inline constexpr std::array<const uint8_t *, 6> FACES HS_PROGMEM_UNIQUE(
    dodecahedral_fold_faces) = {FACE_0.data(), FACE_1.data(), FACE_2.data(),
                                FACE_3.data(), FACE_4.data(), FACE_5.data()};
} // namespace DodecahedralTable

/** @brief Attempts a chamber-validated reflection matrix from a coarse cube cell.
 * @param input Unit direction.
 * @param result Receives the folded direction on success.
 * @param matrices Nonnull immutable group matrix array.
 * @param faces Optional face-major classifiers; null uses the flash maps.
 * @return Whether the candidate lies strictly inside the canonical chamber. */
HS_HOT_INLINE inline bool
try_dodecahedral_table(const math::Vector &input, math::Vector &result,
                       const DodecahedralTable::Matrix *matrices,
                       const uint8_t *faces) {
  const float ax = fabsf(input.x), ay = fabsf(input.y), az = fabsf(input.z);
  float major, u, v;
  int face;
  if (ax >= ay && ax >= az) {
    major = ax;
    face = input.x >= 0.0f ? 0 : 1;
    u = input.x >= 0.0f ? input.z : -input.z;
    v = input.y;
  } else if (ay >= az) {
    major = ay;
    face = input.y >= 0.0f ? 2 : 3;
    u = input.x;
    v = input.y >= 0.0f ? input.z : -input.z;
  } else {
    major = az;
    face = input.z >= 0.0f ? 4 : 5;
    u = input.z >= 0.0f ? input.x : -input.x;
    v = input.y;
  }
  if (!(major > 0.0f))
    return false;
  const float inverse = 1.0f / major;
  constexpr int N = DodecahedralTable::FACE_STEPS;
  const int x =
      std::min(static_cast<int>((u * inverse + 1.0f) * (0.5f * N)), N - 1);
  const int y =
      std::min(static_cast<int>((v * inverse + 1.0f) * (0.5f * N)), N - 1);
  const uint8_t index = faces ? faces[face * N * N + y * N + x]
                              : DodecahedralTable::FACES[face][y * N + x];
  if (index >= DodecahedralTable::GROUP.count)
    return false;
  const auto &m = matrices[index];
  result = {m[0][0] * input.x + m[0][1] * input.y + m[0][2] * input.z,
            m[1][0] * input.x + m[1][1] * input.y + m[1][2] * input.z,
            m[2][0] * input.x + m[2][1] * input.y + m[2][2] * input.z};
  constexpr float MARGIN = DodecahedralTable::WALL_MARGIN;
  return result.x > MARGIN && result.z > MARGIN &&
         math::dot(result, DODECAHEDRAL_MIRRORS[1]) > MARGIN;
}

/** @brief Attempts a validated transform with flash-resident classifiers.
 * @param input Unit direction.
 * @param result Receives the folded direction on success.
 * @param matrices Nonnull immutable group matrix array.
 * @return Whether the candidate lies inside the canonical chamber. */
HS_HOT_INLINE inline bool
try_dodecahedral_table(const math::Vector &input, math::Vector &result,
                       const DodecahedralTable::Matrix *matrices) {
  return try_dodecahedral_table(input, result, matrices, nullptr);
}

/** @brief Attempts the canonical flash-resident reflection matrix.
 * @param input Unit direction.
 * @param result Receives the folded direction on success.
 * @return Whether the candidate lies inside the canonical chamber. */
HS_HOT_INLINE inline bool try_dodecahedral_table(const math::Vector &input,
                                                 math::Vector &result) {
  return try_dodecahedral_table(input, result,
                                DodecahedralTable::GROUP.values.data());
}

/** @brief Folds by a validated matrix, using ordered reflections near walls.
 * @param input Unit direction.
 * @param matrices Nonnull immutable group matrix array.
 * @param faces Optional face-major classifiers; null uses the flash maps.
 * @return The symmetry-equivalent chamber direction. */
HS_HOT_INLINE inline math::Vector
dodecahedral_table_lens(const math::Vector &input,
                        const DodecahedralTable::Matrix *matrices,
                        const uint8_t *faces) {
  math::Vector result;
  return try_dodecahedral_table(input, result, matrices, faces)
             ? result
             : dodecahedral_kaleidoscope_lens(input);
}

/** @brief Folds with explicit matrices and flash-resident classifiers.
 * @param input Unit direction.
 * @param matrices Nonnull immutable group matrix array.
 * @return The symmetry-equivalent chamber direction. */
HS_HOT_INLINE inline math::Vector
dodecahedral_table_lens(const math::Vector &input,
                        const DodecahedralTable::Matrix *matrices) {
  return dodecahedral_table_lens(input, matrices, nullptr);
}

/** @brief Folds with the canonical flash-resident matrices.
 * @param input Unit direction.
 * @return The symmetry-equivalent chamber direction. */
HS_HOT_INLINE inline math::Vector
dodecahedral_table_lens(const math::Vector &input) {
  return dodecahedral_table_lens(input, DodecahedralTable::GROUP.values.data());
}

} // namespace lenses
