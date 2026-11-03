/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file tangent_field.h
 * @brief Conformal generators and polynomial tangent fields on the unit sphere.
 */

#include <array>
#include <concepts>
#include <cstdint>

#include "core/math/lenses.h"
#include "core/math/mobius.h"
#include "core/math/spherical_harmonics.h"
#include "core/math/tangent_field_constants.h"

namespace math::tangent_field {

/** @brief Unit-sphere tangent velocity with uniform speed and derivative bounds.
 * @tparam F Field value type. */
template <typename F>
concept TangentField = requires(const F &field, const Vector &p) {
  { field.velocity(p) } -> std::same_as<Vector>;
  { field.speed_bound() } -> std::convertible_to<float>;
  { field.lipschitz_bound() } -> std::convertible_to<float>;
};

/** @brief Converts a sphere rotation and boost to a traceless complex matrix.
 * @param rotation Right-handed rotation vector.
 * @param boost Potential-gradient boost vector.
 * @return Complex generator matrix. */
inline MobiusGenerator to_mobius_generator(const Vector &rotation,
                                           const Vector &boost) {
  return {{0.5f * boost.y, -0.5f * rotation.y},
          {0.5f * (boost.x + rotation.z), 0.5f * (boost.z - rotation.x)},
          {0.5f * (boost.x - rotation.z), -0.5f * (boost.z + rotation.x)}};
}

/** @brief Rotation plus potential-gradient conformal field. */
struct ConformalGenerator {
  Vector rotation{0.0f, 0.0f, 0.0f}; ///< Right-handed rotation vector.
  Vector boost{0.0f, 0.0f, 0.0f};    ///< Gradient boost toward its axis.

  /** @brief Evaluates the conformal velocity.
   * @param p Unit direction.
   * @return Tangent velocity. */
  HS_FLASH_MEMBER Vector velocity(const Vector &p) const {
    return cross(rotation, p) + boost - p * dot(boost, p);
  }
  /** @brief Uniform speed bound.
   * @return Sum of rotation and boost lengths. */
  float speed_bound() const { return rotation.length() + boost.length(); }
  /** @brief Uniform covariant derivative bound.
   * @return Sum of rotation and boost lengths. */
  float lipschitz_bound() const { return speed_bound(); }
  /** @brief Converts to the complex matrix representation.
   * @return Traceless complex generator. */
  MobiusGenerator matrix() const {
    return to_mobius_generator(rotation, boost);
  }
};

/** @brief Converts a complex matrix to sphere rotation and boost vectors.
 * @param generator Traceless complex generator.
 * @return Equivalent sphere generator. */
inline ConformalGenerator
from_mobius_generator(const MobiusGenerator &generator) {
  const Complex B = generator.beta + generator.gamma.conjugate();
  const Complex W =
      Complex(0.0f, 1.0f) * (generator.beta - generator.gamma.conjugate());
  return {{W.re, -2.0f * generator.alpha.im, W.im},
          {B.re, 2.0f * generator.alpha.re, B.im}};
}

/** @brief Scalar polynomial family. */
enum class TermType : uint8_t {
  MODE,        ///< Real spherical harmonic mode.
  TETRAHEDRAL, ///< Tetrahedral cubic invariant.
  OCTAHEDRAL,  ///< Octahedral quartic invariant.
  DODECAHEDRAL ///< Dodecahedral sextic invariant.
};
/** @brief Tangent assembly from a scalar potential. */
enum class Kind : uint8_t {
  TOROIDAL, ///< Cross product with the outward normal.
  POLOIDAL  ///< Tangent projection of the gradient.
};

/** @brief Authored scalar term and its tangent-field amplitude. */
struct PotentialTerm {
  TermType type = TermType::MODE; ///< Polynomial family.
  Kind kind = Kind::TOROIDAL;     ///< Tangent assembly.
  uint8_t l = 1;                  ///< Harmonic degree in [0, 8].
  uint8_t m = 0;                  ///< Nonnegative harmonic order.
  float amplitude = 0.0f;         ///< Signed peak-speed amplitude.
  float phase = 0.0f;             ///< Harmonic azimuthal phase in radians.
};

/** @brief Extracts unnormalised half-turn axes from the tetrahedral mirrors.
 * @return Three columns of the half-turn matrices plus identity. */
constexpr auto make_tetrahedral_axes() {
  constexpr auto GROUP =
      lenses::DodecahedralTable::make_group(lenses::TETRAHEDRAL_MIRRORS);
  static_assert(GROUP.count == 24);
  std::array<Vector, 3> axes{};
  size_t count = 0;
  for (size_t index = 0; index < GROUP.count; ++index) {
    const auto &matrix = GROUP.values[index];
    const float trace = matrix[0][0] + matrix[1][1] + matrix[2][2];
    if (trace < -1.00002f || trace > -0.99998f)
      continue;
    bool symmetric = true;
    for (int row = 0; row < 3; ++row)
      for (int column = 0; column < 3; ++column) {
        const float delta = matrix[row][column] - matrix[column][row];
        if (delta < -2e-5f || delta > 2e-5f)
          symmetric = false;
      }
    if (!symmetric)
      continue;
    Vector best{0.0f, 0.0f, 0.0f};
    float best_square = 0.0f;
    for (int column = 0; column < 3; ++column) {
      const Vector candidate(matrix[0][column] + (column == 0 ? 1.0f : 0.0f),
                             matrix[1][column] + (column == 1 ? 1.0f : 0.0f),
                             matrix[2][column] + (column == 2 ? 1.0f : 0.0f));
      const float square = dot(candidate, candidate);
      if (square > best_square) {
        best = candidate;
        best_square = square;
      }
    }
    if (count < axes.size())
      axes[count++] = best;
  }
  return std::pair{axes, count};
}

/** @brief Tetrahedral invariant axes derived from the chamber's half-turns. */
inline constexpr auto TETRAHEDRAL_AXES = make_tetrahedral_axes().first;
static_assert(make_tetrahedral_axes().second == 3);

/** @brief Ambient invariant value and gradient. */
struct PotentialSample {
  float value;     ///< Scalar value.
  Vector gradient; ///< Ambient Cartesian gradient.
};

/** @brief Samples a fold-aligned invariant polynomial.
 * @param type Invariant polynomial family.
 * @param p Cartesian point in the fold frame.
 * @return Scalar value and analytic gradient. */
inline PotentialSample invariant_sample(TermType type, const Vector &p) {
  if (type == TermType::TETRAHEDRAL) {
    const float a = dot(TETRAHEDRAL_AXES[0], p);
    const float b = dot(TETRAHEDRAL_AXES[1], p);
    const float c = dot(TETRAHEDRAL_AXES[2], p);
    return {a * b * c, TETRAHEDRAL_AXES[0] * (b * c) +
                           TETRAHEDRAL_AXES[1] * (a * c) +
                           TETRAHEDRAL_AXES[2] * (a * b)};
  }
  const float x2 = p.x * p.x, y2 = p.y * p.y, z2 = p.z * p.z;
  if (type == TermType::OCTAHEDRAL)
    return {x2 * x2 + y2 * y2 + z2 * z2,
            {4.0f * p.x * x2, 4.0f * p.y * y2, 4.0f * p.z * z2}};
  constexpr float PHI2 = PHI * PHI;
  const float a = PHI2 * x2 - y2, b = PHI2 * y2 - z2, c = PHI2 * z2 - x2;
  return {a * b * c,
          {2.0f * p.x * b * (PHI2 * c - a), 2.0f * p.y * c * (PHI2 * a - b),
           2.0f * p.z * a * (PHI2 * b - c)}};
}

/** @brief One polynomial term with cached scale and phase. */
struct PreparedTerm {
  PotentialTerm term{};      ///< Authored term.
  float scale = 0.0f;        ///< Amplitude times peak normalization.
  Complex phase{1.0f, 0.0f}; ///< Prepared azimuthal phase.
  float lipschitz = 0.0f;    ///< Amplitude-weighted covariant bound.

  /** @brief Evaluates the normalized ambient gradient.
   * @param p Unit direction.
   * @return Cartesian gradient including amplitude. */
  Vector gradient(const Vector &p) const {
    return term.type == TermType::MODE ? SHMath::spherical_harmonic_gradient(
                                             term.l, term.m, p, scale, phase)
                                       : sample(p).gradient;
  }

  /** @brief Samples the normalized scalar term.
   * @param p Unit direction.
   * @return Scalar value and ambient gradient including amplitude. */
  PotentialSample sample(const Vector &p) const {
    if (term.type != TermType::MODE) {
      const bool rotated = phase.re != 1.0f || phase.im != 0.0f;
      const Vector local = rotated
                               ? Vector(phase.re * p.x - phase.im * p.z, p.y,
                                        phase.im * p.x + phase.re * p.z)
                               : p;
      const auto raw = invariant_sample(term.type, local);
      const Vector gradient =
          rotated
              ? Vector(phase.re * raw.gradient.x + phase.im * raw.gradient.z,
                       raw.gradient.y,
                       -phase.im * raw.gradient.x + phase.re * raw.gradient.z)
              : raw.gradient;
      return {scale * raw.value, gradient * scale};
    }
    Complex power(1.0f, 0.0f);
    for (int i = 0; i < term.m; ++i)
      power = power * Complex(p.x, p.z);
    return {
        scale * SHMath::reduced_legendre(term.l, term.m, p.y) *
            (phase * power).re,
        SHMath::spherical_harmonic_gradient(term.l, term.m, p, scale, phase)};
  }
};

/** @brief Prepares a scalar term's normalization and azimuthal phase.
 * @param term Authored scalar term.
 * @param drift_phase Accumulated azimuthal drift in radians.
 * @return Cached scalar term. */
inline PreparedTerm prepare_term(const PotentialTerm &term,
                                 float drift_phase = 0.0f) {
  PreparedTerm result;
  result.term = term;
  TermConstants constants;
  if (term.type == TermType::MODE) {
    HS_CHECK(term.l <= 8 && term.m <= term.l,
             "tangent field: invalid harmonic mode");
    constants = MODE_CONSTANTS[term.l * (term.l + 1) / 2 + term.m];
    result.scale = term.amplitude * SHMath::harmonic_scale(term.l, term.m) *
                   constants.scale;
    const float angle = term.phase + term.m * drift_phase;
    result.phase = Complex(cosf(angle), -sinf(angle));
  } else {
    constants = INVARIANT_CONSTANTS[static_cast<unsigned>(term.type) - 1];
    result.scale = term.amplitude * constants.scale;
    if (drift_phase != 0.0f)
      result.phase = Complex(cosf(drift_phase), -sinf(drift_phase));
  }
  result.lipschitz = fabsf(term.amplitude) * constants.lipschitz;
  return result;
}

/** @brief Sum of at most four normalized polynomial tangent terms. */
class PotentialField {
public:
  std::array<PreparedTerm, 4> terms{}; ///< Cached scalar terms.

  /** @brief Prepares terms, a shared drift phase and a local-to-world rotation.
   * @param authored Authored scalar terms.
   * @param drift_phase Accumulated azimuthal phase in radians.
   * @param orientation Unit local-to-world quaternion. */
  explicit PotentialField(const std::array<PotentialTerm, 4> &authored = {},
                          float drift_phase = 0.0f,
                          Quaternion orientation = Quaternion(1.0f, 0.0f, 0.0f,
                                                              0.0f))
      : orientation(orientation),
        oriented(orientation.r != 1.0f ||
                 dot(orientation.v, orientation.v) != 0.0f) {
    for (size_t i = 0; i < terms.size(); ++i) {
      terms[i] = prepare_term(authored[i], drift_phase);
      bound += fabsf(authored[i].amplitude);
      derivative_bound += terms[i].lipschitz;
    }
  }

  /** @brief Evaluates the assembled tangent velocity.
   * @param p Unit world direction.
   * @return Tangent velocity in world coordinates. */
  HS_FLASH_MEMBER Vector velocity(const Vector &p) const {
    const Vector local = oriented ? rotate(p, orientation.conjugate()) : p;
    Vector result(0.0f, 0.0f, 0.0f);
    for (const auto &term : terms) {
      if (term.scale == 0.0f)
        continue;
      const Vector gradient = term.gradient(local);
      result += term.term.kind == Kind::TOROIDAL
                    ? cross(local, gradient)
                    : gradient - local * dot(local, gradient);
    }
    return oriented ? rotate(result, orientation) : result;
  }

  /** @brief Uniform speed bound.
   * @return Sum of absolute amplitudes. */
  float speed_bound() const { return bound; }
  /** @brief Uniform covariant derivative bound.
   * @return Sum of amplitude-weighted term bounds. */
  float lipschitz_bound() const { return derivative_bound; }

private:
  Quaternion orientation;        ///< Local-to-world rotation.
  bool oriented;                 ///< Whether the orientation is nonidentity.
  float bound = 0.0f;            ///< Uniform speed bound.
  float derivative_bound = 0.0f; ///< Uniform derivative bound.
};

/** @brief Conservative strength at which a fixed Euler field step can fold.
 * @param speed Uniform speed bound.
 * @param lipschitz Uniform covariant derivative bound.
 * @return Conservative fold threshold, or infinity for a zero field. */
inline float fold_strength(float speed, float lipschitz) {
  const float denominator = hypotf(lipschitz, sqrtf(2.0f) * speed) + lipschitz;
  return denominator > 0.0f ? 2.0f / denominator
                            : std::numeric_limits<float>::infinity();
}

static_assert(TangentField<PotentialField> && TangentField<ConformalGenerator>);

} // namespace math::tangent_field
