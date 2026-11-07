/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <complex>
#include <vector>

#include "core/math/tangent_field.h"
#include "core/math/spherical.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/vec_test_util.h"

namespace hs_test::tangent_field_tests {
namespace TF = math::tangent_field;

/** @brief Deterministic sphere sample including both poles. */
inline math::Vector direction(int index, int count = 1000) {
  const double y = 1.0 - 2.0 * index / (count - 1.0);
  const double radius = std::sqrt(std::max(0.0, 1.0 - y * y));
  const double angle = index * 2.39996322972865332;
  return {float(radius * std::cos(angle)), float(y),
          float(radius * std::sin(angle))};
}

/** @brief Independent double-precision reduced Legendre recurrence. */
inline std::pair<double, double> legendre_reference(int l, int m, double y) {
  double previous = 0, value = 1, previous_derivative = 0, derivative = 0;
  for (int degree = m; degree < l; ++degree) {
    const double next =
        ((2 * degree + 1) * y * value - (degree + m) * previous) /
        (degree + 1 - m);
    const double next_derivative =
        ((2 * degree + 1) * (value + y * derivative) -
         (degree + m) * previous_derivative) /
        (degree + 1 - m);
    previous = value;
    value = next;
    previous_derivative = derivative;
    derivative = next_derivative;
  }
  return {value, derivative};
}

/** @brief Runtime recurrence reproduces the compiled Legendre table. */
inline void test_legendre_coefficient_table() {
  auto *volatile build = &SHMath::make_legendre_coefficients;
  const auto coefficients = build();
  for (size_t index = 0; index < coefficients.size(); ++index) {
    const float expected = SHMath::LEGENDRE_COEFFICIENTS[index];
    HS_EXPECT_NEAR(coefficients[index], expected, 1e-7f * std::abs(expected));
  }
}

inline void test_horner_values_and_phased_gradients() {
  for (int l = 0; l <= 8; ++l)
    for (int m = 0; m <= l; ++m)
      for (int index = 0; index <= 128; ++index) {
        const float y = -1.0f + index / 64.0f;
        const auto expected = legendre_reference(l, m, y);
        HS_EXPECT_NEAR(SHMath::reduced_legendre(l, m, y), expected.first,
                       2e-5 * std::max(1.0, std::abs(expected.first)));
        HS_EXPECT_NEAR(SHMath::reduced_legendre_derivative(l, m, y),
                       expected.second,
                       4e-5 * std::max(1.0, std::abs(expected.second)));
      }
  for (int l = 1; l <= 8; ++l)
    for (int m = 0; m <= l; ++m)
      for (int index = 1; index < 40; ++index) {
        const auto p = direction(index, 41);
        const float norm = SHMath::harmonic_scale(l, m);
        const math::Complex phase(std::cos(0.73f), -std::sin(0.73f));
        const auto gradient =
            SHMath::spherical_harmonic_gradient(l, m, p, norm, phase);
        const auto value = [=](math::Vector v) {
          const auto q = legendre_reference(l, m, v.y).first;
          return double(norm) * q *
                 (std::complex<double>(phase.re, phase.im) *
                  std::pow(std::complex<double>(v.x, v.z), m))
                     .real();
        };
        const math::Vector axes[] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const float components[] = {gradient.x, gradient.y, gradient.z};
        for (int axis = 0; axis < 3; ++axis) {
          const auto plus = p + axes[axis] * 0.001f;
          const auto minus = p - axes[axis] * 0.001f;
          const double denominator = axis == 0   ? plus.x - minus.x
                                     : axis == 1 ? plus.y - minus.y
                                                 : plus.z - minus.z;
          HS_EXPECT_NEAR(components[axis],
                         (value(plus) - value(minus)) / denominator, 0.003);
        }
      }
}

inline void test_mobius_generator_velocity_and_inverse() {
  float negated_rotation_error = 0.0f;
  for (int index = 0; index < 1000; ++index) {
    const auto p = direction(index);
    const TF::ConformalGenerator generator{
        {0.3f * sinf(index * 0.17f), 0.4f, -0.2f},
        {0.2f, -0.1f, 0.3f * cosf(index * 0.23f)}};
    const auto matrix = generator.matrix();
    const auto recovered = TF::from_mobius_generator(matrix);
    HS_EXPECT_VEC(recovered.rotation, generator.rotation, 1e-6f);
    HS_EXPECT_VEC(recovered.boost, generator.boost, 1e-6f);
    const auto velocity =
        (math::mobius_transform(p, math::mobius_exp(matrix, 0.001f)) -
         math::mobius_transform(p, math::mobius_exp(matrix, -0.001f))) /
        0.002f;
    HS_EXPECT_VEC(velocity, generator.velocity(p), 0.001f);
    negated_rotation_error +=
        (velocity -
         TF::ConformalGenerator{-generator.rotation, generator.boost}.velocity(
             p))
            .length();
  }
  HS_EXPECT_GT(negated_rotation_error, 100.0f);
}

inline void test_mobius_exp_log_projective_roundtrip() {
  for (int index = 0; index < 100; ++index) {
    const TF::ConformalGenerator generator{
        {1.1f * sinf(index * 0.17f), 0.8f, -0.6f},
        {0.7f, -0.4f, 0.9f * cosf(index * 0.23f)}};
    const auto matrix = math::mobius_exp(generator.matrix());
    const auto log = math::mobius_log(matrix);
    const auto recovered = TF::from_mobius_generator(log);
    HS_EXPECT_VEC(recovered.rotation, generator.rotation, 2e-5f);
    HS_EXPECT_VEC(recovered.boost, generator.boost, 2e-5f);
  }
  const math::MobiusParams cases[] = {
      {},
      {{1, 0}, {2, -1}, {0, 0}, {1, 0}},
      {{-1, 0}, {0.001f, 0}, {0, 0}, {-1, 0}},
      {{0, 1}, {0, 0}, {0, 0}, {0, -1}},
      {-1.072f, 0.304f, 0.416f, 0.0f, 0.0f, 0.0f, 0.7071f, 0.0f}};
  for (const auto &matrix : cases) {
    const auto recovered = math::mobius_exp(math::mobius_log(matrix));
    for (int index = 0; index < 100; ++index)
      HS_EXPECT_VEC(math::mobius_transform(direction(index, 100), matrix),
                    math::mobius_transform(direction(index, 100), recovered),
                    1e-5f);
  }
}

inline void test_mobius_exp_domain_corners_are_finite() {
  double maximum_relative_error = 0.0;
  for (unsigned corner = 0; corner < 64; ++corner) {
    const auto sign = [=](int bit) {
      return (corner & (1u << bit)) ? 1.0f : -1.0f;
    };
    const auto matrix = math::mobius_exp(TF::to_mobius_generator(
        {sign(0) * math::TWO_PI_F, sign(1) * math::TWO_PI_F,
         sign(2) * math::TWO_PI_F},
        {sign(3) * 6.0f, sign(4) * 6.0f, sign(5) * 6.0f}));
    const auto complex = [](math::Complex z) {
      return std::complex<double>(z.re, z.im);
    };
    const auto ad = complex(matrix.a) * complex(matrix.d);
    const auto bc = complex(matrix.b) * complex(matrix.c);
    const double relative_error =
        std::abs(ad - bc - 1.0) / (1.0 + std::abs(ad) + std::abs(bc));
    maximum_relative_error = std::max(maximum_relative_error, relative_error);
    for (int index = 0; index < 33; ++index) {
      const auto transformed =
          math::mobius_transform(direction(index, 33), matrix);
      HS_EXPECT_TRUE(std::isfinite(transformed.x) &&
                     std::isfinite(transformed.y) &&
                     std::isfinite(transformed.z));
      HS_EXPECT_NEAR(transformed.length(), 1.0f, 2e-6f);
    }
  }
  HS_EXPECT_LT(maximum_relative_error, 2e-6);
}

/** @brief Constructs a single-term prepared field. */
inline TF::PotentialField single_term(TF::PotentialTerm term) {
  return TF::PotentialField({term, {}, {}, {}});
}

inline void test_constant_potentials_have_zero_bounds() {
  for (auto kind : {TF::Kind::TOROIDAL, TF::Kind::POLOIDAL}) {
    const TF::PotentialTerm constant{TF::TermType::MODE, kind, 0, 0, -4.0f};
    const auto field = single_term(constant);
    HS_EXPECT_EQ(field.speed_bound(), 0.0f);
    HS_EXPECT_EQ(field.lipschitz_bound(), 0.0f);
    HS_EXPECT_TRUE(std::isinf(
        TF::fold_strength(field.speed_bound(), field.lipschitz_bound())));
    HS_EXPECT_VEC(field.velocity(direction(7)), math::Vector(0, 0, 0), 0.0f);
    const TF::PotentialTerm active{TF::TermType::MODE, kind, 1, 0, 0.5f};
    const TF::PotentialField mixed({constant, active, {}, {}});
    const auto active_field = single_term(active);
    HS_EXPECT_EQ(mixed.speed_bound(), 0.5f);
    HS_EXPECT_EQ(mixed.lipschitz_bound(), active_field.lipschitz_bound());
    HS_EXPECT_VEC(mixed.velocity(direction(7)),
                  active_field.velocity(direction(7)), 0.0f);
  }
}

inline void test_potential_divergence_and_curl() {
  const TF::PotentialTerm terms[] = {
      {TF::TermType::MODE, TF::Kind::TOROIDAL, 8, 5, 1.0f, 0.73f},
      {TF::TermType::TETRAHEDRAL, TF::Kind::TOROIDAL, 0, 0, 1.0f},
      {TF::TermType::OCTAHEDRAL, TF::Kind::TOROIDAL, 0, 0, 1.0f},
      {TF::TermType::DODECAHEDRAL, TF::Kind::TOROIDAL, 0, 0, 1.0f}};
  for (auto term : terms)
    for (TF::Kind kind : {TF::Kind::TOROIDAL, TF::Kind::POLOIDAL}) {
      term.kind = kind;
      const auto field = single_term(term);
      for (int index = 1; index < 100; ++index) {
        const auto p = direction(index, 101);
        const auto e = math::cross(p, math::Y_AXIS).normalized();
        const auto f = math::cross(p, e);
        const auto derivative = [&](math::Vector axis) {
          const auto plus = (p * cosf(0.001f) + axis * sinf(0.001f));
          const auto minus = (p * cosf(0.001f) - axis * sinf(0.001f));
          return (field.velocity(plus) - field.velocity(minus)) / 0.002f;
        };
        const auto de = derivative(e), df = derivative(f);
        HS_EXPECT_NEAR(math::dot(p, field.velocity(p)), 0.0f, 2e-6f);
        const float quantity = kind == TF::Kind::TOROIDAL
                                   ? math::dot(e, de) + math::dot(f, df)
                                   : math::dot(f, de) - math::dot(e, df);
        HS_EXPECT_NEAR(quantity, 0.0f, 0.002f);
      }
    }
}

/** @brief Exact sphere exponential. */
inline math::Vector exponential(math::Vector p, math::Vector tangent) {
  const float length = tangent.length();
  return length > 1e-8f ? p * cosf(length) + tangent * (sinf(length) / length)
                        : p;
}

/** @brief Applies a reflection-group matrix. */
inline math::Vector apply(const lenses::DodecahedralTable::Matrix &matrix,
                          math::Vector p) {
  return {matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2] * p.z,
          matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2] * p.z,
          matrix[2][0] * p.x + matrix[2][1] * p.y + matrix[2][2] * p.z};
}

inline void test_invariants_and_vector_fields_follow_fold_groups() {
  const std::array<std::array<math::Vector, 3>, 8> mirrors = {
      lenses::TETRAHEDRAL_MIRRORS,     lenses::OCTAHEDRAL_MIRRORS,
      lenses::DODECAHEDRAL_MIRRORS,    lenses::TRIANGULAR_PRISM_MIRRORS,
      lenses::SQUARE_PRISM_MIRRORS,    lenses::PENTAGONAL_PRISM_MIRRORS,
      lenses::HEXAGONAL_PRISM_MIRRORS, lenses::OCTAGONAL_PRISM_MIRRORS};
  const TF::TermType types[] = {TF::TermType::TETRAHEDRAL,
                                TF::TermType::OCTAHEDRAL,
                                TF::TermType::DODECAHEDRAL, TF::TermType::MODE};
  const int degrees[] = {0, 0, 0, 3, 4, 5, 6, 8};
  for (size_t symmetry = 0; symmetry < mirrors.size(); ++symmetry) {
    const auto group = lenses::DodecahedralTable::make_group(mirrors[symmetry]);
    TF::PotentialTerm term;
    term.type = types[std::min<size_t>(symmetry, 3)];
    term.l = term.m = degrees[symmetry];
    term.amplitude = 1.0f;
    const auto prepared = TF::prepare_term(term);
    const auto toroidal = single_term(term);
    term.kind = TF::Kind::POLOIDAL;
    const auto poloidal = single_term(term);
    const auto image_sample = [&](const TF::PotentialField &field,
                                  math::Vector point, float strength) {
      const auto warped = exponential(point, field.velocity(point) * strength);
      return lenses::polyhedral_kaleidoscope_lens(warped, mirrors[symmetry]);
    };
    float reflected_image_difference = 0.0f;
    for (int index = 0; index < 30; ++index) {
      const auto p = direction(index, 30);
      const auto folded =
          lenses::polyhedral_kaleidoscope_lens(p, mirrors[symmetry]);
      HS_EXPECT_NEAR(prepared.sample(p).value, prepared.sample(folded).value,
                     2e-5f);
      for (size_t element = 0; element < group.count; ++element) {
        const auto &g = group.values[element];
        const auto gp = apply(g, p);
        const float determinant = math::dot(
            apply(g, math::X_AXIS),
            math::cross(apply(g, math::Y_AXIS), apply(g, math::Z_AXIS)));
        HS_EXPECT_VEC(poloidal.velocity(gp), apply(g, poloidal.velocity(p)),
                      4e-5f);
        HS_EXPECT_VEC(toroidal.velocity(gp),
                      apply(g, toroidal.velocity(p)) * determinant, 4e-5f);
        HS_EXPECT_VEC(image_sample(poloidal, gp, 0.05f),
                      image_sample(poloidal, p, 0.05f), 5e-5f);
        HS_EXPECT_VEC(
            image_sample(toroidal, gp, 0.05f),
            image_sample(toroidal, p, determinant > 0.0f ? 0.05f : -0.05f),
            5e-5f);
        if (determinant < 0.0f)
          reflected_image_difference += (image_sample(toroidal, gp, 0.05f) -
                                         image_sample(toroidal, p, 0.05f))
                                            .length();
      }
    }
    HS_EXPECT_GT(reflected_image_difference, 0.01f);
  }
  TF::PotentialTerm azimuthal{TF::TermType::MODE, TF::Kind::POLOIDAL, 6, 6,
                              1.0f};
  const auto prepared = TF::prepare_term(azimuthal);
  const auto azimuthal_poloidal = single_term(azimuthal);
  azimuthal.kind = TF::Kind::TOROIDAL;
  const auto azimuthal_toroidal = single_term(azimuthal);
  const auto image_sample = [](const TF::PotentialField &field,
                               math::Vector point, float strength) {
    const auto warped = exponential(point, field.velocity(point) * strength);
    return lenses::kaleidoscope_lens(warped);
  };
  float reflected_image_difference = 0.0f;
  for (int index = 0; index < 1000; ++index) {
    const auto p = direction(index);
    HS_EXPECT_NEAR(prepared.sample(p).value,
                   prepared.sample(lenses::kaleidoscope_lens(p)).value, 2e-5f);
    const math::Vector reflected(p.x, p.y, -p.z);
    const math::Vector rotated(0.5f * p.x + 0.8660254038f * p.z, p.y,
                               -0.8660254038f * p.x + 0.5f * p.z);
    HS_EXPECT_VEC(image_sample(azimuthal_poloidal, reflected, 0.05f),
                  image_sample(azimuthal_poloidal, p, 0.05f), 5e-5f);
    HS_EXPECT_VEC(image_sample(azimuthal_toroidal, reflected, 0.05f),
                  image_sample(azimuthal_toroidal, p, -0.05f), 5e-5f);
    HS_EXPECT_VEC(image_sample(azimuthal_toroidal, rotated, 0.05f),
                  image_sample(azimuthal_toroidal, p, 0.05f), 5e-5f);
    reflected_image_difference +=
        (image_sample(azimuthal_toroidal, reflected, 0.05f) -
         image_sample(azimuthal_toroidal, p, 0.05f))
            .length();
  }
  HS_EXPECT_GT(reflected_image_difference, 0.01f);
}

inline void test_rossby_haurwitz_matches_spherical_stream_function() {
  for (int r = 1; r <= 7; ++r) {
    double peak = 0.0;
    for (int index = 0; index <= 100000; ++index) {
      const double latitude = index * (math::PI_F / 200000.0);
      const double y = std::sin(latitude), radius = std::cos(latitude);
      const double first = std::pow(radius, r - 1) * (1 - (r + 1) * y * y);
      const double second = r * y * std::pow(radius, r - 1);
      peak = std::max(peak, std::max(std::abs(first), std::abs(second)));
    }
    const TF::PotentialTerm terms[] = {
        {TF::TermType::MODE, TF::Kind::TOROIDAL, 1, 0, -0.3f},
        {TF::TermType::MODE, TF::Kind::TOROIDAL, uint8_t(r + 1), uint8_t(r),
         (r % 2 ? -0.7f : 0.7f)}};
    const TF::PotentialField field({terms[0], terms[1], {}, {}}, 0.37f);
    for (int index = 1; index < 1000; ++index) {
      const auto p = direction(index, 1001);
      const double longitude = std::atan2(p.z, p.x), latitude = std::asin(p.y);
      const double radius = std::cos(latitude), y = std::sin(latitude);
      const double angle = r * (longitude - 0.37f);
      const double d_lat =
          -0.3 * radius / 1.01 + 0.7 / (1.01 * peak) * std::pow(radius, r - 1) *
                                     (1 - (r + 1) * y * y) * std::cos(angle);
      const double d_lon = -0.7 / (1.01 * peak) * r * y *
                           std::pow(radius, r - 1) * std::sin(angle);
      const math::Vector e_lon(-std::sin(longitude), 0, std::cos(longitude));
      const math::Vector e_lat(-y * std::cos(longitude), radius,
                               -y * std::sin(longitude));
      HS_EXPECT_VEC(field.velocity(p),
                    e_lon * float(d_lat) - e_lat * float(d_lon), 5e-5f);
    }
  }
}

inline void test_normalized_terms_have_tight_speed_bounds() {
  const auto measure = [](TF::PotentialTerm term) {
    const auto field = single_term(term);
    float peak = 0.0f;
    for (int lat = 0; lat <= 256; ++lat) {
      const double phi = -math::PI_F / 2 + lat * (math::PI_F / 256);
      for (int lon = 0; lon < 256; ++lon) {
        const double lambda = lon * (math::TWO_PI_F / 256);
        const math::Vector p(std::cos(phi) * std::cos(lambda), std::sin(phi),
                             std::cos(phi) * std::sin(lambda));
        peak = std::max(peak, field.velocity(p).length());
      }
    }
    HS_EXPECT_GE(peak, 0.985f);
    HS_EXPECT_LT(peak, 1.0f);
    HS_EXPECT_EQ(field.speed_bound(), 1.0f);
    if (term.type == TF::TermType::MODE) {
      const float correction =
          TF::MODE_CONSTANTS[term.l * (term.l + 1) / 2 + term.m].scale;
      const float unsold =
          sqrtf(term.l * (term.l + 1) * (2 * term.l + 1) / (4.0f * math::PI_F));
      HS_EXPECT_LE(peak / correction, unsold);
      double reference_peak = 0.0;
      for (int index = 0; index <= 10000; ++index) {
        const double latitude = index * (3.14159265358979323846 / 20000.0);
        const double y = std::sin(latitude), radius = std::cos(latitude);
        const auto q = legendre_reference(term.l, term.m, y);
        const double power = std::pow(radius, std::max(0, int(term.m) - 1));
        const double meridional = q.second * std::pow(radius, term.m + 1) -
                                  term.m * y * q.first * power;
        const double zonal = term.m * q.first * power;
        reference_peak = std::max(
            reference_peak, std::max(std::abs(meridional), std::abs(zonal)));
      }
      const double normalized_peak =
          reference_peak * std::abs(TF::prepare_term(term).scale);
      HS_EXPECT_NEAR(normalized_peak, 1.0 / 1.01, 1e-4);
    }
  };
  for (int l = 1; l <= 8; ++l)
    for (int m = 0; m <= l; ++m)
      measure({TF::TermType::MODE, TF::Kind::TOROIDAL, uint8_t(l), uint8_t(m),
               1.0f});
  for (auto type : {TF::TermType::TETRAHEDRAL, TF::TermType::OCTAHEDRAL,
                    TF::TermType::DODECAHEDRAL})
    measure({type, TF::Kind::POLOIDAL, 0, 0, 1.0f});
}

/** @brief Zonal harmonic peak speed is independent of authored phase and drift. */
inline void test_zonal_potential_phase_preserves_peak_speed() {
  constexpr float COORDINATE = 0.7071067811865475f;
  const math::Vector point(COORDINATE, COORDINATE, 0.0f);
  for (const auto kind : {TF::Kind::TOROIDAL, TF::Kind::POLOIDAL}) {
    for (const float phase : {0.0f, math::PI_F / 2, math::PI_F, -1.3f}) {
      std::array<TF::PotentialTerm, 4> terms{};
      terms[0] = {TF::TermType::MODE, kind, 2, 0, 1.0f, phase};
      for (const float drift : {-2.0f, 0.0f, 0.7f}) {
        const TF::PotentialField field(terms, drift);
        HS_EXPECT_NEAR(field.velocity(point).length(), 1.0f / 1.01f, 1e-5f);
        const math::Vector expected =
            kind == TF::Kind::TOROIDAL
                ? math::Vector(0.0f, 0.0f, 1.0f / 1.01f)
                : math::Vector(-COORDINATE, COORDINATE, 0.0f) / 1.01f;
        HS_EXPECT_VEC(field.velocity(point), expected, 1e-5f);
      }
    }
  }
}

inline void test_fold_advisory_precedes_sampled_jacobian_onset() {
  std::vector<std::array<TF::PotentialTerm, 4>> configurations;
  for (int l = 1; l <= 8; ++l)
    for (int m = 0; m <= l; ++m)
      for (auto kind : {TF::Kind::POLOIDAL, TF::Kind::TOROIDAL})
        configurations.push_back(
            {TF::PotentialTerm{TF::TermType::MODE, kind, uint8_t(l), uint8_t(m),
                               1.0f},
             {},
             {},
             {}});
  for (auto type : {TF::TermType::TETRAHEDRAL, TF::TermType::OCTAHEDRAL,
                    TF::TermType::DODECAHEDRAL})
    for (auto kind : {TF::Kind::POLOIDAL, TF::Kind::TOROIDAL})
      configurations.push_back(
          {TF::PotentialTerm{type, kind, 0, 0, 1.0f}, {}, {}, {}});
  for (int r = 1; r <= 7; ++r)
    configurations.push_back(
        {TF::PotentialTerm{TF::TermType::MODE, TF::Kind::TOROIDAL, 1, 0, -0.3f},
         TF::PotentialTerm{TF::TermType::MODE, TF::Kind::TOROIDAL,
                           uint8_t(r + 1), uint8_t(r), 0.7f},
         {},
         {}});
  for (int index = 0; index < 12; ++index) {
    std::array<TF::PotentialTerm, 4> terms;
    for (int i = 0; i < 4; ++i) {
      const int l = 1 + (index * 3 + i * 7) % 8;
      terms[i] = {TF::TermType::MODE,
                  i % 2 ? TF::Kind::POLOIDAL : TF::Kind::TOROIDAL,
                  uint8_t(l),
                  uint8_t((index + i) % (l + 1)),
                  sinf(index * 1.17f + i * 0.73f),
                  0.31f * i};
    }
    configurations.push_back(terms);
  }
  for (const auto &terms : configurations) {
    const TF::PotentialField field(terms);
    const float advisory =
        TF::fold_strength(field.speed_bound(), field.lipschitz_bound());
    float minimum_at_advisory = 1.0f;
    float onset = std::numeric_limits<float>::infinity();
    for (int index = 0; index < 500; ++index) {
      const auto p = direction(index, 500);
      const auto u = field.velocity(p);
      const float speed = u.length();
      const auto seed = fabsf(p.y) < 0.9f ? math::Y_AXIS : math::X_AXIS;
      const auto e =
          speed > 1e-5f ? u / speed : math::cross(p, seed).normalized();
      const auto f = math::cross(p, e);
      const auto derivative = [&](math::Vector axis) {
        return (field.velocity(exponential(p, axis * 0.0005f)) -
                field.velocity(exponential(p, axis * -0.0005f))) /
               0.001f;
      };
      const auto de = derivative(e), df = derivative(f);
      const float a = math::dot(e, de), b = math::dot(e, df);
      const float c = math::dot(f, de), d = math::dot(f, df);
      const auto jacobian = [&](float strength) {
        const float arc = strength * speed;
        const float sinc = fabsf(arc) > 1e-5f ? sinf(arc) / arc : 1.0f;
        return (1.0f + strength * a) * (cosf(arc) + sinc * strength * d) -
               strength * strength * sinc * b * c;
      };
      minimum_at_advisory = std::min(minimum_at_advisory, jacobian(advisory));
      const float cap = math::PI_F / field.speed_bound();
      if (speed <= 1e-5f) {
        const float discriminant = 0.25f * (a - d) * (a - d) + b * c;
        if (discriminant >= 0.0f) {
          const float root = sqrtf(discriminant);
          for (const float eigenvalue :
               {0.5f * (a + d) - root, 0.5f * (a + d) + root})
            if (eigenvalue < 0.0f && -1.0f / eigenvalue <= cap)
              onset = std::min(onset, -1.0f / eigenvalue);
        }
      }
      float previous = 0.0f;
      for (int scan = 1; scan <= 256; ++scan) {
        const float strength = cap * scan / 256.0f;
        if (jacobian(strength) <= 0.0f) {
          float low = previous, high = strength;
          for (int iteration = 0; iteration < 18; ++iteration) {
            const float middle = 0.5f * (low + high);
            if (jacobian(middle) > 0.0f)
              low = middle;
            else
              high = middle;
          }
          onset = std::min(onset, high);
          break;
        }
        previous = strength;
      }

      if (index % 50 == 0) {
        const float strength = advisory * 1.3f;
        const auto transform = [&](math::Vector point) {
          return exponential(point, field.velocity(point) * strength);
        };
        const auto mapped_de = (transform(exponential(p, e * 0.0005f)) -
                                transform(exponential(p, e * -0.0005f))) /
                               0.001f;
        const auto mapped_df = (transform(exponential(p, f * 0.0005f)) -
                                transform(exponential(p, f * -0.0005f))) /
                               0.001f;
        const float numeric =
            math::dot(transform(p), math::cross(mapped_de, mapped_df));
        HS_EXPECT_NEAR(jacobian(strength), numeric, 0.002f);
      }
    }
    HS_EXPECT_GT(minimum_at_advisory, 0.0f);
    HS_EXPECT_TRUE(std::isfinite(onset));
    HS_EXPECT_LE(advisory, onset);
    if (terms[0].type == TF::TermType::TETRAHEDRAL &&
        terms[0].kind == TF::Kind::POLOIDAL) {
      HS_EXPECT_NEAR(onset, 0.40185f, 0.004f);
      HS_EXPECT_GT((onset - advisory) / onset, 0.05f);
      HS_EXPECT_LT((onset - advisory) / onset, 0.08f);
      std::printf("tetrahedral onset %.6f advisory %.6f\n", onset, advisory);
    }
  }
  const auto tetrahedral =
      single_term({TF::TermType::TETRAHEDRAL, TF::Kind::POLOIDAL, 0, 0, 1.0f});
  const float advisory = TF::fold_strength(tetrahedral.speed_bound(),
                                           tetrahedral.lipschitz_bound());
  HS_EXPECT_NEAR(advisory, 0.3772f, 0.001f);
  HS_EXPECT_LT(advisory, 0.39787f * 1.01f);
}

/** @brief Drift rotates every scalar family consistently in its local frame. */
inline void test_potential_drift_rotates_all_terms() {
  const auto rotate_y = [](const math::Vector &p, double angle) {
    return math::Vector(cos(angle) * p.x + sin(angle) * p.z, p.y,
                        -sin(angle) * p.x + cos(angle) * p.z);
  };
  constexpr TF::PotentialTerm CASES[] = {
      {TF::TermType::MODE, TF::Kind::TOROIDAL, 4, 0, 1.0f, 0.31f},
      {TF::TermType::MODE, TF::Kind::TOROIDAL, 4, 1, 1.0f, 0.31f},
      {TF::TermType::MODE, TF::Kind::TOROIDAL, 4, 2, 1.0f, 0.31f},
      {TF::TermType::TETRAHEDRAL, TF::Kind::TOROIDAL, 0, 0, 1.0f},
      {TF::TermType::OCTAHEDRAL, TF::Kind::TOROIDAL, 0, 0, 1.0f},
      {TF::TermType::DODECAHEDRAL, TF::Kind::TOROIDAL, 0, 0, 1.0f},
  };
  for (auto term : CASES)
    for (TF::Kind kind : {TF::Kind::TOROIDAL, TF::Kind::POLOIDAL}) {
      term.kind = kind;
      std::array<TF::PotentialTerm, 4> terms{};
      terms[0] = term;
      const TF::PotentialField base(terms);
      const auto scalar = TF::prepare_term(term);
      for (const auto orientation :
           {math::Quaternion(1, 0, 0, 0),
            math::make_rotation(math::Y_AXIS,
                                math::Vector::from_spherical(0.7f, 0.9f))})
        for (const float phase : {-0.47f, 1.1f}) {
          const TF::PotentialField field(terms, phase, orientation);
          const TF::PotentialField stationary(terms, 0, orientation);
          const auto drifting_scalar = TF::prepare_term(term, phase);
          float movement = 0;
          for (int index = 0; index < 64; ++index) {
            const auto p = direction(index, 64);
            const auto local = math::rotate(p, orientation.conjugate());
            const auto lookup = rotate_y(local, phase);
            const auto expected = math::rotate(
                rotate_y(base.velocity(lookup), -phase), orientation);
            HS_EXPECT_VEC(field.velocity(p), expected, 4e-5f);
            movement += (field.velocity(p) - stationary.velocity(p)).length();
            const auto sample = drifting_scalar.sample(p);
            const auto original = scalar.sample(rotate_y(p, phase));
            HS_EXPECT_NEAR(sample.value, original.value, 4e-5f);
            HS_EXPECT_VEC(sample.gradient, rotate_y(original.gradient, -phase),
                          4e-5f);
          }
          HS_EXPECT_EQ(field.speed_bound(), base.speed_bound());
          HS_EXPECT_EQ(field.lipschitz_bound(), base.lipschitz_bound());
          if (term.type == TF::TermType::MODE && term.m == 0)
            HS_EXPECT_NEAR(movement, 0.0f, 1e-6f);
          else
            HS_EXPECT_GT(movement, 0.01f);
        }
      terms[0].amplitude = 0;
      HS_EXPECT_VEC(TF::PotentialField(terms, 0.7f).velocity(direction(7)),
                    math::Vector(0, 0, 0), 0.0f);
    }
}

/** @brief Poloidal fields stay tangent to approximate lens directions. */
inline void test_potential_near_unit_tangency() {
  for (const auto type :
       {TF::TermType::MODE, TF::TermType::TETRAHEDRAL, TF::TermType::OCTAHEDRAL,
        TF::TermType::DODECAHEDRAL}) {
    std::array<TF::PotentialTerm, 4> terms{};
    terms[0] = {type, TF::Kind::POLOIDAL, 6, 3, 1.0f, 0.31f};
    const TF::PotentialField field(terms, 0.73f,
                                   math::make_rotation(math::X_AXIS, 0.47f));
    for (float scale : {0.9985f, 1.0f, 1.0015f})
      for (int index = 0; index < 128; ++index) {
        const auto p = direction(index, 128) * scale;
        HS_EXPECT_NEAR(math::dot(p, field.velocity(p)), 0.0f, 3e-6f);
      }
  }
}

/** @brief Velocity reads live prepared terms after public term edits. */
inline void test_potential_velocity_reads_prepared_term_edits() {
  std::array<TF::PotentialTerm, 4> terms{};
  terms[0] = {TF::TermType::MODE, TF::Kind::POLOIDAL, 3, 1, 0.8f, 0.4f};
  TF::PotentialField edited;
  edited.terms[0] = TF::prepare_term(terms[0]);
  const TF::PotentialField expected_poloidal(terms);
  for (int index = 0; index < 64; ++index) {
    const auto point = direction(index, 64);
    HS_EXPECT_TRUE(vec_identical(edited.velocity(point),
                                 expected_poloidal.velocity(point)));
  }
  edited.terms[0].term.kind = TF::Kind::TOROIDAL;
  terms[0].kind = TF::Kind::TOROIDAL;
  const TF::PotentialField expected_toroidal(terms);
  for (int index = 0; index < 64; ++index) {
    const auto point = direction(index, 64);
    HS_EXPECT_TRUE(vec_identical(edited.velocity(point),
                                 expected_toroidal.velocity(point)));
  }
}

inline int run_tangent_field_tests() {
  hs_test::ModuleFixture fixture("tangent_field");
  HS_RUN(test_legendre_coefficient_table);
  HS_RUN(test_horner_values_and_phased_gradients);
  HS_RUN(test_mobius_generator_velocity_and_inverse);
  HS_RUN(test_mobius_exp_log_projective_roundtrip);
  HS_RUN(test_mobius_exp_domain_corners_are_finite);
  HS_RUN(test_potential_divergence_and_curl);
  HS_RUN(test_constant_potentials_have_zero_bounds);
  HS_RUN(test_invariants_and_vector_fields_follow_fold_groups);
  HS_RUN(test_rossby_haurwitz_matches_spherical_stream_function);
  HS_RUN(test_normalized_terms_have_tight_speed_bounds);
  HS_RUN(test_zonal_potential_phase_preserves_peak_speed);
  HS_RUN(test_potential_drift_rotates_all_terms);
  HS_RUN(test_potential_near_unit_tangency);
  HS_RUN(test_potential_velocity_reads_prepared_term_edits);
  HS_RUN(test_fold_advisory_precedes_sampled_jacobian_onset);
  return fixture.result();
}
} // namespace hs_test::tangent_field_tests
