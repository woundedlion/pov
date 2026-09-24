/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstring>

#include "core/math/noise_field.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace noise_field_tests {

inline FastNoiseLite make_noise(int32_t seed) {
  FastNoiseLite noise(seed);
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetFrequency(1.0f);
  return noise;
}

inline void test_noise_field_key_identity() {
  math::NoiseFieldSpec a{math::NoiseDomain::SPHERE_3D,
                         math::NoiseBasis::FBM3,
                         41,
                         2.0f,
                         0.01f,
                         0.25f,
                         math::NoiseChannelLayout::DIRECT_V1};
  math::NoiseFieldSpec b = a;
  b.scale = 7.0f;
  b.rate = -0.01f;
  b.phase = 0.75f;
  HS_EXPECT_TRUE(math::noise_field_key(a) == math::noise_field_key(b));
  b.domain = math::NoiseDomain::PROJECTED_2D;
  HS_EXPECT_FALSE(math::noise_field_key(a) == math::noise_field_key(b));
  b = a;
  b.basis = math::NoiseBasis::RIDGED3;
  HS_EXPECT_FALSE(math::noise_field_key(a) == math::noise_field_key(b));
  b = a;
  b.seed = a.seed + 1;
  HS_EXPECT_FALSE(math::noise_field_key(a) == math::noise_field_key(b));
  b = a;
  b.channel_layout = math::NoiseChannelLayout::CURL_V1;
  HS_EXPECT_FALSE(math::noise_field_key(a) == math::noise_field_key(b));
  b = a;
  b.octave_layout = a.octave_layout + 1;
  HS_EXPECT_FALSE(math::noise_field_key(a) == math::noise_field_key(b));
  b = a;
  b.loop_layout = a.loop_layout + 1;
  HS_EXPECT_FALSE(math::noise_field_key(a) == math::noise_field_key(b));
  b = a;
  b.stencil_layout = a.stencil_layout + 1;
  HS_EXPECT_FALSE(math::noise_field_key(a) == math::noise_field_key(b));
}

inline void test_noise_field_periodic_coordinates() {
  const math::Vector v = math::Vector(0.25f, -0.5f, 0.8291562f).normalized();
  const math::Vector sphere0 = math::noise_sphere_coordinate(v, 3.0f, 0.0f);
  const math::Vector sphere1 = math::noise_sphere_coordinate(v, 3.0f, 1.0f);
  HS_EXPECT_EQ(std::memcmp(&sphere0, &sphere1, sizeof(math::Vector)), 0);
  const math::Complex p(1.25f, -0.75f);
  const math::Vector projected0 =
      math::noise_projected_coordinate(p, 0.5f, 0.0f);
  const math::Vector projected1 =
      math::noise_projected_coordinate(p, 0.5f, 1.0f);
  HS_EXPECT_EQ(std::memcmp(&projected0, &projected1, sizeof(math::Vector)), 0);
}

/**
 * @brief Verifies the hoisted-offset overloads reproduce the phase-taking ones
 *        bit for bit, so lifting the loop point out of a per-pixel walk is not
 *        a numeric change.
 */
inline void test_noise_field_hoisted_loop_offsets() {
  const math::Vector v = math::Vector(-0.6f, 0.3f, 0.7416198f).normalized();
  const math::Complex p(-2.5f, 0.875f);
  for (float phase : {0.0f, 0.125f, 0.4f, 0.9375f, 1.75f, -0.3f}) {
    const math::Vector sphere = math::noise_sphere_coordinate(v, 3.0f, phase);
    const math::Vector sphere_hoisted = math::noise_sphere_coordinate(
        v, 3.0f, math::noise_sphere_loop_offset(phase));
    HS_EXPECT_EQ(std::memcmp(&sphere, &sphere_hoisted, sizeof(math::Vector)),
                 0);
    const math::Vector projected =
        math::noise_projected_coordinate(p, 0.5f, phase);
    const math::Vector projected_hoisted = math::noise_projected_coordinate(
        p, 0.5f, math::noise_projected_loop_offset(phase));
    HS_EXPECT_EQ(
        std::memcmp(&projected, &projected_hoisted, sizeof(math::Vector)), 0);
  }
}

inline void test_noise_field_octave_formulas() {
  const FastNoiseLite noise = make_noise(-317);
  constexpr std::array<math::Vector, 4> POINTS = {
      math::Vector(0.0f, 0.0f, 0.0f), math::Vector(3.25f, -7.5f, 11.0f),
      math::Vector(-31.75f, 0.125f, 2.5f), math::Vector(64.0f, 32.0f, -16.0f)};
  constexpr std::array<std::array<float, 4>, 3> EXPECTED = {{
      {0.0f, 0.111684620f, 0.477153748f, -0.110211231f},
      {0.0f, -0.00533447927f, 0.00685898308f, -0.173915580f},
      {1.0f, 0.734051824f, -0.0769191384f, 0.208394766f},
  }};
  size_t basis_index = 0;
  for (math::NoiseBasis basis :
       {math::NoiseBasis::SIMPLEX, math::NoiseBasis::FBM3,
        math::NoiseBasis::RIDGED3}) {
    for (size_t point = 0; point < POINTS.size(); ++point)
      HS_EXPECT_NEAR(math::sample_noise_octaves(noise, basis, POINTS[point]),
                     EXPECTED[basis_index][point], 2e-6f);
    ++basis_index;
  }
}

inline void test_noise_field_ridged_channel_pairs() {
  const FastNoiseLite noise = make_noise(991);
  const math::Vector q(-3.0f, 8.5f, 29.0f);
  for (size_t channel = 0; channel < 3; ++channel) {
    const float c =
        math::sample_noise_octaves(noise, math::NoiseBasis::RIDGED3,
                                   q + math::NOISE_CHANNEL_OFFSETS[channel]);
    const float d =
        math::sample_noise_octaves(noise, math::NoiseBasis::RIDGED3,
                                   q + math::NOISE_RIDGED_OFFSETS[channel]);
    HS_EXPECT_NEAR(math::sample_noise_vector_channel(
                       noise, math::NoiseBasis::RIDGED3, q, channel),
                   0.5f * (c - d), 1e-7f);
  }
}

inline void test_noise_field_direct_tangent() {
  const FastNoiseLite noise = make_noise(1337);
  constexpr std::array<math::Vector, 6> DIRECTIONS = {
      math::Vector(1.0f, 0.0f, 0.0f), math::Vector(-1.0f, 0.0f, 0.0f),
      math::Vector(0.0f, 1.0f, 0.0f), math::Vector(0.0f, -1.0f, 0.0f),
      math::Vector(0.0f, 0.0f, 1.0f), math::Vector(0.0f, 0.0f, -1.0f)};
  for (const math::Vector &v : DIRECTIONS) {
    const math::Vector q = math::noise_sphere_coordinate(v, 1.0f, 0.375f);
    const math::Vector simplex = math::sample_direct_tangent(
        noise, math::NoiseBasis::SIMPLEX, q, v, 1.0f, 0.0f);
    const math::Vector specialized =
        math::sample_direct_simplex_tangent(noise, q, v);
    const math::Vector quarter_turn = math::sample_direct_tangent(
        noise, math::NoiseBasis::SIMPLEX, q, v, 0.0f, 1.0f);
    const math::Vector expected_turn = math::cross(v, simplex);
    HS_EXPECT_NEAR(quarter_turn.x, expected_turn.x, 1e-6f);
    HS_EXPECT_NEAR(quarter_turn.y, expected_turn.y, 1e-6f);
    HS_EXPECT_NEAR(quarter_turn.z, expected_turn.z, 1e-6f);
    HS_EXPECT_EQ(std::memcmp(&simplex, &specialized, sizeof(math::Vector)), 0);
    for (math::NoiseBasis basis :
         {math::NoiseBasis::FBM3, math::NoiseBasis::RIDGED3}) {
      const math::Vector u0 =
          math::sample_direct_tangent(noise, basis, q, v, 0.0f);
      const math::Vector u1 =
          math::sample_direct_tangent(noise, basis, q, v, 1.0f);
      HS_EXPECT_NEAR(math::dot(u0, v), 0.0f, 1e-6f);
      HS_EXPECT_LE(u0.length(), 1.000001f);
      HS_EXPECT_GT(u0.length(), 0.001f);
      HS_EXPECT_NEAR(u0.x, u1.x, 1e-6f);
      HS_EXPECT_NEAR(u0.y, u1.y, 1e-6f);
      HS_EXPECT_NEAR(u0.z, u1.z, 1e-6f);
    }
  }
}

inline void test_noise_field_tetrahedral_gradient() {
  auto linear = [](const math::Vector &p) {
    return 0.25f * p.x - 0.5f * p.y + p.z;
  };
  constexpr std::array<math::Vector, 3> POINTS = {
      math::Vector(0.0f, 0.0f, 0.0f), math::Vector(1.0f, -2.0f, 3.0f),
      math::Vector(-4.0f, 0.5f, -1.25f)};
  for (const math::Vector &q : POINTS) {
    const math::Vector expected(0.25f, -0.5f, 1.0f);
    const math::Vector actual = math::tetrahedral_gradient(q, linear);
    HS_EXPECT_NEAR(actual.x, expected.x, 2e-4f);
    HS_EXPECT_NEAR(actual.y, expected.y, 2e-4f);
    HS_EXPECT_NEAR(actual.z, expected.z, 2e-4f);
  }
}

inline void test_noise_field_analytic_gradient() {
  constexpr float STEP = 1.0f / 2048.0f;
  constexpr std::array<math::Vector, 5> POINTS = {
      math::Vector(0.0f, 0.0f, 0.0f), math::Vector(0.25f, -0.75f, 1.5f),
      math::Vector(-3.125f, 8.0f, 0.0625f), math::Vector(31.0f, -17.0f, 9.0f),
      math::Vector(-0.499f, 0.501f, -0.001f)};
  for (FastNoiseLite::RotationType3D rotation :
       {FastNoiseLite::RotationType3D_None,
        FastNoiseLite::RotationType3D_ImproveXYPlanes,
        FastNoiseLite::RotationType3D_ImproveXZPlanes}) {
    FastNoiseLite noise = make_noise(-991);
    noise.SetFrequency(0.73f);
    noise.SetRotationType3D(rotation);
    for (const math::Vector &point : POINTS) {
      math::Vector gradient;
      noise.GetNoiseGradientSingle(point.x, point.y, point.z, gradient.x,
                                   gradient.y, gradient.z);
      const auto derivative = [&](const math::Vector &axis) {
        const math::Vector low = point - STEP * axis;
        const math::Vector high = point + STEP * axis;
        return (noise.GetNoiseSingle(high.x, high.y, high.z) -
                noise.GetNoiseSingle(low.x, low.y, low.z)) /
               (2.0f * STEP);
      };
      HS_EXPECT_NEAR(gradient.x, derivative(math::X_AXIS), 8e-3f);
      HS_EXPECT_NEAR(gradient.y, derivative(math::Y_AXIS), 8e-3f);
      HS_EXPECT_NEAR(gradient.z, derivative(math::Z_AXIS), 8e-3f);
    }
  }
}

inline void test_noise_field_simplex_curl_approximation() {
  const FastNoiseLite noise = make_noise(7127);
  float max_error = 0.0f;
  float total_error = 0.0f;
  int samples = 0;
  for (int latitude = -8; latitude <= 8; ++latitude) {
    const float y = latitude / 8.0f;
    const float radius = sqrtf(1.0f - y * y);
    for (int longitude = 0; longitude < 48; ++longitude) {
      const float angle = math::TWO_PI_F * longitude / 48.0f;
      const math::Vector v(radius * cosf(angle), y, radius * sinf(angle));
      const math::Vector q = math::noise_sphere_coordinate(v, 2.0f, 0.125f);
      const math::Vector analytic =
          math::sample_simplex_curl_tangent(noise, q, v);
      const math::Vector reference_gradient =
          math::tetrahedral_gradient(q, [&](const math::Vector &point) {
            return math::sample_noise_octaves(noise, math::NoiseBasis::SIMPLEX,
                                              point);
          });
      math::Vector reference = math::cross(v, reference_gradient);
      const float length = reference.length();
      if (length > 1.0f)
        reference /= length;
      const float error = (analytic - reference).length();
      max_error = std::max(max_error, error);
      total_error += error;
      ++samples;
    }
  }
  std::printf("  simplex curl error: max=%.9g mean=%.9g samples=%d\n",
              max_error, total_error / samples, samples);
  // 816 samples: max 0.176433, mean 0.0279079; margins 2% and 7%.
  HS_EXPECT_LT(max_error, 0.18f);
  HS_EXPECT_LT(total_error / samples, 0.03f);
}

inline void test_noise_field_curl_tangent() {
  const FastNoiseLite noise = make_noise(7127);
  for (math::NoiseBasis basis :
       {math::NoiseBasis::FBM3, math::NoiseBasis::RIDGED3}) {
    float magnitude_sum = 0.0f;
    for (int latitude = -8; latitude <= 8; ++latitude) {
      const float y = latitude / 8.0f;
      const float radius = sqrtf(1.0f - y * y);
      for (int longitude = 0; longitude < 24; ++longitude) {
        const float angle = math::TWO_PI_F * longitude / 24.0f;
        const math::Vector v(radius * cosf(angle), y, radius * sinf(angle));
        const math::Vector q = math::noise_sphere_coordinate(v, 2.0f, 0.125f);
        const math::Vector u = math::sample_curl_tangent(noise, basis, q, v);
        magnitude_sum += u.length();
        HS_EXPECT_TRUE(std::isfinite(u.x) && std::isfinite(u.y) &&
                       std::isfinite(u.z));
        HS_EXPECT_NEAR(math::dot(u, v), 0.0f, 2e-5f);
        HS_EXPECT_LE(u.length(), 1.00001f);
      }
    }
    HS_EXPECT_GT(magnitude_sum, 1.0f);
  }
}

inline void test_sphere_exp_map_and_transport() {
  const math::Vector v(0.0f, 1.0f, 0.0f);
  HS_EXPECT_EQ(math::sphere_exp_map(v, math::Vector()), v);
  const math::Vector tangent(0.25f, 0.0f, -0.1f);
  const math::Vector moved = math::sphere_exp_map(v, tangent);
  HS_EXPECT_NEAR(moved.length(), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::fast_acos(math::dot(v, moved)), tangent.length(),
                 5.1e-5f);
  const math::Vector transported = math::parallel_transport(v, moved, tangent);
  HS_EXPECT_NEAR(math::dot(transported, moved), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(transported.length(), tangent.length(), 1e-6f);
}

inline void test_half_radian_exp_map_approximation() {
  float max_error = 0.0f;
  for (int latitude_step = -16; latitude_step <= 16; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * math::PI_F / 16.0f);
    const float radius = cosf(latitude);
    for (int longitude_step = 0; longitude_step < 64; ++longitude_step) {
      const float longitude = longitude_step * (math::TWO_PI_F / 64.0f);
      const math::Vector v(radius * cosf(longitude), sinf(latitude),
                           radius * sinf(longitude));
      math::Vector tangent = math::cross(v, math::Z_AXIS);
      if (math::dot(tangent, tangent) < 1e-6f)
        tangent = math::cross(v, math::X_AXIS);
      tangent = tangent.normalized();
      for (int distance_step = 0; distance_step <= 64; ++distance_step) {
        const math::Vector displacement =
            (0.5f * distance_step / 64.0f) * tangent;
        const math::Vector exact = math::sphere_exp_map(v, displacement);
        const math::Vector approximate =
            math::sphere_exp_map_half_radian(v, displacement);
        max_error = std::max(
            max_error, std::max(fabsf(exact.x - approximate.x),
                                std::max(fabsf(exact.y - approximate.y),
                                         fabsf(exact.z - approximate.z))));
      }
    }
  }
  HS_EXPECT_LT(max_error, 2e-7f);
}

inline int run_noise_field_tests() {
  hs_test::ModuleFixture fixture("noise_field");
  test_noise_field_key_identity();
  test_noise_field_periodic_coordinates();
  test_noise_field_hoisted_loop_offsets();
  test_noise_field_octave_formulas();
  test_noise_field_ridged_channel_pairs();
  test_noise_field_direct_tangent();
  test_noise_field_tetrahedral_gradient();
  test_noise_field_analytic_gradient();
  test_noise_field_simplex_curl_approximation();
  test_noise_field_curl_tangent();
  test_sphere_exp_map_and_transport();
  test_half_radian_exp_map_approximation();
  return fixture.result();
}

} // namespace noise_field_tests
} // namespace hs_test
