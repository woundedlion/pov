/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Noise field interpolation, derivatives, basis and seed contracts.
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

/** @brief Pins noise field key identity. */
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

/** @brief Pins noise field periodic coordinates. */
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
 * @brief Verifies both time loops trace a NOISE_LOOP_RADIUS circle: the sphere
 *        loop in the lattice xy plane, the projected loop with x == y.
 */
inline void test_noise_field_loop_offset_geometry() {
  const math::Vector sphere0 = math::noise_sphere_loop_offset(0.0f);
  const math::Vector sphere_quarter = math::noise_sphere_loop_offset(0.25f);
  HS_EXPECT_NEAR((sphere_quarter - sphere0).length(),
                 math::NOISE_LOOP_RADIUS * 1.41421356f, 1e-3f);
  const math::Vector projected0 = math::noise_projected_loop_offset(0.0f);
  const math::Vector projected_quarter =
      math::noise_projected_loop_offset(0.25f);
  HS_EXPECT_NEAR((projected_quarter - projected0).length(),
                 math::NOISE_LOOP_RADIUS * 1.41421356f, 1e-3f);
  for (float phase : {0.0f, 0.125f, 0.25f, 0.4f, 0.9375f, -0.3f}) {
    const math::Vector sphere = math::noise_sphere_loop_offset(phase);
    HS_EXPECT_NEAR(sphere.length(), math::NOISE_LOOP_RADIUS, 1e-4f);
    HS_EXPECT_EQ(sphere.z, 0.0f);
    const math::Vector projected = math::noise_projected_loop_offset(phase);
    HS_EXPECT_NEAR(projected.length(), math::NOISE_LOOP_RADIUS, 1e-4f);
    HS_EXPECT_EQ(projected.x, projected.y);
  }
}

/**
 * @brief Verifies the hoisted-offset overloads reproduce the phase-taking ones
 *        bit for bit.
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

/** @brief Pins noise field octave formulas. */
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

/** @brief Pins noise field ridged channel pairs. */
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

/** @brief Pins noise field direct tangent. */
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
    const math::Vector quarter_turns = math::sample_direct_tangent(
        noise, math::NoiseBasis::SIMPLEX, q, v, 0.25f);
    HS_EXPECT_NEAR(quarter_turns.x, quarter_turn.x, 1e-5f);
    HS_EXPECT_NEAR(quarter_turns.y, quarter_turn.y, 1e-5f);
    HS_EXPECT_NEAR(quarter_turns.z, quarter_turn.z, 1e-5f);
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
      const math::Vector u_quarter =
          math::sample_direct_tangent(noise, basis, q, v, 0.25f);
      const math::Vector u0_turned = math::cross(v, u0);
      HS_EXPECT_NEAR(u_quarter.x, u0_turned.x, 1e-5f);
      HS_EXPECT_NEAR(u_quarter.y, u0_turned.y, 1e-5f);
      HS_EXPECT_NEAR(u_quarter.z, u0_turned.z, 1e-5f);
    }
  }
}

/** @brief Pins noise field tetrahedral gradient. */
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

/** @brief Pins noise field analytic gradient. */
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

/** @brief Pins vector noise rotation setter order. */
inline void test_vector_noise_rotation_setter_order() {
  struct Case {
    FastNoiseLite::RotationType3D rotation;
    math::Vector simplex;
    math::Vector untransformed;
  };
  const Case CASES[] = {
      {FastNoiseLite::RotationType3D_None,
       {1.30739737f, -2.88196921f, 0.453307718f},
       {1.20877635f, -2.75325942f, 0.496854395f}},
      {FastNoiseLite::RotationType3D_ImproveXYPlanes,
       {1.25860775f, -2.76022482f, 0.495675534f},
       {1.25860775f, -2.76022482f, 0.495675534f}},
      {FastNoiseLite::RotationType3D_ImproveXZPlanes,
       {1.48228431f, -2.84245944f, 0.170210212f},
       {1.48228431f, -2.84245944f, 0.170210212f}},
  };
  for (const Case &test : CASES) {
    HS_CONTEXT("rotation", static_cast<int>(test.rotation));
    // FASTNOISELITE_ONLY_OPENSIMPLEX2 leaves warp type selecting the 3D transform.
    for (const auto warp : {FastNoiseLite::DomainWarpType_OpenSimplex2,
                            FastNoiseLite::DomainWarpType_BasicGrid}) {
      HS_CONTEXT("warp", static_cast<int>(warp));
      FastNoiseLite first = make_noise(31), second = make_noise(31);
      first.SetDomainWarpType(warp);
      first.SetRotationType3D(test.rotation);
      second.SetRotationType3D(test.rotation);
      second.SetDomainWarpType(warp);
      math::Vector a(1.25f, -2.75f, 0.5f), b = a;
      first.GetVectorNoiseSingle(a.x, a.y, a.z);
      second.GetVectorNoiseSingle(b.x, b.y, b.z);
      HS_EXPECT_EQ(a.x, b.x);
      HS_EXPECT_EQ(a.y, b.y);
      HS_EXPECT_EQ(a.z, b.z);
      const math::Vector &EXPECTED =
          warp == FastNoiseLite::DomainWarpType_OpenSimplex2
              ? test.simplex
              : test.untransformed;
      HS_EXPECT_NEAR(a.x, EXPECTED.x, 2e-6f);
      HS_EXPECT_NEAR(a.y, EXPECTED.y, 2e-6f);
      HS_EXPECT_NEAR(a.z, EXPECTED.z, 2e-6f);
    }
  }
}

/** @brief Pins noise field simplex curl approximation. */
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
      max_error = hs_test::fold_worst(max_error, error);
      total_error += error;
      ++samples;
    }
  }
  std::printf("  simplex curl error: max=%.9g mean=%.9g samples=%d\n",
              max_error, total_error / samples, samples);
  HS_EXPECT_LT(max_error, 0.18f);
  HS_EXPECT_LT(total_error / samples, 0.03f);
}

/** @brief Pins noise field curl tangent. */
inline void test_noise_field_curl_tangent() {
  const FastNoiseLite noise = make_noise(7127);
  for (math::NoiseBasis basis :
       {math::NoiseBasis::FBM3, math::NoiseBasis::RIDGED3}) {
    HS_CONTEXT("curl basis", static_cast<int>(basis));
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
    const float MEAN_MAGNITUDE = magnitude_sum / (17.0f * 24.0f);
    const float EXPECTED_MEAN =
        basis == math::NoiseBasis::FBM3 ? 0.954f : 0.987f;
    HS_EXPECT_NEAR(MEAN_MAGNITUDE, EXPECTED_MEAN, 0.005f);
  }
}

/** @brief Pins sphere exp map and transport. */
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

/** @brief Pins half radian exp map approximation. */
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
        max_error =
            hs_test::fold_worst(max_error, fabsf(exact.x - approximate.x));
        max_error =
            hs_test::fold_worst(max_error, fabsf(exact.y - approximate.y));
        max_error =
            hs_test::fold_worst(max_error, fabsf(exact.z - approximate.z));
      }
    }
  }
  HS_EXPECT_LT(max_error, 2e-7f);
}

inline int run_noise_field_tests() {
  hs_test::ModuleFixture fixture("noise_field");
  test_noise_field_key_identity();
  test_noise_field_periodic_coordinates();
  test_noise_field_loop_offset_geometry();
  test_noise_field_hoisted_loop_offsets();
  test_noise_field_octave_formulas();
  test_noise_field_ridged_channel_pairs();
  test_noise_field_direct_tangent();
  test_noise_field_tetrahedral_gradient();
  test_noise_field_analytic_gradient();
  test_vector_noise_rotation_setter_order();
  test_noise_field_simplex_curl_approximation();
  test_noise_field_curl_tangent();
  test_sphere_exp_map_and_transport();
  test_half_radian_exp_map_approximation();
  return fixture.result();
}

} // namespace noise_field_tests
} // namespace hs_test
