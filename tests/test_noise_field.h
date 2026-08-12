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

inline float reference_basis(const FastNoiseLite &noise, NoiseBasis basis,
                             const Vector &q) {
  const float n0 = noise.GetNoise(q.x, q.y, q.z);
  if (basis == NoiseBasis::SIMPLEX)
    return n0;
  const float n1 = noise.GetNoise(q.x * 2.0f, q.y * 2.0f, q.z * 2.0f);
  const float n2 = noise.GetNoise(q.x * 4.0f, q.y * 4.0f, q.z * 4.0f);
  if (basis == NoiseBasis::FBM3)
    return (4.0f * n0 + 2.0f * n1 + n2) / 7.0f;
  const float r = (4.0f * (1.0f - std::abs(n0)) + 2.0f * (1.0f - std::abs(n1)) +
                   (1.0f - std::abs(n2))) /
                  7.0f;
  return 2.0f * r - 1.0f;
}

inline void test_noise_field_key_identity() {
  NoiseFieldSpec a{
      NoiseDomain::SPHERE_3D,       NoiseBasis::FBM3, 41, 2.0f, 0.01f, 0.25f,
      NoiseChannelLayout::DIRECT_V1};
  NoiseFieldSpec b = a;
  b.scale = 7.0f;
  b.rate = -0.01f;
  b.phase = 0.75f;
  HS_EXPECT_TRUE(noise_field_key(a) == noise_field_key(b));
  b.domain = NoiseDomain::PROJECTED_2D;
  HS_EXPECT_FALSE(noise_field_key(a) == noise_field_key(b));
  b = a;
  b.channel_layout = NoiseChannelLayout::CURL_V1;
  HS_EXPECT_FALSE(noise_field_key(a) == noise_field_key(b));
}

inline void test_noise_field_periodic_coordinates() {
  const Vector v = Vector(0.25f, -0.5f, 0.8291562f).normalized();
  const Vector sphere0 = noise_sphere_coordinate(v, 3.0f, 0.0f);
  const Vector sphere1 = noise_sphere_coordinate(v, 3.0f, 1.0f);
  HS_EXPECT_EQ(std::memcmp(&sphere0, &sphere1, sizeof(Vector)), 0);
  const Complex p(1.25f, -0.75f);
  const Vector projected0 = noise_projected_coordinate(p, 0.5f, 0.0f);
  const Vector projected1 = noise_projected_coordinate(p, 0.5f, 1.0f);
  HS_EXPECT_EQ(std::memcmp(&projected0, &projected1, sizeof(Vector)), 0);
}

inline void test_noise_field_octave_formulas() {
  const FastNoiseLite noise = make_noise(-317);
  constexpr std::array<Vector, 4> POINTS = {
      Vector(0.0f, 0.0f, 0.0f), Vector(3.25f, -7.5f, 11.0f),
      Vector(-31.75f, 0.125f, 2.5f), Vector(64.0f, 32.0f, -16.0f)};
  for (NoiseBasis basis :
       {NoiseBasis::SIMPLEX, NoiseBasis::FBM3, NoiseBasis::RIDGED3})
    for (const Vector &q : POINTS)
      HS_EXPECT_NEAR(sample_noise_octaves(noise, basis, q),
                     reference_basis(noise, basis, q), 1e-7f);
}

inline void test_noise_field_ridged_channel_pairs() {
  const FastNoiseLite noise = make_noise(991);
  const Vector q(-3.0f, 8.5f, 29.0f);
  for (size_t channel = 0; channel < 3; ++channel) {
    const float c = reference_basis(noise, NoiseBasis::RIDGED3,
                                    q + NOISE_CHANNEL_OFFSETS[channel]);
    const float d = reference_basis(noise, NoiseBasis::RIDGED3,
                                    q + NOISE_RIDGED_OFFSETS[channel]);
    HS_EXPECT_NEAR(
        sample_noise_vector_channel(noise, NoiseBasis::RIDGED3, q, channel),
        0.5f * (c - d), 1e-7f);
  }
}

inline void test_noise_field_direct_tangent() {
  const FastNoiseLite noise = make_noise(1337);
  constexpr std::array<Vector, 6> DIRECTIONS = {
      Vector(1.0f, 0.0f, 0.0f), Vector(-1.0f, 0.0f, 0.0f),
      Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, -1.0f, 0.0f),
      Vector(0.0f, 0.0f, 1.0f), Vector(0.0f, 0.0f, -1.0f)};
  for (const Vector &v : DIRECTIONS) {
    const Vector q = noise_sphere_coordinate(v, 1.0f, 0.375f);
    const Vector u0 =
        sample_direct_tangent(noise, NoiseBasis::FBM3, q, v, 0.0f);
    const Vector u1 =
        sample_direct_tangent(noise, NoiseBasis::FBM3, q, v, 1.0f);
    HS_EXPECT_NEAR(dot(u0, v), 0.0f, 1e-6f);
    HS_EXPECT_LE(u0.length(), 1.000001f);
    HS_EXPECT_NEAR(u0.x, u1.x, 1e-6f);
    HS_EXPECT_NEAR(u0.y, u1.y, 1e-6f);
    HS_EXPECT_NEAR(u0.z, u1.z, 1e-6f);
  }
}

inline void test_noise_field_tetrahedral_gradient() {
  auto linear = [](const Vector &p) { return 0.25f * p.x - 0.5f * p.y + p.z; };
  constexpr std::array<Vector, 3> POINTS = {Vector(0.0f, 0.0f, 0.0f),
                                            Vector(1.0f, -2.0f, 3.0f),
                                            Vector(-4.0f, 0.5f, -1.25f)};
  for (const Vector &q : POINTS) {
    const Vector expected(0.25f, -0.5f, 1.0f);
    const Vector actual = tetrahedral_gradient(q, linear);
    HS_EXPECT_NEAR(actual.x, expected.x, 2e-4f);
    HS_EXPECT_NEAR(actual.y, expected.y, 2e-4f);
    HS_EXPECT_NEAR(actual.z, expected.z, 2e-4f);
  }
}

inline void test_noise_field_curl_tangent() {
  const FastNoiseLite noise = make_noise(7127);
  for (int latitude = -8; latitude <= 8; ++latitude) {
    const float y = latitude / 8.0f;
    const float radius = sqrtf(1.0f - y * y);
    for (int longitude = 0; longitude < 24; ++longitude) {
      const float angle = TWO_PI_F * longitude / 24.0f;
      const Vector v(radius * cosf(angle), y, radius * sinf(angle));
      const Vector q = noise_sphere_coordinate(v, 2.0f, 0.125f);
      const Vector u = sample_curl_tangent(noise, NoiseBasis::RIDGED3, q, v);
      HS_EXPECT_TRUE(std::isfinite(u.x) && std::isfinite(u.y) &&
                     std::isfinite(u.z));
      HS_EXPECT_NEAR(dot(u, v), 0.0f, 2e-5f);
      HS_EXPECT_LE(u.length(), 1.00001f);
    }
  }
}

inline void test_sphere_exp_map_and_transport() {
  const Vector v(0.0f, 1.0f, 0.0f);
  HS_EXPECT_TRUE(sphere_exp_map(v, Vector()) == v);
  const Vector tangent(0.25f, 0.0f, -0.1f);
  const Vector moved = sphere_exp_map(v, tangent);
  HS_EXPECT_NEAR(moved.length(), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(fast_acos(dot(v, moved)), tangent.length(), 2e-5f);
  const Vector transported = transport_tangent(v, moved, tangent);
  HS_EXPECT_NEAR(dot(transported, moved), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(transported.length(), tangent.length(), 1e-6f);
}

inline int run_noise_field_tests() {
  hs_test::ModuleFixture fixture("noise_field");
  test_noise_field_key_identity();
  test_noise_field_periodic_coordinates();
  test_noise_field_octave_formulas();
  test_noise_field_ridged_channel_pairs();
  test_noise_field_direct_tangent();
  test_noise_field_tetrahedral_gradient();
  test_noise_field_curl_tangent();
  test_sphere_exp_map_and_transport();
  return fixture.result();
}

} // namespace noise_field_tests
} // namespace hs_test
