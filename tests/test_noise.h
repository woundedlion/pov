/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Golden-reference tests for the noise paths. The existing transformer tests
 * only check divergence/NaN/on-sphere invariants, so a silent change to the
 * FastNoiseLite OpenSimplex2 generator (or to noise_transform's sampling) would
 * still render deterministically and pass. These pin the actual produced values
 * of a fixed sample grid against a golden FNV-1a hash plus explicit samples, so
 * any drift in the noise output fails here.
 *
 * Coverage:
 *   - FastNoiseLite OpenSimplex2 2D/3D sample grids (default seed 1337, fixed
 *     frequency 0.125): golden FNV-1a hash + pinned individual samples.
 *   - noise_transform: pinned displaced output for fixed params and inputs.
 */
#pragma once

#include <cstdint>
#include <cstring>

#include "core/FastNoiseLite.h"
#include "core/animation.h"
#include "core/transformers.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace noise_tests {

/**
 * @brief Folds a float's 4 bytes into a running FNV-1a hash.
 * @param h Running hash, updated in place.
 * @param v The float to fold.
 * @details Hashing the raw IEEE bits makes the golden sensitive to the last ULP,
 *          so even a tiny generator drift changes the digest.
 */
inline void fnv1a_float(uint32_t &h, float v) {
  uint32_t bits;
  std::memcpy(&bits, &v, sizeof(bits));
  for (int s = 0; s < 4; ++s) {
    h ^= (bits >> (8 * s)) & 0xffu;
    h *= 16777619u;
  }
}

/**
 * @brief Golden-hashes a fixed 3D OpenSimplex2 sample grid.
 * @details Seed 1337 (FastNoiseLite default), frequency 0.125, a 5x5x5 grid over
 *          [-3, 3]. The hash is the pinned current output; two explicit samples
 *          guard against a hash collision masking a real change. GetNoise(0,0,0)
 *          is exactly 0 for OpenSimplex2 at the lattice origin.
 */
inline void test_noise3d_golden_grid() {
  FastNoiseLite n; // seed 1337
  n.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  n.SetFrequency(0.125f);

  uint32_t h = 2166136261u;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 5; ++k) {
        float x = static_cast<float>(i) * 1.5f - 3.0f;
        float y = static_cast<float>(j) * 1.5f - 3.0f;
        float z = static_cast<float>(k) * 1.5f - 3.0f;
        fnv1a_float(h, n.GetNoise(x, y, z));
      }
  HS_EXPECT_EQ(h, 0x785a2cd5u);

  HS_EXPECT_NEAR(n.GetNoise(0.0f, 0.0f, 0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(n.GetNoise(1.5f, -1.5f, 0.75f), -0.680702090f, 1e-5f);
}

/**
 * @brief Golden-hashes a fixed 2D OpenSimplex2 sample grid.
 * @details Same generator config, an 8x8 grid over [-3, 3]; backs the 2D
 *          GetNoise path that stereo_noise_warp and the stereo pattern effects
 *          sample.
 */
inline void test_noise2d_golden_grid() {
  FastNoiseLite n;
  n.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  n.SetFrequency(0.125f);

  uint32_t h = 2166136261u;
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j) {
      float x = static_cast<float>(i) * 0.75f - 3.0f;
      float y = static_cast<float>(j) * 0.75f - 3.0f;
      fnv1a_float(h, n.GetNoise(x, y));
    }
  HS_EXPECT_EQ(h, 0x81b99fd4u);

  HS_EXPECT_NEAR(n.GetNoise(0.0f, 0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(n.GetNoise(1.5f, -0.75f), 0.161937639f, 1e-5f);
}

/**
 * @brief Pins noise_transform's displaced output for fixed params and inputs.
 * @details An oracle for the whole noise_transform pipeline (3-channel sampling,
 *          tangent projection, soft-cap, renormalize), not just its on-sphere
 *          invariant: a regression that still produced a unit vector but a
 *          different one would slip past the existing test_transformers checks.
 */
inline void test_noise_transform_golden() {
  NoiseParams p;
  p.amplitude = 0.5f;
  p.scale = 4.0f;
  p.time = 1.0f;
  p.speed = 1.0f;
  p.frequency = 0.125f;
  p.sync();

  struct Case {
    Vector in;
    Vector out;
  };
  const Case cases[] = {
      {Vector(1, 0, 0), Vector(0.999982476f, 0.000889263f, -0.005846164f)},
      {Vector(0, 1, 0), Vector(0.013488254f, 0.999905467f, -0.002678307f)},
      {Vector(0, 0, 1), Vector(0.009780074f, -0.009725943f, 0.999904871f)},
  };
  for (const Case &c : cases) {
    Vector r = noise_transform(c.in, p);
    HS_EXPECT_NEAR(r.x, c.out.x, 1e-5f);
    HS_EXPECT_NEAR(r.y, c.out.y, 1e-5f);
    HS_EXPECT_NEAR(r.z, c.out.z, 1e-5f);
    HS_EXPECT_NEAR(r.length(), 1.0f, 1e-4f);
  }
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every noise-module test and reports the aggregate result.
 * @return 0 on success, non-zero on any failure.
 */
inline int run_noise_tests() {
  hs_test::ModuleFixture fixture("noise");

  test_noise3d_golden_grid();
  test_noise2d_golden_grid();
  test_noise_transform_golden();

  return fixture.result();
}

} // namespace noise_tests
} // namespace hs_test
