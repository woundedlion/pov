/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Golden-reference tests for the noise paths. The existing transformer tests
 * only check divergence/NaN/on-sphere invariants, so a silent change to the
 * FastNoiseLite OpenSimplex2 generator (or to noise_transform's sampling) would
 * still render deterministically and pass. These pin the actual produced values
 * of a fixed sample grid against golden values, so any drift in the noise
 * output fails here with a useful sample index.
 *
 * Coverage:
 *   - FastNoiseLite OpenSimplex2 2D/3D sample grids (default seed 1337, fixed
 *     frequency 0.125): tolerance-checked sample grids.
 *   - noise_transform: pinned displaced output for fixed params and inputs.
 */
#pragma once

#include "core/vendor/FastNoiseLite.h"
#include "core/animation/animation.h"
#include "core/animation/transformer.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace noise_tests {

template <size_t N, typename SampleFn>
inline void expect_noise_grid(const float (&golden)[N], SampleFn sample,
                              const char *label) {
  for (size_t i = 0; i < N; ++i) {
    const float actual = sample(i);
    HS_CONTEXT(label, static_cast<long long>(i));
    HS_EXPECT_NEAR(actual, golden[i], 1e-5f);
  }
}

/**
 * @brief Checks a fixed 3D OpenSimplex2 sample grid.
 * @details Seed 1337 (FastNoiseLite default), frequency 0.125, a 5x5x5 grid over
 *          [-3, 3].
 */
inline void test_noise3d_golden_grid() {
  static constexpr float GOLDEN[] = {
      0.0f,
      0.262275875f,
      0.000182601056f,
      -0.118954547f,
      -0.162505522f,
      -0.268739372f,
      0.00534849893f,
      -0.165682793f,
      -0.324480593f,
      -0.300289452f,
      -0.126105696f,
      0.198956922f,
      0.0510877892f,
      -0.268885821f,
      -0.253398299f,
      -0.279718488f,
      0.226656422f,
      0.235562578f,
      -0.14045912f,
      -0.26499486f,
      -0.455629349f,
      0.0971491933f,
      0.211210907f,
      -0.0184362717f,
      -0.13788563f,
      0.00646347925f,
      4.71958316e-07f,
      -0.333672017f,
      -0.338074178f,
      -0.280233026f,
      -0.0045843292f,
      0.0f,
      -0.479047745f,
      -0.742138505f,
      -0.56183672f,
      0.300476968f,
      0.479041398f,
      0.000537356711f,
      -0.486070722f,
      -0.364384741f,
      0.270055801f,
      0.750464439f,
      0.49426046f,
      -0.00888773147f,
      -0.126702234f,
      0.200929955f,
      0.73918128f,
      0.641688645f,
      0.235991016f,
      0.00400738092f,
      0.126288295f,
      -0.185355574f,
      -0.545353353f,
      -0.492811203f,
      -0.252667874f,
      0.153195083f,
      1.26980722e-05f,
      -0.624002695f,
      -0.965121627f,
      -0.66349566f,
      0.341002256f,
      0.617554426f,
      0.0f,
      -0.623465359f,
      -0.493437111f,
      0.284797311f,
      0.965121627f,
      0.630988359f,
      0.0164001491f,
      -0.082271263f,
      0.252667904f,
      0.940847635f,
      0.802863657f,
      0.351754606f,
      0.253580898f,
      0.399021506f,
      -0.116767988f,
      -0.501583576f,
      -0.416029364f,
      -0.00254672882f,
      0.176378503f,
      0.0137204276f,
      -0.486070722f,
      -0.729165852f,
      -0.329668283f,
      0.0846422315f,
      0.470861167f,
      -0.000537356013f,
      -0.469678491f,
      -0.289182067f,
      -0.0567527749f,
      0.697224557f,
      0.495409817f,
      0.181585968f,
      0.286267698f,
      -0.219366238f,
      0.459591717f,
      0.531462669f,
      0.571771979f,
      0.80513972f,
      0.687077701f,
      0.152628958f,
      -0.252667904f,
      -0.173569068f,
      0.249845952f,
      0.359536022f,
      0.150392592f,
      -0.235340208f,
      -0.263340741f,
      0.27937001f,
      0.0429178178f,
      0.263859242f,
      -0.0502592809f,
      -0.173291832f,
      0.210480481f,
      -0.337730646f,
      0.17010358f,
      0.205373183f,
      0.275570542f,
      0.480365872f,
      -0.731270492f,
      -0.278498858f,
      0.126379594f,
      0.536400378f,
      0.765076876f,
  };
  FastNoiseLite n; // seed 1337
  n.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  n.SetFrequency(0.125f);

  expect_noise_grid(
      GOLDEN,
      [&](size_t index) {
        const int i = static_cast<int>(index / 25);
        const int j = static_cast<int>((index / 5) % 5);
        const int k = static_cast<int>(index % 5);
        const float x = static_cast<float>(i) * 1.5f - 3.0f;
        const float y = static_cast<float>(j) * 1.5f - 3.0f;
        const float z = static_cast<float>(k) * 1.5f - 3.0f;
        return n.GetNoise(x, y, z);
      },
      "OpenSimplex2 3D");
}

/**
 * @brief Checks a fixed 2D OpenSimplex2 sample grid.
 * @details Same generator config, an 8x8 grid over [-3, 3]; backs the 2D
 *          GetNoise path that projected noise consumers sample.
 */
inline void test_noise2d_golden_grid() {
  static constexpr float GOLDEN[] = {
      -0.812879086f, -0.728954434f,  -0.592017591f,  -0.371756405f,
      -0.113768227f, 0.107305534f,   0.220385596f,   0.147303268f,
      -0.893296123f, -0.905770063f,  -0.862219572f,  -0.715895772f,
      -0.464977026f, -0.169023499f,  0.0486510098f,  0.0368131027f,
      -0.776263058f, -0.937507868f,  -0.959209025f,  -0.821530998f,
      -0.52792424f,  -0.160633594f,  0.113930747f,   0.108118348f,
      -0.510652602f, -0.830184221f,  -0.897609711f,  -0.710870326f,
      -0.331727952f, 0.0935799405f,  0.382687628f,   0.348705471f,
      -0.208708704f, -0.617399991f,  -0.688464046f,  -0.432316244f,
      0.0f,          0.432315558f,   0.688701093f,   0.620223224f,
      0.057524804f,  -0.337877125f,  -0.381433129f,  -0.0935799405f,
      0.331727415f,  0.710868835f,   0.897613525f,   0.830047548f,
      0.242577791f,  -0.0939680934f, -0.114093274f,  0.161937639f,
      0.529363453f,  0.821666479f,   0.959208667f,   0.937391639f,
      0.204559118f,  -0.0717258379f, -0.0662616044f, 0.188465014f,
      0.495982766f,  0.73253417f,    0.864474058f,   0.905769944f,
  };
  FastNoiseLite n;
  n.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  n.SetFrequency(0.125f);

  expect_noise_grid(
      GOLDEN,
      [&](size_t index) {
        const int i = static_cast<int>(index / 8);
        const int j = static_cast<int>(index % 8);
        const float x = static_cast<float>(i) * 0.75f - 3.0f;
        const float y = static_cast<float>(j) * 0.75f - 3.0f;
        return n.GetNoise(x, y);
      },
      "OpenSimplex2 2D");
}

/**
 * @brief Pins noise_transform's displaced output for fixed params and inputs.
 * @details An oracle for the whole noise_transform pipeline (3-channel sampling,
 *          tangent projection, soft-cap, renormalize), not just its on-sphere
 *          invariant: a regression that still produced a unit vector but a
 *          different one would slip past the existing test_transformers checks.
 */
inline void test_noise_transform_golden() {
  Animation::NoiseParams p;
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
