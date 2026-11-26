/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Noise test cases. */

#include "tests/feedback_style_fixture.h"
#include "tests/transformers/fixture.h"

namespace hs_test::transformers_tests {

/**
 * @brief Verifies an active noise warp keeps every sample finite and on the
 *        unit sphere, and actually displaces the points.
 * @details Displacement is summed across samples so a single noise
 *          zero-crossing cannot make the test flaky.
 */
inline void test_noise_active_stays_on_sphere() {
  Animation::NoiseParams p;
  p.amplitude = 0.5f;
  p.scale = 4.0f;
  p.time = 1.0f;
  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                                  math::Vector(0.4f, 0.6f, 0.7f).normalized()};
  float total_moved = 0.0f;
  for (const math::Vector &v : samples) {
    math::Vector r = noise_transform(v, p);
    HS_EXPECT_TRUE(finite_vec(r));
    HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
    total_moved += displacement(r, v);
  }
  HS_EXPECT_GT(total_moved, 1e-2f);
}

/** @brief Authored feedback warp sweeps supply finite unit fold directions. */
inline void test_feedback_warp_inputs_preserve_fold_contract() {
  using Style = Feedback::Style;
  const std::array<math::Vector, 3> *mirror_sets[] = {
      &lenses::TETRAHEDRAL_MIRRORS,     &lenses::OCTAHEDRAL_MIRRORS,
      &lenses::DODECAHEDRAL_MIRRORS,    &lenses::TRIANGULAR_PRISM_MIRRORS,
      &lenses::SQUARE_PRISM_MIRRORS,    &lenses::PENTAGONAL_PRISM_MIRRORS,
      &lenses::HEXAGONAL_PRISM_MIRRORS, &lenses::OCTAGONAL_PRISM_MIRRORS};
  const math::Vector axes[] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                               {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
  Animation::NoiseParams noise;
  for (const auto &from : styles_tests::PRESETS)
    for (const auto &to : styles_tests::PRESETS)
      for (float t : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
        Style style;
        style.noise = &noise;
        style.lerp(from.style, to.style, t);
        for (float time : {0.0f, 17.25f, 255.5f}) {
          noise.time = time;
          for (int sample = 0; sample < 30; ++sample) {
            const math::Vector point =
                sample < 6 ? axes[sample]
                           : math::fib_spiral(24, 0.5f, sample - 6);
            const math::Vector warped = style.space_fn(point, style);
            const bool finite = finite_vec(warped);
            HS_EXPECT_TRUE(finite);
            HS_EXPECT_NEAR(math::dot(warped, warped), 1.0f, 2e-4f);
            if (!finite)
              continue;
            const auto check_fold = [](const math::Vector &folded,
                                       const auto &mirrors) {
              HS_EXPECT_TRUE(finite_vec(folded));
              HS_EXPECT_NEAR(math::dot(folded, folded), 1.0f, 2e-4f);
              for (const math::Vector &normal : mirrors)
                HS_EXPECT_GE(math::dot(folded, normal), -2e-6f);
            };
            for (const auto *mirrors : mirror_sets)
              check_fold(lenses::polyhedral_kaleidoscope_lens(warped, *mirrors),
                         *mirrors);
            check_fold(lenses::dodecahedral_kaleidoscope_lens(warped),
                       lenses::DODECAHEDRAL_MIRRORS);
            check_fold(lenses::dodecahedral_table_lens(warped),
                       lenses::DODECAHEDRAL_MIRRORS);
          }
        }
      }
}

/**
 * @brief Verifies the slide cap bounds shipped noise and binds at a large
 *        amplitude.
 */
inline void test_noise_cross_hemisphere_cap() {
  const auto STYLE = Feedback::Style::LooseWormhole();
  const float capped_angle = std::atan(0.5f);

  auto sweep = [&STYLE](float amplitude) {
    const hs_test::ScopedRng saved_rng(20260803);
    Animation::NoiseParams p;
    p.amplitude = amplitude;
    p.frequency = STYLE.frequency;
    p.scale = STYLE.scale;
    p.speed = STYLE.speed;
    p.time = 0.0f;
    p.sync();
    float worst_angle = 0.0f;
    int off_sphere = 0;
    for (int i = 0; i < 4000; ++i) {
      const float r_x = hs::rand_f(-1.0f, 1.0f);
      const float r_y = hs::rand_f(-1.0f, 1.0f);
      const float r_z = hs::rand_f(-1.0f, 1.0f);
      const math::Vector r(r_x, r_y, r_z);
      if (r.length() < 0.1f)
        continue;
      const math::Vector v = r.normalized();
      const math::Vector out = noise_transform(v, p);
      off_sphere += !finite_vec(out) || std::fabs(out.length() - 1.0f) > 1e-3f;
      worst_angle = std::max(worst_angle, math::angle_between(v, out));
    }
    HS_EXPECT_EQ(off_sphere, 0);
    return worst_angle;
  };

  const float worst_high = sweep(STYLE.amplitude);
  HS_EXPECT_LE(worst_high, capped_angle + 1e-4f);
  HS_EXPECT_NEAR(sweep(1000.0f), capped_angle, 1e-4f);

  // The cap stays slack at this lower amplitude.
  const float worst_low = sweep(0.6f);
  HS_EXPECT_LT(worst_low, capped_angle * 0.5f);
}

// ============================================================================
// NoiseParams::set_seed — the seed mirror tracks the generator it stands for
// ============================================================================

/**
 * @brief Pins NoiseParams::seed against the generator it mirrors.
 * @details FastNoiseLite keeps its seed private, so `seed` is the readable
 *          copy; set_seed must write both, and the defaults must agree.
 */
inline void test_noise_params_seed_mirrors_generator() {
  constexpr float SX = 0.3f, SY = 0.7f, SZ = 1.1f;

  // The declared default must be the generator's own, or an unseeded
  // NoiseParams hashes to a seed it is not using.
  Animation::NoiseParams fresh;
  FastNoiseLite reference;
  reference.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  reference.SetSeed(fresh.seed);
  HS_EXPECT_EQ(fresh.noise.GetNoise(SX, SY, SZ),
               reference.GetNoise(SX, SY, SZ));

  Animation::NoiseParams a, b;
  a.set_seed(4242);
  b.set_seed(4242);
  HS_EXPECT_EQ(a.seed, 4242);
  HS_EXPECT_EQ(a.noise.GetNoise(SX, SY, SZ), b.noise.GetNoise(SX, SY, SZ));

  // A changed mirror must actually reach the generator.
  b.set_seed(4243);
  HS_EXPECT_EQ(b.seed, 4243);
  HS_EXPECT_TRUE(a.noise.GetNoise(SX, SY, SZ) != b.noise.GetNoise(SX, SY, SZ));
}

// ============================================================================
// noise_transform
// ============================================================================

/**
 * @brief Verifies zero amplitude short-circuits the noise warp to a no-op.
 */
inline void test_noise_zero_amplitude_is_identity() {
  Animation::NoiseParams p;
  p.amplitude = 0.0f;
  math::Vector v = math::Vector(0.2f, 0.5f, 0.84f).normalized();
  math::Vector r = noise_transform(v, p);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

} // namespace hs_test::transformers_tests
