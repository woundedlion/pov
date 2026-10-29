/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_lut.h
 * @brief Scalar cube storage and source-reference tests. */

/** @brief Compares paired cube-face bakes with independent scalar face walks. */
inline void test_sphere_bakes_match_reference() {
  static std::array<int8_t, Pattern::SphereLutView::SIZE> actual, reference;
#if defined(HS_TEST_FAST_MATH)
  static std::array<float, Pattern::SphereLutView::SIZE> reference_samples;
#endif
  const auto check_bake = [&] {
#if defined(HS_TEST_FAST_MATH)
    // Fast-math reassociation can cross an int8 rounding boundary.
    constexpr float NOISE_TOLERANCE = 2e-5f;
    constexpr float ROUNDING_TOLERANCE = 0.5f + 127.0f * NOISE_TOLERANCE;
    for (size_t i = 0; i < actual.size(); ++i)
      HS_EXPECT_NEAR(static_cast<float>(actual[i]), reference_samples[i],
                     ROUNDING_TOLERANCE);
#else
    HS_EXPECT_EQ(std::memcmp(actual.data(), reference.data(), actual.size()),
                 0);
#endif
  };
  for (int variant = 0; variant < 4; ++variant) {
    FastNoiseLite noise;
    noise.SetSeed(6047 + variant);
    noise.SetNoiseType(variant == 2 ? FastNoiseLite::NoiseType_Perlin
                                    : FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetFrequency(variant == 1 ? 0.7f : 1.0f);
    if (variant == 3)
      noise.SetRotationType3D(FastNoiseLite::RotationType3D_ImproveXYPlanes);
    for (float scale : {0.0625f, 2.0f, 5.0f}) {
      for (float phase : {0.0f, 0.25f, 0.625f}) {
        const math::Vector offset = math::noise_sphere_loop_offset(phase);
        constexpr int N = Pattern::SphereLutView::FACE_STEPS;
        constexpr float STEP = 2.0f / (N - 1);
        for (int face = 0; face < Pattern::SphereLutView::FACE_COUNT; ++face) {
          for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
              const math::Vector d = Pattern::sphere_face_direction(
                  face, -1.0f + STEP * x, -1.0f + STEP * y);
              const math::Vector q = scale * d + offset;
              const float sample =
                  hs::clamp(noise.GetNoiseSingle(q.x, q.y, q.z), -1.0f, 1.0f);
              const int INDEX = face * N * N + y * N + x;
              reference[INDEX] = static_cast<int8_t>(
                  sample * 127.0f + (sample < 0.0f ? -0.5f : 0.5f));
#if defined(HS_TEST_FAST_MATH)
              reference_samples[INDEX] = sample * 127.0f;
#endif
            }
          }
        }
        Pattern::prepare_sphere_lut<Pattern::SphereNoise>(
            actual,
            Pattern::SphereNoise::prepare({.scale = scale}, {phase}, noise));
        check_bake();
        if (variant == 0) {
          Pattern::prepare_sphere_lut<Pattern::SphereNoise>(
              actual,
              Pattern::SphereNoise::prepare({.scale = scale}, {phase}, noise));
          check_bake();
        }
      }
    }
  }
}

/**
 * @brief Verifies the cube-map hue-noise LUT agrees across a shared face edge
 *        and at a three-face corner, where the face tie-break picks different
 *        faces for neighbouring directions.
 */
inline void test_sphere_lut_seamless_across_faces() {
  static std::array<int8_t, Pattern::SphereLutView::SIZE> hue_noise;
  FastNoiseLite noise;
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(1301);
  noise.SetFrequency(1.0f);
  Pattern::prepare_sphere_lut<Pattern::SphereNoise>(
      hue_noise, Pattern::SphereNoise::prepare({.scale = 3.0f}, {0.0f}, noise));
  const Pattern::SphereLutView view{hue_noise.data(), true};

  constexpr float NUDGE = 1e-4f;
  float edge_low = 1.0f;
  float edge_high = -1.0f;
  for (int step = 0; step <= 16; ++step) {
    const float t = -1.0f + 0.125f * step;
    const float x_face = Pattern::sample_sphere_lut(
        view, math::Vector(1.0f, 1.0f - NUDGE, t).normalized());
    const float y_face = Pattern::sample_sphere_lut(
        view, math::Vector(1.0f - NUDGE, 1.0f, t).normalized());
    HS_EXPECT_NEAR(x_face, y_face, 2e-3f);
    edge_low = std::min(edge_low, x_face);
    edge_high = std::max(edge_high, x_face);
  }
  HS_EXPECT_GT(edge_high - edge_low, 1e-2f);

  const float corner_x = Pattern::sample_sphere_lut(
      view, math::Vector(1.0f, 1.0f - NUDGE, 1.0f - NUDGE).normalized());
  const float corner_y = Pattern::sample_sphere_lut(
      view, math::Vector(1.0f - NUDGE, 1.0f, 1.0f - NUDGE).normalized());
  const float corner_z = Pattern::sample_sphere_lut(
      view, math::Vector(1.0f - NUDGE, 1.0f - NUDGE, 1.0f).normalized());
  HS_EXPECT_NEAR(corner_x, corner_y, 2e-3f);
  HS_EXPECT_NEAR(corner_x, corner_z, 2e-3f);
}

inline void test_hue_rotation_lut_clamps_out_of_range_value() {
  static std::array<Pixel, 3 * HueRotationLutView::SIZE> storage;
  storage.fill(Pixel(0, 0, 65535));
  Pixel *const lut = storage.data() + HueRotationLutView::SIZE;
  for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue) {
    lut[hue] = Pixel(65535, 0, 0);
    lut[(HueRotationLutView::VALUE_STEPS - 1) * HueRotationLutView::HUE_STEPS +
        hue] = Pixel(0, 65535, 0);
  }

  const HueRotationLutView view{lut, true};
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, -0.5f, 0.25f), Pixel(65535, 0, 0));
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, 1.5f, 0.25f), Pixel(0, 65535, 0));
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, NAN, 0.25f), Pixel(0, 65535, 0));
  for (float amount : {-4.75f, -1e-7f, 3.5e8f, -3.5e8f, INFINITY, NAN})
    HS_EXPECT_EQ(sample_hue_rotation_lut(view, 0.0f, amount),
                 Pixel(65535, 0, 0));
}

/** @brief Face-coordinate noise taps use the float table when supplied. */
inline void test_sphere_face_linear_table() {
  static std::array<int8_t, Pattern::SphereLutView::SIZE> bytes{};
  static std::array<float, Pattern::SphereLutView::SIZE> linear;
  for (int face = 0; face < Pattern::SphereLutView::FACE_COUNT; ++face)
    for (int y = 0; y < Pattern::SphereLutView::FACE_STEPS; ++y)
      for (int x = 0; x < Pattern::SphereLutView::FACE_STEPS; ++x)
        linear[face * Pattern::SphereLutView::FACE_SIZE +
               y * Pattern::SphereLutView::FACE_STEPS + x] =
            float(face * 10 + x + y);
  const Pattern::SphereLutView VIEW{bytes.data(), true, linear.data()};
  for (int face = 0; face < Pattern::SphereLutView::FACE_COUNT; ++face)
    for (float u : {-1.0f, -0.25f, 1.0f})
      for (float v : {-1.0f, 0.25f, 1.0f}) {
        const float EXPECTED = (face * 10 + (u + v + 2) * 11.5f) / 127.0f;
        HS_EXPECT_NEAR(Pattern::sample_sphere_lut_face(VIEW, face, u, v),
                       EXPECTED, 2e-7f);
      }
}

/** @brief Complete cache identity invalidates source and destination changes. */
inline void test_sphere_bake_cache() {
  std::array<int8_t, Pattern::SphereLutView::SIZE> a, b;
  struct Identity {
    int seed;
    float scale;
    float phase;
    int generation;
    int recipe;
    bool operator==(const Identity &) const = default;
  };
  Identity key{1337, 1.0f, 0.0f, 1, 1};
  Pattern::SphereBakeCache<Identity> cache;
  FastNoiseLite noise;
  noise.SetFrequency(1.0f);
  const auto prepared = Pattern::SphereNoise::prepare({}, {}, noise);
  HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(a, key, prepared));
  a.fill(-128);
  HS_EXPECT_TRUE(!cache.refresh<Pattern::SphereNoise>(a, key, prepared));
  HS_EXPECT_EQ(a.front(), -128);
  for (int component = 0; component < 5; ++component) {
    auto changed = key;
    switch (component) {
    case 0:
      ++changed.seed;
      break;
    case 1:
      changed.scale = 2;
      break;
    case 2:
      changed.phase = .5f;
      break;
    case 3:
      ++changed.generation;
      break;
    case 4:
      ++changed.recipe;
      break;
    }
    HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(a, changed, prepared));
    HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(a, key, prepared));
  }
  HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(b, key, prepared));
  cache.reset();
  HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(b, key, prepared));
}
