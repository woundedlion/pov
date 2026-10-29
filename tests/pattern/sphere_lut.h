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
        constexpr float HALF_STEP = 1.0f / (N - 1);
        for (int face = 0; face < Pattern::SphereLutView::FACE_COUNT; ++face) {
          for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
              const math::Vector d = Pattern::sphere_face_direction(
                  face, (2 * x - (N - 1)) * HALF_STEP,
                  (2 * y - (N - 1)) * HALF_STEP);
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

/** @brief Samples a chosen face at the ratios of a shared edge direction. */
inline float sphere_edge_face_sample(const Pattern::SphereLutView &view,
                                     int face, const math::Vector &direction) {
  float u = 0.0f;
  float v = 0.0f;
  switch (face) {
  case 0:
    u = direction.z;
    v = direction.y;
    break;
  case 1:
    u = -direction.z;
    v = direction.y;
    break;
  case 2:
    u = direction.x;
    v = direction.z;
    break;
  case 3:
    u = direction.x;
    v = -direction.z;
    break;
  case 4:
    u = direction.x;
    v = direction.y;
    break;
  case 5:
    u = -direction.x;
    v = direction.y;
    break;
  }
  return Pattern::sample_sphere_lut_face(view, face, u, v);
}

/** @brief Opposite-index symmetry keeps every cube edge continuous. */
inline void test_sphere_lut_all_edges() {
  using namespace Pattern;
  for (int index = 0; index < SphereLutView::FACE_STEPS; ++index)
    HS_EXPECT_EQ(
        sphere_face_grid_coordinate(index),
        -sphere_face_grid_coordinate(SphereLutView::FACE_STEPS - 1 - index));
  std::array<int8_t, SphereLutView::SIZE> table;
  float maximum = 0.0f;
  for (int seed : {-719, 1337, 6047})
    for (float phase : {0.0f, 0.173f, 0.73f}) {
      FastNoiseLite generator;
      generator.SetSeed(seed);
      generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
      generator.SetFrequency(1.0f);
      for (float scale : {3.0f, 64.0f}) {
        prepare_sphere_lut<SphereNoise>(
            table, SphereNoise::prepare({.scale = scale}, {phase}, generator));
        const SphereLutView VIEW{table.data(), true};
        for (int a : {-1, 1})
          for (int b : {-1, 1})
            for (int index = 0; index <= 200; ++index) {
              const float T = -1.0f + 0.01f * index;
              const math::Vector XY(a, b, T), XZ(a, T, b), YZ(T, a, b);
              const float XY_ERROR =
                  fabsf(sphere_edge_face_sample(VIEW, a > 0 ? 0 : 1, XY) -
                        sphere_edge_face_sample(VIEW, b > 0 ? 2 : 3, XY));
              const float XZ_ERROR =
                  fabsf(sphere_edge_face_sample(VIEW, a > 0 ? 0 : 1, XZ) -
                        sphere_edge_face_sample(VIEW, b > 0 ? 4 : 5, XZ));
              const float YZ_ERROR =
                  fabsf(sphere_edge_face_sample(VIEW, a > 0 ? 2 : 3, YZ) -
                        sphere_edge_face_sample(VIEW, b > 0 ? 4 : 5, YZ));
              maximum = std::max({maximum, XY_ERROR, XZ_ERROR, YZ_ERROR});
              HS_EXPECT_LE(XY_ERROR, 1.0e-5f);
              HS_EXPECT_LE(XZ_ERROR, 1.0e-5f);
              HS_EXPECT_LE(YZ_ERROR, 1.0e-5f);
            }
      }
    }
#if defined(HS_PATTERN_QUALIFICATION_REPORT)
  std::printf("cube_seam_maximum,%.9g\n", maximum);
#else
  (void)maximum;
#endif
}

/** @brief Seeded source-only qualification separates quantization and interpolation. */
inline void test_sphere_noise_bake_qualification() {
  using namespace Pattern;
  constexpr int SAMPLE_COUNT = 4096;
  constexpr float QUANTIZATION_BOUND = 1.0f / 254.0f + 3.0e-6f;
  const std::array<NoiseParams, 9> RECIPES{{
      {},
      {.contrast = 2.0f, .variant = NoiseVariant::OCTAVE_CONTOUR},
      {.contrast = 8.0f, .variant = NoiseVariant::OCTAVE_CONTOUR},
      {.basis = NoiseBasis::FBM3, .variant = NoiseVariant::OCTAVE_CONTOUR},
      {.basis = NoiseBasis::FBM3,
       .contrast = 2.0f,
       .variant = NoiseVariant::OCTAVE_CONTOUR},
      {.basis = NoiseBasis::FBM3,
       .contrast = 8.0f,
       .variant = NoiseVariant::OCTAVE_CONTOUR},
      {.basis = NoiseBasis::RIDGED3, .variant = NoiseVariant::OCTAVE_CONTOUR},
      {.basis = NoiseBasis::RIDGED3,
       .contrast = 2.0f,
       .variant = NoiseVariant::OCTAVE_CONTOUR},
      {.basis = NoiseBasis::RIDGED3,
       .contrast = 8.0f,
       .variant = NoiseVariant::OCTAVE_CONTOUR},
  }};
  std::array<int8_t, SphereLutView::SIZE> table;
  std::array<float, SphereLutView::SIZE> exact;
#if defined(HS_PATTERN_QUALIFICATION_REPORT)
  std::printf(
      "recipe,scale,max_quantization,max_spatial,max_total,mean_total\n");
#endif
  for (size_t recipe = 0; recipe < RECIPES.size(); ++recipe)
    for (float scale :
         {1.0f / 64.0f, 0.0625f, 0.125f, 0.5f, 1.0f, 2.0f, 3.0f, 5.0f, 8.0f}) {
      float max_quantization = 0.0f, max_spatial = 0.0f, max_total = 0.0f;
      double total = 0.0;
      for (int seed : {-719, 1337, 6047})
        for (float phase : {0.0f, 0.173f, 0.73f}) {
          FastNoiseLite generator;
          generator.SetSeed(seed);
          generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
          generator.SetFrequency(1.0f);
          auto params = RECIPES[recipe];
          params.scale = scale;
          const auto PREPARED =
              SphereNoise::prepare(params, {phase}, generator);
          prepare_sphere_lut<SphereNoise>(table, PREPARED);
          for (int face = 0; face < SphereLutView::FACE_COUNT; ++face)
            for (int y = 0; y < SphereLutView::FACE_STEPS; ++y)
              for (int x = 0; x < SphereLutView::FACE_STEPS; ++x) {
                const int INDEX = face * SphereLutView::FACE_SIZE +
                                  y * SphereLutView::FACE_STEPS + x;
                const auto D =
                    sphere_face_coordinate(face, sphere_face_grid_coordinate(x),
                                           sphere_face_grid_coordinate(y));
                const float VALUE = SphereNoise::sample(D, PREPARED);
                exact[INDEX] = VALUE * 127.0f;
                max_quantization = std::max(
                    max_quantization, fabsf(VALUE - table[INDEX] / 127.0f));
              }
          const SphereLutView VIEW{table.data(), true};
          const SphereLutView FLOAT_VIEW{table.data(), true, exact.data()};
          for (int index = 0; index < SAMPLE_COUNT; ++index) {
            const float Z = 1.0f - 2.0f * (index + 0.5f) / SAMPLE_COUNT;
            const float R = sqrtf(1.0f - Z * Z);
            const float ANGLE =
                math::TWO_PI_F * math::wrap_t(index * 0.61803398875f);
            const auto D =
                SphereDirection::from({R * cosf(ANGLE), R * sinf(ANGLE), Z});
            const float DIRECT = SphereNoise::sample(D, PREPARED);
            const float BAKED = sample_sphere_lut(VIEW, D);
            const float SPATIAL = sample_sphere_lut(FLOAT_VIEW, D);
            max_quantization =
                std::max(max_quantization, fabsf(BAKED - SPATIAL));
            max_spatial = std::max(max_spatial, fabsf(DIRECT - SPATIAL));
            max_total = std::max(max_total, fabsf(DIRECT - BAKED));
            total += fabsf(DIRECT - BAKED);
          }
        }
      HS_EXPECT_LE(max_quantization, QUANTIZATION_BOUND);
      constexpr std::array<float, 9> TOTAL_BOUNDS{
          1.5f, 1.75f, 1.9f, 1.2f, 1.55f, 1.85f, 1.6f, 1.85f, 1.96f};
      HS_EXPECT_LE(max_total, TOTAL_BOUNDS[recipe]);
      if (recipe == 0 && scale <= 0.125f)
        HS_EXPECT_LE(max_total, 0.005f);
#if defined(HS_PATTERN_QUALIFICATION_REPORT)
      std::printf("%zu,%.9g,%.9g,%.9g,%.9g,%.9g\n", recipe, scale,
                  max_quantization, max_spatial, max_total,
                  total / (9 * SAMPLE_COUNT));
#else
      (void)max_spatial;
      (void)total;
#endif
    }
}
