/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Scalar cube storage and source-reference tests.

/** @brief Compares cube bakes with direct noise over the same face directions. */
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
              const int index = face * N * N + y * N + x;
              reference[index] = static_cast<int8_t>(
                  sample * 127.0f + (sample < 0.0f ? -0.5f : 0.5f));
#if defined(HS_TEST_FAST_MATH)
              reference_samples[index] = sample * 127.0f;
#endif
            }
          }
        }
        Pattern::prepare_sphere_lut<Pattern::SphereNoise>(
            actual,
            Pattern::SphereNoise::prepare({.scale = scale}, {phase}, noise));
        check_bake();
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
  const Pattern::SphereLutView view{bytes.data(), true, linear.data()};
  for (int face = 0; face < Pattern::SphereLutView::FACE_COUNT; ++face)
    for (float u : {-1.0f, -0.25f, 1.0f})
      for (float v : {-1.0f, 0.25f, 1.0f}) {
        const float expected =
            (face * 10 +
             (u + v + 2) * ((Pattern::SphereLutView::FACE_STEPS - 1) * 0.5f)) /
            127.0f;
        HS_EXPECT_NEAR(Pattern::sample_sphere_lut_face(view, face, u, v),
                       expected, 2e-7f);
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
  HS_EXPECT_FALSE(cache.holds(nullptr, {}));
  FastNoiseLite noise;
  noise.SetFrequency(1.0f);
  const auto prepared = Pattern::SphereNoise::prepare({}, {}, noise);
  a.fill(-128);
  HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(a, key, prepared));
  for (auto value : a)
    HS_EXPECT_NE(value, -128);
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
    a.fill(-128);
    HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(a, changed, prepared));
    for (auto value : a)
      HS_EXPECT_NE(value, -128);
    a.fill(-128);
    HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(a, key, prepared));
    for (auto value : a)
      HS_EXPECT_NE(value, -128);
  }
  b.fill(-128);
  HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(b, key, prepared));
  for (auto value : b)
    HS_EXPECT_NE(value, -128);
  cache.reset();
  b.fill(-128);
  HS_EXPECT_TRUE(cache.refresh<Pattern::SphereNoise>(b, key, prepared));
  for (auto value : b)
    HS_EXPECT_NE(value, -128);
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
        const SphereLutView view{table.data(), true};
        for (int a : {-1, 1})
          for (int b : {-1, 1})
            for (int index = 0; index <= 200; ++index) {
              const float t = -1.0f + 0.01f * index;
              const math::Vector xy(a, b, t), xz(a, t, b), yz(t, a, b);
              const float xy_error =
                  fabsf(sphere_edge_face_sample(view, a > 0 ? 0 : 1, xy) -
                        sphere_edge_face_sample(view, b > 0 ? 2 : 3, xy));
              const float xz_error =
                  fabsf(sphere_edge_face_sample(view, a > 0 ? 0 : 1, xz) -
                        sphere_edge_face_sample(view, b > 0 ? 4 : 5, xz));
              const float yz_error =
                  fabsf(sphere_edge_face_sample(view, a > 0 ? 2 : 3, yz) -
                        sphere_edge_face_sample(view, b > 0 ? 4 : 5, yz));
              maximum = std::max({maximum, xy_error, xz_error, yz_error});
              HS_EXPECT_LE(xy_error, 1.0e-5f);
              HS_EXPECT_LE(xz_error, 1.0e-5f);
              HS_EXPECT_LE(yz_error, 1.0e-5f);
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
  const std::array<NoiseParams, 9> recipes{{
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
  constexpr std::array<float, 9> SCALES{
      1.0f / 64.0f, 0.0625f, 0.125f, 0.5f, 1.0f, 2.0f, 3.0f, 5.0f, 8.0f};
  // Mean baked error at any tested scale.
  constexpr float MEAN_BUDGET = 0.6f;
  std::array<int8_t, SphereLutView::SIZE> table;
  std::array<float, SphereLutView::SIZE> exact;
#if defined(HS_PATTERN_QUALIFICATION_REPORT)
  std::printf(
      "recipe,scale,max_quantization,max_spatial,max_total,mean_total\n");
#endif
  for (size_t recipe = 0; recipe < recipes.size(); ++recipe) {
    double previous_mean = 0.0;
    for (size_t scale_index = 0; scale_index < SCALES.size(); ++scale_index) {
      const float scale = SCALES[scale_index];
      HS_CONTEXT("recipe / scale microunits", recipe,
                 static_cast<int>(scale * 1000000));
      float max_quantization = 0.0f, max_spatial = 0.0f, max_total = 0.0f;
      double total = 0.0;
      for (int seed : {-719, 1337, 6047})
        for (float phase : {0.0f, 0.173f, 0.73f}) {
          FastNoiseLite generator;
          generator.SetSeed(seed);
          generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
          generator.SetFrequency(1.0f);
          auto params = recipes[recipe];
          params.scale = scale;
          const auto prepared =
              SphereNoise::prepare(params, {phase}, generator);
          prepare_sphere_lut<SphereNoise>(table, prepared);
          for (int face = 0; face < SphereLutView::FACE_COUNT; ++face)
            for (int y = 0; y < SphereLutView::FACE_STEPS; ++y)
              for (int x = 0; x < SphereLutView::FACE_STEPS; ++x) {
                const int index = face * SphereLutView::FACE_SIZE +
                                  y * SphereLutView::FACE_STEPS + x;
                const auto d =
                    sphere_face_coordinate(face, sphere_face_grid_coordinate(x),
                                           sphere_face_grid_coordinate(y));
                const float value = SphereNoise::sample(d, prepared);
                exact[index] = value * 127.0f;
                max_quantization = std::max(
                    max_quantization, fabsf(value - table[index] / 127.0f));
              }
          const SphereLutView view{table.data(), true};
          const SphereLutView float_view{table.data(), true, exact.data()};
          for (int index = 0; index < SAMPLE_COUNT; ++index) {
            const float z = 1.0f - 2.0f * (index + 0.5f) / SAMPLE_COUNT;
            const float r = sqrtf(1.0f - z * z);
            const float angle =
                math::TWO_PI_F * math::wrap_t(index * 0.61803398875f);
            const auto d =
                SphereDirection::from({r * cosf(angle), r * sinf(angle), z});
            const float direct = SphereNoise::sample(d, prepared);
            const float baked = sample_sphere_lut(view, d);
            const float spatial = sample_sphere_lut(float_view, d);
            max_quantization =
                std::max(max_quantization, fabsf(baked - spatial));
            max_spatial = std::max(max_spatial, fabsf(direct - spatial));
            max_total = std::max(max_total, fabsf(direct - baked));
            total += fabsf(direct - baked);
          }
        }
      // Quantization dominates at low scale; above it interpolation error
      // grows with feature frequency.
      const double mean = total / (9 * SAMPLE_COUNT);
      HS_EXPECT_LE(mean, MEAN_BUDGET);
      if (scale <= 0.125f)
        HS_EXPECT_LE(mean, QUANTIZATION_BOUND);
      else
        HS_EXPECT_GT(mean, previous_mean);
      previous_mean = mean;
      HS_EXPECT_LE(max_quantization, QUANTIZATION_BOUND);
      constexpr std::array<float, 9> TOTAL_BOUNDS{
          1.5f, 1.75f, 1.9f, 1.2f, 1.55f, 1.85f, 1.6f, 1.85f, 1.96f};
      HS_EXPECT_LE(max_total, TOTAL_BOUNDS[recipe]);
      constexpr std::array<float, 9> LOW_SCALE_BOUNDS{0.005f, 0.010f, 0.040f,
                                                      0.006f, 0.012f, 0.045f,
                                                      0.040f, 0.080f, 0.220f};
      if (scale <= 0.125f) {
        HS_EXPECT_LE(max_total, LOW_SCALE_BOUNDS[recipe]);
        HS_EXPECT_LE(max_spatial, LOW_SCALE_BOUNDS[recipe]);
      }
#if defined(HS_PATTERN_QUALIFICATION_REPORT)
      std::printf("%zu,%.9g,%.9g,%.9g,%.9g,%.9g\n", recipe, scale,
                  max_quantization, max_spatial, max_total,
                  total / (9 * SAMPLE_COUNT));
#else
      (void)max_spatial;
#endif
    }
  }
}
