/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Pullback periodic and burst ripple preparation and sampling tests.
 */

/** @brief Supplied centers preserve the periodic wave's other prepared fields. */
inline void test_pullback_ripple_supplied_center() {
  Pullback::Surface::PeriodicRippleParams params;
  params.center_azimuth = 0.73f;
  params.center_polar = 1.17f;
  params.strength = 0.15f;
  params.thickness = 0.7f;
  params.decay = 0.2f;
  const auto CENTER =
      math::Vector::from_spherical(params.center_azimuth, params.center_polar);
  const auto OTHER = math::Vector(0.2f, 0.8f, -0.3f).normalized();
  for (float cycle : {-0.25f, 0.0f, 0.01f, 0.5f, 1.0f, 1.25f}) {
    const auto FIXED = Pullback::Surface::prepare_ripple(params, cycle);
    const auto SUPPLIED =
        Pullback::Surface::prepare_ripple(params, cycle, CENTER);
    HS_EXPECT_EQ(
        std::memcmp(&FIXED.ripple, &SUPPLIED.ripple, sizeof(FIXED.ripple)), 0);
    auto expected = FIXED;
    expected.ripple.center = OTHER;
    const auto ACTUAL = Pullback::Surface::prepare_ripple(params, cycle, OTHER);
    HS_EXPECT_EQ(
        std::memcmp(&expected.ripple, &ACTUAL.ripple, sizeof(expected.ripple)),
        0);
  }
}

/** @brief Sparse and stale burst slots retain canonical prepared waves. */
inline void test_pullback_ripple_burst_centers() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  params.period = 91.48529f;
  params.strength = 0.15f;
  params.thickness = 0.7f;
  params.decay = 0.2f;
  RippleState state;
  state.seed = 0x12345678;
  state.burst_count = 4;
  state.waves[4] = {5.0f, 1};
  state.waves[5] = {10.0f, 1};
  state.waves[8] = {params.period, 2};
  state.waves[9] = {params.period - 1.0f, 2};
  state.waves[10] = {params.period + 1.0f, 2};
  state.waves[11] = {4.0f, 2};
  state.waves[12] = {1.0f, 99};
  for (bool random : {false, true}) {
    params.random_center = random;
    const auto ACTUAL = prepare_ripples(params, state);
    PreparedRipples reference;
    for (uint32_t lag = 4; lag > 0; --lag) {
      const uint32_t BURST = state.burst_count - lag;
      for (size_t i = 0; i < 4; ++i) {
        const auto &wave = state.waves[(BURST % 4) * 4 + i];
        if (wave.age < 0.0f || wave.age >= params.period || wave.burst != BURST)
          continue;
        auto prepared =
            Pullback::Surface::prepare_ripple(params, wave.age / params.period);
        if (random)
          prepared.ripple.center = ripple_center(state.seed, BURST);
        reference.waves[reference.count++] = prepared;
      }
    }
    HS_EXPECT_EQ(ACTUAL.count, uint8_t{4});
    HS_EXPECT_EQ(ACTUAL.count, reference.count);
    for (uint8_t i = 0; i < ACTUAL.count; ++i)
      HS_EXPECT_EQ(std::memcmp(&ACTUAL.waves[i], &reference.waves[i],
                               sizeof(ACTUAL.waves[i])),
                   0);
    static RadialTables tables;
    prepare_radial_tables(params, state, tables);
    HS_EXPECT_EQ(tables.count, uint8_t{2});
    for (uint8_t i = 0; i < tables.count; ++i) {
      const auto EXPECTED =
          random ? ripple_center(state.seed, i + 1)
                 : Pullback::Surface::prepare_ripple_center(params);
      HS_EXPECT_EQ(tables.groups[i].center, EXPECTED);
    }
  }
}

/** @brief Ripple phase repeats after its configured period. */
inline void test_pullback_periodic_ripple() {
  Pullback::Surface::PeriodicRippleParams clock;
  clock.period = 12.0f;
  HS_EXPECT_EQ(Pullback::Surface::ripple_cycle(0.0f, clock), 0.0f);
  HS_EXPECT_EQ(Pullback::Surface::ripple_cycle(6.0f, clock), 0.5f);
  HS_EXPECT_EQ(Pullback::Surface::ripple_cycle(12.0f, clock), 1.0f);
  Pullback::Surface::PeriodicRippleParams params;
  params.strength = 0.15f;
  params.decay = 0.0f;
  params.thickness = 0.7f;
  params.center_polar = 0.0f;

  const math::Vector midpoint =
      math::Vector::from_spherical(0.0f, 0.5f * math::PI_F);
  const Pullback::SurfaceResult start = Pullback::Surface::periodic_ripple(
      midpoint, Pullback::Surface::prepare_ripple(params, 0.0f), true);
  HS_EXPECT_EQ(start.sphere.x, midpoint.x);
  HS_EXPECT_EQ(start.sphere.y, midpoint.y);
  HS_EXPECT_EQ(start.sphere.z, midpoint.z);
  HS_EXPECT_EQ(start.path_length, 0.0f);

  const Pullback::SurfaceResult crest = Pullback::Surface::periodic_ripple(
      midpoint, Pullback::Surface::prepare_ripple(params, 0.5f), true);
  HS_EXPECT_NEAR(crest.sphere.length(), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(crest.path_length, 0.5f * params.strength, 2e-3f);

  const Pullback::SurfaceResult wrapped = Pullback::Surface::periodic_ripple(
      midpoint, Pullback::Surface::prepare_ripple(params, 1.0f), true);
  HS_EXPECT_EQ(wrapped.sphere.x, start.sphere.x);
  HS_EXPECT_EQ(wrapped.sphere.y, start.sphere.y);
  HS_EXPECT_EQ(wrapped.sphere.z, start.sphere.z);
  HS_EXPECT_EQ(wrapped.path_length, start.path_length);
}

/** @brief Builds one explicitly phased ripple wave.
 * @param center Unit ripple center.
 * @param phase Wave crest distance in radians.
 * @param amplitude Peak displacement in radians.
 * @return Prepared ripple with refreshed rejection bounds. */
inline Pullback::Surface::PreparedRipple
ripple_test_wave(const math::Vector &center, float phase, float amplitude) {
  Animation::RippleWave wave;
  wave.center = center;
  wave.phase = phase;
  wave.amplitude = amplitude;
  wave.decay = 0.0f;
  wave.thickness = 0.7f;
  wave.sync();
  return {wave};
}

/** @brief Applies the individual ripple kernels in their supplied order.
 * @param view Unit sample direction.
 * @param prepared Active waves.
 * @param support_margin Optional closest distance to a wave support boundary.
 * @return Displaced direction. */
inline math::Vector ripple_test_sequential(
    math::Vector view,
    const Pullback::SphereAnimation::PreparedRipples &prepared,
    float *support_margin = nullptr) {
  for (uint8_t i = 0; i < prepared.count; ++i) {
    const auto &wave = prepared.waves[i].ripple;
    if (support_margin && wave.amplitude > 0.001f) {
      const float DISTANCE =
          acosf(hs::clamp(math::dot(view, wave.center), -1.0f, 1.0f));
      for (float edge : {wave.phase - 2.0f * wave.half_width(),
                         wave.phase + 2.0f * wave.half_width()})
        if (edge > 0.0f && edge < math::PI_F)
          *support_margin = fminf(*support_margin, fabsf(DISTANCE - edge));
    }
    view = Pullback::Surface::periodic_ripple(view, prepared.waves[i], false)
               .sphere;
  }
  return view;
}

/** @brief Group rotations preserve the analytic rotation within their domain. */
inline void test_pullback_group_rotation_polynomial() {
  static_assert(RippleKernel::GROUP_ROTATION_POLYNOMIAL_MAX == 0.5f);
  const float OUTSIDE = std::nextafter(0.5f, 1.0f);
  const float ANGLES[] = {-0.9f,  -OUTSIDE, -0.5f, -0.49999f, -0.2f,
                          -0.15f, -0.05f,   0.0f,  0.05f,     0.15f,
                          0.2f,   0.49999f, 0.5f,  OUTSIDE,   0.9f};
  for (int row = 0; row < 13; ++row)
    for (int column = 0; column < 17; ++column) {
      const float Y = -0.99f + 1.98f * row / 12;
      const float ANGLE = math::TWO_PI_F * column / 17;
      const float RADIUS = sqrtf(1.0f - Y * Y);
      const math::Vector START(RADIUS * cosf(ANGLE), Y, RADIUS * sinf(ANGLE));
      const math::Vector CENTER = math::Vector(cosf(ANGLE * 0.37f + 0.2f), 0.3f,
                                               sinf(ANGLE * 0.71f + 0.4f))
                                      .normalized();
      const auto CROSS = math::cross(CENTER, START);
      const float CROSS_SQ = math::dot(CROSS, CROSS);
      if (CROSS_SQ < 1e-4f)
        continue;
      const auto AXIS = CROSS * (1.0f / sqrtf(CROSS_SQ));
      for (float angle : ANGLES) {
        math::Vector output;
        HS_EXPECT_TRUE(
            RippleKernel::rotate_group(START, CENTER, angle, output));
        const auto REFERENCE = RippleKernel::rotate(START, AXIS, angle);
        const auto DELTA = output - REFERENCE;
        HS_EXPECT_LE(math::dot(DELTA, DELTA), 1e-12f);
        HS_EXPECT_NEAR(math::dot(output, output), 1.0f, 1e-6f);
        const double COSINE = std::cos(double(angle));
        const double SINE = std::sin(double(angle));
        const double DOT = double(AXIS.x) * START.x + double(AXIS.y) * START.y +
                           double(AXIS.z) * START.z;
        const double expected[] = {
            START.x * COSINE +
                (double(AXIS.y) * START.z - double(AXIS.z) * START.y) * SINE +
                AXIS.x * DOT * (1.0 - COSINE),
            START.y * COSINE +
                (double(AXIS.z) * START.x - double(AXIS.x) * START.z) * SINE +
                AXIS.y * DOT * (1.0 - COSINE),
            START.z * COSINE +
                (double(AXIS.x) * START.y - double(AXIS.y) * START.x) * SINE +
                AXIS.z * DOT * (1.0 - COSINE)};
        const double ERROR_SQ =
            (output.x - expected[0]) * (output.x - expected[0]) +
            (output.y - expected[1]) * (output.y - expected[1]) +
            (output.z - expected[2]) * (output.z - expected[2]);
        HS_EXPECT_LE(ERROR_SQ, 1e-12);
        if (fabsf(angle) > 0.5f) {
          HS_EXPECT_EQ(output.x, REFERENCE.x);
          HS_EXPECT_EQ(output.y, REFERENCE.y);
          HS_EXPECT_EQ(output.z, REFERENCE.z);
        }
      }
    }
}

/** @brief Empty ripple batches preserve the direction and zero travel. */
inline void test_pullback_grouped_ripple_identity() {
  const Pullback::SphereAnimation::PreparedRipples EMPTY{};
  const math::Vector VIEW = math::Vector(0.3f, 0.4f, 0.5f).normalized();
  for (bool track_path : {false, true}) {
    const auto RESULT =
        Pullback::SphereAnimation::grouped_ripples(VIEW, EMPTY, track_path);
    HS_EXPECT_EQ(RESULT.sphere.x, VIEW.x);
    HS_EXPECT_EQ(RESULT.sphere.y, VIEW.y);
    HS_EXPECT_EQ(RESULT.sphere.z, VIEW.z);
    HS_EXPECT_EQ(RESULT.path_length, 0.0f);
  }
}

/** @brief A crest rotates outward by its amplitude and reports that arc. */
inline void test_pullback_grouped_ripple_analytic_crest() {
  Pullback::SphereAnimation::PreparedRipples prepared{};
  prepared.count = 1;
  prepared.waves[0] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.075f);
  const math::Vector EXPECTED(cosf(0.075f), -sinf(0.075f), 0.0f);
  const auto TRACKED =
      Pullback::SphereAnimation::grouped_ripples(math::X_AXIS, prepared, true);
  const auto UNTRACKED =
      Pullback::SphereAnimation::grouped_ripples(math::X_AXIS, prepared, false);
  HS_EXPECT_NEAR(TRACKED.sphere.x, EXPECTED.x, 1e-6f);
  HS_EXPECT_NEAR(TRACKED.sphere.y, EXPECTED.y, 1e-6f);
  HS_EXPECT_NEAR(TRACKED.sphere.z, EXPECTED.z, 1e-6f);
  HS_EXPECT_NEAR(TRACKED.path_length, 0.075f, 1e-6f);
  HS_EXPECT_GT(TRACKED.path_length, 0.07f);
  HS_EXPECT_EQ(UNTRACKED.sphere.x, TRACKED.sphere.x);
  HS_EXPECT_EQ(UNTRACKED.sphere.y, TRACKED.sphere.y);
  HS_EXPECT_EQ(UNTRACKED.sphere.z, TRACKED.sphere.z);
  HS_EXPECT_EQ(UNTRACKED.path_length, 0.0f);
}

/** @brief Mixed-center batches preserve the launch order of noncommuting waves. */
inline void test_pullback_grouped_ripple_order() {
  Pullback::SphereAnimation::PreparedRipples prepared{};
  prepared.count = 4;
  prepared.waves[0] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.13f);
  prepared.waves[1] = ripple_test_wave(math::Z_AXIS, math::PI_F * 0.5f, 0.11f);
  prepared.waves[2] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.09f);
  prepared.waves[3] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.07f);
  const math::Vector EXPECTED = ripple_test_sequential(math::X_AXIS, prepared);
  const auto RESULT =
      Pullback::SphereAnimation::grouped_ripples(math::X_AXIS, prepared, true);
  HS_EXPECT_NEAR(RESULT.sphere.x, EXPECTED.x, 1e-5f);
  HS_EXPECT_NEAR(RESULT.sphere.y, EXPECTED.y, 1e-5f);
  HS_EXPECT_NEAR(RESULT.sphere.z, EXPECTED.z, 1e-5f);
  std::swap(prepared.waves[0], prepared.waves[3]);
  const math::Vector REVERSED = ripple_test_sequential(math::X_AXIS, prepared);
  HS_EXPECT_GT((REVERSED - EXPECTED).length(), 1e-4f);
}

/** @brief Preset-range batches remain finite, unit length and geometrically faithful. */
inline void test_pullback_grouped_ripple_parameter_grid() {
  using namespace Pullback::SphereAnimation;
  float max_error = 0.0f, max_interior_error = 0.0f;
  unsigned cutoff_exceptions = 0;
  for (float thickness : {0.7f, 1.4f})
    for (float decay : {0.0f, 5.0f})
      for (float period : {30.0f, 143.0f})
        for (bool random_center : {false, true}) {
          BurstRippleParams params;
          params.strength = 0.15f;
          params.thickness = thickness;
          params.decay = decay;
          params.period = period;
          params.burst_rate = 1.0f / 64.0f;
          params.random_center = random_center;
          RippleState state;
          state.seed = 149502542u;
          for (int frame = 0; frame < 256; ++frame) {
            advance_ripple(state, params);
            const auto PREPARED = prepare_ripples(params, state);
            for (int i = 0; i < 24; ++i) {
              const math::Vector VIEW = math::fib_spiral(24, 0.5f, i);
              const auto RESULT = grouped_ripples(VIEW, PREPARED, true);
              const auto UNTRACKED = grouped_ripples(VIEW, PREPARED, false);
              float support_margin = math::PI_F;
              const auto EXPECTED =
                  ripple_test_sequential(VIEW, PREPARED, &support_margin);
              const float ERROR =
                  std::max({fabsf(RESULT.sphere.x - EXPECTED.x),
                            fabsf(RESULT.sphere.y - EXPECTED.y),
                            fabsf(RESULT.sphere.z - EXPECTED.z)});
              max_error = fmaxf(max_error, ERROR);
              if (support_margin > 5e-4f) {
                max_interior_error = fmaxf(max_interior_error, ERROR);
                HS_EXPECT_LE(ERROR, 1e-4f);
              } else
                cutoff_exceptions += ERROR > 1e-4f;
              HS_EXPECT_TRUE(std::isfinite(RESULT.sphere.x));
              HS_EXPECT_TRUE(std::isfinite(RESULT.sphere.y));
              HS_EXPECT_TRUE(std::isfinite(RESULT.sphere.z));
              HS_EXPECT_TRUE(std::isfinite(RESULT.path_length));
              HS_EXPECT_GE(RESULT.path_length, 0.0f);
              HS_EXPECT_NEAR(math::dot(RESULT.sphere, RESULT.sphere), 1.0f,
                             1e-5f);
              HS_EXPECT_LE(ERROR, 0.002f);
              HS_EXPECT_EQ(UNTRACKED.sphere.x, RESULT.sphere.x);
              HS_EXPECT_EQ(UNTRACKED.sphere.y, RESULT.sphere.y);
              HS_EXPECT_EQ(UNTRACKED.sphere.z, RESULT.sphere.z);
              HS_EXPECT_EQ(UNTRACKED.path_length, 0.0f);
            }
          }
        }
  std::printf("  ripple grid max=%g interior=%g cutoff_exceptions=%u\n",
              max_error, max_interior_error, cutoff_exceptions);
}

/** @brief Pole, pole-crossing and large-angle guards preserve sequential geometry. */
inline void test_pullback_grouped_ripple_guards() {
  using namespace Pullback::SphereAnimation;
  PreparedRipples prepared{};
  prepared.count = 3;
  prepared.waves[0] = ripple_test_wave(math::Z_AXIS, 0.5f * math::PI_F, 0.07f);
  prepared.waves[1] = ripple_test_wave(math::Y_AXIS, 0.32310888f, 0.15f);
  prepared.waves[2] = ripple_test_wave(math::Y_AXIS, 0.5f * math::PI_F, 0.2f);
  for (const math::Vector &view :
       {math::Y_AXIS, -math::Y_AXIS, math::X_AXIS, 1.0001f * math::X_AXIS,
        math::Vector(sinf(0.02f), cosf(0.02f), 0.0f)}) {
    const auto EXPECTED = ripple_test_sequential(view, prepared);
    const auto RESULT = grouped_ripples(view, prepared, true);
    HS_EXPECT_NEAR(RESULT.sphere.x, EXPECTED.x, 1e-5f);
    HS_EXPECT_NEAR(RESULT.sphere.y, EXPECTED.y, 1e-5f);
    HS_EXPECT_NEAR(RESULT.sphere.z, EXPECTED.z, 1e-5f);
    HS_EXPECT_TRUE(std::isfinite(RESULT.path_length));
  }
}

/** @brief Tiny-wave travel stays continuous across the radial-basis fallback. */
inline void test_pullback_grouped_ripple_path_continuity() {
  using namespace Pullback::SphereAnimation;
  PreparedRipples prepared{};
  prepared.count = 1;
  const float DISTANCE = asinf(0.01f);
  prepared.waves[0] =
      ripple_test_wave(math::Y_AXIS, DISTANCE + 0.99f * 0.175f, 0.002f);
  const float BELOW = DISTANCE - 1e-6f;
  const float ABOVE = DISTANCE + 1e-6f;
  const auto A = grouped_ripples(math::Vector(sinf(BELOW), cosf(BELOW), 0.0f),
                                 prepared, true);
  const auto B = grouped_ripples(math::Vector(sinf(ABOVE), cosf(ABOVE), 0.0f),
                                 prepared, true);
  HS_EXPECT_GT(A.path_length, 1e-5f);
  HS_EXPECT_GT(B.path_length, 1e-5f);
  HS_EXPECT_NEAR(A.path_length, B.path_length, 2e-6f);
  HS_EXPECT_LT((A.sphere - B.sphere).length(), 1e-5f);
}

/** @brief Baked validity preserves every endpoint-cell admission and transfer bit. */
inline void test_pullback_radial_cell_validity() {
  using namespace Pullback::SphereAnimation;
  static RadialTables tables;
  tables = {};
  HS_EXPECT_FALSE(tables.validity_ready);
  tables.enabled = true;
  tables.count = 1;
  auto &group = tables.groups[0];
  group.previous.fill({.04f, .1f});
  group.current.fill({.08f, .2f});
  const auto CHECK = [&]() {
    prepare_radial_validity(tables);
    HS_EXPECT_TRUE(tables.validity_ready);
    for (size_t cell = 0; cell < RadialGroup::SAMPLES - 1; ++cell) {
      const bool VALID = !(group.previous[cell].path < 0.0f ||
                           group.previous[cell + 1].path < 0.0f ||
                           group.current[cell].path < 0.0f ||
                           group.current[cell + 1].path < 0.0f);
      HS_EXPECT_EQ(bool(tables.valid_cells[0][cell]), VALID);
      const float DISTANCE =
          (float(cell) + .375f) * (math::PI_F / (RadialGroup::SAMPLES - 1));
      for (bool track_path : {false, true}) {
        RadialTransfer expected{.123f, .456f}, actual = expected;
        const bool EXPECTED = RippleKernel::sample_transfer(
            group.previous, group.current, DISTANCE, .625f, track_path,
            expected);
        const bool ACTUAL = RippleKernel::sample_transfer(
            group.previous, group.current, DISTANCE, .625f, track_path, actual,
            tables.valid_cells[0].data());
        HS_EXPECT_EQ(ACTUAL, EXPECTED);
        HS_EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(actual)), 0);
      }
    }
  };
  CHECK();
  for (int endpoint : {0, 1})
    for (size_t node :
         {size_t{0}, size_t{1}, size_t{RadialGroup::SAMPLES - 1} / 2,
          size_t{RadialGroup::SAMPLES - 1}}) {
      auto &curve = endpoint == 0 ? group.previous : group.current;
      const float SAVED = curve[node].path;
      curve[node].path = -1.0f;
      CHECK();
      curve[node].path = SAVED;
    }
  group.previous[0].path = -0.0f;
  group.current[0].path = std::numeric_limits<float>::quiet_NaN();
  CHECK();
  RadialTransfer expected{.123f, .456f}, actual = expected;
  const bool EXPECTED = RippleKernel::sample_transfer(
      group.previous, group.current, math::PI_F, .25f, true, expected);
  const bool ACTUAL = RippleKernel::sample_transfer(
      group.previous, group.current, math::PI_F, .25f, true, actual,
      tables.valid_cells[0].data());
  HS_EXPECT_EQ(ACTUAL, EXPECTED);
  HS_EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(actual)), 0);
}

/** @brief Owned masks match unmasked radial output through supported and rejected updates. */
inline void test_pullback_radial_validity_lifetime() {
  using namespace Pullback::SphereAnimation;
  static RadialTables tables;
  tables = {};
  BurstRippleParams params;
  params.period = 80.0f;
  params.burst_rate = 1.0f / 64.0f;
  RippleState state;
  state.seed = 149502542u;
  uint64_t successes = 0;
  for (int frame = 0; frame < 600; ++frame) {
    if (frame == 173)
      params.strength = .12f;
    if (frame == 279)
      params.period = 143.0f;
    if (frame == 381)
      params.random_center = 0;
    if (frame == 419) {
      auto unsupported = params;
      unsupported.thickness = .1f;
      prepare_radial_tables(unsupported, state, tables, true);
      HS_EXPECT_FALSE(tables.enabled);
      HS_EXPECT_FALSE(tables.validity_ready);
    }
    advance_ripple(state, params);
    prepare_radial_tables(params, state, tables, true);
    HS_EXPECT_TRUE(tables.enabled);
    HS_EXPECT_TRUE(tables.validity_ready);
    for (int probe = 0; probe < 32; ++probe) {
      const float AZIMUTH = math::TWO_PI_F * probe / 32;
      const float Y = -.95f + 1.9f * ((probe * 17 + frame) % 32) / 31;
      const float RADIUS = sqrtf(1.0f - Y * Y);
      const math::Vector VIEW{RADIUS * cosf(AZIMUTH), Y,
                              RADIUS * sinf(AZIMUTH)};
      for (bool track_path : {false, true}) {
        Pullback::SurfaceResult actual, expected;
        const bool ACTUAL =
            try_radial_ripples(VIEW, track_path, &tables, .375f, actual);
        tables.validity_ready = false;
        const bool EXPECTED =
            try_radial_ripples(VIEW, track_path, &tables, .375f, expected);
        tables.validity_ready = true;
        HS_EXPECT_EQ(ACTUAL, EXPECTED);
        HS_EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(actual)), 0);
        successes += ACTUAL;
      }
    }
  }
  HS_EXPECT_GT(successes, uint64_t{10000});
  prepare_radial_tables(params, state, tables);
  HS_EXPECT_FALSE(tables.validity_ready);
  RadialTables manual{};
  HS_EXPECT_FALSE(manual.validity_ready);
  manual.enabled = true;
  manual.count = 1;
  auto &group = manual.groups[0];
  group.center = math::Y_AXIS;
  group.cos_near = 1.0f;
  group.cos_far = -1.0f;
  group.previous.fill({0.0f, .1f});
  group.current.fill({0.0f, .1f});
  Pullback::SurfaceResult result;
  HS_EXPECT_TRUE(try_radial_ripples(math::X_AXIS, true, &manual, .5f, result));
  group.current[(RadialGroup::SAMPLES - 1) / 2].path = -1.0f;
  HS_EXPECT_FALSE(try_radial_ripples(math::X_AXIS, true, &manual, .5f, result));
}

/** @brief Roundoff budget for independently evaluated radial endpoint curves. */
inline float radial_endpoint_roundoff(float strength, float a, float b) {
  constexpr float UNIT = std::numeric_limits<float>::epsilon() * .5f;
  constexpr float ROUNDS =
      8 * Pullback::SphereAnimation::RIPPLE_WAVES_PER_BURST;
  constexpr float GAMMA = ROUNDS * UNIT / (1 - ROUNDS * UNIT);
  return 2 * GAMMA * std::max({strength, fabsf(a), fabsf(b)});
}

/** @brief Cached endpoints preserve reuse bytes and agree with independent bakes. */
inline void test_pullback_radial_ripple_endpoint_reuse() {
  using namespace Pullback::SphereAnimation;
  static RadialTables cached, fresh;
  ScratchScope scratch(scratch_arena_a);
  auto *sentinel = scratch_arena_a.make<uint32_t>(0x12345678u);
  const size_t OFFSET = scratch_arena_a.get_offset();
  const auto CHECK = [&](const BurstRippleParams &params,
                         const RippleState &state) {
    prepare_radial_tables(params, state, cached, true);
    HS_EXPECT_EQ(scratch_arena_a.get_offset(), OFFSET);
    HS_EXPECT_EQ(*sentinel, 0x12345678u);
    fresh.enabled = false;
    prepare_radial_tables(params, state, fresh, true);
    HS_EXPECT_EQ(scratch_arena_a.get_offset(), OFFSET);
    HS_EXPECT_EQ(*sentinel, 0x12345678u);
    HS_EXPECT_EQ(cached.enabled, fresh.enabled);
    HS_EXPECT_EQ(cached.count, fresh.count);
    HS_EXPECT_EQ(cached.validity_ready, fresh.validity_ready);
    if (fresh.validity_ready)
      HS_EXPECT_EQ(std::memcmp(cached.valid_cells.data(),
                               fresh.valid_cells.data(),
                               sizeof(cached.valid_cells)),
                   0);
    for (uint8_t group = 0; group < fresh.count; ++group) {
      const auto &ACTUAL = cached.groups[group];
      const auto &EXPECTED = fresh.groups[group];
      HS_EXPECT_EQ(
          std::memcmp(&ACTUAL.center, &EXPECTED.center, sizeof(ACTUAL.center)),
          0);
      HS_EXPECT_EQ(
          std::memcmp(&ACTUAL.cos_near, &EXPECTED.cos_near, sizeof(float)), 0);
      HS_EXPECT_EQ(
          std::memcmp(&ACTUAL.cos_far, &EXPECTED.cos_far, sizeof(float)), 0);
      HS_EXPECT_TRUE(ACTUAL.current_key.matches(EXPECTED.current_key));
#if HS_TEST_FAST_MATH
      for (size_t sample = 0; sample < ACTUAL.previous.size(); ++sample) {
        const auto &a = ACTUAL.previous[sample];
        const auto &b = EXPECTED.previous[sample];
        HS_EXPECT_EQ(a.path < 0, b.path < 0);
        if (b.path < 0) {
          HS_EXPECT_EQ(std::memcmp(&a, &b, sizeof(a)), 0);
          continue;
        }
        const float ERROR =
            radial_endpoint_roundoff(params.strength, a.path, b.path);
        HS_EXPECT_NEAR(a.angle, b.angle, ERROR);
        HS_EXPECT_NEAR(a.path, b.path, ERROR);
      }
#else
      HS_EXPECT_EQ(std::memcmp(ACTUAL.previous.data(), EXPECTED.previous.data(),
                               sizeof(ACTUAL.previous)),
                   0);
#endif
      HS_EXPECT_EQ(std::memcmp(ACTUAL.current.data(), EXPECTED.current.data(),
                               sizeof(ACTUAL.current)),
                   0);
    }
  };
  for (bool randomize : {false, true}) {
    cached.enabled = false;
    BurstRippleParams params;
    params.period = 143.0f;
    params.burst_rate = 1.0f / 64.0f;
    RippleState state;
    state.seed = 149502542u;
    uint32_t random = 0x719463a5u;
    const auto NEXT = [&]() {
      random ^= random << 13;
      random ^= random >> 17;
      random ^= random << 5;
      return static_cast<float>(random & 0xffffffu) / 16777216.0f;
    };
    for (int frame = 0; frame < 600; ++frame) {
      if (randomize && frame % 29 == 0) {
        params.period = 60.0f + 83.0f * NEXT();
        params.strength = 0.01f + 0.14f * NEXT();
        params.decay = 5.0f * NEXT();
        params.thickness = 0.7f + 0.7f * NEXT();
        params.burst_rate = (1.0f + NEXT()) / 128.0f;
        params.burst_size = static_cast<uint8_t>(4.0f * NEXT());
        params.random_center = static_cast<uint8_t>(NEXT() >= 0.5f);
        params.center_azimuth = math::TWO_PI_F * NEXT();
        params.center_polar = math::PI_F * NEXT();
        state.seed = random;
      }
      if (randomize && frame % 83 == 0)
        for (auto &wave : state.waves)
          wave.age += 0.25f;
      if (randomize && frame % 101 == 0)
        for (int skipped = 0; skipped < 7; ++skipped)
          advance_ripple(state, params);
      advance_ripple(state, params);
      CHECK(params, state);
      if (randomize && frame % 47 == 0) {
        auto unsupported = params;
        unsupported.thickness = 0.1f;
        CHECK(unsupported, state);
        CHECK(params, state);
      }
    }
  }
  cached.enabled = false;
  BurstRippleParams params;
  params.period = 143.0f;
  RippleState reordered;
  reordered.burst_count = 4;
  for (uint32_t burst = 0; burst < 4; ++burst)
    reordered.waves[burst * 4] = {20.0f * (burst + 1), burst};
  CHECK(params, reordered);
  for (uint32_t burst = 0; burst < 4; ++burst)
    reordered.waves[burst * 4].age = 81.0f - 20.0f * burst;
  CHECK(params, reordered);

  for (int control = 0; control < 4; ++control) {
    cached.enabled = false;
    params = {};
    params.period = 143.0f;
    reordered.waves[0] = {40.0f, 0};
    CHECK(params, reordered);
    advance_ripple(reordered, params);
    switch (control) {
    case 0:
      params.period = 142.0f;
      break;
    case 1:
      params.strength = 0.125f;
      break;
    case 2:
      params.decay = 0.2f;
      break;
    case 3:
      params.thickness = 0.8f;
      break;
    }
    CHECK(params, reordered);
  }

  for (bool change_shape : {false, true}) {
    cached.enabled = false;
    params = {};
    params.period = 143.0f;
    reordered = {};
    reordered.burst_count = 1;
    reordered.waves[0] = {40.0f, 0};
    prepare_radial_tables(params, reordered, cached, true);
    HS_EXPECT_EQ(cached.count, uint8_t{1});
    constexpr float SENTINEL = 1234.0f;
    for (size_t sample = 0; sample < cached.groups[0].current.size(); ++sample)
      cached.groups[0].current[sample] = {SENTINEL + float(sample),
                                          SENTINEL + 2 * float(sample)};
    const auto EXPECTED = cached.groups[0].current;
    advance_ripple(reordered, params);
    if (change_shape)
      params.strength = 0.125f;
    prepare_radial_tables(params, reordered, cached, true);
    if (change_shape)
      HS_EXPECT_NE(cached.groups[0].previous[0].angle, SENTINEL);
    else
      HS_EXPECT_EQ(std::memcmp(cached.groups[0].previous.data(),
                               EXPECTED.data(), sizeof(EXPECTED)),
                   0);
  }

  RadialCurveKey key{{80.0f, 0.15f, 0.1f, 0.7f}, {20.0f, 4.0f}, 2};
  HS_EXPECT_TRUE(key.matches(key));
  for (size_t scalar = 0; scalar < key.shape.size(); ++scalar) {
    auto different = key;
    different.shape[scalar] = std::nextafter(key.shape[scalar], 1000.0f);
    HS_EXPECT_FALSE(key.matches(different));
  }
  for (uint8_t wave = 0; wave < key.count; ++wave) {
    auto different = key;
    different.ages[wave] = std::nextafter(key.ages[wave], 1000.0f);
    HS_EXPECT_FALSE(key.matches(different));
  }
  auto different = key;
  --different.count;
  HS_EXPECT_FALSE(key.matches(different));
  different = key;
  different.ages[0] = -0.0f;
  key.ages[0] = 0.0f;
  HS_EXPECT_FALSE(key.matches(different));
}

/** @brief Radial endpoints reproduce an analytic single-wave crest. */
inline void test_pullback_radial_ripple_crest() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  params.period = 80.0f;
  params.decay = 0.0f;
  params.random_center = 0;
  params.center_polar = 0.0f;
  RippleState state;
  state.burst_count = 1;
  state.waves[0] = {40.0f, 0};
  static RadialTables tables;
  prepare_radial_tables(params, state, tables);
  HS_EXPECT_TRUE(tables.enabled);
  HS_EXPECT_EQ(tables.count, 1);
  auto prepared = prepare_ripples(params, state);
  const auto RESULT =
      radial_ripples(math::X_AXIS, prepared, true, &tables, 1.0f);
  // Half-way through its period the wave carries half its strength.
  const float CREST = 0.5f * params.strength;
  HS_EXPECT_NEAR(RESULT.sphere.x, cosf(CREST), 1e-6f);
  HS_EXPECT_NEAR(RESULT.sphere.y, -sinf(CREST), 1e-6f);
  HS_EXPECT_NEAR(RESULT.sphere.z, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(RESULT.path_length, CREST, 1e-6f);
  const auto UNTRACKED =
      radial_ripples(math::X_AXIS, prepared, false, &tables, 1.0f);
  HS_EXPECT_EQ(UNTRACKED.sphere.x, RESULT.sphere.x);
  HS_EXPECT_EQ(UNTRACKED.sphere.y, RESULT.sphere.y);
  HS_EXPECT_EQ(UNTRACKED.path_length, 0.0f);
}

/** @brief Unsupported clocks and narrow wavelets retain grouped evaluation. */
inline void test_pullback_radial_ripple_fallback() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  RippleState state;
  state.burst_count = 1;
  state.waves[0] = {20.0f, 0};
  static RadialTables tables;
  for (bool narrow : {false, true}) {
    params.period = narrow ? 80.0f : 30.0f;
    params.thickness = narrow ? 0.1f : 0.7f;
    prepare_radial_tables(params, state, tables);
    HS_EXPECT_FALSE(tables.enabled);
    auto prepared = prepare_ripples(params, state);
    const math::Vector VIEW = math::Vector(0.4f, 0.5f, 0.6f).normalized();
    const auto EXPECTED = grouped_ripples(VIEW, prepared, true);
    const auto RESULT = radial_ripples(VIEW, prepared, true, &tables, 1.0f);
    HS_EXPECT_EQ(RESULT.sphere.x, EXPECTED.sphere.x);
    HS_EXPECT_EQ(RESULT.sphere.y, EXPECTED.sphere.y);
    HS_EXPECT_EQ(RESULT.sphere.z, EXPECTED.sphere.z);
    HS_EXPECT_EQ(RESULT.path_length, EXPECTED.path_length);
  }
  params = {};
  prepare_radial_tables(params, RippleState{}, tables);
  HS_EXPECT_TRUE(tables.enabled);
  HS_EXPECT_EQ(tables.count, 0);
  PreparedRipples empty;
  const auto IDENTITY =
      radial_ripples(math::X_AXIS, empty, true, &tables, 1.0f);
  HS_EXPECT_EQ(IDENTITY.sphere.x, 1.0f);
  HS_EXPECT_EQ(IDENTITY.path_length, 0.0f);
}

/** @brief Radial rejection resolves the whole original input exactly once. */
inline void test_pullback_radial_ripple_lazy_fallback() {
  using namespace Pullback::SphereAnimation;
  static RadialTables tables;
  tables.enabled = true;
  tables.count = 2;
  for (auto &group : tables.groups) {
    group.cos_near = 1.0f;
    group.cos_far = -1.0f;
    group.previous.fill({0.2f, 0.2f});
    group.current.fill({0.2f, 0.2f});
  }
  tables.groups[0].center = math::Y_AXIS;
  tables.groups[1].center = math::Z_AXIS;
  const Pullback::SurfaceResult EXPECTED{math::Vector(0.2f, 0.3f, 0.4f), 0.75f};
  int calls = 0;
  math::Vector input = math::X_AXIS;
  const auto FALLBACK = [&](const math::Vector &original) {
    ++calls;
    HS_EXPECT_EQ(original, input);
    return EXPECTED;
  };
  const auto CHECK_FALLBACK = [&](const RadialTables *view) {
    calls = 0;
    const auto RESULT = radial_ripples(input, true, view, 0.5f, FALLBACK);
    HS_EXPECT_EQ(calls, 1);
    HS_EXPECT_EQ(RESULT.sphere, EXPECTED.sphere);
    HS_EXPECT_EQ(RESULT.path_length, EXPECTED.path_length);
  };
  CHECK_FALLBACK(nullptr);
  tables.enabled = false;
  CHECK_FALLBACK(&tables);
  tables.enabled = true;
  tables.groups[1].current.fill({0.2f, -1.0f});
  CHECK_FALLBACK(&tables);
  tables.count = 1;
  input = math::Y_AXIS;
  CHECK_FALLBACK(&tables);
  input = 2.0f * math::X_AXIS;
  CHECK_FALLBACK(&tables);
  input = math::X_AXIS;
  calls = 0;
  const auto RESULT = radial_ripples(input, true, &tables, 0.5f, FALLBACK);
  HS_EXPECT_EQ(calls, 0);
  HS_EXPECT_GT((RESULT.sphere - input).length(), 0.1f);
}

template <typename Value> struct RippleTestResource {
  Value value;
  template <auto Key> const Value &get() const { return value; }
};
struct RippleTestClock : Pullback::SphereAnimation::RippleState {
  const Pullback::SphereAnimation::RadialTables *radial_tables = nullptr;
  float radial_fraction = -1.0f;
};
struct RippleTestFrame {
  RippleTestResource<Pullback::SphereAnimation::BurstRippleParams> params;
  RippleTestResource<RippleTestClock> resources;
};
struct RippleTestBinding {
  using FrameState = RippleTestFrame;
};
using RippleTestProvider = Pullback::SurfaceProvider<
    RippleTestBinding, Pullback::SphereAnimation::BurstRippleParams, true>;
struct CustomRippleTestProvider : RippleTestProvider {
  static inline int calls = 0;
  static auto prepare(const FrameState &frame) {
    ++calls;
    auto prepared = RippleTestProvider::prepare(frame);
    prepared.count = 0;
    return prepared;
  }
};
struct GenericRippleTestProvider {
  using Binding = RippleTestBinding;
  using FrameState = RippleTestFrame;
  static auto prepare(const FrameState &frame) {
    return RippleTestProvider::prepare(frame);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

/** @brief Only owned usable radial views defer exact column preparation. */
inline void test_pullback_ripple_lazy_preparation() {
  using namespace Pullback::SphereAnimation;
  using Policy = Pullback::Surface::BurstRipple<RippleTestProvider>;
  using Custom = Pullback::Surface::BurstRipple<CustomRippleTestProvider>;
  using Generic = Pullback::Surface::BurstRipple<GenericRippleTestProvider>;
  static_assert(Policy::LAZY_RADIAL_PREPARE);
  static_assert(!Custom::LAZY_RADIAL_PREPARE);
  static_assert(!Generic::LAZY_RADIAL_PREPARE);
  static_assert(std::is_same_v<Custom::Prepared, PreparedRipples>);
  static_assert(std::is_same_v<Generic::Prepared, PreparedRipples>);
  RippleTestFrame frame;
  auto &params = frame.params.value;
  auto &clock = frame.resources.value;
  params.period = 80.0f;
  params.random_center = 0;
  params.center_polar = 0.0f;
  clock.burst_count = 1;
  clock.waves[0] = {40.0f, 0};
  static RadialTables tables;
  prepare_radial_tables(params, clock, tables);
  clock.radial_tables = &tables;
  clock.radial_fraction = 0.5f;
  const auto PREPARED = Policy::prepare(frame);
  HS_EXPECT_FALSE(PREPARED.has_value());
  const auto EXACT = RippleTestProvider::prepare(frame);
  const auto EXPECTED =
      radial_ripples(math::X_AXIS, EXACT, true, &tables, 0.5f);
  const auto RESULT = Policy::apply(math::X_AXIS, frame, PREPARED);
  HS_EXPECT_EQ(RESULT.sphere, EXPECTED.sphere);
  HS_EXPECT_EQ(RESULT.path_length, EXPECTED.path_length);
  const auto CUSTOM_INTACT = Custom::prepare(frame);
  HS_EXPECT_EQ(Custom::apply(math::X_AXIS, frame, CUSTOM_INTACT).sphere,
               math::X_AXIS);
  tables.groups[0].current.fill({0.0f, -1.0f});
  const auto FALLBACK = Policy::apply(math::X_AXIS, frame, PREPARED);
  const auto GROUPED = grouped_ripples(math::X_AXIS, EXACT, true);
  HS_EXPECT_EQ(FALLBACK.sphere, GROUPED.sphere);
  HS_EXPECT_EQ(FALLBACK.path_length, GROUPED.path_length);
  CustomRippleTestProvider::calls = 0;
  const auto CUSTOM = Custom::prepare(frame);
  HS_EXPECT_EQ(CustomRippleTestProvider::calls, 1);
  HS_EXPECT_EQ(CUSTOM.count, 0);
  HS_EXPECT_EQ(Custom::apply(math::X_AXIS, frame, CUSTOM).sphere, math::X_AXIS);
  HS_EXPECT_EQ(CustomRippleTestProvider::calls, 1);
  HS_EXPECT_EQ(Generic::prepare(frame).count, EXACT.count);
  HS_EXPECT_EQ(Generic::apply(math::X_AXIS, frame, EXACT).sphere,
               GROUPED.sphere);
  {
    clock.radial_fraction = -1.0f;
    const auto ENDPOINT = Policy::prepare(frame);
    HS_EXPECT_TRUE(ENDPOINT.has_value());
    HS_EXPECT_EQ(ENDPOINT->count, EXACT.count);
    HS_EXPECT_EQ(Policy::apply(math::X_AXIS, frame, ENDPOINT).sphere,
                 GROUPED.sphere);
  }
  clock.radial_fraction = 0.5f;
  tables.enabled = false;
  HS_EXPECT_TRUE(Policy::prepare(frame).has_value());
  tables.enabled = true;
  clock.radial_tables = nullptr;
  HS_EXPECT_TRUE(Policy::prepare(frame).has_value());
}
/** @brief Cached ripple geometry stays subpixel across column times. */
inline void test_pullback_radial_ripple_geometry() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  params.period = 91.48529f;
  params.decay = 0.2f;
  params.burst_rate = 1.0f / 128.0f;
  RippleState state;
  state.seed = 149502542u;
  static RadialTables tables;
  for (int frame = 0; frame < 300; ++frame) {
    advance_ripple(state, params);
    prepare_radial_tables(params, state, tables);
    for (float fraction : {0.0f, 0.5f, 1.0f}) {
      auto sampled = state;
      for (auto &wave : sampled.waves)
        wave.age -= 1.0f - fraction;
      auto prepared = prepare_ripples(params, sampled);
      for (int i = 0; i < 32; ++i) {
        const math::Vector VIEW = math::fib_spiral(32, 0.5f, i);
        const auto EXPECTED = grouped_ripples(VIEW, prepared, true);
        const auto RESULT =
            radial_ripples(VIEW, prepared, true, &tables, fraction);
        HS_EXPECT_TRUE(std::isfinite(RESULT.path_length));
        HS_EXPECT_GE(RESULT.path_length, 0.0f);
        HS_EXPECT_NEAR(math::dot(RESULT.sphere, RESULT.sphere), 1.0f, 1e-5f);
        HS_EXPECT_LE((RESULT.sphere - EXPECTED.sphere).length(), 0.005f);
        HS_EXPECT_LE(fabsf(RESULT.path_length - EXPECTED.path_length), 0.02f);
      }
    }
  }
}
