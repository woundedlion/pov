/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Lattice trace: cubic compositor and octet preparation parity
// ============================================================================

/** @brief Palette that ramps hue with depth so layer order shows in colour. */
struct DepthPalette {
  /**
   * @brief Samples the ramp.
   * @param t Palette position in [0, 1].
   * @return Opaque colour at @p t.
   */
  Color4 get(float t) const {
    return {Pixel(uint16_t(1000 + 50000 * t), uint16_t(60000 - 45000 * t),
                  uint16_t(12000 + 30000 * t)),
            1.0f};
  }
};

/** @brief Arena-backed DepthPalette bake. */
struct DepthPaletteFixture {
  alignas(std::max_align_t) uint8_t storage[4096];
  Arena arena{storage, sizeof(storage)};
  BakedPaletteStorage palette;

  DepthPaletteFixture() { palette.bake(arena, DepthPalette{}); }
};

/** @brief Lattice trace settings bound to a baked DepthPalette. */
struct TraceFixture : DepthPaletteFixture {
  SDF::LatticeTrace::CrossingStorage crossings;
  SDF::LatticeTrace::Settings settings;

  TraceFixture() {
    settings.palette = &palette.view();
    settings.crossings = &crossings;
  }
};

/** @brief Axis and spiral view sample count. */
inline constexpr int DIRECTION_SAMPLES = 51;

/**
 * @brief Test view direction: the three axes, then a spiral over the sphere.
 * @param sample Sample index.
 * @return Unit view direction.
 */
inline math::Vector direction(int sample) {
  if (sample == 0)
    return math::X_AXIS;
  if (sample == 1)
    return math::Y_AXIS;
  if (sample == 2)
    return math::Z_AXIS;
  const float z = (sample - 2.5f) / ((DIRECTION_SAMPLES - 3) / 2.0f) - 1.0f;
  const float radius = sqrtf(1.0f - z * z);
  const float phase = sample * 2.39996323f;
  return {radius * cosf(phase), radius * sinf(phase), z};
}

/** @brief Every lattice view sample is finite and unit length. */
inline void test_lattice_view_directions_are_unit_length() {
  for (int i = 0; i < DIRECTION_SAMPLES; ++i) {
    const math::Vector VIEW = direction(i);
    HS_EXPECT_TRUE(std::isfinite(VIEW.x));
    HS_EXPECT_TRUE(std::isfinite(VIEW.y));
    HS_EXPECT_TRUE(std::isfinite(VIEW.z));
    HS_EXPECT_NEAR(VIEW.x * VIEW.x + VIEW.y * VIEW.y + VIEW.z * VIEW.z, 1.0f,
                   2e-6f);
  }
}

/**
 * @brief Checks one composite_crossings specialization against event shading.
 * @tparam SLICE 4D-slice specialization.
 * @tparam SHELLS Specialized shell count, 0 for the generic loop.
 * @param prepared Prepared lattice shading.
 * @param view View direction.
 * @param lit Incremented when the composite has coverage.
 */
template <bool SLICE, uint8_t SHELLS>
void check_cubic(const SDF::Lattice::PreparedShading &prepared,
                 const math::Vector &view, int &lit) {
  const auto composite =
      SDF::Lattice::composite_crossings<SLICE, SHELLS>(view, prepared);
  const auto actual = composite.finish();
  LayerComposite reference;
  SDF::Lattice::Events<SLICE, SHELLS> reference_events(view, prepared.lattice);
  Raycast::trace_events(reference_events, {0, prepared.lattice.far_distance},
                        {}, [&](const Raycast::Contribution &hit) {
                          prepared.appearance.composite(reference, hit.t,
                                                        hit.coverage);
                          return !reference.saturated();
                        });
  const Color4 REFERENCE = reference.finish();
  HS_EXPECT_NEAR(actual.alpha, REFERENCE.alpha, 3e-4f);
  // Independent reciprocal evaluation can cross one final pixel rounding boundary.
  const Pixel ACTUAL_PIXEL = composite.premultiplied();
  const Pixel REFERENCE_PIXEL = reference.premultiplied();
  HS_EXPECT_LE(max_channel_gap(ACTUAL_PIXEL, REFERENCE_PIXEL), 1);
  lit += actual.alpha > 0.0f;
}

/**
 * @brief Verifies the cubic crossing compositor matches Events tracing.
 * @details The generic path runs for every case; each specialization runs on
 * the domain/shell pairs where `HyperLattice` dispatches to the matching
 * `HyperLatticeDetail::Renderer` instantiation.
 */
inline void test_cubic_compositor_matches_event_shading() {
  DepthPaletteFixture fixture;
  SDF::Lattice::CrossingList crossings;
  int lit = 0;
  for (const auto domain :
       {SDF::Lattice::Domain::SPATIAL_3D, SDF::Lattice::Domain::SLICE_4D}) {
    for (const float cell : {.4f, 1.7f}) {
      SDF::Lattice::Settings settings;
      settings.mode = domain;
      settings.cell_size = cell;
      settings.wire_radius = .12f;
      auto embedding = math::Mat4::identity();
      math::rotate_plane(embedding, 0, 1, .31f);
      if (domain == SDF::Lattice::Domain::SLICE_4D)
        math::rotate_plane(embedding, 0, 3, .73f);
      for (const auto shells :
           {SDF::Lattice::ShellCount::ONE, SDF::Lattice::ShellCount::TWO,
            SDF::Lattice::ShellCount::THREE}) {
        settings.shells = shells;
        const SDF::Lattice::PreparedShading prepared{
            SDF::Lattice::prepare(settings, {{.231f, -.437f, .719f, .383f}},
                                  embedding, 7.0f, .018f),
            {1.0f / 7.0f, 0.0f, 2.0f, &fixture.palette.view(), .8f},
            &crossings};
        for (int i = 0; i < DIRECTION_SAMPLES; ++i) {
          const auto view = direction(i);
          check_cubic<false, 0>(prepared, view, lit);
          if (domain == SDF::Lattice::Domain::SLICE_4D &&
              shells == SDF::Lattice::ShellCount::TWO)
            check_cubic<true, 2>(prepared, view, lit);
          else if (domain == SDF::Lattice::Domain::SLICE_4D &&
                   shells == SDF::Lattice::ShellCount::THREE)
            check_cubic<true, 0>(prepared, view, lit);
          else if (shells == SDF::Lattice::ShellCount::TWO)
            check_cubic<false, 2>(prepared, view, lit);
        }
      }
    }
  }
  HS_EXPECT_GT(lit, 100);
}

/** @brief Verifies LatticeTrace octet shading matches world-space octet events. */
inline void test_octet_preparation_matches_world_events() {
  namespace Trace = SDF::LatticeTrace;
  TraceFixture fixture;
  Trace::Settings &settings = fixture.settings;
  settings.pixel_half_angle = .018f;
  int lit = 0;
  for (const auto domain : {Raycast::SamplingDomain::SPATIAL_3D,
                            Raycast::SamplingDomain::SLICE_4D}) {
    settings.domain = domain;
    const bool SLICE = domain == Raycast::SamplingDomain::SLICE_4D;
    settings.center = {{.231f, -.437f, .719f, SLICE ? .383f : 0.0f}};
    settings.embedding = math::Mat4::identity();
    math::rotate_plane(settings.embedding, 0, 1, .31f);
    if (SLICE)
      math::rotate_plane(settings.embedding, 0, 3, .73f);
    for (float cell : {.4f, 1.7f}) {
      settings.cell_size = cell;
      settings.wire_radius = .055f * cell;
      for (float radial : {0.0f, .7f}) {
        settings.radial_start = radial;
        const auto prepared = Trace::prepare(settings);
        HS_EXPECT_TRUE(prepared.valid);
        for (int i = 0; i < DIRECTION_SAMPLES; ++i) {
          const auto view = direction(i);
          const auto ray = prepared.camera.ray(view);
          const auto ambient =
              prepared.camera.embedding.apply({{view.x, view.y, view.z, 0.0f}});
          Raycast::ShadedTrace expected;
          if (SLICE) {
            SDF::OctetEvents4 events(prepared.octet4,
                                     prepared.camera.point4(ray.origin),
                                     ambient, ray.interval, prepared.footprint);
            expected = Raycast::shade_events(
                events, ray.interval, prepared.limits, prepared.appearance);
          } else {
            const Raycast::Ray world{prepared.camera.point3(ray.origin),
                                     {ambient[0], ambient[1], ambient[2]},
                                     ray.interval};
            SDF::OctetEvents events(prepared.octet, world, prepared.footprint);
            expected = Raycast::shade_events(
                events, ray.interval, prepared.limits, prepared.appearance);
          }
          const auto actual = SLICE ? Trace::shade<true>(view, prepared)
                                    : Trace::shade<false>(view, prepared);
          const Pixel premultiplied =
              expected.color.color * expected.color.alpha;
          HS_EXPECT_EQ(actual.status, expected.trace.status);
          HS_EXPECT_LE(max_channel_gap(actual.color, premultiplied), 2);
          lit += actual.color != Pixel{};
        }
      }
    }
  }
  HS_EXPECT_GT(lit, 10);
}

/** @brief Distant octet centres match a metric-aware periodic oracle. */
inline void test_octet_distant_center_periodicity() {
  namespace Trace = SDF::LatticeTrace;
  TraceFixture fixture;
  Trace::Settings &settings = fixture.settings;
  settings.domain = Raycast::SamplingDomain::SLICE_4D;
  settings.cell_size = std::nextafter(1.0f / SDF::OctetFramework4::HALF_CUBE,
                                      std::numeric_limits<float>::infinity());
  settings.center = {{32768.25f, -65536.5f, 98304.75f, -131072.0f}};
  const auto distant = Trace::prepare(settings);
  HS_EXPECT_TRUE(distant.valid);
  float expected[4];
  for (int i = 0; i < 4; ++i) {
    const float NORMALIZED =
        settings.center[i] * distant.octet4_projection.inverse_scale;
    expected[i] = static_cast<float>(
        static_cast<double>(NORMALIZED) -
        2.0 * std::trunc(static_cast<double>(NORMALIZED) / 2.0));
    settings.center[i] = expected[i] / distant.octet4_projection.inverse_scale;
  }
  const auto local = Trace::prepare(settings);
  HS_EXPECT_TRUE(local.valid);
  for (int i = 0; i < 4; ++i) {
    HS_CONTEXT("axis", i);
    HS_EXPECT_NEAR(distant.octet4_projection.origin[i], expected[i], 1e-6f);
    HS_EXPECT_NEAR(distant.octet4_projection.origin[i],
                   local.octet4_projection.origin[i], 1e-6f);
    HS_EXPECT_TRUE(std::abs(distant.octet4_projection.origin[i]) < 2.0f);
  }
}

/** @brief Non-positive and NaN wire radii invalidate octets but not shells. */
inline void test_lattice_prepare_validates_wire_radius_per_geometry() {
  namespace Trace = SDF::LatticeTrace;
  TraceFixture fixture;
  SDF::ShellLayerStorage layers;
  Trace::Settings &settings = fixture.settings;
  settings.shell_layers = &layers;
  for (const auto domain : {Raycast::SamplingDomain::SPATIAL_3D,
                            Raycast::SamplingDomain::SLICE_4D}) {
    settings.domain = domain;
    settings.geometry = Trace::Geometry::OCTET;
    settings.wire_radius = .05f;
    HS_EXPECT_TRUE(Trace::prepare(settings).valid);
    for (float radius :
         {0.0f, -1.0f, std::numeric_limits<float>::quiet_NaN()}) {
      settings.wire_radius = radius;
      settings.geometry = Trace::Geometry::OCTET;
      HS_EXPECT_TRUE(!Trace::prepare(settings).valid);
      settings.geometry = Trace::Geometry::SHELLS;
      HS_EXPECT_TRUE(Trace::prepare(settings).valid);
    }
  }
}

/** @brief Prepares octet trace settings over shared crossing scratch. */
inline SDF::LatticeTrace::Prepared
prepare_octet_trace(SDF::LatticeTrace::Settings settings) {
  static SDF::OctetTrace::CrossingStorage crossings;
  settings.crossings = &crossings;
  return SDF::LatticeTrace::prepare(settings);
}

inline void test_octet_prepared_projection() {
  using Domain = Raycast::SamplingDomain;
  reset_globals();
  namespace Trace = SDF::LatticeTrace;
  DepthPaletteFixture palette;
  Trace::Settings settings;
  settings.palette = &palette.palette.view();
  settings.pixel_half_angle = .5f * math::coarse_pixel_pitch<288, 144>();
  const auto expect_same = [](const Raycast::ShadedTrace &actual,
                              const Raycast::ShadedTrace &expected) {
    HS_EXPECT_EQ(actual.trace.status, expected.trace.status);
    HS_EXPECT_NEAR(actual.color.alpha, expected.color.alpha, 3e-4f);
    HS_EXPECT_NEAR((actual.color.color * actual.color.alpha).r,
                   (expected.color.color * expected.color.alpha).r, 2);
    HS_EXPECT_NEAR((actual.color.color * actual.color.alpha).g,
                   (expected.color.color * expected.color.alpha).g, 2);
    HS_EXPECT_NEAR((actual.color.color * actual.color.alpha).b,
                   (expected.color.color * expected.color.alpha).b, 2);
  };
  const auto expect_premultiplied = [](const SDF::OctetTrace::Sample &actual,
                                       const Raycast::ShadedTrace &expected) {
    HS_EXPECT_EQ(actual.status, expected.trace.status);
    const Pixel PREMULTIPLIED = expected.color.color * expected.color.alpha;
    HS_EXPECT_NEAR(actual.color.r, PREMULTIPLIED.r, 1);
    HS_EXPECT_NEAR(actual.color.g, PREMULTIPLIED.g, 1);
    HS_EXPECT_NEAR(actual.color.b, PREMULTIPLIED.b, 1);
  };
  size_t lit = 0, compared = 0;
  for (Domain domain : {Domain::SPATIAL_3D, Domain::SLICE_4D}) {
    settings.domain = domain;
    settings.center = {
        {.231f, -.437f, .719f, domain == Domain::SLICE_4D ? .383f : 0.0f}};
    for (float angle : {0.0f, .31f, 1.23f}) {
      settings.embedding = math::Mat4::identity();
      math::rotate_plane(settings.embedding, 0, 1, angle);
      math::rotate_plane(settings.embedding, 1, 2, angle * .73f);
      if (domain == Domain::SLICE_4D) {
        math::rotate_plane(settings.embedding, 0, 3, angle * .41f);
        math::rotate_plane(settings.embedding, 2, 3, angle * 1.9f);
      }
      for (float cell : {.4f, 1.7f}) {
        settings.cell_size = cell;
        settings.wire_radius = .055f * cell;
        for (float radial : {0.0f, .7f, 2.3f}) {
          settings.radial_start = radial;
          auto prepared = prepare_octet_trace(settings);
          HS_EXPECT_TRUE(prepared.valid);
          const auto &camera = prepared.camera;
          for (float near : {0.0f, .19f}) {
            prepared.camera.interval.near = near;
            for (int sample = 0; sample < 48; ++sample) {
              const float PHASE = sample * 2.39996323f;
              const float Z = (sample + .5f) / 24.0f - 1.0f;
              const float R = sqrtf(1.0f - Z * Z);
              const math::Vector DIRECTION(R * cosf(PHASE), R * sinf(PHASE), Z);
              const auto RAY = camera.ray(DIRECTION);
              const auto AMBIENT = camera.embedding.apply(
                  {{DIRECTION.x, DIRECTION.y, DIRECTION.z, 0.0f}});
              Raycast::ShadedTrace world;
              Raycast::ShadedTrace projected;
              SDF::OctetTrace::Sample traced;
              if (domain == Domain::SLICE_4D) {
                SDF::OctetEvents4 generic(prepared.octet4,
                                          camera.point4(RAY.origin), AMBIENT,
                                          RAY.interval, prepared.footprint);
                world =
                    Raycast::shade_events(generic, RAY.interval,
                                          prepared.limits, prepared.appearance);
                const SDF::OctetEvents4 events(
                    prepared.octet4, prepared.octet4_projection, DIRECTION,
                    camera.radial_start, near, prepared.footprint);
                auto copy = events;
                projected = Raycast::shade_events(
                    copy, RAY.interval, prepared.limits, prepared.appearance);
                traced = SDF::OctetTrace::trace_events(
                    events, RAY.interval, prepared.limits, prepared.appearance);
              } else {
                const Raycast::Ray WORLD{camera.point3(RAY.origin),
                                         {AMBIENT[0], AMBIENT[1], AMBIENT[2]},
                                         RAY.interval};
                SDF::OctetEvents generic(prepared.octet, WORLD,
                                         prepared.footprint);
                world =
                    Raycast::shade_events(generic, RAY.interval,
                                          prepared.limits, prepared.appearance);
                const SDF::OctetEvents events(prepared.octet_projection,
                                              DIRECTION, camera.radial_start,
                                              near, prepared.footprint);
                auto copy = events;
                projected = Raycast::shade_events(
                    copy, RAY.interval, prepared.limits, prepared.appearance);
                traced = SDF::OctetTrace::trace_events(
                    events, RAY.interval, prepared.limits, prepared.appearance);
              }
              expect_same(projected, world);
              expect_premultiplied(traced, projected);
              const auto SHADED =
                  domain == Domain::SLICE_4D
                      ? Trace::shade<true>(DIRECTION, prepared)
                      : Trace::shade<false>(DIRECTION, prepared);
              HS_EXPECT_EQ(SHADED.status, traced.status);
              // Reassociation and the canonical 4D sum can cross one pixel quantum.
#if defined(HS_TEST_FAST_MATH)
              constexpr bool ALLOW_QUANTIZATION = true;
#else
              const bool ALLOW_QUANTIZATION = domain == Domain::SLICE_4D;
#endif
              if (ALLOW_QUANTIZATION) {
                HS_EXPECT_NEAR(SHADED.color.r, traced.color.r, 1);
                HS_EXPECT_NEAR(SHADED.color.g, traced.color.g, 1);
                HS_EXPECT_NEAR(SHADED.color.b, traced.color.b, 1);
              } else {
                HS_EXPECT_EQ(SHADED.color, traced.color);
              }
              ++compared;
              lit += traced.color != Pixel{};
            }
          }
        }
      }
    }
  }
  HS_EXPECT_GT(lit, compared / 10);
}

/**
 * @brief The canonical-frame 4D octet trace against OctetEvents4 over random
 *        flights, including its plane-share crossing skip.
 */
inline void test_octet_4d_canonical_trace() {
  reset_globals();
  namespace Trace = SDF::LatticeTrace;
  DepthPaletteFixture palette;
  Trace::Settings settings;
  settings.palette = &palette.palette.view();
  settings.pixel_half_angle = .5f * math::coarse_pixel_pitch<288, 144>();
  settings.domain = Raycast::SamplingDomain::SLICE_4D;
  hs::Pcg32 rng(0x0c7e74du);
  const auto uniform = [&rng](float low, float high) {
    return low + (high - low) * hs_test::rand_uniform(rng, 0.0f, 1.0f);
  };
  size_t compared = 0, lit = 0, differing = 0;
  for (int frame = 0; frame < 40; ++frame) {
    settings.cell_size = uniform(.6f, 3.0f);
    settings.wire_radius = uniform(.015f, .09f) * settings.cell_size;
    settings.far_distance = uniform(3.0f, 9.0f);
    settings.aa_strength = uniform(0.0f, 2.0f);
    settings.radial_start = frame % 3 ? 0.0f : uniform(0.0f, 1.5f);
    settings.center = {
        {uniform(-3, 3), uniform(-3, 3), uniform(-3, 3), uniform(-3, 3)}};
    settings.embedding = math::Mat4::identity();
    for (int a = 0; a < 4; ++a)
      for (int b = a + 1; b < 4; ++b)
        math::rotate_plane(settings.embedding, a, b, uniform(0.0f, 6.3f));
    const auto prepared = prepare_octet_trace(settings);
    HS_EXPECT_TRUE(prepared.valid);
    const auto &camera = prepared.camera;
    for (int ray = 0; ray < 400; ++ray) {
      const math::Vector DIRECTION =
          math::Vector{uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)}
              .normalized();
      const SDF::OctetEvents4 events(
          prepared.octet4, prepared.octet4_projection, DIRECTION,
          camera.radial_start, camera.interval.near, prepared.footprint);
      const auto EXPECTED = SDF::OctetTrace::trace_events(
          events, camera.interval, prepared.limits, prepared.appearance);
      const auto ACTUAL = SDF::OctetTrace::trace_4d(
          DIRECTION, camera, prepared.octet4, prepared.octet4_projection,
          prepared.footprint, prepared.limits, prepared.appearance,
          *prepared.crossings);
      HS_EXPECT_EQ(ACTUAL.status, EXPECTED.status);
      const int DELTA = max_channel_gap(ACTUAL.color, EXPECTED.color);
      differing += DELTA > 1;
      HS_EXPECT_LE(DELTA, 16);
      ++compared;
      lit += EXPECTED.color != Pixel{};
    }
  }
  HS_EXPECT_GT(lit, compared / 10);
  // Counted crossings and the fixed-point bound move a far-boundary
  // crossing or a rounding tie in a few rays by a few codes.
  HS_EXPECT_LE(differing, compared / 4000);
}

inline void test_octet_4d_canonical_boundary_directions() {
  reset_globals();
  namespace Trace = SDF::LatticeTrace;
  DepthPaletteFixture palette;
  Trace::Settings settings;
  settings.palette = &palette.palette.view();
  settings.pixel_half_angle = .5f * math::coarse_pixel_pitch<288, 144>();
  settings.domain = Raycast::SamplingDomain::SLICE_4D;
  settings.center = {{.255f, .465f, .645f, .375f}};
  const std::array<math::Vector, 10> DIRECTIONS{
      math::X_AXIS,
      -math::X_AXIS,
      math::Y_AXIS,
      -math::Y_AXIS,
      math::Z_AXIS,
      -math::Z_AXIS,
      math::Vector{1, 1e-5f, 0}.normalized(),
      math::Vector{1, 1, 0}.normalized(),
      math::Vector{-1, 1, 0}.normalized(),
      math::Vector{1, 1, 1}.normalized()};
  size_t lit = 0;
  for (int frame = 0; frame <= 8; ++frame) {
    settings.embedding = math::Mat4::identity();
    math::rotate_plane(settings.embedding, 0, 3, frame * math::PI_F / 16);
    settings.radial_start = frame % 2 ? .5f : 0.0f;
    const auto prepared = prepare_octet_trace(settings);
    HS_EXPECT_TRUE(prepared.valid);
    for (const auto &direction : DIRECTIONS) {
      const auto &camera = prepared.camera;
      const SDF::OctetEvents4 events(
          prepared.octet4, prepared.octet4_projection, direction,
          camera.radial_start, camera.interval.near, prepared.footprint);
      const auto EXPECTED = SDF::OctetTrace::trace_events(
          events, camera.interval, prepared.limits, prepared.appearance);
      const auto ACTUAL = SDF::OctetTrace::trace_4d(
          direction, camera, prepared.octet4, prepared.octet4_projection,
          prepared.footprint, prepared.limits, prepared.appearance,
          *prepared.crossings);
      HS_EXPECT_EQ(ACTUAL.status, EXPECTED.status);
      HS_EXPECT_LE(max_channel_gap(ACTUAL.color, EXPECTED.color), 2);
      lit += EXPECTED.color != Pixel{};
    }
  }
  HS_EXPECT_GT(lit, DIRECTIONS.size() * 9 / 3);
}

/**
 * @brief OctetTrace sorted reference reproduces trace_events over random 3D and
 *        4D flights, under default, candidate-starved and layer-starved
 *        budgets.
 */
inline void test_octet_trace_sorted_matches_events() {
  using Domain = Raycast::SamplingDomain;
  reset_globals();
  namespace Trace = SDF::LatticeTrace;
  DepthPaletteFixture palette;
  Trace::Settings settings;
  settings.palette = &palette.palette.view();
  settings.pixel_half_angle = .5f * math::coarse_pixel_pitch<288, 144>();
  hs::Pcg32 rng(0x50f7ed5u);
  const auto uniform = [&rng](float low, float high) {
    return low + (high - low) * hs_test::rand_uniform(rng, 0.0f, 1.0f);
  };
  SDF::OctetTrace::CrossingStorage storage;
  size_t compared = 0, lit = 0, candidate_exhausted = 0, layer_exhausted = 0;
  const auto compare = [&](const auto &events, Raycast::Interval interval,
                           const Raycast::TraceLimits &limits,
                           const Raycast::Appearance &appearance, int budget) {
    const auto EXPECTED =
        SDF::OctetTrace::trace_events(events, interval, limits, appearance);
    const auto ACTUAL = SDF::OctetTrace::detail::reference::trace_sorted(
        events, interval, limits, appearance, storage);
    HS_EXPECT_EQ(ACTUAL.status, EXPECTED.status);
    HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
    ++compared;
    lit += EXPECTED.color != Pixel{};
    if (EXPECTED.status == Raycast::TraceStatus::BUDGET_EXHAUSTED) {
      if (budget == 1)
        ++candidate_exhausted;
      else if (budget == 2)
        ++layer_exhausted;
    }
  };
  for (Domain domain : {Domain::SPATIAL_3D, Domain::SLICE_4D}) {
    HS_CONTEXT(domain == Domain::SLICE_4D ? "4d" : "3d");
    settings.domain = domain;
    for (int frame = 0; frame < 12; ++frame) {
      settings.cell_size = uniform(.6f, 3.0f);
      settings.wire_radius = uniform(.015f, .09f) * settings.cell_size;
      settings.far_distance = uniform(3.0f, 9.0f);
      settings.radial_start = frame % 3 ? 0.0f : uniform(0.0f, 1.5f);
      settings.center = {{uniform(-3, 3), uniform(-3, 3), uniform(-3, 3),
                          domain == Domain::SLICE_4D ? uniform(-3, 3) : 0.0f}};
      settings.embedding = math::Mat4::identity();
      for (int a = 0; a < 4; ++a)
        for (int b = a + 1; b < 4; ++b)
          if (domain == Domain::SLICE_4D || b < 3)
            math::rotate_plane(settings.embedding, a, b, uniform(0.0f, 6.3f));
      const auto prepared = prepare_octet_trace(settings);
      HS_EXPECT_TRUE(prepared.valid);
      const auto &camera = prepared.camera;
      Raycast::TraceLimits starved_candidates = prepared.limits;
      starved_candidates.max_candidates = 6;
      Raycast::TraceLimits starved_layers = prepared.limits;
      starved_layers.max_layers = 2;
      for (int ray = 0; ray < 100; ++ray) {
        const math::Vector DIRECTION =
            math::Vector{uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)}
                .normalized();
        const Raycast::TraceLimits budgets[] = {
            prepared.limits, starved_candidates, starved_layers};
        for (int budget = 0; budget < 3; ++budget) {
          const auto &limits = budgets[budget];
          if (domain == Domain::SLICE_4D)
            compare(SDF::OctetEvents4(prepared.octet4,
                                      prepared.octet4_projection, DIRECTION,
                                      camera.radial_start, camera.interval.near,
                                      prepared.footprint),
                    camera.interval, limits, prepared.appearance, budget);
          else
            compare(SDF::OctetEvents(prepared.octet_projection, DIRECTION,
                                     camera.radial_start, camera.interval.near,
                                     prepared.footprint),
                    camera.interval, limits, prepared.appearance, budget);
        }
      }
    }
  }
  HS_EXPECT_GT(lit, compared / 10);
  HS_EXPECT_GT(candidate_exhausted, size_t{0});
  HS_EXPECT_GT(layer_exhausted, size_t{0});
}
