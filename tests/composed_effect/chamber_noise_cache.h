/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file chamber_noise_cache.h
 * @brief Temporal chamber-noise ownership and sampling tests. */

/** @brief Chamber patches retain old endpoints and match their complete source identity. */
inline void test_chamber_noise_cache_bake_lifetime() {
  using namespace Pullback;
  Modulation::Config config;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.enabled = true;
  config.seed = 6047;
  config.noise.scale = 2.4029531f;
  config.noise.speed = .00015f;
  Modulation::Resources<false> source;
  auto clock = Modulation::init(config);
  Modulation::ChamberNoiseCache cache;
  auto a = cache.prepare(source.prepare(config, clock), math::Quaternion());
  HS_EXPECT_TRUE(a.patch != nullptr);
  HS_EXPECT_FALSE(source.prepared.cube.active);
  const auto saved_a = a.patch->values;
  for (int y : {0, 7, 23})
    for (int x : {0, 8, 15}) {
      const float U = -.002f + (.385f + .002f) * x / 15;
      const float V = -.002f + (.622f + .002f) * y / 23;
      const auto direction =
          Pattern::SphereDirection::admit(math::Vector(U, 1, V));
      const auto &noise = std::get<0>(source.prepared.source);
      const auto q = direction.value() * config.noise.scale + noise.loop_offset;
      const float expected =
          hs::clamp(noise.noise.GetNoiseSingle(q.x, q.y, q.z), -1.0f, 1.0f);
      HS_EXPECT_NEAR(a.patch->values[y * 16 + x], expected, 1e-5f);
    }
  auto same = cache.prepare(source.prepared, math::Quaternion());
  HS_EXPECT_TRUE(same.patch == a.patch);
  const auto first_source = source.prepared;
  clock.phase.phase_turns = .00015f;
  auto b = cache.prepare(source.prepare(config, clock), math::Quaternion());
  HS_EXPECT_TRUE(a.patch != b.patch);
  HS_EXPECT_TRUE(a.patch->values == saved_a);
  const auto saved_b = b.patch->values;
  clock.phase.phase_turns = .0003f;
  auto c = cache.prepare(source.prepare(config, clock), math::Quaternion());
  HS_EXPECT_TRUE(c.patch == a.patch);
  HS_EXPECT_TRUE(b.patch->values == saved_b);
  HS_EXPECT_FALSE(a.patch->matches(first_source, math::Quaternion()));
  auto changed_projection = math::Quaternion(1.0f, 0.0f, .01f, 0.0f);
  HS_EXPECT_FALSE(c.patch->matches(source.prepared, changed_projection));
  auto unsupported = source.prepared;
  unsupported.config.noise.variant = Pattern::NoiseVariant::OCTAVE_CONTOUR;
  HS_EXPECT_TRUE(cache.prepare(unsupported, math::Quaternion()).patch ==
                 nullptr);
  unsupported = source.prepared;
  unsupported.config.enabled = false;
  HS_EXPECT_TRUE(cache.prepare(unsupported, math::Quaternion()).patch ==
                 nullptr);
}
/** @brief Unqualified controls and discontinuous endpoints never bake a patch. */
inline void test_chamber_noise_cache_admission() {
  using namespace Pullback;
  auto config = AshCloud<96, 20>::initial_params()
                    .template get<"modulation">()
                    .modulation;
  HS_EXPECT_TRUE(Modulation::ChamberNoisePatch::eligible(config));
  Modulation::Resources<false> source;
  auto clock = Modulation::init(config);
  Modulation::ChamberNoiseCache cache;
  const auto first =
      cache.prepare(source.prepare(config, clock), math::Quaternion());
  HS_EXPECT_TRUE(first.patch != nullptr);
  const auto values = first.patch->values;
  const uint8_t SLOT = cache.slot;
  for (int change = 0; change < 7; ++change) {
    auto unsupported = source.prepared;
    if (change == 0)
      unsupported.config.noise.scale = 64.0f;
    if (change == 1)
      unsupported.config.noise.scale = 3.0f;
    if (change == 2)
      unsupported.config.noise.speed = 0.001f;
    if (change == 3)
      unsupported.config.noise.speed = 0.0f;
    if (change == 4)
      ++unsupported.config.seed;
    if (change == 5)
      unsupported.config.version = 2;
    if (change == 6)
      unsupported.config.noise.variant = Pattern::NoiseVariant::OCTAVE_CONTOUR;
    HS_EXPECT_FALSE(
        Modulation::ChamberNoisePatch::eligible(unsupported.config));
    HS_EXPECT_TRUE(cache.prepare(unsupported, math::Quaternion()).patch ==
                   nullptr);
    HS_EXPECT_EQ(cache.slot, SLOT);
    HS_EXPECT_TRUE(first.patch->values == values);
    HS_EXPECT_FALSE(first.patch->matches(unsupported, math::Quaternion()));
  }
  HS_EXPECT_TRUE(cache.prepare(source.prepared, math::Quaternion()).patch ==
                 first.patch);
  clock.phase.phase_turns = 0.00015f;
  HS_EXPECT_TRUE(
      cache.prepare(source.prepare(config, clock), math::Quaternion()).patch !=
      nullptr);
  const uint8_t BEFORE_JUMP = cache.slot;
  const auto old_values = cache.patches[BEFORE_JUMP].values;
  clock.phase.phase_turns = 0.1f;
  source.prepare(config, clock);
  HS_EXPECT_TRUE(cache.prepare(source.prepared, math::Quaternion()).patch ==
                 nullptr);
  HS_EXPECT_TRUE(cache.prepare(source.prepared, math::Quaternion()).patch ==
                 nullptr);
  HS_EXPECT_EQ(cache.slot, BEFORE_JUMP);
  HS_EXPECT_TRUE(cache.patches[BEFORE_JUMP].values == old_values);
  clock.phase.phase_turns += 0.00015f;
  HS_EXPECT_TRUE(
      cache.prepare(source.prepare(config, clock), math::Quaternion()).patch !=
      nullptr);
  clock.phase.phase_turns += 0.00015f;
  const math::Quaternion jumped(cosf(0.05f), 0.0f, sinf(0.05f), 0.0f);
  const uint8_t BEFORE_ROTATION = cache.slot;
  HS_EXPECT_TRUE(cache.prepare(source.prepare(config, clock), jumped).patch ==
                 nullptr);
  HS_EXPECT_EQ(cache.slot, BEFORE_ROTATION);
  HS_EXPECT_TRUE(Modulation::ChamberNoisePatch::continuous(
      0.9999f, math::Quaternion(), 0.0001f, math::Quaternion()));
  HS_EXPECT_TRUE(Modulation::ChamberNoisePatch::continuous(
      0.0f, math::Quaternion(), 0.00015f, math::Quaternion(-1, 0, 0, 0)));
  HS_EXPECT_FALSE(Modulation::ChamberNoisePatch::continuous(
      0.0f, math::Quaternion(), 0.01f, math::Quaternion()));
}

/** @brief Cached noise interpolates both endpoints and falls back outside its chart. */
inline void test_chamber_noise_cache_phase_and_fallback() {
  using namespace Pullback;
  using FX = AshCloud<96, 20>;
  using Provider =
      ModulationProvider<FX::Binding, "modulation", Pattern::SourceKind::NOISE,
                         Pattern::Evaluation::DIRECT, true>;
  static_assert(
      std::is_same_v<Provider::Prepared, Color::ChamberPhaseModulationState>);
  using Grid = ModulationProvider<Millefiori<96, 20>::Binding, "field",
                                  Pattern::SourceKind::SPHERE_SMOOTH_GRID,
                                  Pattern::Evaluation::DIRECT, true>;
  static_assert(std::is_same_v<Grid::Prepared, Color::PhaseModulationState>);
  static_assert(Color::Modulate<Provider>::PHASE_NOISE_APPROXIMATION);
  static_assert(!Color::Modulate<Grid>::PHASE_NOISE_APPROXIMATION);
  reset_globals();
  FX effect;
  effect.init();
  auto from = ComposedFrameWhiteBox::frame(effect);
  auto to = from;
  auto config = from.params.template get<"modulation">().modulation;
  config.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, .25f};
  from.params.template get<"modulation">().modulation = config;
  to.params.template get<"modulation">().modulation = config;
  from.projection_conjugate = math::Quaternion();
  to.projection_conjugate = math::Quaternion(cosf(.001f), 0, sinf(.001f), 0);
  Modulation::Resources<false> previous, current;
  auto clock = Modulation::init(config);
  clock.phase.phase_turns = .9999f;
  previous.prepare(config, clock);
  clock.phase.phase_turns = .0001f;
  current.prepare(config, clock);
  Modulation::ChamberNoiseCache cache;
  const auto a = cache.prepare(previous.prepared, from.projection_conjugate);
  const auto b = cache.prepare(current.prepared, to.projection_conjugate);
  from.resources.template get<"modulation">().modulation = &previous.prepared;
  to.resources.template get<"modulation">().modulation = &current.prepared;
  from.resources.template get<"modulation">().phase_noise_current = a;
  to.resources.template get<"modulation">().phase_noise_current = b;
  auto *ap = const_cast<Modulation::ChamberNoisePatch *>(a.patch);
  auto *bp = const_cast<Modulation::ChamberNoisePatch *>(b.patch);
  for (int y = 0; y < 24; ++y)
    for (int x = 0; x < 16; ++x) {
      ap->values[y * 16 + x] = .01f * x + .02f * y - .3f;
      bp->values[y * 16 + x] = .01f * x + .02f * y - .2f;
    }
  for (float fraction : {0.0f, .25f, .5f, .75f, 1.0f}) {
    auto column = FrameInterval(from, to).at(fraction);
    auto prepared = Provider::prepare(column);
    HS_EXPECT_EQ(prepared.ready, fraction > 0 && fraction < 1);
    if (!prepared.ready)
      continue;
    HS_EXPECT_TRUE(prepared.noise_current == b.patch);
    HS_EXPECT_TRUE(prepared.noise_previous == a.patch);
    const float U = .1915f, V = .31f;
    const auto sphere = math::rotate(math::Vector(U, 1, V).normalized(),
                                     column.projection_conjugate);
    ColorSample carrier{Color4(Pixel(2000, 4000, 6000), .37f),
                        {.5f, .7f, sphere, 2.0f}};
    Color::Modulate<Provider>::apply_phase_in_place(carrier, column, prepared);
    const float X = (U + .002f) / (.385f + .002f) * 15;
    const float Y = (V + .002f) / (.622f + .002f) * 23;
    const float expected = .25f * (.01f * X + .02f * Y - .3f + .1f * fraction);
    HS_EXPECT_NEAR(carrier.pending_hue, math::wrap_t(expected), 2e-6f);
    HS_EXPECT_EQ(carrier.color.alpha, .37f);
    for (const auto direction :
         {math::Vector(0, -1, 0), math::Vector(.8f, 1, .2f).normalized()}) {
      const auto local = math::rotate(direction, column.projection_conjugate);
      ColorSample fallback{Color4(Pixel(2000, 4000, 6000), .37f),
                           {.5f, .7f, local, 2.0f}};
      Color::Modulate<Provider>::apply_phase_in_place(fallback, column,
                                                      prepared);
      const float advance = (.0001f - .9999f) + 1.0f;
      const float phase = math::wrap_t(.9999f + advance * fraction);
      const float angle = math::TWO_PI_F * phase;
      const auto q =
          Pattern::SphereDirection::admit(local).value() * config.noise.scale +
          math::Vector(32 * cosf(angle), 32 * sinf(angle), 0);
      const auto &noise = std::get<0>(current.prepared.source).noise;
      const float value =
          hs::clamp(noise.GetNoiseSingle(q.x, q.y, q.z), -1.0f, 1.0f);
      HS_EXPECT_NEAR(fallback.pending_hue, math::wrap_t(.25f * value), 1e-6f);
    }
  }
  auto column = FrameInterval(from, to).at(.5f);
  column.resources.template get<"modulation">().phase_noise_current = {};
  HS_EXPECT_TRUE(Provider::prepare(column).noise_current == nullptr);
  HS_EXPECT_TRUE(Provider::prepare(column).ready);
  column = FrameInterval(from, to).at(.5f);
  column.resources.template get<"modulation">()
      .phase_noise_current.projection.r += .1f;
  HS_EXPECT_TRUE(Provider::prepare(column).noise_current == nullptr);
  auto changed = to;
  changed.params.template get<"projection">().spin_rate += .01f;
  auto projection_changed = FrameInterval(from, changed).at(.5f);
  HS_EXPECT_TRUE(Provider::prepare(projection_changed).noise_current ==
                 nullptr);
  HS_EXPECT_TRUE(Provider::prepare(projection_changed).ready);
  changed = to;
  changed.params.template get<"modulation">().modulation.seed += 1;
  HS_EXPECT_FALSE(
      Provider::prepare(FrameInterval(from, changed).at(.5f)).ready);
  Modulation::Resources<false> equal_clock_source;
  Modulation::ChamberNoiseCache equal_clock_cache;
  auto equal_clock = Modulation::init(config);
  equal_clock.phase.phase_turns = 0.9999f;
  auto equal_frame = to;
  equal_frame.resources.template get<"modulation">().modulation =
      &equal_clock_source.prepare(config, equal_clock);
  equal_frame.resources.template get<"modulation">().phase_noise_current =
      equal_clock_cache.prepare(equal_clock_source.prepared,
                                equal_frame.projection_conjugate);
  const auto equal_column = FrameInterval(from, equal_frame).at(0.5f);
  HS_EXPECT_EQ(equal_column.temporal.fraction, 1.0f);
  const auto equal_prepared = Provider::prepare(equal_column);
  HS_EXPECT_TRUE(equal_prepared.noise_current != nullptr);
  HS_EXPECT_EQ(equal_prepared.noise_fraction, 0.5f);
  Modulation::Resources<false> far_source;
  auto far_clock = Modulation::init(config);
  far_clock.phase.phase_turns = 0.1f;
  Modulation::ChamberNoiseCache far_cache;
  auto far_frame = to;
  far_frame.resources.template get<"modulation">().modulation =
      &far_source.prepare(config, far_clock);
  far_frame.resources.template get<"modulation">().phase_noise_current =
      far_cache.prepare(far_source.prepared, far_frame.projection_conjugate);
  HS_EXPECT_TRUE(far_frame.resources.template get<"modulation">()
                     .phase_noise_current.patch != nullptr);
  auto far_column = FrameInterval(from, far_frame).at(0.5f);
  HS_EXPECT_TRUE(Provider::prepare(far_column).ready);
  HS_EXPECT_TRUE(Provider::prepare(far_column).noise_current == nullptr);
  far_clock.phase.phase_turns = 0.0001f;
  Modulation::ChamberNoiseCache rotated_cache;
  far_frame.projection_conjugate =
      math::Quaternion(cosf(0.05f), 0, sinf(0.05f), 0);
  far_frame.resources.template get<"modulation">().modulation =
      &far_source.prepare(config, far_clock);
  far_frame.resources.template get<"modulation">().phase_noise_current =
      rotated_cache.prepare(far_source.prepared,
                            far_frame.projection_conjugate);
  HS_EXPECT_TRUE(far_frame.resources.template get<"modulation">()
                     .phase_noise_current.patch != nullptr);
  HS_EXPECT_TRUE(Provider::prepare(FrameInterval(from, far_frame).at(0.5f))
                     .noise_current == nullptr);
  const float saved = current.prepared.state.phase.phase_turns;
  current.prepared.state.phase.phase_turns = .1f;
  HS_EXPECT_TRUE(
      Provider::prepare(FrameInterval(from, to).at(.5f)).noise_current ==
      nullptr);
  current.prepared.state.phase.phase_turns = saved;
}
/** @brief Consecutive draws preserve the retained patch without activating a source cube. */
inline void test_chamber_noise_cache_runtime_history() {
  using FX = AshCloud<96, 20>;
  reset_globals();
  pin_frame_clock(0);
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  const auto first = ComposedFrameWhiteBox::current_frame(effect);
  const auto a =
      first.resources.template get<"modulation">().phase_noise_current;
  HS_EXPECT_TRUE(a.patch != nullptr);
  const auto saved_a = a.patch->values;
  const auto first_source =
      *first.resources.template get<"modulation">().modulation;
  pin_frame_clock(1);
  effect.draw_frame();
  effect.advance_display();
  const auto second = ComposedFrameWhiteBox::current_frame(effect);
  const auto b =
      second.resources.template get<"modulation">().phase_noise_current;
  HS_EXPECT_TRUE(a.patch != b.patch);
  HS_EXPECT_TRUE(a.patch->values == saved_a);
  const auto saved_b = b.patch->values;
  const auto same = ComposedFrameWhiteBox::frame(effect);
  HS_EXPECT_TRUE(
      same.resources.template get<"modulation">().phase_noise_current.patch ==
      b.patch);
  pin_frame_clock(2);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(b.patch->values == saved_b);
  HS_EXPECT_FALSE(a.patch->matches(first_source, a.projection));
  HS_EXPECT_FALSE(ComposedFrameWhiteBox::current_frame(effect)
                      .resources.template get<"modulation">()
                      .modulation->cube.active);
  const auto before_static_a = a.patch->values;
  const auto before_static_b = b.patch->values;
  Spacetime::set_flat_surface(true);
  pin_frame_clock(3);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(ComposedFrameWhiteBox::current_frame(effect)
                     .resources.template get<"modulation">()
                     .phase_noise_current.patch == nullptr);
  HS_EXPECT_TRUE(a.patch->values == before_static_a);
  HS_EXPECT_TRUE(b.patch->values == before_static_b);
  Spacetime::set_flat_surface(false);
}
