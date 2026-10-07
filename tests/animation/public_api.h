/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Public animation APIs
// ============================================================================

/**
 * @brief Verifies a non-repeating RandomTimer fires exactly once, at a frame
 * within the requested inclusive [min, max] delay window.
 */
inline void test_random_timer_fires_within_range() {
  const auto saved_rng = hs::random();
  hs::random().seed(1337);
  Timeline tl;
  struct {
    int fires = 0;
    int fire_frame = -1;
    int frame = 0;
  } st;
  tl.add(0, Animation::RandomTimer({.min = 3, .max = 7}, [&st](Canvas &) {
           st.fires++;
           st.fire_frame = st.frame;
         }));
  for (st.frame = 1; st.frame <= 12; ++st.frame)
    tl.step(fake_canvas());
  HS_EXPECT_EQ(st.fires, 1);
  HS_EXPECT_GE(st.fire_frame, 3);
  HS_EXPECT_LE(st.fire_frame, 7);

  bool saw_minimum = false;
  bool saw_maximum = false;
  for (uint32_t seed = 0; seed < 64; ++seed) {
    hs::random().seed(seed);
    st.fires = 0;
    st.fire_frame = -1;
    tl.add(0, Animation::RandomTimer({.min = 3, .max = 4}, [&st](Canvas &) {
             ++st.fires;
             st.fire_frame = st.frame;
           }));
    for (st.frame = 1; st.frame <= 6; ++st.frame)
      tl.step(fake_canvas());
    HS_EXPECT_EQ(st.fires, 1);
    HS_EXPECT_GE(st.fire_frame, 3);
    HS_EXPECT_LE(st.fire_frame, 4);
    saw_minimum |= st.fire_frame == 3;
    saw_maximum |= st.fire_frame == 4;
  }
  HS_EXPECT_TRUE(saw_minimum);
  HS_EXPECT_TRUE(saw_maximum);
  hs::random() = saved_rng;
}

/**
 * @brief Verifies a one-shot timer ends by completion, not cancellation.
 * @details Timeline's pin-completion guard exempts is_canceled() animations.
 */
inline void test_one_shot_timer_ends_by_completion_not_cancel() {
  Animation::PeriodicTimer periodic(2, [](Canvas &) {}, /*repeat=*/false);
  periodic.step(fake_canvas());
  HS_EXPECT_FALSE(periodic.done());
  periodic.step(fake_canvas()); // t=2: fires and ends itself
  HS_EXPECT_TRUE(periodic.done());
  HS_EXPECT_FALSE(periodic.is_canceled());
  HS_EXPECT_FALSE(periodic.repeats());

  Animation::RandomTimer random({.min = 2, .max = 2}, [](Canvas &) {});
  random.step(fake_canvas());
  HS_EXPECT_FALSE(random.done());
  random.step(fake_canvas());
  HS_EXPECT_TRUE(random.done());
  HS_EXPECT_FALSE(random.is_canceled());
}

/**
 * @brief A repeating animation that ends itself through the protected finish().
 */
class SelfFinishing : public Animation::AnimationBase<SelfFinishing> {
public:
  explicit SelfFinishing(uint32_t at) : AnimationBase(100, true), at(at) {}

  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    if (t >= at)
      finish();
  }

private:
  uint32_t at; /**< Frame at which the animation ends itself. */
};

/**
 * @brief Verifies finish() terminates a repeating animation instead of leaving
 * it done() && repeats().
 * @details Timeline rewinds and re-fires the .then() of anything reporting both.
 */
inline void test_finish_terminates_a_repeating_animation() {
  SelfFinishing anim(2);
  anim.step(fake_canvas());
  HS_EXPECT_FALSE(anim.done());
  anim.step(fake_canvas());
  HS_EXPECT_TRUE(anim.done());
  HS_EXPECT_FALSE(anim.repeats());
  HS_EXPECT_FALSE(anim.is_canceled());

  Timeline tl;
  int fires = 0;
  tl.add(0, SelfFinishing(2).then([&fires]() { fires++; }));
  for (int i = 0; i < 6; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_EQ(fires, 1);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/**
 * @brief Verifies a finished finite parameter animation reports full progress
 * rather than dividing by the zeroed duration.
 */
inline void test_finished_param_animation_progress_is_finite() {
  class SelfFinishingParam
      : public Animation::FiniteParamAnimationBase<SelfFinishingParam> {
  public:
    SelfFinishingParam() : FiniteParamAnimationBase(100, true) {}
    void end_now() { finish(); }
    float progress() const { return normalized_progress(); }
  } anim;
  anim.end_now();
  const float progress = anim.progress();
  HS_EXPECT_TRUE(std::isfinite(progress));
  HS_EXPECT_NEAR(progress, 1.0f, 1e-6f);
}

/**
 * @brief Verifies PeriodicTimer::set_period reschedules the next trigger from
 * now (t + new_period), not from the original schedule.
 */
inline void test_periodic_timer_set_period_reschedules_from_now() {
  struct {
    int fire_frame = -1;
    int fires = 0;
    int frame = 0;
  } st;
  Animation::PeriodicTimer timer(
      5,
      [&st](Canvas &) {
        st.fires++;
        st.fire_frame = st.frame;
      },
      /*repeat=*/true);
  st.frame = 1;
  timer.step(fake_canvas()); // t=1, no trigger (next=5)
  timer.set_period(3);       // reschedule: next = 1 + 3 = 4
  for (st.frame = 2; st.frame <= 4; ++st.frame)
    timer.step(fake_canvas());
  HS_EXPECT_EQ(st.fires, 1);
  HS_EXPECT_EQ(st.fire_frame, 4);
}

/**
 * @brief Verifies PeriodicTimer::set_period called every frame with an
 * unchanged period still lets the timer fire.
 * @details Period 4 over 8 frames must fire at t=4 and t=8; non-positive periods
 * clamp to one frame and keep firing under repeated set_period calls.
 */
inline void test_periodic_timer_set_period_unchanged_does_not_defer() {
  for (int period : {0, -1, std::numeric_limits<int>::min()}) {
    int calls = 0;
    Animation::PeriodicTimer clamped(period, [&](Canvas &) { ++calls; }, true);
    for (int frame = 0; frame < 3; ++frame) {
      clamped.step(fake_canvas());
      HS_EXPECT_EQ(calls, frame + 1);
    }
    clamped.set_period(10);
    clamped.step(fake_canvas());
    HS_EXPECT_EQ(calls, 3);
    clamped.set_period(period);
    clamped.step(fake_canvas());
    HS_EXPECT_EQ(calls, 4);
  }
  struct {
    int fires = 0;
  } st;
  Animation::PeriodicTimer timer(
      4, [&st](Canvas &) { st.fires++; }, /*repeat=*/true);
  for (int frame = 1; frame <= 8; ++frame) {
    timer.set_period(4); // live-slider write with no change
    timer.step(fake_canvas());
    HS_EXPECT_EQ(st.fires, frame / 4);
  }
  HS_EXPECT_EQ(st.fires, 2);
}

/** @brief Verifies invalid ring and line counts keep MobiusFlow finite with a unit product. */
inline void test_mobiusflow_degenerate_inputs_remain_finite() {
  const float NAN_VALUE = std::numeric_limits<float>::quiet_NaN();
  const float INF_VALUE = std::numeric_limits<float>::infinity();
  for (float rings : {-1.0f, NAN_VALUE, INF_VALUE}) {
    for (float lines : {0.0f, NAN_VALUE, INF_VALUE}) {
      math::MobiusParams params;
      Animation::MobiusFlow flow(params, rings, lines, 8, false);
      flow.step(fake_canvas());
      HS_EXPECT_TRUE(std::isfinite(params.a.re));
      HS_EXPECT_TRUE(std::isfinite(params.a.im));
      HS_EXPECT_TRUE(std::isfinite(params.d.re));
      HS_EXPECT_TRUE(std::isfinite(params.d.im));
      HS_EXPECT_NEAR(params.a.re * params.d.re - params.a.im * params.d.im,
                     1.0f, 1e-4f);
      HS_EXPECT_NEAR(params.a.re * params.d.im + params.a.im * params.d.re,
                     0.0f, 1e-4f);
    }
  }
}

/**
 * @brief Verifies MobiusFlow::step keeps the transform's a·d product at unity
 * (a and d are reciprocal: d = 1/a) while actually moving the parameters.
 */
inline void test_mobiusflow_step_preserves_unit_product() {
  math::MobiusParams params;
  const float rings = 2.0f, lines = 4.0f;
  const int duration = 8;
  Animation::MobiusFlow flow(params, rings, lines, duration, /*repeat=*/false);
  for (int i = 0; i < duration; ++i) {
    flow.step(fake_canvas());
    const float re = params.a.re * params.d.re - params.a.im * params.d.im;
    const float im = params.a.re * params.d.im + params.a.im * params.d.re;
    HS_EXPECT_NEAR(re, 1.0f, 1e-4f);
    HS_EXPECT_NEAR(im, 0.0f, 1e-4f);
  }
  HS_EXPECT_GT(std::abs(params.a.im), 1e-3f); // moved off the identity
}

/**
 * @brief Verifies a registered emitter runs on every ParticleSystem::step and
 * can spawn into the pool.
 */
inline void test_particle_system_emitter_dispatch() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  ps.init(arena);
  int calls = 0;
  ps.add_emitter([&](Animation::ParticleSystem<32, 4> &sys) {
    calls++;
    sys.spawn(math::Vector(1, 0, 0), math::Vector(0, 0, 0), 0);
  });
  ps.step(fake_canvas());
  HS_EXPECT_EQ(calls, 1);
  HS_EXPECT_GT(static_cast<int>(ps.active()), 0);
  ps.step(fake_canvas());
  HS_EXPECT_EQ(calls, 2); // runs every frame
}

/**
 * @brief Verifies Motion::set_duration rescales the elapsed frame count instead
 * of completing the motion when the new duration is below the current position.
 */
inline void test_motion_set_duration_below_position_rescales() {
  using Ori = math::Orientation<16>;
  ProceduralPath path;
  path.f = [](float t) {
    float a = 2.0f * math::PI_F * t;
    return math::Vector(std::cos(a), std::sin(a), 0.0f);
  };
  Ori o;
  Animation::Motion<288, 16> motion(o, path, 60, /*repeat=*/true);

  const math::Vector probe = math::X_AXIS;
  for (int i = 0; i < 40; ++i)
    motion.step(fake_canvas());
  HS_EXPECT_FALSE(motion.done());
  const math::Vector before = o.orient(probe);

  // 40 frames into a 60-frame loop; 20 frames is behind that position.
  motion.set_duration(20);
  HS_EXPECT_FALSE(motion.done());

  // Reanchoring keeps the first relative step incremental.
  motion.step(fake_canvas());
  HS_EXPECT_LT(math::angle_between(before, o.orient(probe)), 0.5f);
  for (int i = 1; i < 6; ++i)
    motion.step(fake_canvas());
  HS_EXPECT_FALSE(motion.done());
  motion.step(fake_canvas());
  HS_EXPECT_TRUE(motion.done());
}

/** @brief Verifies Progress pause behavior and eased output bounds. */
inline void test_progress_pause_and_eased_bounds() {
  bool paused = true;
  std::vector<float> values;
  Animation::Progress progress([&](float t) { values.push_back(t); }, 4,
                               [](float t) { return 2.0f * t; },
                               {.paused = &paused});
  progress.step(fake_canvas());
  HS_EXPECT_TRUE(values.empty());
  paused = false;
  progress.step(fake_canvas());
  HS_EXPECT_EQ(values.size(), size_t{1});
  HS_EXPECT_EQ(values.back(), 0.5f);
  paused = true;
  progress.step(fake_canvas());
  HS_EXPECT_EQ(values.size(), size_t{1});
  paused = false;
  for (int i = 0; i < 3; ++i)
    progress.step(fake_canvas());
  HS_EXPECT_EQ(values.size(), size_t{4});
  HS_EXPECT_EQ(values[1], 1.0f);
  HS_EXPECT_EQ(values[2], 1.0f);
  HS_EXPECT_EQ(values[3], 1.0f);
  HS_EXPECT_TRUE(progress.done());
}

/** @brief Verifies recorded trail orientations remain independent of the body's live rotation. */
inline void test_trail_body_records_independent_orientation_history() {
  Animation::TrailBody<2, 2> body;
  HS_EXPECT_VEC(body.v, math::Y_AXIS, 0.0f);
  HS_EXPECT_EQ(body.trail.length(), size_t{0});
  body.trail.record(body.orientation);
  body.orientation.set(math::make_rotation(math::Z_AXIS, math::PI_F * 0.5f));
  HS_EXPECT_VEC(body.trail.get(0).orient(body.v), math::Y_AXIS, 1e-6f);
  body.trail.record(body.orientation);
  body.orientation.set(math::Quaternion());
  body.trail.record(body.orientation);
  HS_EXPECT_EQ(body.trail.length(), size_t{2});
  HS_EXPECT_VEC(body.trail.get(0).orient(body.v), -math::X_AXIS, 1e-5f);
  HS_EXPECT_VEC(body.trail.get(1).orient(body.v), math::Y_AXIS, 1e-6f);
}
