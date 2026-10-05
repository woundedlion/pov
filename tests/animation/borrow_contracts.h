/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// ============================================================================
// Borrow-contract guards (compile-time)
// ----------------------------------------------------------------------------
// Motion/Lerp/MobiusFlow store non-owning borrows of effect-owned
// state. These static_asserts lock the contract: an effect-owned lvalue is
// accepted, a temporary (which would dangle) is rejected.
// ============================================================================
namespace borrow_guard {
/** @brief Orientation alias used by the borrow-contract static_asserts. */
using Ori = math::Orientation<16>;

static_assert(std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                      ProceduralPath &, int>,
              "Motion must accept an lvalue (effect-owned) path");
static_assert(!std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                       ProceduralPath &&, int>,
              "Motion must REJECT a temporary path (would dangle)");

static_assert(
    std::is_constructible_v<Animation::Lerp, Lerpable &, const Lerpable &,
                            const Lerpable &, int, EasingFn>,
    "Lerp must accept lvalue (effect-owned) start/target");
static_assert(
    !std::is_constructible_v<Animation::Lerp, Lerpable &, const Lerpable &&,
                             const Lerpable &, int, EasingFn>,
    "Lerp must REJECT a temporary start (would dangle)");
static_assert(
    !std::is_constructible_v<Animation::Lerp, Lerpable &, const Lerpable &,
                             const Lerpable &&, int, EasingFn>,
    "Lerp must REJECT a temporary target (would dangle)");

static_assert(
    std::is_constructible_v<Animation::ColorWipe, GenerativePalette &,
                            const GenerativePalette::Snapshot &,
                            const GenerativePalette::Snapshot &, int, EasingFn>,
    "ColorWipe must accept effect-owned snapshots");
static_assert(
    !std::is_constructible_v<Animation::ColorWipe, GenerativePalette &,
                             GenerativePalette::Snapshot &&,
                             const GenerativePalette::Snapshot &, int,
                             EasingFn>,
    "ColorWipe must REJECT a temporary start snapshot (would dangle)");
static_assert(
    !std::is_constructible_v<Animation::ColorWipe, GenerativePalette &,
                             const GenerativePalette::Snapshot &,
                             GenerativePalette::Snapshot &&, int, EasingFn>,
    "ColorWipe must REJECT a temporary target snapshot (would dangle)");

static_assert(
    std::is_constructible_v<Animation::MobiusFlow, math::MobiusParams &,
                            const float &, const float &, int>,
    "MobiusFlow must accept lvalue (effect-owned) scalars");
static_assert(
    !std::is_constructible_v<Animation::MobiusFlow, math::MobiusParams &,
                             const float &&, const float &, int>,
    "MobiusFlow must REJECT a temporary num_rings (would dangle)");
static_assert(
    !std::is_constructible_v<Animation::MobiusFlow, math::MobiusParams &,
                             const float &, const float &&, int>,
    "MobiusFlow must REJECT a temporary num_lines (would dangle)");
} // namespace borrow_guard

/**
 * @brief Verifies rotation_substeps returns a tight ceil with each sub-interval
 * within MAX.
 */
inline void test_rotation_substeps_shared_and_tight() {
  constexpr float MAX = 0.1f;
  // Always at least 1, even for a sub-threshold angle.
  HS_EXPECT_EQ(Animation::rotation_substeps(0.0f, MAX), 1);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 0.5f, MAX), 1);
  // Tight ceil: N*MAX needs exactly N subdivisions.
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX, MAX), 1);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.0f, MAX), 3);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.2f, MAX), 4);
  // Every sub-interval stays within MAX.
  for (float a = 0.0f; a < 2.0f; a += 0.013f) {
    int n = Animation::rotation_substeps(a, MAX);
    HS_EXPECT_GE(n, 1);
    HS_EXPECT_LE(a / n, MAX + 1e-6f);
  }
}

/**
 * @brief Verifies Rotation::step does not discard sub-MIN_STEP_ANGLE increments.
 * @details A rotation slow enough that each frame's delta is below MIN_STEP_ANGLE
 * (1e-4 rad) must still accumulate those deltas so the orientation actually
 * turns over many frames. ease_linear is linear here, so the per-frame delta is
 * total_angle / duration.
 */
inline void test_rotation_accumulates_subthreshold_deltas() {
  using Ori = math::Orientation<16>;
  Ori o; // identity
  // 0.05 rad over 1000 frames => 5e-5 rad/frame, half of MIN_STEP_ANGLE, so every
  // frame's raw delta is below the early-out threshold.
  Animation::Rotation<288, 16> rot(o, math::Z_AXIS, 0.05f, 1000,
                                   math::ease_linear);
  for (int i = 0; i < 20; ++i)
    rot.step(fake_canvas());
  // A Z rotation sends +X toward +Y.
  math::Vector v = o.orient(math::X_AXIS, o.length() - 1);
  HS_EXPECT_GT(v.y, 5e-4f);
  HS_EXPECT_NEAR(v.x, 1.0f, 1e-3f);
}

/**
 * @brief Verifies a repeating Rotation lands its full sweep every cycle.
 * @details An easing with zero slope at t=1 leaves a sub-MIN_STEP_ANGLE residual on
 * the final frame. That frame has no successor to accumulate into, so dropping
 * it slips the residual (~1e-4 rad here) once per cycle.
 */
inline void test_rotation_applies_final_frame_residual() {
  using Ori = math::Orientation<16>;
  constexpr float ANGLE = 0.2f;
  constexpr int DURATION = 1000;
  constexpr int CYCLES = 10;
  Ori o; // identity
  Animation::Rotation<288, 16> rot(o, math::Z_AXIS, ANGLE, DURATION,
                                   math::ease_in_out_sin);
  for (int c = 0; c < CYCLES; ++c) {
    for (int i = 0; i < DURATION; ++i)
      rot.step(fake_canvas());
    rot.rewind();
  }
  // A Z rotation sends +X to (cos, sin) of the accumulated angle.
  math::Vector v = o.orient(math::X_AXIS, o.length() - 1);
  HS_EXPECT_NEAR(std::atan2(v.y, v.x), ANGLE * CYCLES, 2e-4f);
}

/**
 * @brief Verifies two animations sharing one Orientation COMPOSE their
 * sub-frame motion-blur history within a frame instead of clobbering it.
 * @details A per-animation collapse would discard the first animation's
 * freshly-built sub-frame trail. The decisive signature is the OLDEST sub-frame
 * (index 0): with composition it still reflects the pre-frame orientation
 * (identity here). This test touches the global Timeline; Rotation::step never
 * dereferences the canvas.
 */
inline void test_timeline_shared_orientation_composes_motion_blur() {
  using Ori = math::Orientation<16>;
  Ori o; // identity, single frame
  Timeline tl;
  // Two quarter-turn rotations about the same axis, each completing in one frame.
  tl.add(0, Animation::Rotation<288, 16>(o, math::Z_AXIS, math::PI_F / 2, 1,
                                         math::ease_linear));
  tl.add(0, Animation::Rotation<288, 16>(o, math::Z_AXIS, math::PI_F / 2, 1,
                                         math::ease_linear));
  tl.step(fake_canvas());

  HS_EXPECT_GE(o.length(), 2);
  // Oldest sub-frame is the pre-frame orientation (identity): +X stays +X.
  math::Vector oldest = o.orient(math::X_AXIS, 0);
  HS_EXPECT_NEAR(oldest.x, 1.0f, 1e-3f);
  HS_EXPECT_NEAR(oldest.y, 0.0f, 1e-3f);
  // Newest reflects both rotations (a half turn): +X -> -X.
  math::Vector newest = o.orient(math::X_AXIS, o.length() - 1);
  HS_EXPECT_NEAR(newest.x, -1.0f, 1e-3f);
}

/**
 * @brief Verifies every distinct Orientation is still collapsed exactly once
 * when the distinct count exceeds step()'s id cache and the pass falls back to
 * rescanning the earlier events.
 * @details Each Orientation carries a two-frame history into the frame, so a
 * collapse leaves the +90 pose as the oldest sub-frame (+X -> +Y) while a
 * skipped one leaves the identity it started at (+X stays +X).
 */
inline void test_timeline_collapse_past_id_cache() {
  struct CountedRotation : Animation::Rotation<288, 16> {
    int &collapses;

    CountedRotation(math::Orientation<16> &orientation, int &collapses)
        : Animation::Rotation<288, 16>(orientation, math::Z_AXIS,
                                       math::PI_F / 2, 1, math::ease_linear),
          collapses(collapses) {}

    void collapse_orientation() override {
      ++collapses;
      Animation::Rotation<288, 16>::collapse_orientation();
    }
  };

  constexpr int N = Timeline::MAX_COLLAPSE_IDS + 2;
  math::Orientation<16> orientations[N];
  int collapses[N] = {};
  Timeline tl;
  for (int i = 0; i < N; ++i) {
    orientations[i].push(math::make_rotation(math::Z_AXIS, math::PI_F / 2));
    tl.add(0, CountedRotation(orientations[i], collapses[i]));
  }
  tl.add(0, CountedRotation(orientations[N - 1], collapses[N - 1]));
  tl.step(fake_canvas());

  for (int i = 0; i < N; ++i) {
    math::Vector oldest = orientations[i].orient(math::X_AXIS, 0);
    HS_EXPECT_NEAR(oldest.y, 1.0f, 1e-3f);
    HS_EXPECT_EQ(collapses[i], 1);
  }
}

/**
 * @brief Verifies Timeline schedules events by start frame and removes
 * completed one-shots so later events can run after earlier ones finish.
 * @details An event added with in_frames > 0 stays dormant until t reaches its
 * start, then steps. Uses the global Timeline; each Timeline is scoped so the
 * live-guard balances, and Transition::step never dereferences the canvas.
 */
inline void test_timeline_sequences_events_by_start_frame() {
  Timeline tl;
  float a = 0.0f, b = 0.0f;
  tl.add(0,
         Animation::Transition(a, 10.0f, 2, math::ease_linear)); // starts now
  tl.add(3, Animation::Transition(b, 20.0f, 2,
                                  math::ease_linear)); // delayed 3 frames

  tl.step(fake_canvas()); // t=1
  HS_EXPECT_GT(a, 0.0f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f);

  tl.step(fake_canvas()); // t=2: a completes and is removed
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f); // still dormant (start=3)
  HS_EXPECT_EQ(tl.event_count(), 1);

  tl.step(fake_canvas()); // t=3: b starts
  HS_EXPECT_GT(b, 0.0f);
  HS_EXPECT_LT(b, 20.0f);

  tl.step(fake_canvas()); // t=4: b completes
  HS_EXPECT_NEAR(b, 20.0f, 1e-3f);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/** @brief Verifies event-level pause freezes delays, steps, and callbacks. */
inline void test_timeline_pausable_event_uses_active_time() {
  Timeline tl;
  bool paused = true;
  float value = 0.0f;
  float ambient = 0.0f;
  int completions = 0;
  tl.add_pausable(
      3, Animation::Transition(value, 9.0f, 3, math::ease_linear).then([&]() {
        ++completions;
      }),
      &paused);
  tl.add(0, Animation::Transition(ambient, 1.0f, 2, math::ease_linear));

  for (int i = 0; i < 5; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(ambient, 1.0f, 1e-6f);
  HS_EXPECT_EQ(completions, 0);

  paused = false;
  tl.step(fake_canvas());
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 0.0f, 1e-6f);
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 3.0f, 1e-6f);

  paused = true;
  for (int i = 0; i < 5; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 3.0f, 1e-6f);
  HS_EXPECT_EQ(completions, 0);

  paused = false;
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 6.0f, 1e-6f);
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 9.0f, 1e-6f);
  HS_EXPECT_EQ(completions, 1);
}

/**
 * @brief Verifies the latest representable start frame remains schedulable.
 */
inline void test_timeline_accepts_maximum_start_frame() {
  Timeline tl;
  global_timeline_t = UINT32_MAX - 2;
  float value = 0.0f;
  tl.add(2, Animation::Transition(value, 1.0f, 1, math::ease_linear));

  HS_EXPECT_EQ(global_timeline_events[0].start, UINT32_MAX);
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 0.0f, 1e-6f);
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 1.0f, 1e-6f);
  HS_EXPECT_EQ(global_timeline_t, UINT32_MAX);
}

/**
 * @brief Verifies a repeating animation rewinds at the end of each cycle
 * instead of being removed, replaying the curve.
 * @details Mutation writes f(easing(t/duration)) each step, so after completion
 * and rewind the next step drops back to the mid-cycle value; a non-rewinding
 * timer would clamp at 1.
 */
inline void test_timeline_repeating_animation_rewinds_each_cycle() {
  Timeline tl;
  float v = -1.0f;
  tl.add(0, Animation::Mutation(
                v, [](float e) { return e; }, 2, math::ease_linear,
                /*repeat=*/true));

  tl.step(fake_canvas()); // t=1 -> v = eased(0.5) = 0.5
  HS_EXPECT_NEAR(v, 0.5f, 1e-3f);
  tl.step(fake_canvas()); // t=2 -> v = eased(1.0) = 1.0, done -> rewind
  HS_EXPECT_NEAR(v, 1.0f, 1e-3f);
  tl.step(fake_canvas()); // rewound to t=1 -> v = 0.5
  HS_EXPECT_NEAR(v, 0.5f, 1e-3f);

  HS_EXPECT_EQ(tl.event_count(), 1);
}

/**
 * @brief Verifies a canceled repeating animation is removed, not rewound and
 * replayed forever.
 * @details cancel() makes done() permanently true; repeats() must drop on
 * cancel so Timeline routes it through the removal branch instead of keeping it
 * as a per-frame, callback-firing zombie.
 */
inline void test_timeline_cancel_removes_repeating_animation() {
  Timeline tl;
  float v = -1.0f;
  auto *h = tl.add_get(0,
                       Animation::Mutation(
                           v, [](float e) { return e; }, 2, math::ease_linear,
                           /*repeat=*/true),
                       Timeline::Pin::PINNED);
  tl.step(fake_canvas());
  tl.step(fake_canvas()); // completes a cycle, rewinds, and stays (repeating)
  HS_EXPECT_EQ(tl.event_count(), 1);

  const float before_cancel = v;
  h->cancel();
  tl.step(fake_canvas()); // canceled: done() && !repeats() -> removed
  HS_EXPECT_EQ(tl.event_count(), 0);
  HS_EXPECT_EQ(v, before_cancel);

  // A removed event stops stepping: v stays frozen instead of oscillating.
  const float v_frozen = v;
  for (int i = 0; i < 6; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_NEAR(v, v_frozen, 1e-6f);
}

/**
 * @brief Verifies a repeating animation canceled from inside its own .then()
 * fires the callback exactly once and is removed in that frame.
 * @details The repeating branch rewinds and fires the per-cycle callback; a
 * cancel taken there leaves the animation done() and non-repeating, so keeping
 * the event would route it through the removal branch on the next frame and
 * fire .then() again.
 */
inline void test_timeline_repeating_canceled_in_callback_fires_then_once() {
  Timeline tl;
  float v = -1.0f;
  struct {
    int thens = 0;
    Animation::Mutation *anim = nullptr;
  } st; // one capture keeps the callback inside Fn's inplace budget
  st.anim = tl.add_get(0,
                       Animation::Mutation(
                           v, [](float e) { return e; }, 2, math::ease_linear,
                           /*repeat=*/true)
                           .then([&st]() {
                             st.thens++;
                             st.anim->cancel();
                           }),
                       Timeline::Pin::UNPINNED);
  for (int i = 0; i < 6; ++i)
    tl.step(fake_canvas()); // cycles complete at t=2, 4, 6 without the cancel
  HS_EXPECT_EQ(st.thens, 1);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/** @brief Cancellation suppresses due timer callbacks and paused drawing. */
inline void test_timeline_cancel_suppresses_step_side_effects() {
  {
    Timeline tl;
    int fires = 0;
    auto *timer = tl.add_get(0,
                             Animation::PeriodicTimer(
                                 3, [&](Canvas &) { ++fires; }, true),
                             Timeline::Pin::UNPINNED);
    tl.step(fake_canvas());
    tl.step(fake_canvas());
    timer->cancel();
    tl.step(fake_canvas());
    HS_EXPECT_EQ(fires, 0);
    HS_EXPECT_EQ(tl.event_count(), 0);
  }
  {
    struct PausedProbe : Animation::AnimationBase<PausedProbe> {
      int *draws;
      explicit PausedProbe(int &count) : draws(&count) {}
      void step_paused(Canvas &) override { ++*draws; }
    };
    Timeline tl;
    int draws = 0;
    bool paused = true;
    auto *probe =
        tl.add_get(0, PausedProbe(draws), Timeline::Pin::UNPINNED, &paused);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(draws, 1);
    probe->cancel();
    tl.step(fake_canvas());
    HS_EXPECT_EQ(draws, 1);
    HS_EXPECT_EQ(tl.event_count(), 0);
  }
}

/**
 * @brief Verifies cancel() fires the animation's .then() as the event is
 * removed.
 * @details cancel() reaches Timeline's removal branch through done(), so the
 * post callback runs there. TransformerPool::spawn_pinned reclaims its pool slot
 * from exactly this path: a pinned animation is infinite or repeating, so
 * cancellation is its only route to the callback.
 */
inline void test_timeline_cancel_fires_post_callback() {
  Timeline tl;
  float v = -1.0f;
  int thens = 0;
  auto *h = tl.add_get(0,
                       Animation::Mutation(
                           v, [](float e) { return e; }, 8, math::ease_linear)
                           .then([&]() { thens++; }),
                       Timeline::Pin::UNPINNED);
  tl.step(fake_canvas());
  HS_EXPECT_EQ(thens, 0);
  HS_EXPECT_EQ(tl.event_count(), 1);

  h->cancel();
  tl.step(fake_canvas());
  HS_EXPECT_EQ(tl.event_count(), 0);
  HS_EXPECT_EQ(thens, 1);
}

/**
 * @brief Cancellation before the start frame removes paused and unpaused events.
 */
inline void test_timeline_cancel_before_start() {
  for (bool paused : {false, true}) {
    Timeline tl;
    float value = -1.0f;
    int callbacks = 0;
    auto *event = tl.add_get(
        100,
        Animation::Mutation(
            value, [](float e) { return e; }, 4, math::ease_linear, true)
            .then([&] { ++callbacks; }),
        Timeline::Pin::PINNED, &paused);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(tl.event_count(), 1);
    event->cancel();
    tl.step(fake_canvas());
    HS_EXPECT_EQ(tl.event_count(), 0);
    HS_EXPECT_EQ(value, -1.0f);
    HS_EXPECT_EQ(callbacks, 1);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(callbacks, 1);
  }
}

/** @brief Cancellation of a started, paused event fires its callback once. */
inline void test_timeline_cancel_while_paused_removes_event() {
  Timeline tl;
  bool paused = false;
  float v = -1.0f;
  int thens = 0;
  auto *h = tl.add_get(0,
                       Animation::Mutation(
                           v, [](float e) { return e; }, 4, math::ease_linear,
                           /*repeat=*/true)
                           .then([&]() { thens++; }),
                       Timeline::Pin::PINNED, &paused);
  tl.step(fake_canvas()); // starts while unpaused
  HS_EXPECT_EQ(tl.event_count(), 1);

  paused = true;
  tl.step(fake_canvas()); // paused hold keeps the event
  HS_EXPECT_EQ(tl.event_count(), 1);
  HS_EXPECT_EQ(thens, 0);

  h->cancel();
  tl.step(fake_canvas()); // canceled while paused: removed, .then() fires
  HS_EXPECT_EQ(tl.event_count(), 0);
  HS_EXPECT_EQ(thens, 1);

  // Removal is final: no further callbacks from later paused frames.
  for (int i = 0; i < 3; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_EQ(thens, 1);
}

/**
 * @brief Verifies step() compacts the event array when a non-repeating event is
 * removed, and relocated survivors keep stepping from their new positions.
 * @details Later survivors are relocated (move_into) into the freed slots. The
 * decisive check is that the originally-LAST event (relocated furthest) still
 * reaches its own target.
 */
inline void test_timeline_compaction_preserves_later_events() {
  Timeline tl;
  float a = 0.0f, b = 0.0f, c = 0.0f;
  tl.add(0, Animation::Transition(a, 10.0f, 1,
                                  math::ease_linear)); // completes at t=1
  tl.add(0, Animation::Transition(b, 100.0f, 5,
                                  math::ease_linear)); // in-flight survivor
  tl.add(0, Animation::Transition(c, 200.0f, 5,
                                  math::ease_linear)); // in-flight survivor
  HS_EXPECT_EQ(tl.event_count(), 3);

  tl.step(fake_canvas()); // t=1: a done+removed; b,c step once and shift down
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_GT(b, 0.0f);
  HS_EXPECT_GT(c, 0.0f);
  HS_EXPECT_EQ(tl.event_count(), 2);
  float b_after1 = b, c_after1 = c;

  for (int i = 0; i < 4; ++i)
    tl.step(fake_canvas()); // t=2..5: finish the relocated survivors
  HS_EXPECT_GT(b, b_after1);
  HS_EXPECT_GT(c, c_after1);
  HS_EXPECT_NEAR(b, 100.0f, 1e-2f);
  HS_EXPECT_NEAR(c, 200.0f, 1e-2f);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/**
 * @brief Verifies .then() fires on completion and a callback may schedule a
 * follow-up on the same Timeline mid-step.
 * @details step() appends such events past the active snapshot, then gap-fills
 * them into the freed slots (the add-during-callback path). The follow-up is
 * added this frame but only runs on the next.
 */
inline void test_timeline_then_chains_follow_up_event() {
  Timeline tl;
  float a = 0.0f, b = 0.0f;
  // 'a' completes in one frame; its .then() schedules 'b' to start immediately.
  tl.add(0, Animation::Transition(a, 10.0f, 1, math::ease_linear).then([&]() {
    tl.add(0, Animation::Transition(b, 20.0f, 1, math::ease_linear));
  }));

  tl.step(fake_canvas()); // t=1: a completes -> callback adds b (gap-filled in)
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f); // b added this frame, not yet run
  HS_EXPECT_EQ(tl.event_count(), 1);

  tl.step(fake_canvas()); // t=2: the chained event runs and completes
  HS_EXPECT_NEAR(b, 20.0f, 1e-3f);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/**
 * @brief Verifies a repeating timer routes its hook through post_callback() so
 * an attached .then() fires on every trigger.
 * @details A repeating timer (duration=-1, self-resetting) never reaches
 * done(), so the Timeline never fires its per-cycle .then(); the timer must
 * fire the hook itself to match the then() contract.
 */
inline void test_repeating_timer_fires_then_each_cycle() {
  Timeline tl;
  int triggers = 0, thens = 0;
  tl.add(0, Animation::PeriodicTimer(
                3, [&](Canvas &) { triggers++; }, /*repeat=*/true)
                .then([&]() { thens++; }));
  for (int i = 0; i < 9; ++i)
    tl.step(fake_canvas()); // triggers at t=3,6,9
  HS_EXPECT_EQ(triggers, 3);
  HS_EXPECT_EQ(thens, 3);
}

/**
 * @brief Verifies a repeating timer canceled from inside its own callback fires
 * .then() exactly once and is removed.
 * @details The timer fires the per-cycle hook itself, and Timeline fires it
 * again from the removal branch; both run in the trigger frame unless the timer
 * reads repeats() (which drops on cancel) instead of the raw repeat flag.
 */
inline void test_repeating_timer_canceled_in_callback_fires_then_once() {
  Timeline tl;
  int thens = 0;
  struct {
    int triggers = 0;
    Animation::PeriodicTimer *timer = nullptr;
  } st; // one capture keeps the callback inside TimerFn's inplace budget
  st.timer = tl.add_get(0,
                        Animation::PeriodicTimer(
                            3,
                            [&st](Canvas &) {
                              st.triggers++;
                              st.timer->cancel();
                            },
                            /*repeat=*/true)
                            .then([&]() { thens++; }),
                        Timeline::Pin::UNPINNED);
  for (int i = 0; i < 9; ++i)
    tl.step(fake_canvas()); // would trigger at t=3,6,9 without the cancel
  HS_EXPECT_EQ(st.triggers, 1);
  HS_EXPECT_EQ(thens, 1);
  HS_EXPECT_EQ(global_timeline_num_events, 0);
}

/** @brief Verifies cancellation within a timer's completion callback does not re-enter it. */
inline void test_timer_then_self_cancellation_completes_once() {
  for (bool random : {false, true}) {
    Timeline tl;
    struct {
      int triggers = 0;
      int completions = 0;
      Animation::PeriodicTimer *periodic = nullptr;
      Animation::RandomTimer *random = nullptr;
    } state;
    auto trigger = [&state](Canvas &) { ++state.triggers; };
    auto complete = [&state]() {
      ++state.completions;
      if (state.periodic)
        state.periodic->cancel();
      else
        state.random->cancel();
    };
    if (random)
      state.random = tl.add_get(
          0,
          Animation::RandomTimer({.min = 1, .max = 1, .repeat = true}, trigger)
              .then(complete),
          Timeline::Pin::PINNED);
    else
      state.periodic = tl.add_get(
          0, Animation::PeriodicTimer(1, trigger, true).then(complete),
          Timeline::Pin::PINNED);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(state.triggers, 1);
    HS_EXPECT_EQ(state.completions, 1);
    HS_EXPECT_EQ(tl.event_count(), 0);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(state.completions, 1);
  }
}

/**
 * @brief Verifies clear() destroys all events and leaves the timeline reusable,
 * without rewinding the global frame cursor.
 * @details This is the in-place reset the singleton offers in lieu of
 * reassignment. The cursor is shared with every consumer deriving a phase from
 * frame(), so a runtime clear() must leave it running.
 */
inline void test_timeline_clear_destroys_events_keeping_frame() {
  Timeline tl;
  float a = 0.0f;
  tl.add(0, Animation::Transition(a, 10.0f, 5, math::ease_linear));
  tl.step(fake_canvas());
  HS_EXPECT_EQ(tl.event_count(), 1);
  HS_EXPECT_EQ(global_timeline_t, 1u);

  tl.clear();
  HS_EXPECT_EQ(tl.event_count(), 0);
  HS_EXPECT_EQ(global_timeline_t, 1u); // cursor keeps running

  // Reusable after clear().
  float b = 0.0f;
  tl.add(0, Animation::Transition(b, 5.0f, 1, math::ease_linear));
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(b, 5.0f, 1e-3f);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/**
 * @brief Verifies construction/destruction tears down a pinned event and
 * rewinds the frame cursor, where the public clear() would trap.
 * @details The pin guard covers the runtime API only; the instance boundary is
 * unguarded because no retained add_get() handle can outlive it. The trapping
 * half is death case "timeline_clear_pinned".
 */
inline void test_timeline_instance_boundary_reclaims_pinned_event() {
  {
    Timeline tl;
    tl.add_get(0,
               Animation::PeriodicTimer(
                   1, [](Canvas &) {}, /*repeat=*/true),
               Timeline::Pin::PINNED);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(tl.event_count(), 1);
    HS_EXPECT_EQ(global_timeline_t, 1u);
  }
  HS_EXPECT_EQ(global_timeline_num_events, 0);
  HS_EXPECT_EQ(global_timeline_t, 0u); // cursor rewound at the boundary

  // A fresh instance starts from a clean, unpinned table.
  Timeline tl;
  HS_EXPECT_EQ(tl.event_count(), 0);
  tl.clear(); // no pinned event survived, so the guard passes
}

/**
 * @brief Verifies add() is bounded by MAX_EVENTS (the shared global array
 * size): once full, a further add is rejected rather than overrunning, and the
 * once-per-episode drop log re-arms when the table drains.
 */
inline void test_timeline_full_guard_rejects_overflow() {
  Timeline tl;
  float sink = 0.0f;
  HS_EXPECT_EQ(tl.remaining(), Timeline::MAX_EVENTS);
  // Fill to capacity (these events share one float).
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  HS_EXPECT_EQ(tl.event_count(), Timeline::MAX_EVENTS);
  HS_EXPECT_EQ(tl.remaining(), 0);

  // The overflow event has its own target so a silent enqueue would show up.
  float rejected = 0.0f;
  const uint32_t dropped_before = Timeline::dropped_events();
  tl.add(0, Animation::Transition(rejected, 42.0f, 10,
                                  math::ease_linear)); // past full
  HS_EXPECT_EQ(tl.event_count(), Timeline::MAX_EVENTS);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 1);

  tl.step(fake_canvas());
  HS_EXPECT_NEAR(rejected, 0.0f, 1e-6f); // never ran
  HS_EXPECT_TRUE(global_timeline_drop_logged);

  // Draining the table ends the saturation episode: the next overflow logs
  // again rather than riding the first episode's flag.
  for (int i = 0; i < 64 && tl.event_count() > 0; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_EQ(tl.event_count(), 0);
  HS_EXPECT_FALSE(global_timeline_drop_logged);
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  tl.add(0, Animation::Transition(rejected, 42.0f, 10, math::ease_linear));
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 2);
  HS_EXPECT_TRUE(global_timeline_drop_logged);
}

/** @brief Clear hook that counts its invocations through its own ctx. */
inline void count_clear_hook(void *ctx) { ++*static_cast<int *>(ctx); }

/**
 * @brief Verifies remove_clear_hook() unregisters by ctx, leaves the surviving
 *        hooks registered, and ignores a ctx that was never added.
 * @details The removal backfills the hole with the last entry, so dropping the
 * first of several must not drop the one moved into its place — the pattern
 * TransformerPool's destructor relies on when several pools share a Timeline.
 */
inline void test_timeline_remove_clear_hook_unregisters_by_ctx() {
  Timeline tl;
  int first = 0;
  int second = 0;
  tl.add_clear_hook(&first, count_clear_hook);
  tl.add_clear_hook(&second, count_clear_hook);

  tl.clear();
  HS_EXPECT_EQ(first, 1);
  HS_EXPECT_EQ(second, 1);

  // Removing the first entry backfills it with the second.
  tl.remove_clear_hook(&first);
  tl.clear();
  HS_EXPECT_EQ(first, 1);
  HS_EXPECT_EQ(second, 2);

  // An unregistered ctx is a no-op, not a removal of whatever is left.
  int absent = 0;
  tl.remove_clear_hook(&absent);
  tl.clear();
  HS_EXPECT_EQ(absent, 0);
  HS_EXPECT_EQ(second, 3);

  tl.remove_clear_hook(&second);
  tl.clear();
  HS_EXPECT_EQ(first, 1);
  HS_EXPECT_EQ(second, 3);
}

/**
 * @brief Verifies Orientation::upsample SLERP-interpolates the recorded
 * sub-frames up to a target count (preserving endpoints) and collapse()
 * discards all but the newest.
 * @details These are the two primitives multi-animation motion blur is built
 * on.
 */
inline void test_orientation_upsample_then_collapse() {
  math::Orientation<8> o; // identity, 1 frame
  o.push(math::make_rotation(
      math::Z_AXIS, math::PI_F / 2)); // 2 frames: identity, +90 about Z
  HS_EXPECT_EQ(o.length(), 2);

  o.upsample(5);
  HS_EXPECT_EQ(o.length(), 5);

  // Endpoints preserved: frame 0 ~ identity (+X stays +X); frame 4 ~ the +90
  // rotation about Z (+X -> +Y).
  math::Vector f0 = o.orient(math::X_AXIS, 0);
  HS_EXPECT_NEAR(f0.x, 1.0f, 1e-3f);
  math::Vector f4 = o.orient(math::X_AXIS, 4);
  HS_EXPECT_NEAR(f4.y, 1.0f, 1e-3f);

  // SLERP is monotone: +X decreases across the interpolated frames.
  float prevx = 2.0f;
  for (int i = 0; i < 5; ++i) {
    float x = o.orient(math::X_AXIS, i).x;
    HS_EXPECT_LE(x, prevx + 1e-4f);
    prevx = x;
  }

  o.collapse();
  HS_EXPECT_EQ(o.length(), 1);
  math::Vector c = o.orient(math::X_AXIS, 0);
  HS_EXPECT_NEAR(c.y, 1.0f, 1e-3f);
}

/** @brief Internal-angle allowance for float drift after 600 cycles, radians. */
constexpr float MOTION_WARP_TOL = 1e-4f;

/**
 * @brief Verifies a repeating Motion does not drift across many cycles.
 * @details A repeating Motion advances its Orientation by relative deltas taken
 * between consecutive path frames, where each frame is a pure function of the
 * path parameter (point + tangent). Because the frame depends only on the phase,
 * the per-cycle product of deltas telescopes — there is no accumulating
 * quaternion chain to warp the traced curve. The decisive, precession-immune
 * signature is the set of rotation-INVARIANT internal angles between heads
 * sampled at fixed phases within a cycle: a rigid drift (holonomy) leaves them
 * unchanged, so any growth is genuine warp. A late cycle is compared against the
 * ideal Lissajous internal angles within accumulated float drift.
 */
inline void test_motion_repeating_does_not_drift() {
  using Ori = math::Orientation<16>;
  constexpr int duration = 40;
  ProceduralPath path;
  // lissajous(.,.,.,0) == +Y, so the identity-start orientation places the head
  // on the path at phase 0.
  path.f = [](float t) {
    return math::lissajous(1.06f, 1.06f, 0.0f, t * 5.909f);
  };

  Ori o; // identity, single frame; orient(+Y) starts on the path
  const math::Vector node_v = math::Y_AXIS;

  Timeline tl;
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));

  math::Vector late_heads[duration + 1]; // indexed by phase 1..duration
  const int cycles = 600;
  for (int c = 0; c < cycles; ++c) {
    for (int fr = 1; fr <= duration; ++fr) {
      tl.step(fake_canvas());
      if (c == cycles - 1)
        late_heads[fr] = o.orient(node_v);
    }
  }

  // Each phase's internal angle to the phase-1 anchor must match the ideal
  // Lissajous internal angle; a rigid precession leaves these untouched, so any
  // growth is genuine warp. Interior phases only (the boundary frame rewinds).
  const int anchor = 1;
  const math::Vector ideal_anchor = path.f((float)anchor / duration);
  float max_error = 0.0f;
  for (int fr = 2; fr < duration; ++fr) {
    const math::Vector ideal_fr = path.f((float)fr / duration);
    const float ERROR =
        fabsf(math::angle_between(late_heads[anchor], late_heads[fr]) -
              math::angle_between(ideal_anchor, ideal_fr));
    max_error = hs_test::fold_worst(max_error, ERROR);
  }
  std::printf("  late-cycle motion angle drift: %.9g rad\n", max_error);
  HS_EXPECT_LT(max_error, MOTION_WARP_TOL);
}

/** @brief Verifies reanchor preserves orientation continuity when the live path changes. */
inline void test_motion_reanchor_after_path_swap() {
  const auto step_after_swap = [](bool reanchor) {
    ProceduralPath path;
    path.f = [](float t) {
      const float ANGLE = math::TWO_PI_F * t;
      return math::Vector(cosf(ANGLE), sinf(ANGLE), 0);
    };
    math::Orientation<16> orientation;
    Animation::Motion<288, 16> motion(orientation, path, 1000, true);
    for (int frame = 0; frame < 5; ++frame)
      motion.step(fake_canvas());
    const math::Vector BEFORE = orientation.orient(math::X_AXIS);
    path.f = [](float t) {
      const float ANGLE = math::TWO_PI_F * t + math::PI_F * 0.5f;
      return math::Vector(cosf(ANGLE), sinf(ANGLE), 0);
    };
    if (reanchor)
      motion.reanchor();
    motion.step(fake_canvas());
    return math::angle_between(BEFORE, orientation.orient(math::X_AXIS));
  };
  HS_EXPECT_LT(step_after_swap(true), 0.02f);
  HS_EXPECT_GT(step_after_swap(false), 1.0f);
}

/**
 * @brief A co-driver sharing a repeating Motion's Orientation survives the
 * repeat seam.
 * @details Motion re-seats via a relative delta; the co-driver's accumulated
 * rotation persists across the seam. With a CLOSED path Motion's per-cycle
 * contribution telescopes to identity, so the only thing that should move the
 * shared orientation at a seam is the co-driver's own small step — never a
 * large snap-back. Assert the probe's per-frame angular step stays bounded
 * across many seams while its cumulative travel is large (so the co-driver is
 * provably active, not a no-op).
 */
inline void test_motion_codriven_survives_repeat_seam() {
  using Ori = math::Orientation<16>;
  const int duration = 30;
  ProceduralPath path;
  // A closed great circle: path(0) == path(1) with matching tangent, so Motion's
  // per-cycle delta product is identity.
  path.f = [](float t) {
    float a = 2.0f * math::PI_F * t;
    return math::Vector(std::cos(a), std::sin(a), 0.0f);
  };

  Ori o; // identity
  Timeline tl;
  // Repeating Motion + a repeating co-driver rotation about Y, both driving `o`.
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));
  tl.add(0, Animation::Rotation<288, 16>(o, math::Y_AXIS, 2.0f * math::PI_F,
                                         duration, math::ease_linear,
                                         /*repeat=*/true));

  const math::Vector probe = math::Z_AXIS;
  math::Vector prev = o.orient(probe);
  float max_step = 0.0f;
  float total_travel = 0.0f;
  const int cycles = 8;
  for (int c = 0; c < cycles; ++c) {
    for (int fr = 1; fr <= duration; ++fr) {
      tl.step(fake_canvas());
      math::Vector cur = o.orient(probe);
      float step = math::angle_between(prev, cur);
      max_step = std::max(max_step, step);
      total_travel += step;
      prev = cur;
    }
  }

  // No single frame (seams included) snaps the orientation; 0.8 clears the
  // largest legitimate one-frame step but is far below a multi-radian snap.
  HS_EXPECT_LT(max_step, 0.8f);
  HS_EXPECT_GT(total_travel, 5.0f);
}
