/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Timeline scheduling, cancellation and storage lifecycle.

/**
 * @brief Verifies two animations sharing one Orientation COMPOSE their
 * sub-frame motion-blur history within a frame instead of clobbering it.
 * @details With composition the oldest sub-frame (index 0) still reflects the
 * pre-frame orientation (identity here). Uses the global Timeline;
 * Rotation::step never dereferences the canvas.
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

/**
 * @brief Verifies event-level pause freezes delays, steps, and callbacks.
 * @details A one-frame event ahead of the pausable one completes first, so the
 * pause gate must survive the compaction that relocates the pausable event.
 */
inline void test_timeline_pausable_event_uses_active_time() {
  Timeline tl;
  bool paused = true;
  float lead = 0.0f;
  float value = 0.0f;
  float ambient = 0.0f;
  int completions = 0;
  tl.add(0, Animation::Transition(lead, 1.0f, 1, math::ease_linear));
  tl.add_pausable(
      3, Animation::Transition(value, 9.0f, 3, math::ease_linear).then([&]() {
        ++completions;
      }),
      &paused);
  tl.add(0, Animation::Transition(ambient, 1.0f, 2, math::ease_linear));

  for (int i = 0; i < 5; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_NEAR(lead, 1.0f, 1e-6f);
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

  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 0.0f, 1e-6f);
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 1.0f, 1e-6f);
  HS_EXPECT_EQ(global_timeline_t, UINT32_MAX);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/**
 * @brief Verifies a repeating animation rewinds at the end of each cycle
 * instead of being removed, replaying the curve.
 * @details Mutation writes f(easing(t/duration)) each step, so after completion
 * and rewind the next step drops back to the mid-cycle value.
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
 * @details cancel() makes done() permanently true and repeats() false, so
 * Timeline routes the event through the removal branch.
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
 * @details A cancel taken in the repeating branch leaves the animation done()
 * and non-repeating; .then() must not fire again from the removal branch.
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
  tl.step(fake_canvas());
  tl.step(fake_canvas()); // first cycle completes; the callback cancels
  HS_EXPECT_EQ(st.thens, 1);
  HS_EXPECT_EQ(tl.event_count(), 0);
  for (int i = 0; i < 4; ++i)
    tl.step(fake_canvas()); // cycles complete at t=4, 6 without the cancel
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
 * @brief Verifies a paused event redraws on the first step at or after its
 * start frame; delays 0 and 1 both start on the next step().
 */
inline void test_timeline_paused_event_redraws_from_start_frame() {
  struct PausedProbe : Animation::AnimationBase<PausedProbe> {
    int *draws;
    explicit PausedProbe(int &count) : draws(&count) {}
    void step_paused(Canvas &) override { ++*draws; }
  };
  for (int delay : {0, 1}) {
    Timeline tl;
    int draws = 0;
    bool paused = true;
    tl.add_pausable(delay, PausedProbe(draws), &paused);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(draws, 1);
    tl.step(fake_canvas());
    HS_EXPECT_EQ(draws, 2);
  }
}

/**
 * @brief Verifies cancel() fires the animation's .then() as the event is
 * removed.
 * @details cancel() reaches Timeline's removal branch through done(), so the
 * post callback runs there.
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
 * @details Cancellation completes through Timeline without rescheduling or
 * consuming another random delay draw.
 */
inline void test_repeating_timer_canceled_in_callback_fires_then_once() {
  const auto SAVED_RNG = hs::random();
  hs::random().seed(1337);
  Timeline tl;
  int thens = 0;
  struct {
    int triggers = 0;
    Animation::RandomTimer *timer = nullptr;
  } st; // one capture keeps the callback inside TimerFn's inplace budget
  st.timer =
      tl.add_get(0,
                 Animation::RandomTimer({.min = 3, .max = 3, .repeat = true},
                                        [&st](Canvas &) {
                                          st.triggers++;
                                          st.timer->cancel();
                                        })
                     .then([&]() { thens++; }),
                 Timeline::Pin::UNPINNED);
  auto expected_rng = hs::random();
  for (int i = 0; i < 9; ++i)
    tl.step(fake_canvas()); // would trigger at t=3,6,9 without the cancel
  HS_EXPECT_EQ(st.triggers, 1);
  HS_EXPECT_EQ(thens, 1);
  HS_EXPECT_EQ(global_timeline_num_events, 0);
  HS_EXPECT_EQ(hs::random()(), expected_rng());
  hs::random() = SAVED_RNG;
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
 * @details The frame cursor is shared, so a runtime clear() leaves it running.
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
 * @details The pin guard covers the runtime API only. The trapping half is
 * death case "timeline_clear_pinned".
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
 * @brief Runs @p fn with fd 1 redirected into a pipe and counts the timeline
 * drop-log lines it writes.
 * @return The line count, or -1 if the capture could not be established.
 * @details The captured output must fit the pipe buffer or the write blocks.
 */
template <typename Fn> int count_timeline_drop_logs(Fn &&fn) {
  int fds[2] = {-1, -1};
  std::fflush(stdout);
  if (fd_pipe(fds) != 0)
    return -1;
  const int saved_out = fd_dup(1);
  if (saved_out < 0) {
    fd_close(fds[0]);
    fd_close(fds[1]);
    return -1;
  }
  fd_dup2(fds[1], 1);
  fn();
  std::fflush(stdout);
  fd_dup2(saved_out, 1);
  fd_close(saved_out);
  fd_close(fds[1]);
  char buf[1024];
  size_t n = 0;
  for (;;) {
    const long got = fd_read(fds[0], buf + n, sizeof(buf) - 1 - n);
    if (got <= 0)
      break;
    n += static_cast<size_t>(got);
  }
  buf[n] = '\0';
  fd_close(fds[0]);
  int lines = 0;
  for (const char *p = buf;
       (p = std::strstr(p, "Timeline full, failed to add animation!")); ++p)
    ++lines;
  return lines;
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
  HS_EXPECT_EQ(count_timeline_drop_logs([&] {
                 for (int i = 0; i < 2; ++i)
                   tl.add(0, Animation::Transition(rejected, 42.0f, 10,
                                                   math::ease_linear));
               }),
               1);
  HS_EXPECT_EQ(tl.event_count(), Timeline::MAX_EVENTS);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 2);

  tl.step(fake_canvas());
  HS_EXPECT_NEAR(rejected, 0.0f, 1e-6f); // never ran

  // Draining the table ends the saturation episode: the next overflow logs
  // again rather than riding the first episode's flag.
  for (int i = 0; i < 64 && tl.event_count() > 0; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_EQ(tl.event_count(), 0);
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  HS_EXPECT_EQ(count_timeline_drop_logs([&] {
                 tl.add(0, Animation::Transition(rejected, 42.0f, 10,
                                                 math::ease_linear));
               }),
               1);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 3);
}

/** @brief Clear hook that counts its invocations through its own ctx. */
inline void count_clear_hook(void *ctx) { ++*static_cast<int *>(ctx); }

/**
 * @brief Verifies remove_clear_hook() unregisters by ctx, leaves the surviving
 *        hooks registered, and ignores a ctx that was never added.
 * @details The removal backfills the hole with the last entry, so dropping the
 * first of several must not drop the one moved into its place.
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
