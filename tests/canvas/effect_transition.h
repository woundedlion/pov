/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Effect transition controller state and adapter outcomes.

struct TransitionAdapter : hs::EffectTransitionAdapter {
  hs::EffectTransitionStatus preflight_status = hs::EffectTransitionStatus::OK;
  hs::EffectTransitionStatus construct_status = hs::EffectTransitionStatus::OK;
  hs::EffectTransitionStatus handoff_status = hs::EffectTransitionStatus::OK;
  hs::EffectTransitionStatus prepare_status = hs::EffectTransitionStatus::OK;
  hs::EffectTransitionStatus restore_status = hs::EffectTransitionStatus::OK;
  hs::EffectTransitionStatus restore_prepare_status =
      hs::EffectTransitionStatus::OK;
  hs::EffectRestoreCapability capability =
      hs::EffectRestoreCapability::EXACT_STATE;
  std::vector<float> envelopes;
  bool fenced = false;
  bool destroyed = false;
  bool outgoing_published = false;
  bool incoming_published = false;
  bool restored_published = false;
  bool committed = false;
  std::string_view committed_effect_id;
  bool failsafe = false;
  hs::EffectTransitionStatus failsafe_reason = hs::EffectTransitionStatus::OK;
  int preflights = 0;
  int restores = 0;
  int discards = 0;
  hs::EffectHandoffState imported_handoff;

  hs::EffectTransitionStatus
  preflight(const hs::EffectTransitionRequest &,
            hs::EffectRestoreToken &outgoing) override {
    ++preflights;
    outgoing.effect_id = "outgoing";
    outgoing.capability = capability;
    return preflight_status;
  }
  void set_output_envelope(float value) override { envelopes.push_back(value); }
  bool presentation_complete() const override { return fenced; }
  void publish_outgoing_frame() override {
    HS_EXPECT_FALSE(envelopes.empty());
    if (envelopes.empty())
      return;
    HS_EXPECT_EQ(envelopes.back(), 0.0f);
    outgoing_published = true;
  }
  void destroy_outgoing() override {
    HS_EXPECT_TRUE(outgoing_published);
    destroyed = true;
  }
  hs::EffectTransitionStatus
  construct_incoming(const hs::EffectTransitionRequest &) override {
    return construct_status;
  }
  hs::EffectTransitionStatus
  import_handoff(const hs::EffectHandoffState &state) override {
    imported_handoff = state;
    return handoff_status;
  }
  hs::EffectTransitionStatus prepare_incoming_frame() override {
    return prepare_status;
  }
  void publish_incoming_frame() override { incoming_published = true; }
  void commit_identity(const hs::EffectTransitionRequest &request) override {
    committed = true;
    committed_effect_id = request.effect_id;
  }
  void destroy_incoming() override { ++discards; }
  hs::EffectTransitionStatus
  restore_outgoing(const hs::EffectRestoreToken &) override {
    ++restores;
    return restore_status;
  }
  hs::EffectTransitionStatus prepare_restore_frame() override {
    return restore_prepare_status;
  }
  void publish_restore_frame() override { restored_published = true; }
  void enter_clear_failsafe(hs::EffectTransitionStatus reason) override {
    failsafe_reason = reason;
    failsafe = true;
  }
};

// ============================================================================
// Construction
// ============================================================================

/**
 * @brief Verifies a freshly constructed Effect reports its dimensions, starts
 * fully black with a free display buffer, and has its clip set to the whole
 * canvas.
 */
inline void test_construction_dims_and_clear() {
  TestEffect fx(96, 20);
  HS_EXPECT_EQ(fx.width(), 96);
  HS_EXPECT_EQ(fx.height(), 20);
  HS_EXPECT_PIXEL(fx.get_pixel(10, 5), 0, 0, 0);
  HS_EXPECT_TRUE(fx.buffer_free());
  HS_EXPECT_TRUE(fx.clip().is_full());
  HS_EXPECT_EQ(fx.clip().y_end, 20);
  HS_EXPECT_EQ(fx.clip().x_end, 96);
}

/**
 * @brief Verifies the smallest and largest supported dimensions are accepted.
 */
inline void test_construction_dimension_boundaries() {
  {
    TestEffect fx(1, 1);
    HS_EXPECT_EQ(fx.width(), 1);
    HS_EXPECT_EQ(fx.height(), 1);
    HS_EXPECT_PIXEL(fx.get_pixel(0, 0), 0, 0, 0);
  }
  {
    TestEffect fx(MAX_W, MAX_H);
    HS_EXPECT_EQ(fx.width(), MAX_W);
    HS_EXPECT_EQ(fx.height(), MAX_H);
    HS_EXPECT_PIXEL(fx.get_pixel(MAX_W - 1, MAX_H - 1), 0, 0, 0);
  }
}

inline void test_output_envelope_endpoints() {
  TestEffect fx(8, 8);
  const Pixel source{65535, 32768, 1};
  fx.set_output_envelope(0.0f);
  HS_EXPECT_PIXEL(fx.apply_output_envelope(source), 0, 0, 0);
  fx.set_output_envelope(1.0f);
  HS_EXPECT_PIXEL(fx.apply_output_envelope(source), 65535, 32768, 1);
  fx.set_output_envelope(0.5f);
  HS_EXPECT_EQ(fx.output_envelope_u16(), uint16_t(32768));
  HS_EXPECT_PIXEL(fx.apply_output_envelope(source), 32768, 16384, 1);
}

inline void test_effect_transition_fenced_commit() {
  TransitionAdapter adapter;
  hs::EffectTransitionController controller(adapter);
  const hs::EffectTransitionRequest request{
      "kaleidoscope-smooth", "coupled-grid", hs::EffectTransitionOrigin::MANUAL,
      2};
  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::OK);
  controller.tick();
  const size_t ENVELOPES_BEFORE_REFUSAL = adapter.envelopes.size();
  adapter.preflight_status = hs::EffectTransitionStatus::RESOURCE_REJECTED;
  const hs::EffectTransitionRequest REFUSED{
      "alien-brain", "alien-brain", hs::EffectTransitionOrigin::MANUAL, 1};
  HS_EXPECT_EQ(controller.request(REFUSED),
               hs::EffectTransitionStatus::RESOURCE_REJECTED);
  HS_EXPECT_EQ(adapter.envelopes.size(), ENVELOPES_BEFORE_REFUSAL);
  adapter.preflight_status = hs::EffectTransitionStatus::OK;
  controller.tick();
  controller.tick();
  HS_EXPECT_GE(adapter.envelopes.size(), size_t(4));
  if (adapter.envelopes.size() < 4)
    return;
  HS_EXPECT_EQ(adapter.envelopes[1], 1.0f);
  HS_EXPECT_EQ(adapter.envelopes[2], 0.5f);
  HS_EXPECT_EQ(adapter.envelopes[3], 0.0f);
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::CLEAR_PRESENTED);
  controller.tick();
  HS_EXPECT_TRUE(adapter.outgoing_published);
  HS_EXPECT_FALSE(adapter.destroyed);
  adapter.fenced = true;
  controller.tick();
  HS_EXPECT_TRUE(adapter.destroyed);
  controller.tick();
  controller.tick();
  controller.tick();
  HS_EXPECT_TRUE(adapter.incoming_published);
  HS_EXPECT_FALSE(adapter.committed);
  controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::COMMIT_READY);
  controller.tick();
  HS_EXPECT_TRUE(adapter.committed);
  HS_EXPECT_EQ(adapter.committed_effect_id,
               std::string_view(request.effect_id));
  controller.tick();
  controller.tick();
  controller.tick();
  HS_EXPECT_GE(adapter.envelopes.size(), size_t(3));
  if (adapter.envelopes.size() < 3)
    return;
  HS_EXPECT_EQ(adapter.envelopes[adapter.envelopes.size() - 3], 0.0f);
  HS_EXPECT_EQ(adapter.envelopes[adapter.envelopes.size() - 2], 0.5f);
  HS_EXPECT_EQ(adapter.envelopes.back(), 1.0f);
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::STEADY_IN);
  const size_t SECOND_START = adapter.envelopes.size();
  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::OK);
  HS_EXPECT_EQ(adapter.envelopes.back(), 1.0f);
  controller.tick();
  controller.tick();
  controller.tick();
  HS_EXPECT_EQ(adapter.envelopes.size(), SECOND_START + 4);
  HS_EXPECT_EQ(adapter.envelopes[SECOND_START + 1], 1.0f);
  HS_EXPECT_EQ(adapter.envelopes[SECOND_START + 2], 0.5f);
  HS_EXPECT_EQ(adapter.envelopes[SECOND_START + 3], 0.0f);
}

/** @brief Hidden and restored frames remain fenced and receive the latest handoff. */
inline void test_effect_transition_fences_hold() {
  for (bool restoring : {false, true}) {
    for (bool replacement : {false, true}) {
      TransitionAdapter adapter;
      adapter.fenced = true;
      if (restoring)
        adapter.prepare_status =
            hs::EffectTransitionStatus::FIRST_FRAME_REJECTED;
      hs::EffectTransitionController controller(adapter);
      const hs::EffectTransitionRequest REQUEST{
          "alien-brain", "alien-brain", hs::EffectTransitionOrigin::MANUAL, 1};
      const hs::EffectHandoffState FIRST{2, 1.25f, 2.5f, 17};
      const hs::EffectHandoffState SECOND{3, 3.25f, 4.5f, 29};
      HS_EXPECT_EQ(controller.request(REQUEST, FIRST),
                   hs::EffectTransitionStatus::OK);
      if (replacement) {
        controller.tick();
        HS_EXPECT_EQ(controller.request(REQUEST, SECOND),
                     hs::EffectTransitionStatus::OK);
      }
      const auto TARGET =
          restoring ? hs::EffectTransitionState::RESTORE_FRAME_READY
                    : hs::EffectTransitionState::HIDDEN_FRAME_PRESENTED;
      for (int step = 0; step < 20 && controller.current_state() != TARGET;
           ++step)
        controller.tick();
      HS_EXPECT_EQ(controller.current_state(), TARGET);
      const auto &expected = replacement ? SECOND : FIRST;
      HS_EXPECT_EQ(adapter.imported_handoff.schema_version,
                   expected.schema_version);
      HS_EXPECT_EQ(adapter.imported_handoff.projection_clock,
                   expected.projection_clock);
      HS_EXPECT_EQ(adapter.imported_handoff.palette_clock,
                   expected.palette_clock);
      HS_EXPECT_EQ(adapter.imported_handoff.choreography_position,
                   expected.choreography_position);
      adapter.fenced = false;
      const size_t ENVELOPES = adapter.envelopes.size();
      for (int step = 0; step < 4; ++step) {
        controller.tick();
        HS_EXPECT_EQ(controller.current_state(), TARGET);
        HS_EXPECT_FALSE(adapter.committed);
        HS_EXPECT_EQ(adapter.envelopes.size(), ENVELOPES);
      }
      adapter.fenced = true;
      controller.tick();
      HS_EXPECT_EQ(controller.current_state(),
                   restoring ? hs::EffectTransitionState::FADING_BACK
                             : hs::EffectTransitionState::COMMIT_READY);
    }
  }
}

inline void test_effect_transition_replacement_pause_and_restore() {
  TransitionAdapter adapter;
  hs::EffectTransitionController controller(adapter);
  hs::EffectTransitionRequest first{"alien-brain", "alien-brain",
                                    hs::EffectTransitionOrigin::AUTOMATIC, 2};
  HS_EXPECT_EQ(controller.request(first), hs::EffectTransitionStatus::OK);
  controller.tick();
  const size_t pre_replacement_envelopes = adapter.envelopes.size();
  first.effect_id = "kaleidoscope-hex-soft";
  HS_EXPECT_EQ(controller.request(first), hs::EffectTransitionStatus::OK);
  // Replacing the destination mid-fade writes no envelope of its own: the fade
  // carries on down from where it stood instead of snapping back to full.
  HS_EXPECT_EQ(adapter.envelopes.size(), pre_replacement_envelopes);
  HS_EXPECT_EQ(adapter.preflights, 2);
  controller.set_paused(true);
  const size_t envelope_count = adapter.envelopes.size();
  controller.tick();
  HS_EXPECT_EQ(adapter.envelopes.size(), envelope_count);
  controller.set_paused(false);
  controller.tick();
  HS_EXPECT_FALSE(adapter.envelopes.empty());
  if (adapter.envelopes.empty())
    return;
  HS_EXPECT_EQ(adapter.envelopes.back(), 0.5f);
  controller.tick();
  controller.tick();
  adapter.fenced = true;
  controller.tick();
  adapter.construct_status = hs::EffectTransitionStatus::RESOURCE_REJECTED;
  controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::RESTORING_OUT);
  controller.tick();
  HS_EXPECT_TRUE(adapter.restored_published);
  controller.tick();
  controller.tick();
  controller.tick();
  controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::STEADY_OUT);

  TransitionAdapter none;
  none.capability = hs::EffectRestoreCapability::NONE;
  hs::EffectTransitionController failing(none);
  HS_EXPECT_EQ(failing.request(first), hs::EffectTransitionStatus::OK);
  failing.tick();
  failing.tick();
  failing.tick();
  none.fenced = true;
  failing.tick();
  none.construct_status = hs::EffectTransitionStatus::RESOURCE_REJECTED;
  failing.tick();
  failing.tick();
  HS_EXPECT_TRUE(none.failsafe);
  HS_EXPECT_EQ(failing.current_state(),
               hs::EffectTransitionState::CLEAR_FAILSAFE);
  HS_EXPECT_FALSE(none.envelopes.empty());
  if (none.envelopes.empty())
    return;
  HS_EXPECT_EQ(none.envelopes.back(), 0.0f);

  none.construct_status = hs::EffectTransitionStatus::OK;
  HS_EXPECT_EQ(failing.request(first), hs::EffectTransitionStatus::OK);
  HS_EXPECT_EQ(failing.current_state(),
               hs::EffectTransitionState::CONSTRUCTING);
  HS_EXPECT_FALSE(none.envelopes.empty());
  if (none.envelopes.empty())
    return;
  HS_EXPECT_EQ(none.envelopes.back(), 0.0f);
  for (int i = 0; i < 8; ++i)
    failing.tick();
  HS_EXPECT_TRUE(none.committed);
  HS_EXPECT_EQ(failing.current_state(), hs::EffectTransitionState::STEADY_IN);
  HS_EXPECT_FALSE(none.envelopes.empty());
  if (none.envelopes.empty())
    return;
  HS_EXPECT_EQ(none.envelopes.back(), 1.0f);
}

inline void test_effect_transition_shorter_replacement_preserves_progress() {
  TransitionAdapter adapter;
  hs::EffectTransitionController controller(adapter);
  const hs::EffectTransitionRequest initial{
      "alien-brain", "alien-brain", hs::EffectTransitionOrigin::AUTOMATIC, 60};
  HS_EXPECT_EQ(controller.request(initial), hs::EffectTransitionStatus::OK);
  for (int i = 0; i < 55; ++i)
    controller.tick();
  HS_EXPECT_FALSE(adapter.envelopes.empty());
  if (adapter.envelopes.empty())
    return;
  HS_EXPECT_NEAR(adapter.envelopes.back(), 0.1f, 1e-6f);

  const hs::EffectTransitionRequest replacement{
      "kaleidoscope-hex-soft", "kaleidoscope-hex-soft",
      hs::EffectTransitionOrigin::AUTOMATIC, 1};
  HS_EXPECT_EQ(controller.request(replacement), hs::EffectTransitionStatus::OK);
  controller.tick();
  HS_EXPECT_FALSE(adapter.envelopes.empty());
  if (adapter.envelopes.empty())
    return;
  HS_EXPECT_EQ(adapter.envelopes.back(), 0.0f);
}

inline void test_effect_transition_failsafe_retry() {
  TransitionAdapter adapter;
  adapter.capability = hs::EffectRestoreCapability::NONE;
  adapter.fenced = true;
  adapter.construct_status = hs::EffectTransitionStatus::RESOURCE_REJECTED;
  hs::EffectTransitionController controller(adapter);
  const hs::EffectTransitionRequest request{
      "alien-brain", "alien-brain", hs::EffectTransitionOrigin::AUTOMATIC, 1};
  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::OK);
  for (int i = 0; i < 5; ++i)
    controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::CLEAR_FAILSAFE);
  HS_EXPECT_EQ(controller.failure(),
               hs::EffectTransitionStatus::RESTORE_REJECTED);
  HS_EXPECT_EQ(adapter.restores, 0);
  HS_EXPECT_EQ(adapter.discards, 1);

  // a preflight advertising a restorable outgoing cannot resurrect one
  adapter.capability = hs::EffectRestoreCapability::EXACT_STATE;
  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::OK);
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::CONSTRUCTING);
  controller.tick();
  controller.tick();
  HS_EXPECT_EQ(adapter.restores, 0);
  HS_EXPECT_EQ(adapter.discards, 2);
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::CLEAR_FAILSAFE);

  adapter.construct_status = hs::EffectTransitionStatus::OK;
  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::OK);
  for (int i = 0; i < 7; ++i)
    controller.tick();
  HS_EXPECT_TRUE(adapter.committed);
  HS_EXPECT_TRUE(adapter.incoming_published);
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::STEADY_IN);
  HS_EXPECT_FALSE(adapter.envelopes.empty());
  if (adapter.envelopes.empty())
    return;
  HS_EXPECT_EQ(adapter.envelopes.back(), 1.0f);
}

inline void test_effect_transition_refusal_branches() {
  TransitionAdapter adapter;
  hs::EffectTransitionController controller(adapter);
  const hs::EffectTransitionRequest request{
      "alien-brain", "alien-brain", hs::EffectTransitionOrigin::AUTOMATIC, 1};

  hs::EffectTransitionRequest malformed = request;
  malformed.effect_id = "";
  HS_EXPECT_EQ(controller.request(malformed),
               hs::EffectTransitionStatus::UNAVAILABLE);
  malformed = request;
  malformed.fade_ticks = 0;
  HS_EXPECT_EQ(controller.request(malformed),
               hs::EffectTransitionStatus::UNAVAILABLE);
  HS_EXPECT_EQ(adapter.preflights, 0);

  adapter.preflight_status = hs::EffectTransitionStatus::INVALID_RESTORE;
  HS_EXPECT_EQ(controller.request(request),
               hs::EffectTransitionStatus::INVALID_RESTORE);
  HS_EXPECT_EQ(adapter.preflights, 1);
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::STEADY_OUT);
  HS_EXPECT_TRUE(adapter.envelopes.empty());
  adapter.preflight_status = hs::EffectTransitionStatus::OK;

  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::OK);
  controller.tick();
  controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::CLEAR_PRESENTED);
  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::BUSY);
  HS_EXPECT_EQ(adapter.preflights, 2);

  adapter.fenced = true;
  adapter.handoff_status = hs::EffectTransitionStatus::INVALID_HANDOFF;
  controller.tick();
  controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::PREPARING_FIRST_FRAME);
  controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::RESTORING_OUT);
  HS_EXPECT_EQ(controller.failure(),
               hs::EffectTransitionStatus::INVALID_HANDOFF);
  HS_EXPECT_EQ(adapter.discards, 1);
  HS_EXPECT_FALSE(adapter.incoming_published);
  controller.tick();
  HS_EXPECT_EQ(adapter.restores, 1);
  HS_EXPECT_TRUE(adapter.restored_published);
  controller.tick();
  controller.tick();
  controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::STEADY_OUT);
  HS_EXPECT_FALSE(adapter.envelopes.empty());
  if (adapter.envelopes.empty())
    return;
  HS_EXPECT_EQ(adapter.envelopes.back(), 1.0f);
  HS_EXPECT_FALSE(adapter.committed);
  HS_EXPECT_EQ(controller.failure(),
               hs::EffectTransitionStatus::INVALID_HANDOFF);

  adapter.handoff_status = hs::EffectTransitionStatus::OK;
  adapter.prepare_status = hs::EffectTransitionStatus::FIRST_FRAME_REJECTED;
  adapter.restore_status = hs::EffectTransitionStatus::INVALID_RESTORE;
  adapter.restored_published = false;
  HS_EXPECT_EQ(controller.request(request), hs::EffectTransitionStatus::OK);
  HS_EXPECT_EQ(controller.failure(), hs::EffectTransitionStatus::OK);
  for (int i = 0; i < 5; ++i)
    controller.tick();
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::RESTORING_OUT);
  HS_EXPECT_EQ(controller.failure(),
               hs::EffectTransitionStatus::FIRST_FRAME_REJECTED);
  HS_EXPECT_EQ(adapter.discards, 2);
  controller.tick();
  HS_EXPECT_EQ(adapter.restores, 2);
  HS_EXPECT_FALSE(adapter.restored_published);
  HS_EXPECT_TRUE(adapter.failsafe);
  HS_EXPECT_EQ(adapter.failsafe_reason, controller.failure());
  HS_EXPECT_EQ(controller.current_state(),
               hs::EffectTransitionState::CLEAR_FAILSAFE);
  HS_EXPECT_EQ(controller.failure(),
               hs::EffectTransitionStatus::INVALID_RESTORE);
  HS_EXPECT_FALSE(adapter.envelopes.empty());
  if (adapter.envelopes.empty())
    return;
  HS_EXPECT_EQ(adapter.envelopes.back(), 0.0f);

  TransitionAdapter frame;
  frame.fenced = true;
  frame.prepare_status = hs::EffectTransitionStatus::FIRST_FRAME_REJECTED;
  frame.restore_prepare_status = hs::EffectTransitionStatus::RESTORE_REJECTED;
  hs::EffectTransitionController second(frame);
  HS_EXPECT_EQ(second.request(request), hs::EffectTransitionStatus::OK);
  for (int i = 0; i < 6; ++i)
    second.tick();
  HS_EXPECT_EQ(frame.discards, 1);
  HS_EXPECT_EQ(frame.restores, 1);
  HS_EXPECT_FALSE(frame.restored_published);
  HS_EXPECT_TRUE(frame.failsafe);
  HS_EXPECT_EQ(frame.failsafe_reason, second.failure());
  HS_EXPECT_EQ(second.current_state(),
               hs::EffectTransitionState::CLEAR_FAILSAFE);
  HS_EXPECT_EQ(second.failure(), hs::EffectTransitionStatus::RESTORE_REJECTED);
  HS_EXPECT_FALSE(frame.envelopes.empty());
  if (frame.envelopes.empty())
    return;
  HS_EXPECT_EQ(frame.envelopes.back(), 0.0f);
}
