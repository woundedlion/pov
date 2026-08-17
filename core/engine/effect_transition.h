/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <string_view>

namespace hs {

enum class EffectTransitionState : uint8_t {
  STEADY_OUT,
  FADING_OUT,
  CLEAR_PRESENTED,
  CONSTRUCTING,
  PREPARING_FIRST_FRAME,
  PUBLISHING_HIDDEN_FRAME,
  HIDDEN_FRAME_PRESENTED,
  COMMIT_READY,
  FADING_IN,
  STEADY_IN,
  RESTORING_OUT,
  RESTORE_FRAME_READY,
  FADING_BACK,
  CLEAR_FAILSAFE,
};

enum class EffectTransitionOrigin : uint8_t {
  MANUAL,
  AUTOMATIC,
  SYNCHRONIZED,
  RESTORE,
  AUTHORING,
};

enum class EffectRestoreCapability : uint8_t { EXACT_STATE, RESET_ONLY, NONE };

enum class EffectTransitionStatus : uint8_t {
  OK,
  BUSY,
  UNAVAILABLE,
  INVALID_RESTORE,
  INVALID_HANDOFF,
  RESOURCE_REJECTED,
  FIRST_FRAME_REJECTED,
  RESTORE_REJECTED,
};

struct EffectTransitionRequest {
  std::string_view effect_id;
  std::string_view preset_id;
  EffectTransitionOrigin origin = EffectTransitionOrigin::MANUAL;
  uint16_t fade_ticks = 1;
};

struct EffectRestoreToken {
  static constexpr uint16_t SCHEMA_VERSION = 1;
  uint16_t schema_version = SCHEMA_VERSION;
  std::string_view effect_id;
  uint64_t seed_identity = 0;
  uint32_t visit_position = 0;
  EffectTransitionOrigin origin = EffectTransitionOrigin::RESTORE;
  EffectRestoreCapability capability = EffectRestoreCapability::NONE;
  bool animations_paused = false;
};

struct EffectHandoffState {
  static constexpr uint16_t SCHEMA_VERSION = 1;
  uint16_t schema_version = SCHEMA_VERSION;
  float projection_clock = 0.0f;
  float palette_clock = 0.0f;
  uint32_t choreography_position = 0;
};

class EffectTransitionAdapter {
public:
  virtual ~EffectTransitionAdapter() = default;
  virtual EffectTransitionStatus
  preflight(const EffectTransitionRequest &request,
            EffectRestoreToken &outgoing) = 0;
  virtual void set_output_envelope(float value) = 0;
  virtual bool presentation_complete() const = 0;
  virtual void destroy_outgoing() = 0;
  virtual EffectTransitionStatus
  construct_incoming(const EffectTransitionRequest &request) = 0;
  virtual EffectTransitionStatus
  import_handoff(const EffectHandoffState &handoff) = 0;
  virtual EffectTransitionStatus prepare_incoming_frame() = 0;
  virtual void publish_incoming_frame() = 0;
  virtual void commit_identity(const EffectTransitionRequest &request) = 0;
  virtual EffectTransitionStatus
  restore_outgoing(const EffectRestoreToken &token) = 0;
  virtual EffectTransitionStatus prepare_restore_frame() = 0;
  virtual void publish_restore_frame() = 0;
  virtual void enter_clear_failsafe(EffectTransitionStatus reason) = 0;
};

class EffectTransitionController {
public:
  explicit EffectTransitionController(EffectTransitionAdapter &adapter)
      : adapter(adapter) {}

  EffectTransitionStatus request(const EffectTransitionRequest &next,
                                 const EffectHandoffState &next_handoff = {}) {
    if (next.effect_id.empty() || next.fade_ticks == 0)
      return EffectTransitionStatus::UNAVAILABLE;
    if (state != EffectTransitionState::STEADY_OUT &&
        state != EffectTransitionState::STEADY_IN &&
        state != EffectTransitionState::FADING_OUT &&
        state != EffectTransitionState::CLEAR_FAILSAFE)
      return EffectTransitionStatus::BUSY;

    const bool cleared = state == EffectTransitionState::CLEAR_FAILSAFE;
    EffectRestoreToken next_restore;
    const EffectTransitionStatus status = adapter.preflight(next, next_restore);
    if (status != EffectTransitionStatus::OK)
      return status;

    destination = next;
    handoff = next_handoff;
    restore = next_restore;
    evaluation = 0;
    last_failure = EffectTransitionStatus::OK;
    if (cleared) {
      // post-teardown: output is already clear and no outgoing remains
      restore.capability = EffectRestoreCapability::NONE;
      state = EffectTransitionState::CONSTRUCTING;
      return EffectTransitionStatus::OK;
    }
    state = EffectTransitionState::FADING_OUT;
    adapter.set_output_envelope(1.0f);
    return EffectTransitionStatus::OK;
  }

  void set_paused(bool value) { paused = value; }
  bool is_paused() const { return paused; }
  EffectTransitionState current_state() const { return state; }
  EffectTransitionStatus failure() const { return last_failure; }

  void tick() {
    if (paused || state == EffectTransitionState::STEADY_OUT ||
        state == EffectTransitionState::STEADY_IN ||
        state == EffectTransitionState::CLEAR_FAILSAFE)
      return;

    switch (state) {
    case EffectTransitionState::FADING_OUT:
      fade_out();
      break;
    case EffectTransitionState::CLEAR_PRESENTED:
      if (adapter.presentation_complete()) {
        adapter.destroy_outgoing();
        state = EffectTransitionState::CONSTRUCTING;
      }
      break;
    case EffectTransitionState::CONSTRUCTING:
      advance_or_restore(adapter.construct_incoming(destination),
                         EffectTransitionState::PREPARING_FIRST_FRAME);
      break;
    case EffectTransitionState::PREPARING_FIRST_FRAME: {
      EffectTransitionStatus status = adapter.import_handoff(handoff);
      if (status == EffectTransitionStatus::OK)
        status = adapter.prepare_incoming_frame();
      advance_or_restore(status,
                         EffectTransitionState::PUBLISHING_HIDDEN_FRAME);
      break;
    }
    case EffectTransitionState::PUBLISHING_HIDDEN_FRAME:
      adapter.set_output_envelope(0.0f);
      adapter.publish_incoming_frame();
      state = EffectTransitionState::HIDDEN_FRAME_PRESENTED;
      break;
    case EffectTransitionState::HIDDEN_FRAME_PRESENTED:
      if (adapter.presentation_complete())
        state = EffectTransitionState::COMMIT_READY;
      break;
    case EffectTransitionState::COMMIT_READY:
      adapter.commit_identity(destination);
      evaluation = 0;
      state = EffectTransitionState::FADING_IN;
      break;
    case EffectTransitionState::FADING_IN:
      fade_in(EffectTransitionState::STEADY_IN);
      break;
    case EffectTransitionState::RESTORING_OUT:
      restore_step();
      break;
    case EffectTransitionState::RESTORE_FRAME_READY:
      if (adapter.presentation_complete()) {
        evaluation = 0;
        state = EffectTransitionState::FADING_BACK;
      }
      break;
    case EffectTransitionState::FADING_BACK:
      fade_in(EffectTransitionState::STEADY_OUT);
      break;
    default:
      break;
    }
  }

private:
  void fade_out() {
    const float progress = static_cast<float>(evaluation) /
                           static_cast<float>(destination.fade_ticks);
    adapter.set_output_envelope(
        evaluation == destination.fade_ticks ? 0.0f : 1.0f - progress);
    if (evaluation == destination.fade_ticks) {
      state = EffectTransitionState::CLEAR_PRESENTED;
      return;
    }
    ++evaluation;
  }

  void fade_in(EffectTransitionState complete) {
    const float progress = static_cast<float>(evaluation) /
                           static_cast<float>(destination.fade_ticks);
    adapter.set_output_envelope(
        evaluation == destination.fade_ticks ? 1.0f : progress);
    if (evaluation == destination.fade_ticks) {
      state = complete;
      return;
    }
    ++evaluation;
  }

  void advance_or_restore(EffectTransitionStatus status,
                          EffectTransitionState success) {
    if (status == EffectTransitionStatus::OK) {
      state = success;
      return;
    }
    last_failure = status;
    state = EffectTransitionState::RESTORING_OUT;
  }

  void restore_step() {
    if (restore.capability == EffectRestoreCapability::NONE) {
      fail_safe(EffectTransitionStatus::RESTORE_REJECTED);
      return;
    }
    EffectTransitionStatus status = adapter.restore_outgoing(restore);
    if (status == EffectTransitionStatus::OK)
      status = adapter.prepare_restore_frame();
    if (status != EffectTransitionStatus::OK) {
      fail_safe(status);
      return;
    }
    adapter.set_output_envelope(0.0f);
    adapter.publish_restore_frame();
    state = EffectTransitionState::RESTORE_FRAME_READY;
  }

  void fail_safe(EffectTransitionStatus status) {
    last_failure = status;
    adapter.set_output_envelope(0.0f);
    adapter.enter_clear_failsafe(status);
    state = EffectTransitionState::CLEAR_FAILSAFE;
  }

  EffectTransitionAdapter &adapter;
  EffectTransitionState state = EffectTransitionState::STEADY_OUT;
  EffectTransitionRequest destination{};
  EffectRestoreToken restore{};
  EffectHandoffState handoff{};
  EffectTransitionStatus last_failure = EffectTransitionStatus::OK;
  uint16_t evaluation = 0;
  bool paused = false;
};

} // namespace hs
