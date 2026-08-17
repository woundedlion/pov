/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file effect_transition.h
 * @brief Fenced effect-to-effect transition: the state machine
 *        (EffectTransitionController) and the host-side operations it drives
 *        through EffectTransitionAdapter.
 */

#include <cstdint>
#include <string_view>

namespace hs {

/**
 * @brief Where the controller stands in a transition.
 * @details Every edge below is taken by one EffectTransitionController::tick(),
 *   except the two request() edges, which are noted as such. STEADY_OUT and
 *   STEADY_IN are the resting states; the controller treats them alike, and
 *   they differ only in naming which side of a completed transition is showing.
 *   CLEAR_FAILSAFE is terminal for tick(); only request() leaves it.
 *
 *   Forward path:   STEADY_OUT/STEADY_IN --request()--> FADING_OUT ->
 *   CLEAR_PRESENTED -> CONSTRUCTING -> PREPARING_FIRST_FRAME ->
 *   PUBLISHING_HIDDEN_FRAME -> HIDDEN_FRAME_PRESENTED -> COMMIT_READY ->
 *   FADING_IN -> STEADY_IN.
 *
 *   Rollback: any adapter failure in CONSTRUCTING or PREPARING_FIRST_FRAME
 *   diverts to RESTORING_OUT -> RESTORE_FRAME_READY -> FADING_BACK ->
 *   STEADY_OUT. A rollback that itself fails, or one with no restorable
 *   outgoing, goes to CLEAR_FAILSAFE, whose only exit is
 *   CLEAR_FAILSAFE --request()--> CONSTRUCTING (output is already dark and
 *   the outgoing is gone, so the fade-out and clear steps are skipped).
 */
enum class EffectTransitionState : uint8_t {
  STEADY_OUT, /**< Resting on the pre-transition effect. */
  FADING_OUT, /**< Envelope stepping 1 -> 0 over the request's fade_ticks. */
  CLEAR_PRESENTED, /**< Dark frame published; awaiting presentation_complete()
                        before destroy_outgoing(). */
  CONSTRUCTING,    /**< construct_incoming() runs this tick. */
  PREPARING_FIRST_FRAME, /**< import_handoff() then prepare_incoming_frame(). */
  PUBLISHING_HIDDEN_FRAME, /**< Envelope forced to 0, then
                                publish_incoming_frame(). */
  HIDDEN_FRAME_PRESENTED,  /**< Awaiting presentation_complete() on that hidden
                                frame. */
  COMMIT_READY,  /**< commit_identity() runs this tick; the swap is now
                     irreversible. */
  FADING_IN,     /**< Envelope stepping 0 -> 1 over the request's fade_ticks. */
  STEADY_IN,     /**< Resting on the effect this transition installed. */
  RESTORING_OUT, /**< restore_outgoing() then prepare_restore_frame(). */
  RESTORE_FRAME_READY, /**< Awaiting presentation_complete() on the restored
                            frame. */
  FADING_BACK, /**< Envelope stepping 0 -> 1 back onto the outgoing effect. */
  CLEAR_FAILSAFE, /**< Output dark, nothing installed; awaiting a request(). */
};

/**
 * @brief Who asked for a transition.
 * @details Carried verbatim in the request and the restore token for the
 *   adapter's own bookkeeping; the controller never branches on it.
 */
enum class EffectTransitionOrigin : uint8_t {
  MANUAL,
  AUTOMATIC,
  SYNCHRONIZED,
  RESTORE,
  AUTHORING,
};

/**
 * @brief How completely the outgoing effect can be brought back on a rollback.
 * @details The controller only tests for NONE, which sends RESTORING_OUT
 *   straight to the fail-safe; the other two differ only to the adapter.
 */
enum class EffectRestoreCapability : uint8_t {
  EXACT_STATE, /**< Restorable to the state the token captured. */
  RESET_ONLY,  /**< Restorable, but from its initial state. */
  NONE         /**< Not restorable; a rollback must fail safe instead. */
};

/**
 * @brief Outcome of a request() or of one adapter step.
 * @details The controller tests only OK versus not-OK; the specific failure is
 *   what request() returns to the caller and what
 *   EffectTransitionController::failure() reports afterwards. BUSY and
 *   UNAVAILABLE are raised by request() itself, the rest by the adapter.
 */
enum class EffectTransitionStatus : uint8_t {
  OK,          /**< Step accepted. */
  BUSY,        /**< A transition is already in flight past its fade-out. */
  UNAVAILABLE, /**< Request malformed: empty effect_id or zero fade_ticks. */
  INVALID_RESTORE,      /**< Restore token unusable. */
  INVALID_HANDOFF,      /**< Handoff state unusable. */
  RESOURCE_REJECTED,    /**< Incoming effect could not be constructed. */
  FIRST_FRAME_REJECTED, /**< Incoming effect could not render its first frame. */
  RESTORE_REJECTED,     /**< Rollback failed or had no restorable outgoing. */
};

/** @brief What a caller asks the controller to transition to. */
struct EffectTransitionRequest {
  std::string_view effect_id; /**< Destination effect; empty is rejected. */
  std::string_view preset_id; /**< Destination preset, for the adapter. */
  EffectTransitionOrigin origin = EffectTransitionOrigin::MANUAL;
  uint16_t fade_ticks = 1; /**< tick()s each of the fade-out and fade-in spans;
                                0 is rejected. */
};

/**
 * @brief Outgoing-effect snapshot preflight() fills in, for a possible rollback.
 * @details schema_version is for the adapter that reads the token back; the
 *   controller never inspects it.
 */
struct EffectRestoreToken {
  static constexpr uint16_t SCHEMA_VERSION = 1;
  uint16_t schema_version = SCHEMA_VERSION;
  std::string_view effect_id;
  uint64_t seed_identity = 0;
  uint32_t visit_position = 0;
  EffectTransitionOrigin origin = EffectTransitionOrigin::RESTORE;
  EffectRestoreCapability capability =
      EffectRestoreCapability::NONE; /**< The one field the controller reads. */
  bool animations_paused = false;
};

/**
 * @brief Continuity the incoming effect inherits from the outgoing one.
 * @details schema_version is for the importing adapter; the controller never
 *   inspects it.
 */
struct EffectHandoffState {
  static constexpr uint16_t SCHEMA_VERSION = 1;
  uint16_t schema_version = SCHEMA_VERSION;
  float projection_clock = 0.0f;
  float palette_clock = 0.0f;
  uint32_t choreography_position = 0;
};

/**
 * @brief Host-side operations EffectTransitionController drives.
 * @details One call per state edge, in the order EffectTransitionState
 *   documents. The status-returning methods are the failure points: returning
 *   anything but OK diverts the controller (to the caller for preflight(), to
 *   RESTORING_OUT for the construct/prepare steps, to CLEAR_FAILSAFE for the
 *   restore steps). Only OK-ness is tested, so the enumerator named on each
 *   method is a convention that makes failure() legible, not a checked one.
 */
class EffectTransitionAdapter {
public:
  virtual ~EffectTransitionAdapter() = default;
  /**
   * @brief Vets a request and snapshots the outgoing effect before anything
   *        is torn down.
   * @param request The destination being requested.
   * @param outgoing Restore token to fill in; its capability decides whether a
   *        later rollback is possible at all.
   * @return OK to accept; UNAVAILABLE, RESOURCE_REJECTED or BUSY to refuse,
   *         returned verbatim from request() with no state change.
   */
  virtual EffectTransitionStatus
  preflight(const EffectTransitionRequest &request,
            EffectRestoreToken &outgoing) = 0;
  /**
   * @brief Sets the output brightness envelope.
   * @param value Envelope in [0,1]; 0 is fully dark, 1 fully lit.
   */
  virtual void set_output_envelope(float value) = 0;
  /**
   * @brief Whether the last published frame has actually reached the display.
   * @return True once the frame is out; the fence the tear-down and commit
   *         steps wait on.
   */
  virtual bool presentation_complete() const = 0;
  /** @brief Destroys the outgoing effect. Reached only behind a dark frame. */
  virtual void destroy_outgoing() = 0;
  /**
   * @brief Builds the incoming effect.
   * @param request The destination this transition committed to.
   * @return OK, or RESOURCE_REJECTED to roll back.
   */
  virtual EffectTransitionStatus
  construct_incoming(const EffectTransitionRequest &request) = 0;
  /**
   * @brief Hands the outgoing effect's continuity state to the incoming one.
   * @param handoff Clocks and choreography position to adopt.
   * @return OK, or INVALID_HANDOFF to roll back.
   */
  virtual EffectTransitionStatus
  import_handoff(const EffectHandoffState &handoff) = 0;
  /**
   * @brief Renders the incoming effect's first frame, still hidden.
   * @return OK, or FIRST_FRAME_REJECTED to roll back.
   */
  virtual EffectTransitionStatus prepare_incoming_frame() = 0;
  /** @brief Publishes that first frame, with the envelope already at 0. */
  virtual void publish_incoming_frame() = 0;
  /**
   * @brief Adopts the destination as the current identity.
   * @param request The destination now installed.
   * @details The point of no return: past it the controller has no rollback.
   */
  virtual void commit_identity(const EffectTransitionRequest &request) = 0;
  /**
   * @brief Rebuilds the outgoing effect from its restore token.
   * @param token The snapshot preflight() took.
   * @return OK, or INVALID_RESTORE to fail safe.
   */
  virtual EffectTransitionStatus
  restore_outgoing(const EffectRestoreToken &token) = 0;
  /**
   * @brief Renders the restored effect's frame, still hidden.
   * @return OK, or RESTORE_REJECTED to fail safe.
   */
  virtual EffectTransitionStatus prepare_restore_frame() = 0;
  /** @brief Publishes that frame, with the envelope already at 0. */
  virtual void publish_restore_frame() = 0;
  /**
   * @brief Tears down to a dark, effect-less output.
   * @param reason The failure that forced it; also what failure() reports.
   * @details Last resort — the envelope is already 0 and neither the incoming
   *   nor the outgoing effect can be shown.
   */
  virtual void enter_clear_failsafe(EffectTransitionStatus reason) = 0;
};

/**
 * @brief Drives EffectTransitionState, calling one adapter operation per edge.
 * @details Holds no effect and renders nothing: request() arms a destination
 *   and tick() advances one edge per call. See EffectTransitionState for the
 *   graph and EffectTransitionAdapter for the per-edge contract.
 */
class EffectTransitionController {
public:
  /**
   * @brief Binds a controller to the adapter it drives.
   * @param adapter Host operations; retained by reference and must outlive the
   *        controller.
   */
  explicit EffectTransitionController(EffectTransitionAdapter &adapter)
      : adapter(adapter) {}

  /**
   * @brief Arms a destination, starting or re-aiming a transition.
   * @param next Destination effect, preset and fade length.
   * @param next_handoff Continuity state the incoming effect will import.
   * @return OK once armed; UNAVAILABLE for an empty effect_id or zero
   *         fade_ticks; BUSY from any state past FADING_OUT; otherwise
   *         whatever the adapter's preflight() refused with.
   * @details Accepted only from STEADY_OUT, STEADY_IN, FADING_OUT and
   *   CLEAR_FAILSAFE. Re-requesting during FADING_OUT re-aims the in-flight
   *   fade: the elapsed fraction is remapped onto the new fade_ticks so the
   *   envelope carries on down rather than snapping back to full. A rejected
   *   request changes nothing.
   */
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
    const bool refading = state == EffectTransitionState::FADING_OUT;
    // Fraction of an in-flight fade-out already spent, remapped onto the new
    // request's tick count so the envelope carries on down instead of jumping
    // back to full.
    const float spent = refading
                            ? static_cast<float>(evaluation) /
                                  static_cast<float>(destination.fade_ticks)
                            : 0.0f;
    EffectRestoreToken next_restore;
    const EffectTransitionStatus status = adapter.preflight(next, next_restore);
    if (status != EffectTransitionStatus::OK)
      return status;

    destination = next;
    handoff = next_handoff;
    restore = next_restore;
    evaluation =
        static_cast<uint16_t>(spent * static_cast<float>(next.fade_ticks));
    last_failure = EffectTransitionStatus::OK;
    if (cleared) {
      // post-teardown: output is already clear and no outgoing remains
      restore.capability = EffectRestoreCapability::NONE;
      state = EffectTransitionState::CONSTRUCTING;
      return EffectTransitionStatus::OK;
    }
    state = EffectTransitionState::FADING_OUT;
    if (!refading)
      adapter.set_output_envelope(1.0f);
    return EffectTransitionStatus::OK;
  }

  /**
   * @brief Freezes and unfreezes tick().
   * @param value True to make tick() a no-op, holding the current state and
   *        envelope; request() stays available.
   */
  void set_paused(bool value) { paused = value; }
  /**
   * @brief Whether tick() is frozen.
   * @return The current pause flag.
   */
  bool is_paused() const { return paused; }
  /**
   * @brief Where the transition stands.
   * @return The current state.
   */
  EffectTransitionState current_state() const { return state; }
  /**
   * @brief Why the last transition diverted.
   * @return The status that forced the most recent rollback or fail-safe; OK
   *         once a request is accepted, and never reset by tick().
   */
  EffectTransitionStatus failure() const { return last_failure; }

  /**
   * @brief Advances the state machine by one edge.
   * @details A no-op while paused and in the three states that wait on
   *   request() (STEADY_OUT, STEADY_IN, CLEAR_FAILSAFE). Otherwise it takes
   *   exactly one EffectTransitionState edge; a fade state consumes one of its
   *   fade_ticks per call, and a fence state re-polls until
   *   EffectTransitionAdapter::presentation_complete().
   */
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
