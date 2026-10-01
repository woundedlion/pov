/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file choreography.h
 * @brief ChoreographedEffect: preset choreography and schema-versioned parameter
 *        snapshots over one Params type.
 */

#include <cmath>
#include <cstdint>
#include <variant>

#include "animation/animation.h"
#include "control/presets.h"
#include "platform/platform.h"
#include "math/easing.h"
#include "render/canvas.h"

/**
 * @brief Effect base owning one parameter set and its preset choreography:
 *        automatic transitions by each preset's departure policy, manual preset
 *        snaps and schema-versioned parameter snapshots.
 * @details Curiously recurring: `Derived` supplies `PARAMETER_SCHEMA_VERSION`,
 * `PRESET_DWELL_FRAMES` and `valid_params(params)`, plus its presets as rows
 * carrying their parameters and departure policy (`PresetEntry`): a `PRESETS`
 * table or a static `preset(index)`, and/or `PRESET_IDS` naming them. A member
 * `preset_params(index)` may derive a preset's live parameters from its row,
 * for an effect that patches runtime state into each preset; the row still
 * names the departure. The effect calls `begin_choreography()` once from
 * init() and `step_choreography()` every frame; each preset dwells for
 * `PRESET_DWELL_FRAMES`, then departs to the next over its own policy's frames,
 * and a single preset compiles the countdown out. Transition hooks:
 * `blend_params(progress)` writes the interpolated parameters of a
 * Segue::Preset::Lerp departure, `set_preset_opacity(value)` receives a
 * Segue::Preset::Fade departure's opacity, and shadowing `adopt_params(target)`
 * / `transition_armed(target)` keeps state derived from the parameters
 * consistent across snaps and crossfade arming. `initial_params()` overrides
 * the first preset as the startup default. A `Derived` keeping its hooks
 * non-public befriends this base.
 * @tparam Derived The effect class deriving from this base.
 * @tparam ParamsT The effect's parameter-set type.
 */
template <typename Derived, typename ParamsT>
class ChoreographedEffect : public Effect {
public:
  using Params = ParamsT;

  /** @brief A parameter set tagged with the schema version that produced it. */
  struct ParameterSnapshot {
    /** `Derived::PARAMETER_SCHEMA_VERSION` at capture time. */
    uint32_t schema_version;
    Params params; /**< The captured parameter values. */
  };

  /**
   * @brief Captures the live parameters with the effect's schema version.
   * @return A snapshot restore_parameters() accepts while the effect's schema
   *         version is unchanged.
   */
  ParameterSnapshot serialize_parameters() const {
    return {Derived::PARAMETER_SCHEMA_VERSION, params};
  }

  /**
   * @brief Adopts a snapshot's parameters if it is admissible.
   * @details An effect bumps `PARAMETER_SCHEMA_VERSION` whenever its `Params`
   * layout changes, so a snapshot taken under a different layout is rejected
   * rather than reinterpreted. On success any in-flight preset transition is
   * cancelled and the preset dwell restarts. Cancelling an unadopted fade
   * restores its departing preset index; other cancellations retain the index.
   * On rejection nothing is touched.
   * @param snapshot Snapshot to adopt.
   * @return True when the snapshot's schema version matches and its parameters
   *         are valid, false otherwise.
   */
  bool restore_parameters(const ParameterSnapshot &snapshot) {
    if (snapshot.schema_version != Derived::PARAMETER_SCHEMA_VERSION ||
        !Derived::valid_params(snapshot.params))
      return false;
    end_transition();
    derived().adopt_params(snapshot.params);
    preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
    return true;
  }

  /** @brief Authored preset count, for the effect registry: the same
      count preset_row() indexes. */
  static consteval size_t authored_preset_count() { return preset_count_of(); }

  /** @brief The departure policy of the preset at @p index; a single preset
      with no row snaps. */
  static constexpr Segue::Preset::Departure preset_departure(size_t index) {
    if constexpr (requires { Derived::preset(index); })
      return Derived::preset(index).segue;
    else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS[index].segue;
    else
      return Segue::Preset::Snap{};
  }

  /** @brief Frames one pass through every preset spans: each preset's dwell
      and each departure but the last preset's. */
  static consteval size_t cycle_frames() {
    size_t frames = preset_count_of() * Derived::PRESET_DWELL_FRAMES;
    for (size_t index = 0; index + 1 < preset_count_of(); ++index)
      frames += Segue::Preset::frames(preset_departure(index));
    return frames;
  }

protected:
  /** @brief An automatic preset transition's endpoints and progress. */
  struct Transition {
    Params from{};               /**< Parameters the transition departs from. */
    Params to{};                 /**< Parameters it lands on. */
    bool active = false;         /**< False when no transition is in flight. */
    uint16_t elapsed_frames = 0; /**< Unpaused steps elapsed. */
    uint16_t frames = 0;         /**< Steps the departure spans. */
    bool fades = false;          /**< A Fade departure rather than a Lerp. */
    bool adopted = false;        /**< A fade has adopted its target. */
    uint16_t serial = 0;
    size_t from_index = 0; /**< Preset before an unadopted fade. */
  };

  HS_COLD_MEMBER ChoreographedEffect(int W, int H, EffectConfig cfg = {})
      : Effect(W, H, cfg) {}

  /** @brief Adopts a snap target; shadow to re-derive dependent state. */
  void adopt_params(const Params &target) { params = target; }

  /** @brief Notified when a crossfade arms; shadow to capture endpoint state
      the blend hook interpolates alongside the parameters. */
  void transition_armed(const Params &) {}

  /**
   * @brief Overrides the dwell remaining before the next automatic transition.
   * @param frames Frames to hold; 0 lets the next frame start a transition.
   */
  HS_COLD_MEMBER void hold_initial_preset(uint16_t frames) {
    preset_dwell_remaining = frames;
  }

  /**
   * @brief Adopts a preset through the departing preset's policy.
   * @details A manual or synchronized change snaps regardless of policy, since
   * a user driving the control expects the preset it names immediately. An
   * AUTOMATIC change follows the departing preset's policy: Segue::Preset::Snap adopts
   * immediately, Segue::Preset::Lerp arms a crossfade from the live
   * parameters, and Segue::Preset::Fade dims to black, adopts in the dark and
   * brightens. A transition the timeline has no slot for restarts the dwell,
   * so the next attempt is a dwell away rather than on the following frame.
   * @param change The requested preset move.
   * @return False if an automatic transition cannot be scheduled.
   */
  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    static_assert(departures_are_supported(),
                  "a preset's departure outlasts PRESET_DWELL_FRAMES, or its "
                  "effect lacks the departure's hook (blend_params for Lerp, "
                  "set_preset_opacity for Fade)");
    const Params target = preset_target(change.to);
    const Segue::Preset::Departure DEPARTURE = preset_departure(change.from);
    const uint16_t FRAMES = Segue::Preset::frames(DEPARTURE);
    if (change.origin == PresetChangeOrigin::AUTOMATIC && FRAMES > 0) {
      if (transition.active)
        return false;
      const auto *LERP = std::get_if<Segue::Preset::Lerp>(&DEPARTURE);
      const bool *paused = !LERP || LERP->pausable ? &anims_paused : nullptr;
      const uint16_t TRANSITION_SERIAL = transition.serial + 1;
      if (timeline.add_get(0,
                           Animation::Progress(
                               [this, TRANSITION_SERIAL](float progress) {
                                 if (transition.serial == TRANSITION_SERIAL)
                                   run_transition(progress);
                               },
                               FRAMES, LERP ? LERP->easing : math::ease_linear),
                           Timeline::Pin::UNPINNED, paused) == nullptr) {
        preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
        return false;
      }
      transition = {params,     target,          true,  0,
                    FRAMES,     LERP == nullptr, false, TRANSITION_SERIAL,
                    change.from};
      if (LERP)
        derived().transition_armed(target);
      return true;
    }
    end_transition(false);
    derived().adopt_params(target);
    preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
    return true;
  }

  /**
   * @brief Ends an in-flight transition when the user takes a parameter over.
   * @details A crossfade rewrites the whole parameter set every frame, so a
   * transition left running would overwrite the write that just landed; a fade
   * returns to full opacity. A manual edit restarts the preset dwell.
   */
  HS_COLD_MEMBER void parameter_written() override {
    end_transition();
    preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
  }

  /**
   * @brief Configures the preset controller from `PRESET_IDS` or the preset
   *        rows. Call once from init().
   */
  HS_COLD_MEMBER void begin_choreography() {
    configure_presets(preset_count_of());
  }

  /// Retires the preset dwell and starts the next automatic preset transition.
  /// @details Call every frame. Pause suppresses preset selection, so no new
  /// transition begins while paused, and a paused fade shows full opacity.
  HS_COLD_MEMBER void step_choreography() {
    if constexpr (preset_count_of() == 1)
      return;
    else {
      if (anims_paused) {
        if (transition.active && transition.fades)
          set_opacity(1.0f);
        return;
      }
      if (transition.active)
        return;
      if (preset_dwell_remaining > 0 && --preset_dwell_remaining > 0)
        return;
      if (advance_preset()) {
        if (!transition.active || !transition.fades)
          log_preset();
      } else {
        preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
      }
    }
  }

  /**
   * @brief Advances the in-flight transition one step.
   * @param progress Eased transition progress in [0, 1].
   * @details A transition cancelled by a manual preset, an edit or a snapshot
   * restore keeps stepping but writes nothing. A fade holds the departing
   * parameters at falling opacity, adopts the target at the dark midpoint,
   * and rises back to full.
   */
  HS_COLD_MEMBER void run_transition(float progress) {
    if (!transition.active)
      return;
    const bool COMPLETE = ++transition.elapsed_frames >= transition.frames;
    const float PROGRESS = COMPLETE ? 1.0f : progress;
    if (transition.fades) {
      if (PROGRESS >= 0.5f && !transition.adopted) {
        derived().adopt_params(transition.to);
        transition.adopted = true;
        log_preset();
      }
      set_opacity(fabsf(1.0f - 2.0f * PROGRESS));
    } else if constexpr (BLENDS) {
      derived().blend_params(PROGRESS);
    }
    if (COMPLETE) {
      transition.active = false;
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
    }
  }

  /** Live parameters; the registered sliders write straight into these. */
  Params params = initial_params_of();
  Transition transition;
  Timeline timeline;

private:
  Derived &derived() { return static_cast<Derived &>(*this); }

  /** @brief Whether the effect interpolates its parameters. */
  static constexpr bool BLENDS =
      requires(Derived &effect) { effect.blend_params(0.0f); };

  /** @brief Feeds a fade's opacity to an effect that takes one. */
  void set_opacity(float value) {
    if constexpr (requires { derived().set_preset_opacity(value); })
      derived().set_preset_opacity(value);
  }

  /** @brief Marks the preset whose parameters now drive the frame, so a
      profile attributes a fade's dimming frames to the departing preset. */
  void log_preset() {
#ifdef HS_PROFILE_ENABLE
    hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex() + 1),
            static_cast<unsigned>(getPresetCount()));
#endif
  }

  /** @brief Cancels any in-flight transition, restoring a fade's opacity.
   * @param restore_preset Restore an unadopted fade's index and notify observers;
   *        false when a replacement preset will commit its own notification.
   */
  void end_transition(bool restore_preset = true) {
    const bool FADING = transition.active && transition.fades;
    const bool RESTORE = FADING && !transition.adopted && restore_preset &&
                         preset_index != transition.from_index;
    const PresetChange change{preset_index, transition.from_index,
                              PresetChangeOrigin::AUTOMATIC};
    transition.active = false;
    if (RESTORE)
      preset_index = change.to;
    if (FADING)
      set_opacity(1.0f);
    if (RESTORE)
      preset_changed(change);
  }

  /** @brief Preset count: `PRESET_IDS` when present, else the `PRESETS`
      table, else one. */
  static consteval size_t preset_count_of() {
    if constexpr (requires { Derived::PRESET_IDS; }) {
      if constexpr (requires { Derived::PRESETS; })
        static_assert(Derived::PRESET_IDS.size() == Derived::PRESETS.size(),
                      "PRESET_IDS and PRESETS must name the same presets: "
                      "preset_row() indexes the table with a PRESET_IDS "
                      "index");
      return Derived::PRESET_IDS.size();
    } else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS.size();
    else
      return 1;
  }

  /** @brief The row of the preset at @p index: the static `preset(index)`,
      else the `PRESETS` table, else a single preset of the startup params. */
  static auto preset_row(size_t index) {
    if constexpr (requires { Derived::preset(index); })
      return Derived::preset(index);
    else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS[index];
    else {
      static_assert(preset_count_of() <= 1,
                    "an effect with several presets must define a PRESETS "
                    "table or preset(index)");
      return PresetEntry<Params>{initial_params_of(), Segue::Preset::Snap{}};
    }
  }

  /** @brief Whether every departure fits inside the dwell, so a cancelled
      crossfade's blend never outlives its own transition into the next one,
      and the effect takes each departure's hook: blend_params for a Lerp,
      set_preset_opacity for a Fade. */
  static consteval bool departures_are_supported() {
    constexpr bool TAKES_OPACITY =
        requires(Derived &effect) { effect.set_preset_opacity(1.0f); };
    for (size_t index = 0; index < preset_count_of(); ++index) {
      const auto DEPARTURE = preset_departure(index);
      if (Segue::Preset::frames(DEPARTURE) >= Derived::PRESET_DWELL_FRAMES)
        return false;
      if (std::holds_alternative<Segue::Preset::Fade>(DEPARTURE) &&
          !TAKES_OPACITY)
        return false;
      if (std::holds_alternative<Segue::Preset::Lerp>(DEPARTURE) && !BLENDS)
        return false;
    }
    return true;
  }

  /** @brief Startup parameters: `Derived::initial_params()` when declared,
      else the first preset row's parameters. */
  static Params initial_params_of() {
    if constexpr (requires { Derived::initial_params(); })
      return Derived::initial_params();
    else if constexpr (requires { Derived::preset(size_t{0}); })
      return Derived::preset(0).params;
    else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS[0].params;
    else
      return Params{};
  }

  /** @brief Params for the preset at @p index: the member `preset_params`
      hook when the effect patches runtime state in, else the row's own. */
  Params preset_target(size_t index) {
    if constexpr (requires { derived().preset_params(index); })
      return derived().preset_params(index);
    else
      return preset_row(index).params;
  }

  uint16_t preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
};
