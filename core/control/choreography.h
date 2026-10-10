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
#include "control/fields.h"
#include "platform/platform.h"
#include "math/easing.h"
#include "render/canvas.h"

/**
 * @brief Base class for an effect driven by one parameter struct that cycles
 *        through a list of presets.
 * @details Each preset holds for `PRESET_DWELL_FRAMES`, then moves to the next
 * by its departure policy: Snap switches at once, Lerp interpolates the
 * parameters, and Fade dims out, switches at the midpoint and dims back in.
 * Manual preset selection always snaps. The live parameters can be saved and
 * restored as snapshots tagged with `PARAMETER_SCHEMA_VERSION`.
 *
 * `Derived` must provide:
 * - `PARAMETER_SCHEMA_VERSION` and `PRESET_DWELL_FRAMES`.
 * - parameter_fields() (defaults to `Params::FIELDS`) or valid_params(params).
 * - Its presets as `PresetEntry` rows, in a `PRESETS` array or a static
 *   `preset(index)`, optionally named by `PRESET_IDS`. With none of these the
 *   effect has a single preset.
 *
 * `Derived` calls `begin_choreography()` once from init() and
 * `step_choreography()` every frame.
 *
 * Optional members `Derived` may define (befriend this base if non-public):
 * - `blend_params(progress)`: writes a Lerp's interpolated parameters; the
 *   default interpolates the fields parameter_fields() describes.
 * - `set_preset_opacity(value)`: receives a Fade's opacity; required for Fade.
 * - `adopt_params(target)`: installs new parameters; define it to update state
 *   computed from them.
 * - `transition_armed(target)`: called when a Lerp starts.
 * - `finish_blend(target)`: used instead of adopt_params() at the end of a Lerp
 *   and at a Fade's midpoint.
 * - `preset_params(index)`: computes a preset's parameters at runtime instead
 *   of taking them from its row.
 * - `initial_params()`: startup parameters, instead of the first preset's.
 * @tparam Derived The effect class deriving from this base.
 * @tparam ParamsT The effect's parameter struct.
 */
template <typename Derived, typename ParamsT>
class ChoreographedEffect : public Effect {
public:
  using Params = ParamsT; ///< The effect's parameter aggregate.

  /**
   * @brief Uses the parameter aggregate's static descriptions when present.
   * @return `Params::FIELDS`.
   */
  static constexpr auto parameter_fields()
    requires requires { Params::FIELDS; }
  {
    return Params::FIELDS;
  }

  /**
   * @brief Checks the derived effect's described parameter ranges.
   * @param value Parameters to check.
   * @return True when every described field is in range.
   */
  static constexpr bool valid_params(const Params &value)
    requires requires { Derived::parameter_fields(); }
  {
    return Control::valid_fields(value, Derived::parameter_fields());
  }

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
   * @brief Replaces the live parameters with a snapshot's.
   * @details Rejects, changing nothing, a snapshot with a different
   * `PARAMETER_SCHEMA_VERSION` or invalid parameters. On success any
   * transition in progress is cancelled and the dwell restarts; a cancelled
   * fade that had not reached its midpoint returns to the preset it was
   * leaving.
   * @param snapshot Snapshot to restore.
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

#if HS_ENABLE_EFFECT_CONTROL_API
  /**
   * @brief Selects and pauses one authored preset for profiling.
   * @param index Authored preset index.
   */
  void profile_select_preset(size_t index) {
    HS_CHECK(index < authored_preset_count(),
             "profile preset index out of range");
    HS_CHECK(this->selectPreset(index), "profile preset selection failed");
    if constexpr (requires { derived().after_profile_preset(); })
      derived().after_profile_preset();
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(authored_preset_count()));
  }
#endif

  /**
   * @brief Number of presets the effect defines.
   * @return Authored preset count.
   */
  static consteval size_t authored_preset_count() { return preset_count_of(); }

  /** @brief How the preset at @p index moves to the next one; Snap for an
      effect with no preset rows.
   * @param index Preset index.
   * @return The preset's departure.
   */
  static constexpr Segue::Preset::Departure preset_departure(size_t index) {
    if constexpr (requires { Derived::preset(index); })
      return Derived::preset(index).segue;
    else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS[index].segue;
    else
      return Segue::Preset::Snap{};
  }

  /** @brief Frames to play every preset once: all dwells plus every
      departure except the last preset's.
   * @return Frame count of one full preset cycle.
   */
  static consteval size_t cycle_frames() {
    size_t frames = preset_count_of() * Derived::PRESET_DWELL_FRAMES;
    for (size_t index = 0; index + 1 < preset_count_of(); ++index)
      frames += Segue::Preset::frames(preset_departure(index));
    return frames;
  }

protected:
  /** @return The preset being left while a fade is before its midpoint, else
   *          the committed preset. */
  size_t displayed_preset_index() const override {
    return transition.active && transition.fades && !transition.adopted
               ? transition.from_index
               : preset_index;
  }

  /** @brief Registers the derived effect's typed parameter descriptions. */
  HS_COLD_MEMBER void register_described_params() {
    Control::register_fields(
        params, Derived::parameter_fields(),
        [this](const char *name, auto *target, const auto &spec) {
          this->register_param(name, target, spec);
        });
  }

  /**
   * @brief Interpolates the derived effect's described parameter members.
   * @param progress Transition progress in [0, 1].
   */
  HS_COLD_MEMBER void blend_params(float progress)
    requires requires { Derived::parameter_fields(); }
  {
    Control::interpolate_fields(params, transition.from, transition.to,
                                progress, Derived::parameter_fields());
  }

  /** @brief An automatic preset transition's endpoints and progress. */
  struct Transition {
    Params from{};               /**< Parameters the transition departs from. */
    Params to{};                 /**< Parameters it lands on. */
    bool active = false;         /**< False when no transition is in flight. */
    uint16_t elapsed_frames = 0; /**< Unpaused steps elapsed. */
    uint16_t frames = 0;         /**< Steps the departure spans. */
    bool fades = false;          /**< A Fade departure rather than a Lerp. */
    bool adopted = false;        /**< A fade has passed its midpoint. */
    uint16_t serial = 0; ///< Increments per armed transition; stale steps skip.
    size_t from_index = 0; /**< Preset the transition is leaving. */
  };

  /**
   * @brief Constructs the effect base.
   * @param W Canvas width in pixels.
   * @param H Canvas height in pixels.
   * @param cfg Effect configuration.
   */
  HS_COLD_MEMBER ChoreographedEffect(int W, int H, EffectConfig cfg = {})
      : Effect(W, H, cfg) {}

  /**
   * @brief Installs @p target as the live parameters.
   * @param target Parameters to install.
   */
  void adopt_params(const Params &target) { params = target; }

  /** @brief Called when a Lerp transition to @p target starts. */
  void transition_armed(const Params &) {}

  /**
   * @brief Sets the frames left before the next automatic transition.
   * @param frames Frames to hold; 0 starts a transition on the next step.
   */
  HS_COLD_MEMBER void hold_initial_preset(uint16_t frames) {
    preset_dwell_remaining = frames;
  }

  /**
   * @brief Switches to the preset @p change names.
   * @details A manual or synchronized change snaps. An automatic change follows
   * the departing preset's policy; if a transition is already running or the
   * timeline is full it is refused, and a full timeline restarts the dwell.
   * @param change The requested preset change.
   * @return False if the change was refused.
   */
  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    if constexpr (requires { Derived::parameter_fields(); })
      static_assert(Control::curves_supported(Derived::parameter_fields()),
                    "discrete fields require MIDPOINT or SNAP curves");
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
   * @brief Cancels any running transition after a manual parameter edit.
   * @details Restores full opacity and restarts the dwell.
   */
  HS_COLD_MEMBER void parameter_written() override {
    end_transition();
    preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
  }

  /**
   * @brief Sets the preset count and registers a timeline-clear hook that
   *        cancels transitions. Call once from init().
   */
  HS_COLD_MEMBER void begin_choreography() {
    configure_presets(preset_count_of());
    timeline.add_clear_hook(this, [](void *context) {
      auto &effect = *static_cast<ChoreographedEffect *>(context);
      effect.end_transition();
      effect.preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
      effect.set_opacity(1.0f);
    });
  }

  /// Counts down the dwell and starts the next preset's transition when it
  /// expires. Call every frame. While paused it does nothing, except that a
  /// fade in progress shows at full opacity.
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
   * @brief Advances the running transition one step.
   * @param progress Eased transition progress in [0, 1].
   * @details Does nothing once the transition has been cancelled.
   */
  HS_COLD_MEMBER void run_transition(float progress) {
    if (!transition.active)
      return;
    const bool COMPLETE = ++transition.elapsed_frames >= transition.frames;
    const float PROGRESS = COMPLETE ? 1.0f : progress;
    if (transition.fades) {
      if (PROGRESS >= 0.5f && !transition.adopted) {
        if constexpr (requires { derived().finish_blend(transition.to); })
          derived().finish_blend(transition.to);
        else
          derived().adopt_params(transition.to);
        transition.adopted = true;
        log_preset();
      }
      set_opacity(fabsf(1.0f - 2.0f * PROGRESS));
    } else if constexpr (BLENDS) {
      derived().blend_params(PROGRESS);
    }
    if (COMPLETE) {
      if (!transition.fades) {
        if constexpr (requires { derived().finish_blend(transition.to); })
          derived().finish_blend(transition.to);
        else
          derived().adopt_params(transition.to);
      }
      transition.active = false;
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
    }
  }

  /** Live parameters; the registered sliders write straight into these. */
  Params params = initial_params_of();
  Transition transition; ///< The in-flight automatic preset transition.
  Timeline timeline;     ///< Runs the effect's animations.

private:
  Derived &derived() { return static_cast<Derived &>(*this); }

  /** @brief Whether `Derived` has blend_params(). */
  static constexpr bool BLENDS =
      requires(Derived &effect) { effect.blend_params(0.0f); };

  /** @brief Forwards a fade's opacity to set_preset_opacity() if defined. */
  void set_opacity(float value) {
    if constexpr (requires { derived().set_preset_opacity(value); })
      derived().set_preset_opacity(value);
  }

  /** @brief Logs the current preset number (profiling builds only). */
  void log_preset() {
#ifdef HS_PROFILE_ENABLE
    hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex() + 1),
            static_cast<unsigned>(getPresetCount()));
#endif
  }

  /** @brief Cancels any running transition and restores a fade's opacity.
   * @param restore_preset If a fade has not reached its midpoint, return to the
   *        preset it was leaving and call preset_changed().
   */
  HS_COLD_MEMBER void end_transition(bool restore_preset = true) {
    const bool FADING = transition.active && transition.fades;
    const bool RESTORE = FADING && !transition.adopted && restore_preset &&
                         preset_index != transition.from_index;
    const PresetChange change{preset_index, transition.from_index,
                              PresetChangeOrigin::RESTORED};
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

  /** @brief Whether every departure fits inside the dwell and the effect takes
      each departure's hook: blend_params for a Lerp, set_preset_opacity for a
      Fade. */
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

  /** @brief Parameters for the preset at @p index: `Derived::preset_params()`
      if defined, else the preset row's. */
  Params preset_target(size_t index) {
    if constexpr (requires { derived().preset_params(index); })
      return derived().preset_params(index);
    else
      return preset_row(index).params;
  }

  uint16_t preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
};
