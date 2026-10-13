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

#include "core/animation/animation.h"
#include "core/control/presets.h"
#include "core/control/choreography_hooks.h"
#include "core/control/fields.h"
#include "core/platform/platform.h"
#include "core/math/easing.h"
#include "core/render/canvas.h"

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
 * - Parameter descriptions through the `PARAMETER_FIELDS` hook binding or
 *   `Params::FIELDS`, or validation through the `VALID_PARAMS` hook binding.
 * - Its presets as `PresetEntry` rows in a `PRESETS` array, or through the
 *   `PRESET_ROW` hook binding with `PRESET_IDS` naming its rows. `PRESETS` may
 *   also supply `PRESET_IDS`. With none of these the effect has a single preset.
 *
 * `Derived` calls `begin_choreography()` once from init() and
 * `step_choreography()` every frame.
 *
 * `choreography_hooks()` returns a `Control::HookSet` owned by `HookOwner`.
 * Bind customizations through `Control::HookBinding`; a method without a binding
 * is not a customization. Signatures and declaration ownership are checked at
 * construction. An empty hook list uses the base defaults.
 *
 * Optional hook bindings:
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
 * @tparam HookOwner Type declaring the hook list; defaults to Derived.
 */
template <typename Derived, typename ParamsT, typename HookOwner = Derived>
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
  static constexpr bool valid_params(const Params &value) {
    return valid_params_of(value);
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
        !valid_params_of(snapshot.params))
      return false;
    end_transition();
    adopt_target(snapshot.params);
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
    if constexpr (has_hook<Control::ChoreographyHook::AFTER_PROFILE_PRESET>())
      invoke_hook<Control::ChoreographyHook::AFTER_PROFILE_PRESET>();
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
    if constexpr (has_hook<Control::ChoreographyHook::PRESET_ROW>())
      return decltype(Derived::choreography_hooks())::template invoke_static<
                 Control::ChoreographyHook::PRESET_ROW>(index)
          .segue;
    else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS[index].segue;
    else
      return Segue::Preset::Snap{};
  }

  /** @brief Conservative show budget: all dwells plus every departure
      except the last preset's. Completing a Lerp or Fade also consumes the
      next dwell's first frame, leaving one slack frame per included departure.
   * @return Conservative frame budget through the final dwell.
   */
  static consteval size_t cycle_frames() {
    size_t frames = preset_count_of() * Derived::PRESET_DWELL_FRAMES;
    for (size_t index = 0; index + 1 < preset_count_of(); ++index)
      frames += Segue::Preset::frames(preset_departure(index));
    return frames;
  }

protected:
  /** @brief Whether an explicit hook binding is present.
   * @tparam Kind Customization point.
   * @return True for a bound hook. */
  template <Control::ChoreographyHook Kind> static consteval bool has_hook() {
    return decltype(Derived::choreography_hooks())::template HAS<Kind>;
  }

  /** @brief Calls an explicitly bound hook.
   * @tparam Kind Customization point.
   * @tparam Args Argument types.
   * @param args Hook arguments.
   * @return Hook result. */
  template <Control::ChoreographyHook Kind, typename... Args>
  constexpr decltype(auto) invoke_hook(Args &&...args) {
    return decltype(Derived::choreography_hooks())::template invoke<Kind>(
        derived(), std::forward<Args>(args)...);
  }

  /** @brief Whether hook-provided or embedded parameter descriptions exist.
   * @return True when fields are available. */
  static consteval bool has_parameter_fields() {
    return has_hook<Control::ChoreographyHook::PARAMETER_FIELDS>() ||
           requires { Params::FIELDS; };
  }

  /** @brief Resolves parameter descriptions from the bound hook or Params.
   * @return The parameter field descriptions. */
  static constexpr auto parameter_fields_of() {
    if constexpr (has_hook<Control::ChoreographyHook::PARAMETER_FIELDS>())
      return decltype(Derived::choreography_hooks())::template invoke_static<
          Control::ChoreographyHook::PARAMETER_FIELDS>();
    else
      return Params::FIELDS;
  }

  /** @brief Validates parameters through the bound hook or field descriptions.
   * @param value Parameter state to check.
   * @return True when the state is valid. */
  static constexpr bool valid_params_of(const Params &value) {
    if constexpr (has_hook<Control::ChoreographyHook::VALID_PARAMS>())
      return decltype(Derived::choreography_hooks())::template invoke_static<
          Control::ChoreographyHook::VALID_PARAMS>(value);
    else
      return Control::valid_fields(value, parameter_fields_of());
  }

  /** @brief Adopts a target through its hook or default assignment.
   * @param target Target parameter state. */
  void adopt_target(const Params &target) {
    if constexpr (has_hook<Control::ChoreographyHook::ADOPT_PARAMS>())
      invoke_hook<Control::ChoreographyHook::ADOPT_PARAMS>(target);
    else
      adopt_params(target);
  }

  /** @brief Finishes a transition through its hook or target adoption.
   * @param target Final parameter state. */
  void finish_target(const Params &target) {
    if constexpr (has_hook<Control::ChoreographyHook::FINISH_BLEND>())
      invoke_hook<Control::ChoreographyHook::FINISH_BLEND>(target);
    else
      adopt_target(target);
  }

  /** @brief Rejects bindings to inherited choreography defaults.
   * @tparam Kind Customization point.
   * @tparam Default Inherited default function pointer. */
  template <Control::ChoreographyHook Kind, auto Default>
  static consteval void reject_default_binding() {
    if constexpr (has_hook<Kind>()) {
      constexpr auto FUNCTION =
          decltype(Derived::choreography_hooks())::template function<Kind>();
      if constexpr (std::is_same_v<std::remove_cv_t<decltype(FUNCTION)>,
                                   decltype(Default)>)
        static_assert(FUNCTION != Default,
                      "choreography hooks: bound inherited default");
    }
  }

  /** @brief Checks optional bindings against inherited defaults. */
  static consteval void validate_hook_defaults() {
    reject_default_binding<Control::ChoreographyHook::ADOPT_PARAMS,
                           &ChoreographedEffect::adopt_params>();
    reject_default_binding<Control::ChoreographyHook::TRANSITION_ARMED,
                           &ChoreographedEffect::transition_armed>();
    reject_default_binding<Control::ChoreographyHook::VALID_PARAMS,
                           &ChoreographedEffect::valid_params>();
    if constexpr (has_parameter_fields()) {
      reject_default_binding<Control::ChoreographyHook::BLEND_PARAMS,
                             &ChoreographedEffect::blend_params>();

      if constexpr (requires { Params::FIELDS; })
        reject_default_binding<Control::ChoreographyHook::PARAMETER_FIELDS,
                               &ChoreographedEffect::parameter_fields>();
    }
  }

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
        params, parameter_fields_of(),
        [this](const char *name, auto *target, const auto &spec) {
          this->register_param(name, target, spec);
        });
  }

  /**
   * @brief Interpolates the derived effect's described parameter members.
   * @param progress Transition progress in [0, 1].
   */
  HS_COLD_MEMBER void blend_params(float progress) {
    Control::interpolate_fields(params, transition.from, transition.to,
                                progress, parameter_fields_of());
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
      : Effect(W, H, cfg) {
    static_assert(
        std::is_same_v<typename decltype(Derived::choreography_hooks())::Owner,
                       HookOwner>,
        "choreography hooks: declaration owner mismatch");
    static_assert(decltype(Derived::choreography_hooks())::template valid<
                      Params, PresetEntry<Params>>(),
                  "choreography hooks: wrong hook signature");
    validate_hook_defaults();
  }

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
    if constexpr (has_parameter_fields())
      static_assert(Control::curves_supported(parameter_fields_of()),
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
      const uint16_t TRANSITION_SERIAL = transition.serial + 1;
      if (timeline.add_get(0,
                           Animation::Progress(
                               [this, TRANSITION_SERIAL](float progress) {
                                 if (transition.serial == TRANSITION_SERIAL)
                                   run_transition(progress);
                               },
                               FRAMES, LERP ? LERP->easing : math::ease_linear),
                           Timeline::Pin::UNPINNED, &anims_paused) == nullptr) {
        preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
        return false;
      }
      transition = {params,     target,          true,  0,
                    FRAMES,     LERP == nullptr, false, TRANSITION_SERIAL,
                    change.from};
      if (LERP) {
        if constexpr (has_hook<Control::ChoreographyHook::TRANSITION_ARMED>())
          invoke_hook<Control::ChoreographyHook::TRANSITION_ARMED>(target);
      }
      return true;
    }
    end_transition(false);
    adopt_target(target);
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
        finish_target(transition.to);
        transition.adopted = true;
        log_preset();
      }
      set_opacity(fabsf(1.0f - 2.0f * PROGRESS));
    } else if constexpr (BLENDS) {
      if constexpr (has_hook<Control::ChoreographyHook::BLEND_PARAMS>())
        invoke_hook<Control::ChoreographyHook::BLEND_PARAMS>(PROGRESS);
      else
        blend_params(PROGRESS);
    }
    if (COMPLETE) {
      if (!transition.fades) {
        finish_target(transition.to);
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

  /** @brief Whether a bound blend or described-field interpolation is available. */
  static constexpr bool BLENDS =
      has_hook<Control::ChoreographyHook::BLEND_PARAMS>() ||
      has_parameter_fields();

  /** @brief Forwards a fade's opacity through its bound hook. */
  void set_opacity(float value) {
    if constexpr (has_hook<Control::ChoreographyHook::PRESET_OPACITY>())
      invoke_hook<Control::ChoreographyHook::PRESET_OPACITY>(value);
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
    static_assert(
        !has_hook<Control::ChoreographyHook::PRESET_ROW>() || requires {
          Derived::PRESET_IDS;
        }, "preset(index) requires PRESET_IDS naming its rows");
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

  /** @brief The row of the preset at @p index: the `PRESET_ROW` hook with
      `PRESET_IDS`, else the `PRESETS` table, else a single startup preset. */
  static auto preset_row(size_t index) {
    if constexpr (has_hook<Control::ChoreographyHook::PRESET_ROW>())
      return decltype(Derived::choreography_hooks())::template invoke_static<
          Control::ChoreographyHook::PRESET_ROW>(index);
    else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS[index];
    else {
      static_assert(preset_count_of() <= 1,
                    "an effect with several presets must define a PRESETS "
                    "table or a PRESET_ROW hook with PRESET_IDS");
      return PresetEntry<Params>{initial_params_of(), Segue::Preset::Snap{}};
    }
  }

  /** @brief Whether every departure fits inside the dwell and the effect takes
      each departure's hook: blend_params for a Lerp, set_preset_opacity for a
      Fade. */
  static consteval bool departures_are_supported() {
    constexpr bool TAKES_OPACITY =
        has_hook<Control::ChoreographyHook::PRESET_OPACITY>();
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

  /** @brief Startup parameters: the `INITIAL_PARAMS` hook when bound,
      else the first preset row's parameters. */
  static Params initial_params_of() {
    if constexpr (has_hook<Control::ChoreographyHook::INITIAL_PARAMS>())
      return decltype(Derived::choreography_hooks())::template invoke_static<
          Control::ChoreographyHook::INITIAL_PARAMS>();
    else if constexpr (has_hook<Control::ChoreographyHook::PRESET_ROW>())
      return decltype(Derived::choreography_hooks())::template invoke_static<
                 Control::ChoreographyHook::PRESET_ROW>(size_t{0})
          .params;
    else if constexpr (requires { Derived::PRESETS; })
      return Derived::PRESETS[0].params;
    else
      return Params{};
  }

  /** @brief Parameters for the preset at @p index: `Derived::preset_params()`
      if defined, else the preset row's. */
  Params preset_target(size_t index) {
    if constexpr (has_hook<Control::ChoreographyHook::PRESET_PARAMS>())
      return invoke_hook<Control::ChoreographyHook::PRESET_PARAMS>(index);
    else
      return preset_row(index).params;
  }

  uint16_t preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
};
