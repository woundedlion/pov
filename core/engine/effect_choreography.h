/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file effect_choreography.h
 * @brief ChoreographedEffect: preset choreography, schema-versioned parameter
 *        snapshots and field-table slider registration over one Params type.
 */

#include <cstdint>

#include "animation/animation.h"
#include "engine/platform.h"
#include "math/easing.h"
#include "render/canvas.h"

/**
 * @brief Effect base owning one parameter set and its preset choreography:
 *        automatic dwell/crossfade transitions, manual preset snaps and
 *        schema-versioned parameter snapshots.
 * @details Curiously recurring: `Derived` supplies `PRESET_IDS`,
 * `PARAMETER_SCHEMA_VERSION`, `PRESET_DWELL_FRAMES`, `initial_params()`,
 * `preset_params(index)` and `valid_params(params)`, plus the transition
 * hooks: `blend_params(progress)` writes the interpolated parameters while a
 * transition is in flight, and shadowing `adopt_params(target)` /
 * `transition_armed(target)` keeps state derived from the parameters
 * consistent across snaps and crossfade arming. `register_fields()`
 * additionally requires a `field_gate_open(gate)` predicate. A `Derived`
 * keeping its hooks non-public befriends this base.
 * @tparam Derived The effect class deriving from this base.
 * @tparam ParamsT The effect's parameter-set type.
 */
template <typename Derived, typename ParamsT>
class ChoreographedEffect : public Effect {
public:
  using Params = ParamsT;
  /** Frames an automatic preset transition spans. */
  static constexpr uint16_t TRANSITION_DURATION = 480;

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
   * cancelled and the preset dwell restarts; on rejection nothing is touched.
   * @param snapshot Snapshot to adopt.
   * @return True when the snapshot's schema version matches and its parameters
   *         are valid, false otherwise.
   */
  bool restore_parameters(const ParameterSnapshot &snapshot) {
    if (snapshot.schema_version != Derived::PARAMETER_SCHEMA_VERSION ||
        !Derived::valid_params(snapshot.params))
      return false;
    transition.active = false;
    derived().adopt_params(snapshot.params);
    preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
    return true;
  }

protected:
  /** @brief An automatic preset transition's endpoints. */
  struct Transition {
    Params from{};       /**< Parameters the transition departs from. */
    Params to{};         /**< Parameters it lands on. */
    bool active = false; /**< False when no transition is in flight. */
  };

  /**
   * @brief Timeline-lerp subject for a preset transition.
   * @details Animation::Lerp drives this through the transition's stored
   * endpoints each frame. A lerp whose transition was cancelled by a manual
   * preset or a snapshot restore keeps stepping but writes nothing.
   * Deliberately never pause-gated: a started transition always runs to its
   * authored endpoint.
   */
  struct PresetBlend {
    ChoreographedEffect *runtime;

    HS_COLD_MEMBER void lerp(const PresetBlend &, const PresetBlend &,
                             float progress) {
      runtime->run_blend(progress);
    }
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
   * @brief Adopts a preset, crossfading only when choreography asked for it.
   * @details An AUTOMATIC change arms a TRANSITION_DURATION crossfade from the
   * live parameters; a manual or synchronized change snaps, since a user
   * driving the control expects the preset it names immediately.
   * @param change The requested preset move.
   * @return Always true: the runtime never vetoes a preset.
   */
  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const Params target = Derived::preset_params(change.to);
    if (change.origin != PresetChangeOrigin::AUTOMATIC) {
      transition.active = false;
      derived().adopt_params(target);
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
      return true;
    }
    transition = {params, target, true};
    derived().transition_armed(target);
    timeline.add(0, Animation::Lerp(preset_blend, preset_blend, preset_blend,
                                    TRANSITION_DURATION, ease_in_out_sin));
    return true;
  }

  /// Retires the preset dwell and starts the next automatic preset transition.
  /// @details Pause suppresses preset selection, so no new transition begins
  /// while paused.
  HS_COLD_MEMBER void begin_automatic_transition() {
    if constexpr (Derived::PRESET_IDS.size() == 1)
      return;
    if (anims_paused || transition.active)
      return;
    if (preset_dwell_remaining > 0 && --preset_dwell_remaining > 0)
      return;
    if (advancePreset()) {
#ifdef HS_PROFILE_ENABLE
      hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex() + 1),
              static_cast<unsigned>(getPresetCount()));
#endif
    }
  }

  /** @brief Registers a slider for every named field in the family's table. */
  template <typename T> HS_COLD_MEMBER void register_fields(T &family) {
    for (const auto &field : T::FIELDS)
      if (field.name != nullptr && Derived::field_gate_open(field.gate))
        register_animated_param(field.name, &(family.*(field.member)),
                                field.min, field.max);
  }

  /** Live parameters; the registered sliders write straight into these. */
  Params params = Derived::initial_params();
  Transition transition;
  PresetBlend preset_blend{this};
  Timeline timeline;

private:
  Derived &derived() { return static_cast<Derived &>(*this); }

  HS_COLD_MEMBER void run_blend(float progress) {
    if (!transition.active)
      return;
    derived().blend_params(progress);
    if (progress >= 1.0f) {
      transition.active = false;
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
    }
  }

  uint16_t preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
};
