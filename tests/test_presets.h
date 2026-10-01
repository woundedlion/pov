/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Preset transitions, admission saturation, re-arm rejection, and fade cancellation.
 */
#pragma once

#include <array>
#include <vector>

#include "core/control/choreography.h"
#include "core/control/params.h"
#include "core/control/presets.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace presets_tests {

/**
 * @brief Minimal stand-in payload for exercising the preset container.
 * @details Avoids depending on any real preset struct; `id` doubles as an
 *          identity marker in assertions, while `value` checks float copying.
 */
struct DummyParams {
  int id;
  float value;
};

/** @brief Range predicate every entry of the fixture satisfies. */
constexpr bool id_below_four(const DummyParams &d) { return d.id < 4; }
/** @brief Range predicate the fixture's last entry fails. */
constexpr bool id_below_three(const DummyParams &d) { return d.id < 3; }
/** @brief Range predicate the fixture's first entry fails. */
constexpr bool id_above_one(const DummyParams &d) { return d.id > 1; }

/** @brief The fixture's entries, as a constant expression. */
constexpr std::array<PresetEntry<DummyParams>, 3> CONST_ENTRIES{{
    {DummyParams{1, 1.5f}},
    {DummyParams{2, 2.5f}},
    {DummyParams{3, 3.5f}},
}};

static_assert(all_presets_in_ranges(CONST_ENTRIES, id_below_four));
static_assert(!all_presets_in_ranges(CONST_ENTRIES, id_below_three));

/**
 * @brief Verifies all_presets_in_ranges() folds the predicate over every entry.
 * @details The static_asserts above cover the constant-expression use the helper
 *          exists for; these calls pin that a failure at either end of the table
 *          is reported.
 */
inline void test_all_presets_in_ranges_folds_predicate() {
  HS_EXPECT_FALSE(all_presets_in_ranges(CONST_ENTRIES, id_above_one));
}

// --- apply_if_changed -------------------------------------------------------

/**
 * @brief Verifies apply_if_changed invokes the callable only when the value changes.
 * @details The callable fires only when the incoming value differs from the latched
 *          `last`, then `last` is updated — the live-slider debounce idiom. The test
 *          covers no-change (no call), change (one call, latched), and repeat (no call).
 */
inline void test_apply_if_changed() {
  int last = 5;
  int applied = -1;
  int call_count = 0;

  apply_if_changed(5, last, [&](int v) {
    applied = v;
    ++call_count;
  });
  HS_EXPECT_EQ(call_count, 0);
  HS_EXPECT_EQ(last, 5);

  apply_if_changed(8, last, [&](int v) {
    applied = v;
    ++call_count;
  });
  HS_EXPECT_EQ(call_count, 1);
  HS_EXPECT_EQ(applied, 8);
  HS_EXPECT_EQ(last, 8);

  apply_if_changed(8, last, [&](int v) {
    applied = v;
    ++call_count;
  });
  HS_EXPECT_EQ(call_count, 1);
  HS_EXPECT_EQ(last, 8);
}

/** @brief Params whose default differs from preset 0, so the startup source is
    observable. */
struct BootParams {
  float value = -1.0f;
};

/**
 * @brief Effect declaring PRESET_IDS and a static preset(index), no
 *        initial_params().
 * @details Preset 0 carries a value the struct default cannot produce, so the
 * parameters the base boots with name which resolver supplied them.
 */
struct PresetZeroBootEffect
    : public ChoreographedEffect<PresetZeroBootEffect, BootParams> {
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"first",
                                                              "second"};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 60;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  static constexpr PresetEntry<BootParams> preset(size_t index) {
    return {{index == 0 ? 7.0f : 9.0f}, Segue::Preset::Snap{}};
  }
  static constexpr bool valid_params(const BootParams &value) {
    return value.value >= 0.0f;
  }

  PresetZeroBootEffect() : ChoreographedEffect(8, 8) {}
  void draw_frame() override {}

  float boot_value() const { return params.value; }
};

/**
 * @brief Verifies a PRESET_IDS-shaped effect boots at preset(0).
 * @details The base reports preset 0 from construction, so starting at the
 * struct defaults instead would render parameters no preset names while
 * claiming to be on the first one.
 */
inline void test_preset_zero_supplies_startup_params() {
  hs_test::reset_globals();
  PresetZeroBootEffect effect;
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.boot_value(),
               PresetZeroBootEffect::preset(0).params.value);
  HS_EXPECT_NE(effect.boot_value(), BootParams{}.value);
}

/** @brief Params for the dwell-override fixture; presets carry distinct
    values so a snap is observable. */
struct HoldParams {
  float value = 0.0f;
};

/**
 * @brief Two-preset snapping effect exposing the choreography's dwell controls.
 * @details Segue::Preset::Snap keeps the advance synchronous, so a preset
 * index change lands on the frame the dwell retires with no crossfade in
 * between.
 */
struct HoldEffect : public ChoreographedEffect<HoldEffect, HoldParams> {
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"first",
                                                              "second"};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 40;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  static constexpr PresetEntry<HoldParams> preset(size_t index) {
    return {{index == 0 ? 1.0f : 2.0f}, Segue::Preset::Snap{}};
  }
  static constexpr bool valid_params(const HoldParams &) { return true; }

  HoldEffect() : ChoreographedEffect(8, 8) {}
  void draw_frame() override {}

  void arm() { begin_choreography(); }
  void hold(uint16_t frames) { hold_initial_preset(frames); }
  void tick() { step_choreography(); }
  float value() const { return params.value; }
};

/**
 * @brief Verifies hold_initial_preset() replaces the dwell for the first
 *        transition only, and that zero holds nothing.
 * @details The override exists so an effect can stagger its first preset move
 * off the shared cadence; a hold that leaked into later moves would retune the
 * whole choreography instead of its opening frame.
 */
inline void test_hold_initial_preset_overrides_first_dwell() {
  hs_test::reset_globals();
  HoldEffect effect;
  effect.arm();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{2});
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});

  effect.hold(3);
  effect.tick();
  effect.tick();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  effect.tick();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_EQ(effect.value(), HoldEffect::preset(1).params.value);

  // The advance restored the authored dwell, so the next move is a full
  // PRESET_DWELL_FRAMES away rather than another three frames.
  for (uint16_t f = 1; f < HoldEffect::PRESET_DWELL_FRAMES; ++f)
    effect.tick();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  effect.tick();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});

  // Zero holds nothing: the next frame starts the transition.
  effect.hold(0);
  effect.tick();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
}

/** @brief Preset fixture with saturated animation admission. */
struct SaturatedPresetEffect
    : ChoreographedEffect<SaturatedPresetEffect, HoldParams> {
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"first",
                                                              "second"};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 40;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  SaturatedPresetEffect() : ChoreographedEffect(8, 8) {}
  static constexpr PresetEntry<HoldParams> preset(size_t index) {
    return {{index == 0 ? 1.0f : 2.0f},
            Segue::Preset::Lerp{4, math::ease_linear}};
  }
  static constexpr bool valid_params(const HoldParams &) { return true; }
  void draw_frame() override {}
  void arm() { begin_choreography(); }
  void tick() { step_choreography(); }
  bool attempt() { return advance_preset(); }
  float value() const { return params.value; }
  bool blending() const { return transition.active; }
  void blend_params(float t) {
    params.value = transition.from.value +
                   (transition.to.value - transition.from.value) * t;
  }
  void saturate() {
    for (size_t i = 0; i < Timeline::MAX_EVENTS; ++i)
      timeline.add(10000, Animation::PeriodicTimer(10000, [](Canvas &) {}));
  }
  void clear_events() { timeline.clear(); }
  void cancel() { parameter_written(); }
  void step_events(Canvas &canvas) { timeline.step(canvas); }
  uint16_t elapsed() const { return transition.elapsed_frames; }
};

/** @brief Pins preset saturation veto restarts dwell. */
inline void test_preset_saturation_veto_restarts_dwell() {
  hs_test::reset_globals();
  SaturatedPresetEffect effect;
  effect.arm();
  for (int i = 0; i < 10; ++i)
    effect.tick();
  effect.saturate();
  HS_EXPECT_FALSE(effect.attempt());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.value(), 1.0f);
  HS_EXPECT_FALSE(effect.blending());
  effect.clear_events();
  for (int i = 1; i < SaturatedPresetEffect::PRESET_DWELL_FRAMES; ++i)
    effect.tick();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  effect.tick();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_TRUE(effect.blending());
}

/** @brief Preset fixture for fade cancellation. */
struct FadePresetEffect : ChoreographedEffect<FadePresetEffect, HoldParams> {
  using Origin = PresetChangeOrigin;
  struct Notification {
    PresetChange change;
    size_t visible_index;
    bool active;
    float opacity;
  };
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"first",
                                                              "second"};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 40;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  FadePresetEffect() : ChoreographedEffect(8, 8) {}
  static constexpr PresetEntry<HoldParams> preset(size_t index) {
    return {{index == 0 ? 1.0f : 2.0f}, Segue::Preset::Fade{4}};
  }
  static constexpr bool valid_params(const HoldParams &) { return true; }
  void draw_frame() override {}
  void arm() { begin_choreography(); }
  bool attempt() { return advance_preset(); }
  void clear_timeline() { timeline.clear(); }
  bool active() const { return transition.active; }
  void progress(float value) { run_transition(value); }
  void cancel() { parameter_written(); }
  void edit(float value) {
    params.value = value;
    parameter_written();
  }
  void set_preset_opacity(float value) { opacity = value; }
  void preset_changed(const PresetChange &change) override {
    notifications.push_back(
        {change, getPresetIndex(), transition.active, opacity});
    if (cancel_from_notification)
      cancel();
  }
  float value() const { return params.value; }
  float opacity = 1.0f;
  bool cancel_from_notification = false;
  std::vector<Notification> notifications;
};

inline void test_timeline_clear_releases_fade() {
  hs_test::reset_globals();
  FadePresetEffect effect;
  effect.arm();
  HS_EXPECT_TRUE(effect.attempt());
  effect.progress(0.25f);
  HS_EXPECT_TRUE(effect.active());
  HS_EXPECT_EQ(effect.opacity, 0.5f);
  const size_t notifications = effect.notifications.size();
  effect.clear_timeline();
  HS_EXPECT_FALSE(effect.active());
  HS_EXPECT_EQ(effect.opacity, 1.0f);
  HS_EXPECT_EQ(effect.notifications.size(), notifications);
  HS_EXPECT_TRUE(effect.attempt());
  HS_EXPECT_TRUE(effect.active());
}

/** @brief Pins cancelled fade notifies committed index. */
inline void test_cancelled_fade_notifies_committed_index() {
  enum class Cancel { EDIT, RESTORE, MANUAL_FROM, MANUAL_TO, SYNCHRONIZED };
  for (const bool adopted : {false, true}) {
    for (const Cancel action :
         {Cancel::EDIT, Cancel::RESTORE, Cancel::MANUAL_FROM, Cancel::MANUAL_TO,
          Cancel::SYNCHRONIZED}) {
      hs_test::reset_globals();
      FadePresetEffect effect;
      effect.arm();
      HS_EXPECT_TRUE(effect.attempt());
      effect.progress(adopted ? 0.5f : 0.25f);
      HS_EXPECT_EQ(effect.notifications.size(), size_t{1});
      const bool REPLACEMENT = action == Cancel::MANUAL_FROM ||
                               action == Cancel::MANUAL_TO ||
                               action == Cancel::SYNCHRONIZED;
      const size_t INDEX =
          REPLACEMENT ? size_t(action == Cancel::MANUAL_TO) : size_t(adopted);
      float expected_value = float(INDEX + 1);
      switch (action) {
      case Cancel::EDIT:
        effect.cancel_from_notification = true;
        effect.edit(3.0f);
        expected_value = 3.0f;
        break;
      case Cancel::RESTORE:
        HS_EXPECT_TRUE(effect.restore_parameters({1, {4.0f}}));
        expected_value = 4.0f;
        break;
      case Cancel::MANUAL_FROM:
      case Cancel::MANUAL_TO:
        HS_EXPECT_TRUE(effect.selectPreset(INDEX));
        break;
      case Cancel::SYNCHRONIZED:
        HS_EXPECT_TRUE(effect.synchronizePreset(INDEX));
        break;
      }
      const size_t COUNT = REPLACEMENT || !adopted ? 2 : 1;
      HS_EXPECT_EQ(effect.notifications.size(), COUNT);
      if (COUNT == 2 && effect.notifications.size() == COUNT) {
        const auto &notification = effect.notifications.back();
        HS_EXPECT_EQ(notification.change.from, size_t{1});
        HS_EXPECT_EQ(notification.change.to, INDEX);
        const auto ORIGIN = action == Cancel::SYNCHRONIZED
                                ? FadePresetEffect::Origin::SYNCHRONIZED
                            : REPLACEMENT ? FadePresetEffect::Origin::MANUAL
                                          : FadePresetEffect::Origin::AUTOMATIC;
        HS_EXPECT_TRUE(notification.change.origin == ORIGIN);
        HS_EXPECT_EQ(notification.visible_index, INDEX);
        HS_EXPECT_FALSE(notification.active);
        HS_EXPECT_EQ(notification.opacity, 1.0f);
      }
      effect.cancel();
      effect.progress(1.0f);
      HS_EXPECT_EQ(effect.notifications.size(), COUNT);
      HS_EXPECT_EQ(effect.getPresetIndex(), INDEX);
      HS_EXPECT_EQ(effect.value(), expected_value);
      HS_EXPECT_EQ(effect.opacity, 1.0f);
    }
  }
}

/** @brief Pins cancelled fade names visible preset. */
inline void test_cancelled_fade_names_visible_preset() {
  hs_test::reset_globals();
  FadePresetEffect effect;
  effect.arm();
  HS_EXPECT_TRUE(effect.attempt());
  effect.progress(0.25f);
  effect.cancel();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.value(), 1.0f);
  HS_EXPECT_EQ(effect.opacity, 1.0f);
  HS_EXPECT_TRUE(effect.attempt());
  effect.progress(0.5f);
  effect.cancel();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_EQ(effect.value(), 2.0f);
}

/** @brief Pins cancelled crossfade event cannot step replacement. */
inline void test_cancelled_crossfade_event_cannot_step_replacement() {
  hs_test::reset_globals();
  SaturatedPresetEffect effect;
  Canvas canvas(effect);
  effect.arm();
  HS_EXPECT_TRUE(effect.attempt());
  effect.step_events(canvas);
  effect.cancel();
  HS_EXPECT_TRUE(effect.attempt());
  effect.step_events(canvas);
  HS_EXPECT_EQ(effect.elapsed(), uint16_t{1});
}

/** @brief Pins preset crossfade rejects rearming. */
inline void test_preset_crossfade_rejects_rearming() {
  hs_test::reset_globals();
  SaturatedPresetEffect effect;
  effect.arm();
  HS_EXPECT_TRUE(effect.attempt());
  HS_EXPECT_FALSE(effect.attempt());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_TRUE(effect.blending());
}

/**
 * @brief Runs all preset-container test cases.
 * @return The module's failure count, as reported by end_module().
 */
inline int run_presets_tests() {
  hs_test::ModuleFixture fixture("presets");

  test_all_presets_in_ranges_folds_predicate();
  test_apply_if_changed();
  test_preset_zero_supplies_startup_params();
  test_hold_initial_preset_overrides_first_dwell();
  test_preset_saturation_veto_restarts_dwell();
  test_preset_crossfade_rejects_rearming();
  test_cancelled_fade_names_visible_preset();
  test_timeline_clear_releases_fade();
  test_cancelled_fade_notifies_committed_index();
  test_cancelled_crossfade_event_cannot_step_replacement();

  return fixture.result();
}

} // namespace presets_tests
} // namespace hs_test
