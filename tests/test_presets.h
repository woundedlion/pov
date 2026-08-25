/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/control/presets.h — the PresetEntry table row and the
 * all_presets_in_ranges compile-time range check — plus core/control/params.h's
 * apply_if_changed live-value change gate.
 *
 * Self-contained header. run_presets_tests() returns the module failure count.
 */
#pragma once

#include <array>

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
  HS_EXPECT_TRUE(all_presets_in_ranges(CONST_ENTRIES, id_below_four));
  HS_EXPECT_FALSE(all_presets_in_ranges(CONST_ENTRIES, id_below_three));
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
 * @brief Effect declaring PRESET_IDS and a static preset_params, no
 *        initial_params().
 * @details Preset 0 carries a value the struct default cannot produce, so the
 * parameters the base boots with name which resolver supplied them.
 */
struct PresetZeroBootEffect
    : public ChoreographedEffect<PresetZeroBootEffect, BootParams> {
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"first",
                                                              "second"};
  static constexpr Segue::Preset::Snap PRESET_SEGUE{};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 60;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  static constexpr BootParams preset_params(size_t index) {
    return {index == 0 ? 7.0f : 9.0f};
  }
  static constexpr bool valid_params(const BootParams &value) {
    return value.value >= 0.0f;
  }

  PresetZeroBootEffect() : ChoreographedEffect(8, 8) {}
  void draw_frame() override {}

  float boot_value() const { return params.value; }
};

/**
 * @brief Verifies a PRESET_IDS-shaped effect boots at preset_params(0).
 * @details The base reports preset 0 from construction, so starting at the
 * struct defaults instead would render parameters no preset names while
 * claiming to be on the first one.
 */
inline void test_preset_zero_supplies_startup_params() {
  hs_test::reset_globals();
  PresetZeroBootEffect effect;
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.boot_value(),
               PresetZeroBootEffect::preset_params(0).value);
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
  static constexpr Segue::Preset::Snap PRESET_SEGUE{};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 40;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  static constexpr HoldParams preset_params(size_t index) {
    return {index == 0 ? 1.0f : 2.0f};
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
  HS_EXPECT_EQ(effect.value(), HoldEffect::preset_params(1).value);

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

  return fixture.result();
}

} // namespace presets_tests
} // namespace hs_test
