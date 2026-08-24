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

/**
 * @brief Runs all preset-container test cases.
 * @return The module's failure count, as reported by end_module().
 */
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

inline int run_presets_tests() {
  hs_test::ModuleFixture fixture("presets");

  test_all_presets_in_ranges_folds_predicate();
  test_apply_if_changed();

  return fixture.result();
}

} // namespace presets_tests
} // namespace hs_test
