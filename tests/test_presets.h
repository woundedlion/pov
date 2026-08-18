/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/control/presets.h — the Presets<> preset-cycle container
 * (get/next/prev/apply/prev_get, wraparound, and CTAD deduction).
 *
 * Self-contained header. run_presets_tests() returns the module failure count.
 */
#pragma once

#include <array>

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

/**
 * @brief Builds a 3-entry Presets fixture with ids 1..3.
 * @return A Presets<DummyParams, 3> seeded with the common test fixture.
 */
inline auto make_presets() {
  return Presets<DummyParams, 3>{std::array<PresetEntry<DummyParams>, 3>{{
      {DummyParams{1, 1.5f}},
      {DummyParams{2, 2.5f}},
      {DummyParams{3, 3.5f}},
  }}};
}

/**
 * @brief Verifies a freshly built container points at the first entry.
 * @details Checks current_idx and prev_idx both start at 0 and the entry
 *          count is 3.
 */
inline void test_initial_state() {
  auto p = make_presets();
  HS_EXPECT_EQ(p.get().id, 1);
  HS_EXPECT_NEAR(p.get().value, 1.5f, 1e-6f);
  HS_EXPECT_EQ(static_cast<int>(p.current_index()), 0);
  HS_EXPECT_EQ(static_cast<int>(p.prev_index()), 0);
  HS_EXPECT_EQ(static_cast<int>(p.get_entries().size()), 3);
}

/**
 * @brief Verifies next() advances forward, wraps, and tracks the prev entry.
 * @details Confirms next() steps through ids 1->2->3->1 and that prev_get()
 *          always points at the entry that was current before the call.
 */
inline void test_next_cycles_forward_and_tracks_prev() {
  auto p = make_presets();

  p.next();
  HS_EXPECT_EQ(p.get().id, 2);
  HS_EXPECT_EQ(p.prev_get().id, 1); // prev points at the pre-next entry

  p.next();
  HS_EXPECT_EQ(p.get().id, 3);
  HS_EXPECT_EQ(p.prev_get().id, 2);

  p.next();
  HS_EXPECT_EQ(p.get().id, 1);
  HS_EXPECT_EQ(p.prev_get().id, 3);
}

/**
 * @brief Verifies prev() steps backward and wraps from the first to the last.
 */
inline void test_prev_cycles_backward() {
  auto p = make_presets();

  p.prev();
  HS_EXPECT_EQ(p.get().id, 3);
  HS_EXPECT_EQ(p.prev_get().id, 1);

  p.prev();
  HS_EXPECT_EQ(p.get().id, 2);
  HS_EXPECT_EQ(p.prev_get().id, 3);
}

/**
 * @brief Verifies apply() overwrites the target with a copy of the current entry.
 */
inline void test_apply_copies_current() {
  auto p = make_presets();
  p.next(); // now on id 2

  DummyParams target{0, 0.0f};
  p.apply(target);
  HS_EXPECT_EQ(target.id, 2);
  HS_EXPECT_NEAR(target.value, 2.5f, 1e-6f);
}

/**
 * @brief Verifies that with a single entry next()/prev() keep current_idx at 0.
 */
inline void test_single_entry_wraps_in_place() {
  Presets<DummyParams, 1> p{std::array<PresetEntry<DummyParams>, 1>{{
      {DummyParams{42, 0.0f}},
  }}};
  p.next();
  HS_EXPECT_EQ(p.get().id, 42);
  HS_EXPECT_EQ(static_cast<int>(p.current_index()), 0);
  p.prev();
  HS_EXPECT_EQ(p.get().id, 42);
  HS_EXPECT_EQ(static_cast<int>(p.current_index()), 0);
}

/**
 * @brief Verifies select() jumps to an index, records the outgoing entry, and
 *        rejects an out-of-range index without moving.
 */
inline void test_select_jumps_and_rejects_out_of_range() {
  auto p = make_presets();

  HS_EXPECT_TRUE(p.select(2));
  HS_EXPECT_EQ(p.get().id, 3);
  HS_EXPECT_EQ(p.prev_get().id, 1);
  HS_EXPECT_EQ(static_cast<int>(p.current_index()), 2);

  HS_EXPECT_FALSE(p.select(3));
  HS_EXPECT_EQ(p.get().id, 3);
  HS_EXPECT_EQ(p.prev_get().id, 1);
  HS_EXPECT_EQ(static_cast<int>(p.current_index()), 2);

  // Re-selecting the current index still records it as the outgoing entry.
  HS_EXPECT_TRUE(p.select(2));
  HS_EXPECT_EQ(p.prev_get().id, 3);
}

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
 * @brief Verifies class template argument deduction infers the entry type and count.
 * @details Constructs Presets from an array without spelling out <DummyParams, 2>
 *          and confirms the deduced size and forward cycling.
 */
inline void test_ctad_deduces_size() {
  // No explicit <DummyParams, 2>: the deduction guide infers it from the array.
  Presets ctad{std::array<PresetEntry<DummyParams>, 2>{{
      {DummyParams{7, 0.0f}},
      {DummyParams{8, 0.0f}},
  }}};
  HS_EXPECT_EQ(static_cast<int>(ctad.get_entries().size()), 2);
  ctad.next();
  HS_EXPECT_EQ(ctad.get().id, 8);
}

/**
 * @brief Runs all preset-container test cases.
 * @return The module's failure count, as reported by end_module().
 */
inline int run_presets_tests() {
  hs_test::ModuleFixture fixture("presets");

  test_initial_state();
  test_next_cycles_forward_and_tracks_prev();
  test_prev_cycles_backward();
  test_apply_copies_current();
  test_single_entry_wraps_in_place();
  test_select_jumps_and_rejects_out_of_range();
  test_all_presets_in_ranges_folds_predicate();
  test_ctad_deduces_size();

  return fixture.result();
}

} // namespace presets_tests
} // namespace hs_test
