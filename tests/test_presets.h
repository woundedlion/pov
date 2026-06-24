/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/presets.h — the Presets<> preset-cycle container
 * (get/next/prev/apply/prev_get, wraparound, and CTAD deduction).
 *
 * Self-contained header. run_presets_tests() returns the module failure count.
 */
#pragma once

#include <array>

#include "core/presets.h"
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
  HS_EXPECT_EQ(p.current_index(), 0);
  HS_EXPECT_EQ(p.prev_index(), 0);
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
  HS_EXPECT_EQ(p.current_index(), 0);
  p.prev();
  HS_EXPECT_EQ(p.get().id, 42);
  HS_EXPECT_EQ(p.current_index(), 0);
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
  auto scope = hs_test::begin_module("presets");

  test_initial_state();
  test_next_cycles_forward_and_tracks_prev();
  test_prev_cycles_backward();
  test_apply_copies_current();
  test_single_entry_wraps_in_place();
  test_ctad_deduces_size();

  return hs_test::end_module(scope);
}

} // namespace presets_tests
} // namespace hs_test
