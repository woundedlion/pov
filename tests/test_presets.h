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

// Minimal stand-in payload so the container is exercised without depending on
// any real preset struct; id doubles as an identity marker in assertions.
struct DummyParams {
  int id;
  float value;
};

// Builds a 3-entry Presets with ids 1..3 for the common test fixture.
inline auto make_presets() {
  return Presets<DummyParams, 3>{std::array<PresetEntry<DummyParams>, 3>{{
      {DummyParams{1, 1.5f}},
      {DummyParams{2, 2.5f}},
      {DummyParams{3, 3.5f}},
  }}};
}

// A freshly built container points at the first entry, with prev mirroring it.
inline void test_initial_state() {
  auto p = make_presets();
  HS_EXPECT_EQ(p.get().id, 1);
  HS_EXPECT_NEAR(p.get().value, 1.5f, 1e-6f);
  HS_EXPECT_EQ(p.current_idx, 0);
  HS_EXPECT_EQ(p.prev_idx, 0);
  HS_EXPECT_EQ(static_cast<int>(p.get_entries().size()), 3);
}

// next() advances forward, wraps past the last entry, and leaves prev_get()
// pointing at the entry that was current before the call.
inline void test_next_cycles_forward_and_tracks_prev() {
  auto p = make_presets();

  p.next();
  HS_EXPECT_EQ(p.get().id, 2);
  HS_EXPECT_EQ(p.prev_get().id, 1); // prev_idx points at the pre-next entry

  p.next();
  HS_EXPECT_EQ(p.get().id, 3);
  HS_EXPECT_EQ(p.prev_get().id, 2);

  // Wraps back to the first entry.
  p.next();
  HS_EXPECT_EQ(p.get().id, 1);
  HS_EXPECT_EQ(p.prev_get().id, 3);
}

// prev() steps backward and wraps from the first entry to the last.
inline void test_prev_cycles_backward() {
  auto p = make_presets();

  // prev() from index 0 wraps to the last entry.
  p.prev();
  HS_EXPECT_EQ(p.get().id, 3);
  HS_EXPECT_EQ(p.prev_get().id, 1);

  p.prev();
  HS_EXPECT_EQ(p.get().id, 2);
  HS_EXPECT_EQ(p.prev_get().id, 3);
}

// apply() overwrites the caller's target with a copy of the current entry.
inline void test_apply_copies_current() {
  auto p = make_presets();
  p.next(); // now on id 2

  DummyParams target{0, 0.0f};
  p.apply(target);
  HS_EXPECT_EQ(target.id, 2);
  HS_EXPECT_NEAR(target.value, 2.5f, 1e-6f);
}

// With a single entry, next()/prev() are no-ops that keep current_idx at 0.
inline void test_single_entry_wraps_in_place() {
  Presets<DummyParams, 1> p{std::array<PresetEntry<DummyParams>, 1>{{
      {DummyParams{42, 0.0f}},
  }}};
  p.next();
  HS_EXPECT_EQ(p.get().id, 42);
  HS_EXPECT_EQ(p.current_idx, 0);
  p.prev();
  HS_EXPECT_EQ(p.get().id, 42);
  HS_EXPECT_EQ(p.current_idx, 0);
}

// Class template argument deduction infers the entry type and count from the
// array, so Presets can be constructed without spelling out <DummyParams, 2>.
inline void test_ctad_deduces_size() {
  // The deduction guide deduces <DummyParams, 2> from the array's element count.
  Presets ctad{std::array<PresetEntry<DummyParams>, 2>{{
      {DummyParams{7, 0.0f}},
      {DummyParams{8, 0.0f}},
  }}};
  HS_EXPECT_EQ(static_cast<int>(ctad.get_entries().size()), 2);
  ctad.next();
  HS_EXPECT_EQ(ctad.get().id, 8);
}

// Runs all preset-container cases; returns the module's failure count.
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
