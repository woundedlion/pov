/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Per-module fixture for the canonical process-global state. Tests mutate
 * process-wide singletons (the arena split, the shared Timeline event array and
 * frame cursor, the RNG); reset_globals() restores them to a known baseline so
 * module independence does not rest on every test hand-restoring what it
 * touched. ModuleFixture pairs the harness scope with that reset on entry.
 */
#pragma once

#include "core/animation/animation.h"
#include "core/engine/memory.h"
#include "core/engine/platform.h"
#include "tests/test_harness.h"

namespace hs_test {

/**
 * @brief Resets the canonical process-global state to a known baseline.
 * @details Restores the default arena split, clears the shared Timeline (events
 * and frame cursor), and reseeds the RNG to the device-matching seed. No
 * Timeline may be live at the call site: the temporary clears global event state
 * through the singleton's clear() the same way the per-test sites do.
 */
inline void reset_globals() {
  configure_arenas_default();
  Timeline().clear();
  hs::random().seed(1337u);
}

/**
 * @brief Module scope that resets canonical global state on entry.
 * @details Constructed at the top of a module's run_*_tests(); resets the shared
 * singletons before any test runs and forwards begin_module/end_module so the
 * module still reports its pass/fail delta. Use result() as the module return.
 */
struct ModuleFixture {
  ModuleScope scope; /**< Underlying harness scope for the pass/fail delta. */

  /**
   * @brief Resets globals, then opens the named module scope.
   * @param name Module name echoed in the header and footer.
   */
  explicit ModuleFixture(const char *name) : scope((reset_globals(),
                                                    begin_module(name))) {}

  /**
   * @brief Closes the module scope and returns its failure count.
   * @return The module's failure count (delta since construction).
   */
  int result() const { return end_module(scope); }
};

} // namespace hs_test
