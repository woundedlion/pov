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
#include "core/render/canvas.h"
#include "tests/test_harness.h"

#include <cstdlib>

namespace hs_test {

/**
 * @brief Default per-effect frame count for every roster sweep.
 * @details Kept small so the local pre-commit gate stays quick. 8 frames never
 * reaches the long, cyclic code paths (effect morph cycles, particle/trail
 * wraps, arena compaction, and the effect-lifecycle transitions — RingShower
 * slot reuse, Thrusters fire/FIFO expiry, ShapeShifter's 48-frame cut), so CI
 * sets HS_SMOKE_FRAMES=120 (.github/workflows/ci.yml) to exercise them on every
 * push/PR. The effects smoke/determinism passes and the arena/stack budget
 * gates all resolve their window through smoke_frames(), so the budgets are
 * measured over the same window the correctness sweep covers.
 */
constexpr int DEFAULT_SMOKE_FRAMES = 8;

/**
 * @brief Resolves the per-effect frame count from the environment.
 * @return HS_SMOKE_FRAMES if set to a positive int, else DEFAULT_SMOKE_FRAMES.
 */
inline int smoke_frames() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (const char *e = std::getenv("HS_SMOKE_FRAMES")) {
#pragma clang diagnostic pop
    int n = std::atoi(e);
    if (n > 0)
      return n;
  }
  return DEFAULT_SMOKE_FRAMES;
}

/**
 * @brief Concrete Effect that draws nothing, for tests that only need a Canvas.
 * @details A fresh effect's buffer_free() is true, so a Canvas built over it
 * does not spin in its constructor. Shows no background, so the canvas starts
 * black and holds only what the test explicitly plots.
 */
struct StubEffect : public Effect {
  /**
   * @brief Constructs the stub at the given canvas resolution.
   * @param w Canvas width in pixels.
   * @param h Canvas height in pixels.
   */
  StubEffect(int w, int h) : Effect(w, h) {}
  /**
   * @brief Per-frame draw hook; a no-op for every user of this fixture.
   */
  void draw_frame() override {}
};

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
  explicit ModuleFixture(const char *name)
      : scope((reset_globals(), begin_module(name))) {}

  /**
   * @brief Closes the module scope and returns its failure count.
   * @return The module's failure count (delta since construction).
   */
  int result() const { return end_module(scope); }
};

} // namespace hs_test
