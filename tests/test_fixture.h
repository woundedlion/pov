/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Per-module fixture that resets process-global state (arena split, shared
 * Timeline, RNG, mock clock) to a known baseline.
 */
#pragma once

#include "core/animation/animation.h"
#include "core/memory.h"
#include "core/platform/platform.h"
#include "core/render/canvas.h"
#include "core/render/render_policy.h"
#include "tests/test_harness.h"

#include <cmath>
#include <cstdlib>

namespace hs_test {

/**
 * @brief Default per-effect frame count returned by smoke_frames().
 * @details Overridden by HS_SMOKE_FRAMES; CI requires at least
 * CI_MIN_SMOKE_FRAMES.
 */
constexpr int DEFAULT_SMOKE_FRAMES = 8;
constexpr int CI_MIN_SMOKE_FRAMES = 120;

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

/** @brief Rejects a shallow roster window when running under CI. */
inline bool require_ci_smoke_frames() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *ci = std::getenv("CI");
  const char *frames = std::getenv("HS_SMOKE_FRAMES");
#pragma clang diagnostic pop
  if (!ci || ci[0] == '\0' ||
      (frames && std::atoi(frames) >= CI_MIN_SMOKE_FRAMES))
    return true;
  std::fprintf(stderr,
               "CI=on but HS_SMOKE_FRAMES is unset or below %d — "
               "a shallower window skips frame-cyclic paths and arms no "
               "preset transition. Set HS_SMOKE_FRAMES in the workflow "
               "step's env.\n",
               CI_MIN_SMOKE_FRAMES);
  return false;
}

/**
 * @brief Per-frame clock advance in milliseconds for the roster sweeps (~30fps).
 * @details Fixed cadence, so clock-driven animation reproduces exactly.
 */
constexpr unsigned long FRAME_MS = 33;
/**
 * @brief Per-frame clock advance in microseconds, paired with FRAME_MS.
 */
constexpr unsigned long FRAME_US = 33000;

/**
 * @brief Pins the mock clock to frame @p f of the canonical sweep cadence.
 * @param f Zero-based frame index; f = 0 is the pre-init epoch.
 */
inline void pin_frame_clock(int f) {
  hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                    static_cast<unsigned long>(f) * FRAME_US);
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
   * @brief Per-frame draw hook; a no-op.
   */
  void draw_frame() override {}
};

/**
 * @brief Resets the canonical process-global state to a known baseline.
 * @details No Timeline may be live at the call site.
 */
inline void reset_globals() {
  configure_arenas_default();
  Timeline().clear();
  hs::random().seed(1337u);
  hs::clear_mock_time();
  Render::pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
  HS_SCAN_METRIC(hs::g_scan_metrics.reset());
}

/**
 * @brief Draws one uniform float in [lo, hi] from a locally seeded generator.
 * @param rng Generator private to the test, so the draw stream is reproducible
 * without disturbing the process-wide hs::random().
 * @param lo Lower bound (inclusive).
 * @param hi Upper bound (inclusive through rounding).
 * @return A float in [lo, hi], bit-identical on every platform.
 * @details Draw one component at a time: argument evaluation order is
 * unspecified, so a multi-argument call reorders the stream per compiler.
 */
inline float rand_uniform(hs::Pcg32 &rng, float lo, float hi) {
  return lo + hs::random_to_unit(rng(), hs::Pcg32::max()) * (hi - lo);
}

/**
 * @brief Draws a cube-sampled, normalized random unit vector (not uniform).
 * @return A unit Vector; draws inside a 0.1-radius ball are rejected and
 *         redrawn, so the normalize is well conditioned.
 */
inline math::Vector rand_unit() {
  for (;;) {
    const float rx = hs::rand_f(-1, 1);
    const float ry = hs::rand_f(-1, 1);
    const float rz = hs::rand_f(-1, 1);
    math::Vector r(rx, ry, rz);
    if (r.length() > 0.1f)
      return r.normalized();
  }
}

/**
 * @brief Module scope that resets canonical global state on entry.
 * @details Construct at the top of a module's run_*_tests() and return
 * result().
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

/**
 * @brief Expects every parameter finite and inside its registered range, with
 * option params set to one of their option values.
 * @param effect Effect whose current parameter values are checked.
 * @param label Context label for failures.
 */
template <typename FX>
inline void expect_params_in_range(const FX &effect, const char *label) {
  HS_CONTEXT(label);
  for (const auto &def : effect.getParameters()) {
    HS_CONTEXT(def.name);
    const float v = def.get();
    HS_EXPECT_TRUE(std::isfinite(v));
    HS_EXPECT_GE(v, def.min);
    HS_EXPECT_LE(v, def.max);
    if (def.option_count > 0) {
      HS_EXPECT_EQ(v, std::floor(v));
      if (def.option_values) {
        bool listed = false;
        for (int i = 0; i < def.option_count; ++i)
          listed |= v == static_cast<float>(def.option_values[i]);
        HS_EXPECT_TRUE(listed);
      } else {
        HS_EXPECT_LT(v, static_cast<float>(def.option_count));
      }
    }
  }
}

} // namespace hs_test
