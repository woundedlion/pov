/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Roster-wide effect sweeps: smoke, cross-run determinism, clip-clear parity
 * and paused rendering.
 *
 * These passes use no fixed IEEE reference values. Assertions are same-binary
 * comparisons or coarse properties (lit, moving, aliased). This is the effects
 * coverage the -ffast-math -fno-finite-math-only axis can run, which is the flag
 * pair both shipping targets build with (platformio.ini, CMakeLists.txt). The
 * white-box block in tests/test_effects.h checks against fixed references and
 * stays excluded from that axis.
 */
#pragma once

#include "tests/test_effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace effects_smoke_tests {

using namespace hs_test::effects_tests;

constexpr int PAUSED_FRAMES = 4;

template <template <int, int> class E, int W = SMALL_W, int H = SMALL_H>
inline void paused_render_one(const char *name) {
  reset_effect_globals();

  const int frames = PAUSED_FRAMES;
  E<W, H> effect;
  effect.setAnimationsPaused(true);
  effect.init();
  HS_EXPECT_TRUE(effect.animations_paused());
  for (int f = 0; f < frames; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }

  const uint64_t acc = frame_energy<W, H>(effect);

  if (acc == 0)
    std::printf("  PAUSED-BLANK %-20s produced no lit pixel over %d paused "
                "frames @ %dx%d\n",
                name, frames, W, H);
  HS_EXPECT(acc > 0, "effect must produce non-black output while paused");
}

inline void test_every_effect_renders_while_paused() {
  std::printf("  -- paused render, %d frames --\n", PAUSED_FRAMES);
#define HS_PAUSED_ONE(name) paused_render_one<name>(#name);
  HS_EFFECT_LIST(HS_PAUSED_ONE)
#undef HS_PAUSED_ONE
}

/**
 * @brief Module entry point for the roster-wide effect sweeps.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Runs the smoke and determinism passes over every registered effect at
 * the small-aspect resolution, then the clip-clear parity sweep. The FULL tier
 * (HS_EFFECTS_FULL=1; CI on every master push and PR) adds the same smoke and determinism
 * passes at the 288x144 production resolution, which are the bulk of the cost —
 * full-frame software raster over 41,472-pixel frames.
 */
inline int run_effects_smoke_tests() {
  hs_test::ModuleFixture fixture("effects_smoke");

  if (effects_full_suite()) {
    // Full production-resolution roster passes (288x144): smoke, then cross-run
    // determinism under the injected clock.
#define HS_SMOKE_ONE(name) smoke_one<name>(#name);
    HS_EFFECT_LIST(HS_SMOKE_ONE)
#undef HS_SMOKE_ONE
#define HS_DET_ONE(name) determinism_one<name>(#name);
    HS_EFFECT_LIST(HS_DET_ONE)
#undef HS_DET_ONE
  } else {
    std::printf("  [TIER] production-resolution smoke and determinism omitted; "
                "set HS_EFFECTS_FULL=1\n");
  }

  // Small-aspect <96,20> roster passes — always run; this is the QUICK tier's
  // core and the only place that specialization runs under native asserts
  // (see SMALL_W/SMALL_H).
  std::printf("  -- small-aspect resolution %dx%d --\n", SMALL_W, SMALL_H);
#define HS_SMOKE_ONE_SMALL(name) smoke_one<name, SMALL_W, SMALL_H>(#name);
  HS_EFFECT_LIST(HS_SMOKE_ONE_SMALL)
#undef HS_SMOKE_ONE_SMALL
#define HS_DET_ONE_SMALL(name) determinism_one<name, SMALL_W, SMALL_H>(#name);
  HS_EFFECT_LIST(HS_DET_ONE_SMALL)
#undef HS_DET_ONE_SMALL

  // Every effect that does not read outside its display band clip-clears, so the
  // roster is swept rather than the one effect the optimization started on.
  std::printf("  -- clip-clear display parity --\n");
#define HS_CLIP_PARITY_ONE(name)                                               \
  clip_clear_parity_one<name, SMALL_W, SMALL_H>(#name);
  HS_EFFECT_LIST(HS_CLIP_PARITY_ONE)
#undef HS_CLIP_PARITY_ONE

  test_every_effect_renders_while_paused();
  if (effects_full_suite()) {
#define HS_PAUSED_FULL(name)                                                   \
  paused_render_one<name, DEFAULT_W, DEFAULT_H>(#name);
    HS_EFFECT_LIST(HS_PAUSED_FULL)
#undef HS_PAUSED_FULL
  } else {
    std::printf(
        "  [TIER] production-resolution paused renders omitted; set HS_EFFECTS_FULL=1\n");
  }

  return fixture.result();
}

} // namespace effects_smoke_tests
} // namespace hs_test
