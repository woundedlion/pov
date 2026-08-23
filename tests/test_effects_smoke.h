/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Roster-wide effect sweeps: smoke, cross-run determinism, clip-clear parity.
 *
 * Separate from the effects module because these three passes carry no IEEE
 * reference value — every assertion compares two renders produced by the SAME
 * binary, so the float configuration cancels. That makes this the one effects
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

/**
 * @brief Module entry point for the roster-wide effect sweeps.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Runs the smoke and determinism passes over every registered effect at
 * the small-aspect resolution, then the clip-clear parity sweep. The FULL tier
 * (HS_EFFECTS_FULL=1; CI on every master push) adds the same smoke and determinism
 * passes at the 288x144 production resolution, which are the bulk of the cost —
 * full-frame software raster over 41,472-pixel frames.
 */
inline int run_effects_smoke_tests() {
  hs_test::ModuleFixture fixture("effects_smoke");

  // SSOT anti-drift guard, mirroring the WASM startup check
  // (targets/wasm/engine_bindings.h): the self-registering effect count (each
  // header's
  // REGISTER_EFFECT) must equal the static HS_EFFECT_LIST roster, or an effect
  // present in one and missing from the other silently drops smoke coverage
  // below. Active because run_tests enables the effect registry.
  HS_EXPECT_EQ(EffectRegistry::entries().size(),
               static_cast<size_t>(HS_EFFECT_COUNT));

  if (effects_full_suite()) {
    // Full production-resolution roster passes (288x144): smoke, then cross-run
    // determinism under the injected clock.
#define HS_SMOKE_ONE(name) smoke_one<name>(#name);
    HS_EFFECT_LIST(HS_SMOKE_ONE)
#undef HS_SMOKE_ONE
#define HS_DET_ONE(name) determinism_one<name>(#name);
    HS_EFFECT_LIST(HS_DET_ONE)
#undef HS_DET_ONE
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

  return fixture.result();
}

} // namespace effects_smoke_tests
} // namespace hs_test
