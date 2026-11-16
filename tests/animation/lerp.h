/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Type-erased Lerp animation tests.
 */
#pragma once

#include "core/animation/animation.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/animation/canvas_fixture.h"
#include "tests/animation/lerp_fixture.h"

namespace hs_test::animation_tests {

// ============================================================================
// Lerp (type-erased subject.lerp(start, target, t))
// ============================================================================

/**
 * @brief Verifies Lerp drives its subject from start to target over `duration`
 * frames, landing on the target and reporting done().
 */
inline void test_lerp_drives_subject_to_target() {
  Lerpable subject, start, target;
  start.value = 0.0f;
  target.value = 200.0f;
  const int DURATION = 10;
  Animation::Lerp l(subject, start, target, DURATION, math::ease_linear);
  for (int i = 0; i < DURATION; ++i) {
    l.step(fake_canvas());
  }
  HS_EXPECT_NEAR(subject.value, 200.0f, 1e-3f);
  HS_EXPECT_TRUE(l.done());
}

/**
 * @brief Verifies at the halfway eased progress the subject sits exactly
 * between start and target.
 */
inline void test_lerp_midpoint() {
  Lerpable subject, start, target;
  start.value = 10.0f;
  target.value = 20.0f;
  const int DURATION = 4;
  Animation::Lerp l(subject, start, target, DURATION, math::ease_linear);
  l.step(fake_canvas()); // t=1 -> progress 0.25 -> 12.5
  l.step(fake_canvas()); // t=2 -> progress 0.50 -> 15.0
  HS_EXPECT_NEAR(subject.value, 15.0f, 1e-3f);
}

} // namespace hs_test::animation_tests
