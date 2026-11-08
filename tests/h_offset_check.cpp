/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Engine built with -DHS_TEST_H_OFFSET=3, so the bottom physical row stops
// short of the south pole; renders an effect roster at that offset and checks
// Y-clip renormalization at that edge against an energy-conservation oracle.
// Separate TU: offset-3 and offset-0 PhiLUT<H>/TrigLUT<W,H> clash under ODR.
#include "core/engine/engine.h"
#include "effects/Comets.h"
#include "effects/RingSpin.h"
#include "tests/pixel_test_util.h"
#include "tests/test_h_offset_renorm.h"
#include <memory>

template <template <int, int> class E> void render_offset_effect() {
  hs_test::reset_globals();
  auto effect = std::make_unique<E<96, 20>>();
  effect->init();
  for (int frame = 0; frame < hs_test::smoke_frames(); ++frame) {
    HS_CONTEXT("frame", frame);
    hs_test::pin_frame_clock(frame);
    effect->draw_frame();
    effect->advance_display();
    if (frame > 0)
      HS_EXPECT_GT((hs_test::frame_energy<96, 20>(*effect)), uint64_t{0});
  }
}

int main() {
  if (!hs_test::require_ci_smoke_frames())
    return 1;
  int failed = 0;
  {
    hs_test::ModuleFixture fixture("offset3 full-resolution roster");
    render_offset_effect<Comets>();
    render_offset_effect<RingSpin>();
    failed |= fixture.result();
  }
  failed |= hs_test::h_offset_renorm::run_h_offset_renorm_tests();
  return failed ? 1 : 0;
}
