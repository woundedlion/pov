/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Engine built with -DHS_TEST_H_OFFSET=3, so the bottom physical row stops
// short of the south pole; checks Y-clip renormalization at that edge against
// an energy-conservation oracle.
// Separate TU: offset-3 and offset-0 PhiLUT<H>/TrigLUT<W,H> clash under ODR.
#include "core/engine/engine.h"
#include "tests/test_h_offset_renorm.h"

#ifdef HS_OFFSET_FULL_RESOLUTION
#include "effects/Comets.h"
#include "effects/RingSpin.h"
#include "tests/pixel_test_util.h"
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
#endif

int main() {
#ifdef HS_OFFSET_FULL_RESOLUTION
  hs_test::ModuleFixture fixture("offset3 full-resolution roster");
  render_offset_effect<Comets>();
  render_offset_effect<RingSpin>();
  if (fixture.result())
    return 1;
#endif
  return hs_test::h_offset_renorm::run_h_offset_renorm_tests() ? 1 : 0;
}
