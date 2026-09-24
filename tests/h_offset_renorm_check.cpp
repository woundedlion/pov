/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Recompiles the engine with -DHS_TEST_H_OFFSET=3 (set in tests/CMakeLists.txt)
// so the south-pole Y-clip renormalization (only active when hs::H_OFFSET > 0,
// which a normal host build never sets) runs against an energy-conservation
// oracle. Its own TU because offset-3 and offset-0 instantiations of
// PhiLUT<H>/TrigLUT<W,H> would clash under ODR.
#include "core/engine/engine.h"
#include "tests/test_h_offset_renorm.h"

#ifdef HS_OFFSET_DEVICE_RESOLUTION
#include "effects/Comets.h"
#include "effects/RingSpin.h"
#include <memory>

template <template <int, int> class E> void render_offset_effect() {
  hs_test::reset_globals();
  auto effect = std::make_unique<E<96, 20>>();
  effect->init();
  for (int frame = 0; frame < hs_test::smoke_frames(); ++frame) {
    HS_CONTEXT("frame", frame);
    uint64_t energy = 0;
    hs_test::pin_frame_clock(frame);
    effect->draw_frame();
    effect->advance_display();
    for (int y = 0; y < 20; ++y)
      for (int x = 0; x < 96; ++x) {
        const Pixel &pixel = effect->get_pixel(x, y);
        energy += pixel.r + pixel.g + pixel.b;
      }
    if (frame > 0)
      HS_EXPECT_GT(energy, uint64_t{0});
  }
}
#endif

int main() {
#ifdef HS_OFFSET_DEVICE_RESOLUTION
  hs_test::ModuleFixture fixture("offset3 device roster");
  render_offset_effect<Comets>();
  render_offset_effect<RingSpin>();
  if (fixture.result())
    return 1;
#endif
  return hs_test::h_offset_renorm::run_h_offset_renorm_tests() ? 1 : 0;
}
