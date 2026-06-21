// The south-pole Y-clip renormalization in the rasterizer (Screen::AntiAlias's
// bilinear-tap fold, plus the H_VIRT handling threaded through geometry/scan/
// plot/sdf) only does real work when hs::H_OFFSET > 0 — which a normal host
// build never sets, so the device's most numerically subtle hot-loop path ships
// with zero automated coverage on the one platform (Teensy) that has no debugger
// and no console. This standalone executable recompiles the SAME engine code
// with -DHS_TEST_H_OFFSET=3 (set in tests/CMakeLists.txt), so the whole pipeline
// is built with the hardware offset and the renorm runs against an energy-
// conservation oracle. Same recompile-under-device-config tactic as
// fastmath_clamp_check; it lives in its own TU because an offset-3 and an
// offset-0 instantiation of PhiLUT<H>/TrigLUT<W,H> would clash under ODR.
#include <cstdio>

#include "core/effects_engine.h"
#include "tests/test_h_offset_renorm.h"

int main() {
  const int failed = hs_test::h_offset_renorm::run_h_offset_renorm_tests();
  std::printf("=== h_offset_renorm: %d passed, %d failed "
              "(HS_TEST_H_OFFSET=3) ===\n",
              hs_test::stats().passed, hs_test::stats().failed);
  return failed ? 1 : 0;
}
