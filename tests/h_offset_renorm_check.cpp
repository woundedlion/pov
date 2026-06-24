// Recompiles the engine with -DHS_TEST_H_OFFSET=3 (set in tests/CMakeLists.txt)
// so the south-pole Y-clip renormalization (only active when hs::H_OFFSET > 0,
// which a normal host build never sets) runs against an energy-conservation
// oracle. Its own TU because offset-3 and offset-0 instantiations of
// PhiLUT<H>/TrigLUT<W,H> would clash under ODR.
#include <cstdio>

#include "core/engine.h"
#include "tests/test_h_offset_renorm.h"

int main() {
  const int failed = hs_test::h_offset_renorm::run_h_offset_renorm_tests();
  std::printf("=== h_offset_renorm: %d passed, %d failed "
              "(HS_TEST_H_OFFSET=3) ===\n",
              hs_test::stats().passed, hs_test::stats().failed);
  return failed ? 1 : 0;
}
