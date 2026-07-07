// Recompiles the engine with -DHS_TEST_H_OFFSET=3 (set in tests/CMakeLists.txt)
// so the south-pole Y-clip renormalization (only active when hs::H_OFFSET > 0,
// which a normal host build never sets) runs against an energy-conservation
// oracle. Its own TU because offset-3 and offset-0 instantiations of
// PhiLUT<H>/TrigLUT<W,H> would clash under ODR.
#include <cstdio>

#include "core/engine/engine.h"
#include "tests/test_h_offset_renorm.h"

int main() {
  const int failed = hs_test::h_offset_renorm::run_h_offset_renorm_tests();
  const int total = hs_test::stats().passed + hs_test::stats().failed;
  std::printf("=== h_offset_renorm: %d passed, %d failed "
              "(HS_TEST_H_OFFSET=3) ===\n",
              hs_test::stats().passed, hs_test::stats().failed);
  // Floor against silent drift: bump when adding assertions.
  constexpr int MIN_ASSERTIONS = 108;
  if (total < MIN_ASSERTIONS) {
    std::printf("=== h_offset_renorm: only %d assertions ran, expected >= %d "
                "(a check was dropped) ===\n",
                total, MIN_ASSERTIONS);
    return 1;
  }
  return failed ? 1 : 0;
}
