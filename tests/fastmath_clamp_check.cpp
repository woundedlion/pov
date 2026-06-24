// Compiles the same NaN-clamp assertions as run_tests but under the WASM release
// math flags (-O3 -ffast-math -fno-finite-math-only, set in tests/CMakeLists.txt)
// so the hs::clamp NaN->hi contract is verified against the codegen that ships,
// not only the default-IEEE-semantics build.
#include <cstdio>

#include "core/engine.h"
#include "tests/test_color.h"

int main() {
  using namespace hs_test::color_tests;

  test_blend_alpha_clamps_before_cast();
  test_pixel16_scale_clamps_before_cast();
  test_gradient_get_clamps_out_of_range();
  test_mobius_longitude_singularity_saturates_to_endpoint();

  const int failed = hs_test::stats().failed;
  std::printf("=== fastmath_clamp: %d passed, %d failed (-ffast-math "
              "-fno-finite-math-only) ===\n",
              hs_test::stats().passed, failed);
  return failed ? 1 : 0;
}
