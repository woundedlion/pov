// Finding 369: the hs::clamp NaN->hi contract that the engine's float->int
// domain safety depends on is regression-tested in run_tests under the native
// toolchain's DEFAULT IEEE semantics (no math flags), where the contract holds
// trivially — the at-risk codegen is never exercised. The WASM release build
// runs under -O3 -ffast-math -fno-finite-math-only, and the whole safety story
// rests on -fno-finite-math-only preserving NaN/Inf so the saturating clamp
// guard survives. This standalone executable compiles the SAME NaN-clamp
// assertions under those exact release math flags (set in tests/CMakeLists.txt)
// so the contract is verified against the codegen that actually ships, not only
// the default-semantics build. Paired with the compile-time #error guard in
// core/platform.h, which fails the build if -fno-finite-math-only is ever lost.
#include <cstdio>

#include "core/engine.h"
#include "tests/test_color.h"

int main() {
  using namespace hs_test::color_tests;

  // The clamp-before-cast and singularity-saturation assertions, run under the
  // real release math flags. Each verifies a possibly-NaN/Inf value folds to a
  // saturating endpoint instead of invoking float->int cast UB.
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
