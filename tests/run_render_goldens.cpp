/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "tests/test_plot_parity.h"
int main() {
  hs_test::ModuleFixture fixture("render_goldens");
  hs_test::plot_parity_tests::test_platform_fastmath_golden_reference();
  hs_test::plot_parity_tests::
      test_render_fingerprint_rejects_channel_and_spatial_changes();
  return (fixture.result() != 0 ||
          (hs_test::stats().skipped > 0 && hs_test::skips_are_errors()))
             ? 1
             : 0;
}
