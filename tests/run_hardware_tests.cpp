/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_dma_core.h"
#include "tests/test_hd107s_frame.h"
#include "tests/test_dma_controller.h"
#include "tests/test_pov_segmented.h"
#include "tests/test_pov_single.h"
#include "tests/test_pov_sync.h"

namespace hs_test::module_runners {

int run_dma_core_tests() {
  return hs_test::dma_core_tests::run_dma_core_tests();
}
int run_hd107s_tests() { return hs_test::hd107s_tests::run_hd107s_tests(); }
int run_dma_controller_tests() {
  return hs_test::dma_controller_tests::run_dma_controller_tests();
}
int run_pov_segmented_tests() {
  return hs_test::pov_segmented_tests::run_pov_segmented_tests();
}
int run_pov_single_tests() {
  return hs_test::pov_single_tests::run_pov_single_tests();
}
int run_pov_sync_tests() {
  return hs_test::pov_sync_tests::run_pov_sync_tests();
}

} // namespace hs_test::module_runners
