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

int main() {
  return hs_test::h_offset_renorm::run_h_offset_renorm_tests() ? 1 : 0;
}
