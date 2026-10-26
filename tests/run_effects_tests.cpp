/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_effects.h"
#include "tests/test_mindsplatter.h"
#include "tests/test_effects_smoke.h"
#include "tests/test_effect_factory.h"
#include "tests/test_kaleidoscope_melt.h"
#include "tests/test_kaleidoscope_smooth.h"
#include "tests/test_composed_effect.h"
#include "tests/test_shapeshifter_oracle.h"
#include "tests/test_shapeshifter_tiles.h"

namespace hs_test::module_runners {

int run_effects_tests() { return hs_test::effects_tests::run_effects_tests(); }
int run_mindsplatter_tests() {
  return hs_test::mindsplatter_tests::run_mindsplatter_tests();
}
int run_effects_smoke_tests() {
  return hs_test::effects_smoke_tests::run_effects_smoke_tests();
}
int run_effect_factory_tests() {
  return hs_test::effect_factory_tests::run_effect_factory_tests();
}
int run_kaleidoscope_melt_tests() {
  return hs_test::kaleidoscope_melt_tests::run_kaleidoscope_melt_tests();
}
int run_kaleidoscope_smooth_tests() {
  return hs_test::kaleidoscope_smooth_tests::run_kaleidoscope_smooth_tests();
}
int run_composed_effect_tests() {
  return hs_test::composed_effect_tests::run_composed_effect_tests();
}
int run_shapeshifter_oracle_tests() {
  return hs_test::shapeshifter_oracle_tests::run_shapeshifter_oracle_tests();
}
int run_shapeshifter_tiles_tests() {
  return hs_test::shapeshifter_tiles_tests::run_shapeshifter_tiles_tests();
}

bool effects_full_suite() {
  return hs_test::effects_tests::effects_full_suite();
}

} // namespace hs_test::module_runners
