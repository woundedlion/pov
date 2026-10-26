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

HS_FORWARD_MODULE(effects_tests, run_effects_tests)
HS_FORWARD_MODULE(mindsplatter_tests, run_mindsplatter_tests)
HS_FORWARD_MODULE(effects_smoke_tests, run_effects_smoke_tests)
HS_FORWARD_MODULE(effect_factory_tests, run_effect_factory_tests)
HS_FORWARD_MODULE(kaleidoscope_melt_tests, run_kaleidoscope_melt_tests)
HS_FORWARD_MODULE(kaleidoscope_smooth_tests, run_kaleidoscope_smooth_tests)
HS_FORWARD_MODULE(composed_effect_tests, run_composed_effect_tests)
HS_FORWARD_MODULE(shapeshifter_oracle_tests, run_shapeshifter_oracle_tests)
HS_FORWARD_MODULE(shapeshifter_tiles_tests, run_shapeshifter_tiles_tests)

bool effects_full_suite() {
  return hs_test::effects_tests::effects_full_suite();
}

} // namespace hs_test::module_runners
