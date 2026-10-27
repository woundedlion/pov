/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_param_marshal.h"
#include "tests/test_effect_factory.h"
#include "tests/test_kaleidoscope_melt.h"
#include "tests/test_kaleidoscope_smooth.h"
#include "tests/test_composed_effect.h"

namespace hs_test::module_runners {

HS_FORWARD_MODULE(param_marshal_tests, run_param_marshal_tests)
HS_FORWARD_MODULE(effect_factory_tests, run_effect_factory_tests)
HS_FORWARD_MODULE(kaleidoscope_melt_tests, run_kaleidoscope_melt_tests)
HS_FORWARD_MODULE(kaleidoscope_smooth_tests, run_kaleidoscope_smooth_tests)
HS_FORWARD_MODULE(composed_effect_tests, run_composed_effect_tests)

} // namespace hs_test::module_runners
