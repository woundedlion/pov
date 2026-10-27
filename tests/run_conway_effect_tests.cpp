/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "tests/vec_test_util.h"

#include "tests/test_conway_morph.h"
#include "tests/test_conway_continuity.h"
#include "tests/test_conway_soak.h"
#include "tests/test_opchain_probe.h"
#include "tests/test_opchain_arena_survey.h"

namespace hs_test::module_runners {

HS_FORWARD_MODULE(conway_morph_tests, run_conway_morph_tests)
HS_FORWARD_MODULE(conway_continuity_tests, run_conway_continuity_tests)
HS_FORWARD_MODULE(conway_soak_tests, run_conway_soak_tests)
HS_FORWARD_MODULE(opchain_probe_tests, run_opchain_probe_tests)
HS_FORWARD_MODULE(opchain_arena_survey_tests, run_opchain_arena_survey_tests)

} // namespace hs_test::module_runners
