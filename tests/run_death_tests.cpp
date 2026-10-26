/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_death.h"

namespace hs_test::module_runners {

int run_death_tests() { return hs_test::death_tests::run_death_tests(); }

void set_death_self_exe(const char *path) { death_tests::set_self_exe(path); }

void run_death_child_case(const char *name) {
  death_tests::run_child_case(name);
}

} // namespace hs_test::module_runners
