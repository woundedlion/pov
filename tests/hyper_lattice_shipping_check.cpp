/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/engine/engine.h"
#include "tests/test_hyper_lattice.h"

int main() {
  return hs_test::hyper_lattice_tests::run_hyper_lattice_tests() ? 1 : 0;
}
