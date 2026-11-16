/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "tests/death/context.h"

namespace hs_test::death_tests {

/** @file
 * @brief Callables death cases.
 */

// Callables death cases.

/**
 * @brief Death case: calling an empty (default-constructed) Fn must trap.
 */
inline void case_empty_fn_call() {
  Fn<int(int), 16> f;
  if (opaque(false))
    f = [](int x) { return x; };
  int v = f(opaque(7)); // empty invoke -> check_fail -> trap
  (void)v;
}

/**
 * @brief Death case: invoking an empty FunctionRef must trap.
 * @details Concepts surface — the empty state's thunk diverges through
 *          function_ref_empty_call rather than calling through a null
 *          context.
 */
inline void case_empty_function_ref_call() {
  FunctionRef<int(int)> f;
  if (opaque(false))
    f = [](int x) { return x; };
  int v = f(opaque(7)); // empty invoke -> check_fail -> trap
  (void)v;
}

} // namespace hs_test::death_tests
