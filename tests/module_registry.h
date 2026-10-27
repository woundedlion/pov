/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Module registry contract for native test executables.
 */
#pragma once
#include <span>

namespace hs_test {
/** @brief Named module entry point and CI depth policy. */
struct TestModule {
  const char *name;  /**< Module name accepted by the CLI. */
  int (*run)();      /**< Returns the module failure count. */
  bool effects_tier; /**< Requires the effects CI tier levers. */
};

/** @brief Executable-owned modules and optional death-harness hooks. */
struct TestSuite {
  std::span<const TestModule> modules; /**< Ordered module roster. */
  void (*set_death_self_exe)(
      const char *); /**< Records the death-child executable. */
  void (*run_death_child_case)(const char *); /**< Runs a named death child. */
};

/** @brief Returns the executable's module registry. */
const TestSuite &test_suite();
} // namespace hs_test
