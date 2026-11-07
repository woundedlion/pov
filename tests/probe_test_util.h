/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Captured-failure probes that restore the surrounding test harness. */
#pragma once
#include "tests/test_harness.h"
#include "tests/fd_capture_util.h"

namespace hs_test {
/** @brief Captures a probe's failures and diagnostics, restoring harness state. */
inline void expect_probe_failures(auto &&probe, int expected,
                                  const char *diagnostic) {
  const auto saved = stats();
  const auto saved_budget = fail_print_budget();
  fail_print_budget() = {};
  int failures = 0;
  const auto output = capture_stdout([&] {
    probe();
    failures = stats().failed - saved.failed;
  });
  stats() = saved;
  fail_print_budget() = saved_budget;
  HS_EXPECT_TRUE(output.has_value());
  HS_EXPECT_EQ(failures, expected);
  if (output && diagnostic)
    HS_EXPECT_NE(output->find(diagnostic), std::string::npos);
}

} // namespace hs_test
