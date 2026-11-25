/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Fail-fast traps and cross-process cold/warm effect determinism.
 * @details Each case runs in a re-executed child. The parent requires the
 * native trap status and a breadcrumb naming the exact guard. */

#include "tests/death/process.h"

namespace hs_test::death_tests {

/**
 * @brief Parent entry point for the death module.
 * @return The module's failure count.
 * @details Spawn-checks the harness, then runs each case of its slice in a
 *          child and asserts each died by the exact trap status. The whole
 *          module and shard 0 also run compare_determinism_folds(). A shard
 *          writes a record instead of reporting coverage; the census merges the
 *          records and reports coverage without spawning children.
 */
inline int run_death_tests() {
  hs_test::ModuleFixture fixture("death");

  ShardSpec shard;
  if (!shard_spec_from_env(shard)) {
    report_unrunnable("HS_DEATH_SHARD is malformed or HS_DEATH_SHARD_DIR is "
                      "unset",
                      0);
    return fixture.result();
  }

  int n;
  const Case *cs = all_cases(n);
  std::vector<int> lines(static_cast<size_t>(n));
  int duplicates = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < i; ++j)
      if (std::strcmp(cs[i].name, cs[j].name) == 0)
        ++duplicates;
  HS_EXPECT_EQ(duplicates, 0);
  if (duplicates != 0)
    return fixture.result();

  if (shard.role == ShardSpec::Role::CENSUS) {
    merge_shard_records(shard, cs, n, lines.data());
    report_guard_coverage(cs, n, lines.data());
    return fixture.result();
  }
  check_shader_chain_colorize_companion();
  check_islamicstars_landing_handoff_companion();
  check_reaction_lattice_arena_companion();
  check_hyperlattice_trace_companion();
  check_catalog_field_curve_companion();
  if (shard.role == ShardSpec::Role::SHARD)
    std::remove(shard_record_path(shard.dir, shard.index).c_str());

  if (!self_exe() || self_exe()[0] == '\0') {
    report_unrunnable("no argv[0] to re-exec", 0);
    return fixture.result();
  }

  HS_EXPECT_EQ(spawn_child("__timeout_check__", 100), -1);

  int fields[2]{};
  HS_EXPECT_TRUE(parse_record_row("0\t12\r\n", fields, 2));
  HS_EXPECT_EQ(fields[1], 12);
  HS_EXPECT_FALSE(parse_record_row("0+12", fields, 2));
  HS_EXPECT_FALSE(parse_record_row("4294967296 12", fields, 2));

  // Control: a child given an unknown case must exit cleanly.
  int control = spawn_child("__spawn_check__");
  if (!child_exited_clean(control)) {
    report_unrunnable("cannot re-exec self", control);
    set_case_env("");
    return fixture.result();
  }

  // The sentinel proves the trap is observable independently of real cases.
  const int probe = spawn_child(SHAPE_PROBE_CASE);
  if (!child_trapped(probe)) {
    report_unrunnable("trap sentinel did not trap; trap status is unobservable",
                      probe);
    set_case_env("");
    return fixture.result();
  }

  HS_EXPECT_EQ(
      breadcrumb_names_guard(
          "HS_CHECK failed: C:\\tree\\core\\render\\sdf\\shapes.h:12: (false) probe\n",
          "core/render/sdf/shapes.h", "(false) probe"),
      12);
  HS_EXPECT_FALSE(breadcrumb_names_guard(
      "HS_CHECK failed: core/other/shapes.h:12: (false) probe\n",
      "core/render/sdf/shapes.h", "(false) probe"));
  HS_EXPECT_FALSE(breadcrumb_names_guard(
      "HS_CHECK failed: notcore/render/sdf/shapes.h:12: (false) probe\n",
      "core/render/sdf/shapes.h", "(false) probe"));

  HS_EXPECT_FALSE(breadcrumb_names_guard(
      "HS_CHECK failed: core/render/sdf/shapes.h:12: (true) other\n",
      "core/render/sdf/shapes.h", "(false) probe"));
  HS_EXPECT_FALSE(breadcrumb_names_guard(
      "HS_CHECK failed: core/render/sdf/shapes.h:12: (false) pro\n",
      "core/render/sdf/shapes.h", "(false) probe"));

  // The same sentinel proves the capture channel.
  if (!breadcrumb_names_guard(child_output(), "tests/death/registry.h",
                              "(false) death-harness trap-shape probe")) {
    report_unrunnable("cannot capture child output; which guard fired is "
                      "unverifiable",
                      0);
    set_case_env("");
    return fixture.result();
  }

#if defined(_WIN32)
  const int literal_exit = spawn_child("__literal_trap_exit__");
  HS_EXPECT_EQ(literal_exit, TRAP_STATUS);
  HS_EXPECT_FALSE(child_trapped(literal_exit));
#endif

  if (shard.role == ShardSpec::Role::WHOLE || shard.index == 0)
    compare_determinism_folds();

#if defined(_WIN32)
  DWORD handles_before = 0;
  HS_EXPECT_TRUE(GetProcessHandleCount(GetCurrentProcess(), &handles_before));
#endif

  std::vector<bool> ran(static_cast<size_t>(n), false);
  for (int i = 0; i < n; ++i) {
    if (!case_enabled(cs[i]) || i % shard.count != shard.index)
      continue;
    ran[static_cast<size_t>(i)] = true;
    int rc = spawn_child(cs[i].name);
    bool trapped = child_trapped(rc);
    // The child must die at this case's guard, not at another trap.
    int line = breadcrumb_names_guard(child_output(), cs[i].guard_file,
                                      cs[i].guard_text);
    bool at_guard =
        line != 0 && (!cs[i].output_detail ||
                      std::strstr(child_output(), cs[i].output_detail));
    lines[i] = trapped ? line : 0;
    HS_EXPECT_TRUE(trapped);
    HS_EXPECT_TRUE(at_guard);
    std::printf("  [%s] trap fires: %-26s (child rc=%d)\n",
                trapped && at_guard ? "ok" : "FAIL", cs[i].name, rc);
    if (!(trapped && at_guard))
      std::printf("      expected %s: %s\n      child logged: %s\n",
                  cs[i].guard_file, cs[i].guard_text, child_output());
  }

#if defined(_WIN32)
  DWORD handles_after = 0;
  HS_EXPECT_TRUE(GetProcessHandleCount(GetCurrentProcess(), &handles_after));
  HS_EXPECT_EQ(handles_after, handles_before);
#endif

  std::remove(child_capture_path());
  set_case_env(""); // leave the env clean for anything that runs after us
  if (shard.role == ShardSpec::Role::SHARD)
    write_shard_record(shard, n, ran, lines);
  else
    report_guard_coverage(cs, n, lines.data());
  return fixture.result();
}

} // namespace hs_test::death_tests
