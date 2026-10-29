/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Shared native module CLI and suite policy.
 */
#include "tests/module_registry.h"
#include "tests/test_fixture.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

using hs_test::test_suite;

/**
 * @brief Prints the roster's module names, one indented per line.
 * @param out Destination stream (e.g. stdout for --list, stderr on error).
 */
static void print_modules(std::FILE *out) {
  for (const hs_test::TestModule &m : test_suite().modules)
    std::fprintf(out, "  %s\n", m.name);
}

/**
 * @brief Reports whether this invocation runs an effects module.
 * @param argc Argument count as passed to main.
 * @param argv Argument vector as passed to main; names the modules to run.
 * @return True for an unfiltered run, or a filtered run naming an
 * effects-tier module.
 */
static bool runs_effects(int argc, char **argv) {
  if (argc <= 1) {
    for (const hs_test::TestModule &module : test_suite().modules)
      if (module.effects_tier)
        return true;
    return false;
  }
  for (int i = 1; i < argc; ++i)
    for (const hs_test::TestModule &module : test_suite().modules)
      if (module.effects_tier && std::strcmp(argv[i], module.name) == 0)
        return true;
  return false;
}

/**
 * @brief Watchdog bound every CI leg has to set, in microseconds.
 * @details A buffer_free() timeout traps the process instead of failing a
 * test, so CI raises the shipping 2 s bound.
 */
constexpr unsigned long CI_MIN_BUFFER_FREE_WATCHDOG_US = 30000000UL;

/**
 * @brief Verifies the environment carries the depth levers a CI run must set.
 * @param effects_invocation Whether this invocation runs an effects module.
 * @return 0 when every lever this invocation needs is set, else 1.
 * @details Unset levers default to the shallow local tier, so under CI an
 * unset or blank lever is a failure. HS_EFFECTS_FULL and
 * HS_REQUIRE_EFFECTS_FULL are required only when an effects module runs.
 */
static int check_ci_levers(bool effects_invocation) {
  if (!hs_test::runs_in_ci())
    return 0;

  int missing = 0;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *effects_full = std::getenv("HS_EFFECTS_FULL");
  const char *require_effects_full = std::getenv("HS_REQUIRE_EFFECTS_FULL");
  const char *watchdog = std::getenv("HS_BUFFER_FREE_WATCHDOG_US");
  const char *skip_policy = std::getenv("HS_SKIPS_ARE_ERRORS");
#pragma clang diagnostic pop
  if (!skip_policy || skip_policy[0] == '\0') {
    std::fprintf(stderr,
                 "run_tests: CI=on but HS_SKIPS_ARE_ERRORS carries no explicit "
                 "value — a case that retires itself asserts nothing and still "
                 "reports green. Set it in the workflow step's env.\n");
    ++missing;
  }
  if (!hs_test::require_ci_smoke_frames())
    ++missing;
  if (!watchdog ||
      std::strtoul(watchdog, nullptr, 10) < CI_MIN_BUFFER_FREE_WATCHDOG_US) {
    std::fprintf(stderr,
                 "run_tests: CI=on but HS_BUFFER_FREE_WATCHDOG_US is unset or "
                 "below %lu — the Canvas buffer_free() spin traps rather than "
                 "failing, so a runner deschedule kills the shard instead of "
                 "reporting a test. Set HS_BUFFER_FREE_WATCHDOG_US in the "
                 "workflow step's env.\n",
                 CI_MIN_BUFFER_FREE_WATCHDOG_US);
    ++missing;
  }
  for (const char *name :
       {"HS_RING_TEMPORAL_SCENARIO", "HS_SHAPESHIFTER_SURFACE_CENSUS",
        "HS_SURFACE_LOBE_DIAGNOSTICS", "HS_SHAPESHIFTER_PRESET_HASHES",
        "HS_SHAPESHIFTER_IMAGES", "HS_SHAPESHIFTER_SURFACE_ONLY",
        "HS_SURFACE_DARK_FRAMES"}) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    const char *value = std::getenv(name);
#pragma clang diagnostic pop
    if (value != nullptr) {
      std::fprintf(
          stderr, "run_tests: CI=on but %s narrows the run because it is set\n",
          name);
      ++missing;
    }
  }
  if (effects_invocation) {
    for (const char *name :
         {"HS_RING_TEMPORAL_ONLY", "HS_MOBIUSRINGS_ONLY", "HS_DYNAMO_ONLY",
          "HS_ISLAMICSTARS_ONLY", "HS_MESH_CAPTURE_ONLY"}) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
      const char *value = std::getenv(name);
#pragma clang diagnostic pop
      if (value && value[0] != '\0') {
        std::fprintf(stderr,
                     "run_tests: CI=on but %s narrows the run because "
                     "it is set\n",
                     name);
        ++missing;
      }
    }
    const bool tier_declared = effects_full && effects_full[0] != '\0';
    const bool policy_declared =
        require_effects_full && require_effects_full[0] != '\0';
    if (!tier_declared || !policy_declared) {
      std::fprintf(
          stderr,
          "run_tests: CI=on and an effects module is running, but "
          "HS_EFFECTS_FULL and HS_REQUIRE_EFFECTS_FULL are not both set to an "
          "explicit value — the effects tier would be chosen by omission. Set "
          "both in the workflow step's env.\n");
      ++missing;
    } else if (std::atoi(require_effects_full) > 0 &&
               !hs_test::effects_full_suite()) {
      std::fprintf(
          stderr,
          "run_tests: CI=on and HS_REQUIRE_EFFECTS_FULL is on, but "
          "the effects modules would run the QUICK tier, dropping the "
          "288x144 roster passes and the FULL-tier white-box cases. Set "
          "HS_EFFECTS_FULL=1 in the workflow step's env.\n");
      ++missing;
    }
  }
  return missing ? 1 : 0;
}

/**
 * @brief Verifies an expected module set on argv matches the roster exactly.
 * @param argc Count of trailing names (argv past --check-modules).
 * @param argv The expected module names.
 * @return 0 if the set matches test_suite().modules exactly, 3 on any divergence.
 * @details Every test_suite().modules name must appear in argv, and every argv name must be
 * a roster module.
 */
static int check_modules(int argc, char **argv) {
  int mismatches = 0;
  for (const hs_test::TestModule &m : test_suite().modules) {
    bool found = false;
    for (int i = 0; i < argc; ++i)
      if (std::strcmp(m.name, argv[i]) == 0) {
        found = true;
        break;
      }
    if (!found) {
      std::fprintf(stderr,
                   "run_tests: roster module '%s' missing from expected "
                   "list\n",
                   m.name);
      ++mismatches;
    }
  }
  for (int i = 0; i < argc; ++i) {
    bool found = false;
    for (const hs_test::TestModule &m : test_suite().modules)
      if (std::strcmp(m.name, argv[i]) == 0) {
        found = true;
        break;
      }
    if (!found) {
      std::fprintf(stderr,
                   "run_tests: expected name '%s' is not a roster module\n",
                   argv[i]);
      ++mismatches;
    }
  }
  return mismatches ? 3 : 0;
}

/**
 * @brief Test-suite entry point.
 * @param argc Argument count from the C runtime.
 * @param argv Argument vector; argv[0] is the self path used by death tests,
 * remaining args (if any) name the modules to run, or
 * --list/-h/--help/--check-modules.
 * @return 0 on success, 1 if any test failed, a case was skipped under
 * HS_SKIPS_ARE_ERRORS, or a required CI depth lever is missing, 2 on an unknown
 * module name or an invalid death-child invocation, 3 on a --check-modules divergence.
 * @details Dispatches marked death-harness children, else runs the full
 * roster or only the modules named on argv.
 */
int main(int argc, char **argv) {
  std::setvbuf(stdout, nullptr, _IONBF, 0);

  if (test_suite().set_death_self_exe)
    test_suite().set_death_self_exe((argc > 0) ? argv[0] : nullptr);

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *death_child = std::getenv("HS_DEATH_CHILD");
  const char *death_case = std::getenv("HS_DEATH_CASE");
#pragma clang diagnostic pop
  const bool IS_DEATH_CHILD =
      death_child && std::strcmp(death_child, "harness") == 0;
  if (IS_DEATH_CHILD && (argc > 1 || !death_case || death_case[0] == '\0')) {
    std::fprintf(stderr, "run_tests: invalid death-child invocation\n");
    return 2;
  }

  if (argc > 1) {
    if (std::strcmp(argv[1], "--list") == 0 ||
        std::strcmp(argv[1], "-h") == 0 ||
        std::strcmp(argv[1], "--help") == 0) {
      std::printf("usage: run_tests [module...]\nmodules:\n");
      print_modules(stdout);
      return 0;
    }
    if (std::strcmp(argv[1], "--check-modules") == 0)
      return check_modules(argc - 2, argv + 2);
  }

  if (check_ci_levers(!IS_DEATH_CHILD && runs_effects(argc, argv)))
    return 1;

  if (IS_DEATH_CHILD) {
    if (!test_suite().run_death_child_case) {
      std::fprintf(stderr,
                   "run_tests: death cases are unavailable in this binary\n");
      return 2;
    }
    test_suite().run_death_child_case(death_case);
    return 0;
  }

  int failures = 0;
  if (argc > 1) {
    // Filtered run: execute only the named modules; an unknown name exits 2.
    for (int i = 1; i < argc; ++i) {
      const hs_test::TestModule *match = nullptr;
      for (const hs_test::TestModule &m : test_suite().modules) {
        if (std::strcmp(m.name, argv[i]) == 0) {
          match = &m;
          break;
        }
      }
      if (!match) {
        std::fprintf(stderr, "run_tests: unknown module '%s'\navailable:\n",
                     argv[i]);
        print_modules(stderr);
        return 2;
      }
      failures += match->run();
    }
  } else {
    for (const hs_test::TestModule &m : test_suite().modules)
      failures += m.run();
  }
  const int skipped = hs_test::stats().skipped;
  if (skipped > 0) {
    std::printf("=== suite: %d case(s) SKIPPED ===\n", skipped);
    if (hs_test::skips_are_errors()) {
      std::fprintf(stderr,
                   "run_tests: HS_SKIPS_ARE_ERRORS is on and %d case(s) "
                   "retired themselves — this leg is declared to assert every "
                   "case it compiles.\n",
                   skipped);
      ++failures;
    }
  }
  // Collapse to 0/1: a POSIX exit status is 8 bits, so a raw count can wrap.
  return (failures || hs_test::stats().failed > 0) ? 1 : 0;
}
