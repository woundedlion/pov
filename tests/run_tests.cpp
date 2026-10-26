/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#if defined(_WIN32)
#include <crtdbg.h>
#endif

#include "tests/test_fixture.h"

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

/**
 * @brief Test-module name, entry point and effects-tier participation.
 * @details An unfiltered run executes every module in array order; names on
 * argv run only those modules, in the order given.
 */
struct TestModule {
  const char *name; /**< Short module name matched against argv. */
  int (*run)();     /**< Entry point; returns the module's failure count. */
  bool
      effects_tier; /**< Reads effects_full_suite(); selects CI tier-lever checks. */
};

// Expands into MODULES[]; mirrored by
// _hs_test_modules in tests/CMakeLists.txt.
#define HS_TEST_MODULE_LIST(X)                                                 \
  X("3dmath", run_3dmath_tests, false)                                         \
  X("concepts", run_concepts_tests, false)                                     \
  X("memory", run_memory_tests, false)                                         \
  X("spatial", run_spatial_tests, false)                                       \
  X("scb", run_static_circular_buffer_tests, false)                            \
  X("ray", run_ray_tests, false)                                               \
  X("sdf", run_sdf_tests, false)                                               \
  X("conway", run_conway_tests, false)                                         \
  X("conway_morph", run_conway_morph_tests, false)                             \
  X("conway_continuity", run_conway_continuity_tests, false)                   \
  X("partition_seam", run_partition_seam_tests, false)                         \
  X("conway_soak", run_conway_soak_tests, false)                               \
  X("opchain_probe", run_opchain_probe_tests, false)                           \
  X("opchain_arena_survey", run_opchain_arena_survey_tests, false)             \
  X("hankin", run_hankin_tests, false)                                         \
  X("ray_demonstrators", run_ray_demonstrator_tests, false)                    \
  X("sdf_patterns", run_sdf_pattern_tests, false)                              \
  X("ray_events", run_ray_event_tests, false)                                  \
  X("hyper_lattice", run_hyper_lattice_tests, false)                           \
  X("geometry", run_geometry_tests, false)                                     \
  X("spherical_field", run_spherical_field_tests, false)                       \
  X("spherical_harmonics", run_spherical_harmonics_tests, false)               \
  X("mesh", run_mesh_tests, false)                                             \
  X("solids", run_solids_tests, false)                                         \
  X("reaction_graph", run_reaction_graph_tests, false)                         \
  X("color", run_color_tests, false)                                           \
  X("palettes", run_palettes_tests, false)                                     \
  X("easing_waves", run_easing_waves_tests, false)                             \
  X("interpolate", run_interpolate_tests, false)                               \
  X("platform", run_platform_tests, false)                                     \
  X("profiling", run_profiling_tests, false)                                   \
  X("pullback", run_pullback_tests, false)                                     \
  X("shader_chain", run_shader_chain_tests, false)                             \
  X("filter", run_filter_tests, false)                                         \
  X("plot_parity", run_plot_parity_tests, false)                               \
  X("plot_scan", run_plot_scan_tests, false)                                   \
  X("canvas", run_canvas_tests, false)                                         \
  X("scan", run_scan_tests, false)                                             \
  X("mesh_raster", run_mesh_raster_tests, false)                               \
  X("transformers", run_transformers_tests, false)                             \
  X("noise", run_noise_tests, false)                                           \
  X("noise_field", run_noise_field_tests, false)                               \
  X("projections", run_projections_tests, false)                               \
  X("animation", run_animation_tests, false)                                   \
  X("effects", run_effects_tests, true)                                        \
  X("mindsplatter", run_mindsplatter_tests, true)                              \
  X("effects_smoke", run_effects_smoke_tests, true)                            \
  X("effect_factory", run_effect_factory_tests, true)                          \
  X("kaleidoscope_melt", run_kaleidoscope_melt_tests, false)                   \
  X("kaleidoscope_smooth", run_kaleidoscope_smooth_tests, false)               \
  X("composed_effect", run_composed_effect_tests, false)                       \
  X("shapeshifter_oracle", run_shapeshifter_oracle_tests, false)               \
  X("shapeshifter_tiles", run_shapeshifter_tiles_tests, false)                 \
  X("dma_core", run_dma_core_tests, false)                                     \
  X("hd107s", run_hd107s_tests, false)                                         \
  X("dma_controller", run_dma_controller_tests, false)                         \
  X("pov_segmented", run_pov_segmented_tests, false)                           \
  X("pov_single", run_pov_single_tests, false)                                 \
  X("pov_sync", run_pov_sync_tests, false)                                     \
  X("param_marshal", run_param_marshal_tests, false)                           \
  X("wasm_predicates", run_wasm_predicates_tests, false)                       \
  X("led", run_led_tests, false)                                               \
  X("presets", run_presets_tests, false)                                       \
  X("styles", run_styles_tests, false)                                         \
  X("shading", run_shading_tests, false)                                       \
  X("death", run_death_tests, false)

namespace hs_test::module_runners {
void set_death_self_exe(const char *path);
void run_death_child_case(const char *name);
bool effects_full_suite();
#define HS_TEST_MODULE_ENTRY(name, fn, effects_tier) int fn();
HS_TEST_MODULE_LIST(HS_TEST_MODULE_ENTRY)
#undef HS_TEST_MODULE_ENTRY
} // namespace hs_test::module_runners

#define HS_TEST_MODULE_ENTRY(name, fn, effects_tier)                           \
  {name, hs_test::module_runners::fn, effects_tier},
static const TestModule MODULES[] = {HS_TEST_MODULE_LIST(HS_TEST_MODULE_ENTRY)};
#undef HS_TEST_MODULE_ENTRY

/**
 * @brief Prints the roster's module names, one indented per line.
 * @param out Destination stream (e.g. stdout for --list, stderr on error).
 */
static void print_modules(std::FILE *out) {
  for (const TestModule &m : MODULES)
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
  if (argc <= 1)
    return true;
  for (int i = 1; i < argc; ++i)
    for (const TestModule &module : MODULES)
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

/** @brief Reports whether a skipped case must fail the run. */
static bool skips_are_errors() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *lever = std::getenv("HS_SKIPS_ARE_ERRORS");
#pragma clang diagnostic pop
  return lever && std::atoi(lever) != 0;
}

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
  if (effects_invocation) {
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
               !hs_test::module_runners::effects_full_suite()) {
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
 * @return 0 if the set matches MODULES exactly, 3 on any divergence.
 * @details Every MODULES name must appear in argv, and every argv name must be
 * a roster module.
 */
static int check_modules(int argc, char **argv) {
  int mismatches = 0;
  for (const TestModule &m : MODULES) {
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
    for (const TestModule &m : MODULES)
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
#if defined(_WIN32)
  SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
  _set_error_mode(_OUT_TO_STDERR);
  _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
#endif
  // Unbuffered stdout so progress survives a trap/abort in a death-case child.
  std::setvbuf(stdout, nullptr, _IONBF, 0);

  hs_test::module_runners::set_death_self_exe((argc > 0) ? argv[0] : nullptr);

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
    hs_test::module_runners::run_death_child_case(death_case);
    return 0;
  }

  int failures = 0;
  if (argc > 1) {
    // Filtered run: execute only the named modules; an unknown name exits 2.
    for (int i = 1; i < argc; ++i) {
      const TestModule *match = nullptr;
      for (const TestModule &m : MODULES) {
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
    for (const TestModule &m : MODULES)
      failures += m.run();
  }
  const int skipped = hs_test::stats().skipped;
  if (skipped > 0) {
    std::printf("=== suite: %d case(s) SKIPPED ===\n", skipped);
    if (skips_are_errors()) {
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
