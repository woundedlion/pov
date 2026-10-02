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

// Include orientation.h and the engine barrel before filter.h and the effects,
// as a real target does.
#include "core/engine/engine.h"

// Each module needs an include, an HS_TEST_MODULE_LIST row, and a matching
// _hs_test_modules entry in tests/CMakeLists.txt.
#include "tests/test_3dmath.h"
#include "tests/test_concepts.h"
#include "tests/test_memory.h"
#include "tests/test_spatial.h"
#include "tests/test_static_circular_buffer.h"
#include "tests/test_sdf.h"
#include "tests/test_conway.h"
#include "tests/test_conway_morph.h"
#include "tests/test_conway_continuity.h"
#include "tests/test_partition_seam.h"
#include "tests/test_conway_soak.h"
#include "tests/test_opchain_probe.h"
#include "tests/test_opchain_arena_survey.h"
#include "tests/test_hankin.h"
#include "tests/test_hyper_lattice.h"
#include "tests/test_ray_events.h"
#include "tests/test_ray_demonstrators.h"
#include "tests/test_sdf_patterns.h"
#include "tests/test_geometry.h"
#include "tests/test_spherical_field.h"
#include "tests/test_spherical_harmonics.h"
#include "tests/test_mesh.h"
#include "tests/test_solids.h"
#include "tests/test_reaction_graph.h"
#include "tests/test_color.h"
#include "tests/test_palettes.h"
#include "tests/test_easing_waves.h"
#include "tests/test_interpolate.h"
#include "tests/test_platform.h"
#include "tests/test_profiling.h"
#include "tests/test_pullback.h"
#include "tests/test_shader_chain.h"
#include "tests/test_filter.h"
#include "tests/test_plot_scan.h"
#include "tests/test_canvas.h"
#include "tests/test_scan.h"
#include "tests/test_ray.h"
#include "tests/test_mesh_raster.h"
#include "tests/test_transformers.h"
#include "tests/test_noise.h"
#include "tests/test_noise_field.h"
#include "tests/test_projections.h"
#include "tests/test_animation.h"
#include "tests/test_effects.h"
#include "tests/test_mindsplatter.h"
#include "tests/test_effects_smoke.h"
#include "tests/test_effect_factory.h"
#include "tests/test_lattice_melt.h"
#include "tests/test_kaleidoscope_smooth.h"
#include "tests/test_composed_effect.h"
#include "tests/test_shapeshifter_oracle.h"
#include "tests/test_shapeshifter_tiles.h"
#include "tests/test_dma_core.h"
#include "tests/test_hd107s_frame.h"
#include "tests/test_dma_controller.h"
#include "tests/test_pov_segmented.h"
#include "tests/test_pov_single.h"
#include "tests/test_pov_sync.h"
#include "tests/test_param_marshal.h"
#include "tests/test_wasm_predicates.h"
#include "tests/test_led.h"
#include "tests/test_presets.h"
#include "tests/test_styles.h"
#include "tests/test_shading.h"
#include "tests/test_death.h"

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
 * @brief One entry in the test-module roster: a short name plus its entry
 * point.
 * @details An unfiltered run executes every module in array order; passing one
 * or more names on argv runs ONLY those modules, in the order given — the
 * iteration-speed counterpart to the HS_DEATH_CASE single-case dispatch below.
 * Modules own their fixtures. self_exe() is set in main, and hs_test::stats()
 * accumulates their results for the process exit status. Any subset can run
 * in isolation.
 */
struct TestModule {
  const char *name; /**< Short module name matched against argv. */
  int (*run)();     /**< Entry point; returns the module's failure count. */
  bool effects_tier;
};

// Expands into MODULES[] and HS_TEST_MODULE_COUNT; mirrored by
// _hs_test_modules in tests/CMakeLists.txt.
#define HS_TEST_MODULE_LIST(X)                                                 \
  X("3dmath", hs_test::math3d_tests::run_3dmath_tests, false)                  \
  X("concepts", hs_test::concepts_tests::run_concepts_tests, false)            \
  X("memory", hs_test::memory_tests::run_memory_tests, false)                  \
  X("spatial", hs_test::spatial_tests::run_spatial_tests, false)               \
  X("scb", hs_test::scb_tests::run_static_circular_buffer_tests, false)        \
  X("ray", hs_test::ray_tests::run_ray_tests, false)                           \
  X("sdf", hs_test::sdf_tests::run_sdf_tests, false)                           \
  X("conway", hs_test::conway_tests::run_conway_tests, false)                  \
  X("conway_morph", hs_test::conway_morph_tests::run_conway_morph_tests,       \
    false)                                                                     \
  X("conway_continuity",                                                       \
    hs_test::conway_continuity_tests::run_conway_continuity_tests, false)      \
  X("partition_seam", hs_test::partition_seam_tests::run_partition_seam_tests, \
    false)                                                                     \
  X("conway_soak", hs_test::conway_soak_tests::run_conway_soak_tests, false)   \
  X("opchain_probe", hs_test::opchain_probe_tests::run_opchain_probe_tests,    \
    false)                                                                     \
  X("opchain_arena_survey",                                                    \
    hs_test::opchain_arena_survey_tests::run_opchain_arena_survey_tests,       \
    false)                                                                     \
  X("hankin", hs_test::hankin_tests::run_hankin_tests, false)                  \
  X("ray_demonstrators",                                                       \
    hs_test::ray_demonstrator_tests::run_ray_demonstrator_tests, false)        \
  X("sdf_patterns", hs_test::sdf_pattern_tests::run_sdf_pattern_tests, false)  \
  X("ray_events", hs_test::ray_event_tests::run_ray_event_tests, false)        \
  X("hyper_lattice", hs_test::hyper_lattice_tests::run_hyper_lattice_tests,    \
    false)                                                                     \
  X("geometry", hs_test::geometry_tests::run_geometry_tests, false)            \
  X("spherical_field",                                                         \
    hs_test::spherical_field_tests::run_spherical_field_tests, false)          \
  X("spherical_harmonics",                                                     \
    hs_test::spherical_harmonics_tests::run_spherical_harmonics_tests, false)  \
  X("mesh", hs_test::mesh_tests::run_mesh_tests, false)                        \
  X("solids", hs_test::solids_tests::run_solids_tests, false)                  \
  X("reaction_graph", hs_test::reaction_graph_tests::run_reaction_graph_tests, \
    false)                                                                     \
  X("color", hs_test::color_tests::run_color_tests, false)                     \
  X("palettes", hs_test::palettes_tests::run_palettes_tests, false)            \
  X("easing_waves", hs_test::easing_waves_tests::run_easing_waves_tests,       \
    false)                                                                     \
  X("interpolate", hs_test::interpolate_tests::run_interpolate_tests, false)   \
  X("platform", hs_test::platform_tests::run_platform_tests, false)            \
  X("profiling", hs_test::profiling_tests::run_profiling_tests, false)         \
  X("pullback", hs_test::pullback_tests::run_pullback_tests, false)            \
  X("shader_chain", hs_test::shader_chain_tests::run_shader_chain_tests,       \
    false)                                                                     \
  X("filter", hs_test::filter_tests::run_filter_tests, false)                  \
  X("plot_scan", hs_test::plot_scan_tests::run_plot_scan_tests, false)         \
  X("canvas", hs_test::canvas_tests::run_canvas_tests, false)                  \
  X("scan", hs_test::scan_tests::run_scan_tests, false)                        \
  X("mesh_raster", hs_test::mesh_raster_tests::run_mesh_raster_tests, false)   \
  X("transformers", hs_test::transformers_tests::run_transformers_tests,       \
    false)                                                                     \
  X("noise", hs_test::noise_tests::run_noise_tests, false)                     \
  X("noise_field", hs_test::noise_field_tests::run_noise_field_tests, false)   \
  X("projections", hs_test::projections_tests::run_projections_tests, false)   \
  X("animation", hs_test::animation_tests::run_animation_tests, false)         \
  X("effects", hs_test::effects_tests::run_effects_tests, true)                \
  X("mindsplatter", hs_test::mindsplatter_tests::run_mindsplatter_tests, true) \
  X("effects_smoke", hs_test::effects_smoke_tests::run_effects_smoke_tests,    \
    true)                                                                      \
  X("effect_factory", hs_test::effect_factory_tests::run_effect_factory_tests, \
    true)                                                                      \
  X("lattice_melt", hs_test::lattice_melt_tests::run_lattice_melt_tests,       \
    false)                                                                     \
  X("kaleidoscope_smooth",                                                     \
    hs_test::kaleidoscope_smooth_tests::run_kaleidoscope_smooth_tests, false)  \
  X("composed_effect",                                                         \
    hs_test::composed_effect_tests::run_composed_effect_tests, false)          \
  X("shapeshifter_oracle",                                                     \
    hs_test::shapeshifter_oracle_tests::run_shapeshifter_oracle_tests, false)  \
  X("shapeshifter_tiles",                                                      \
    hs_test::shapeshifter_tiles_tests::run_shapeshifter_tiles_tests, false)    \
  X("dma_core", hs_test::dma_core_tests::run_dma_core_tests, false)            \
  X("hd107s", hs_test::hd107s_tests::run_hd107s_tests, false)                  \
  X("dma_controller", hs_test::dma_controller_tests::run_dma_controller_tests, \
    false)                                                                     \
  X("pov_segmented", hs_test::pov_segmented_tests::run_pov_segmented_tests,    \
    false)                                                                     \
  X("pov_single", hs_test::pov_single_tests::run_pov_single_tests, false)      \
  X("pov_sync", hs_test::pov_sync_tests::run_pov_sync_tests, false)            \
  X("param_marshal", hs_test::param_marshal_tests::run_param_marshal_tests,    \
    false)                                                                     \
  X("wasm_predicates",                                                         \
    hs_test::wasm_predicates_tests::run_wasm_predicates_tests, false)          \
  X("led", hs_test::led_tests::run_led_tests, false)                           \
  X("presets", hs_test::presets_tests::run_presets_tests, false)               \
  X("styles", hs_test::styles_tests::run_styles_tests, false)                  \
  X("shading", hs_test::shading_tests::run_shading_tests, false)               \
  X("death", hs_test::death_tests::run_death_tests, false)

#define HS_TEST_MODULE_ENTRY(name, fn, effects_tier) {name, fn, effects_tier},
static const TestModule MODULES[] = {HS_TEST_MODULE_LIST(HS_TEST_MODULE_ENTRY)};
#undef HS_TEST_MODULE_ENTRY

#define HS_TEST_MODULE_COUNT_ADD(name, fn, effects_tier) +1
constexpr int HS_TEST_MODULE_COUNT =
    0 HS_TEST_MODULE_LIST(HS_TEST_MODULE_COUNT_ADD);
#undef HS_TEST_MODULE_COUNT_ADD
static_assert(
    sizeof(MODULES) / sizeof(MODULES[0]) == HS_TEST_MODULE_COUNT,
    "MODULES and HS_TEST_MODULE_COUNT disagree: both derive from "
    "HS_TEST_MODULE_LIST, so this fires only if the list is malformed.");

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
 * @return True for an unfiltered run, or a filtered run naming effects,
 * effects_smoke, effect_factory or mindsplatter.
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
 * @details The Canvas constructor's buffer_free() spin trips a trap, not a
 * failed assertion, so a shared runner that deschedules the test thread past
 * the shipping 2 s bound kills the whole shard. CI raises it; 30 s still
 * catches a genuinely stalled display hand-off.
 */
constexpr unsigned long CI_MIN_BUFFER_FREE_WATCHDOG_US = 30000000UL;

/** @brief Reports whether CI-only runner depth levers must be enforced. */
static bool runs_in_ci() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *ci = std::getenv("CI");
#pragma clang diagnostic pop
  return ci && ci[0] != '\0';
}

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
 * @details Every lever defaults to the shallow local tier when unset, so a
 * workflow step that stops exporting one drops the deep smoke window or the
 * full-resolution roster passes while still reporting green. Under CI that is
 * a failure, the same stance the death harness takes on a suite it cannot run.
 * The effects tier is scored only when an effects module actually runs, so jobs
 * selecting other modules are untouched; when one does run, both
 * HS_EFFECTS_FULL (which tier) and HS_REQUIRE_EFFECTS_FULL (whether FULL is
 * mandatory for this leg) must carry an explicit value. An absent or blank key
 * is a step that lost its declaration, not a vote for QUICK. The watchdog
 * lever is scored on every CI invocation: unlike the depth levers its absence
 * traps the process instead of reporting a test. HS_SKIPS_ARE_ERRORS is
 * scored on every invocation too: a skipped case asserts nothing and moves no
 * counter a green run shows, so whether this leg tolerates one has to be
 * declared rather than defaulted.
 */
static int check_ci_levers(bool effects_invocation) {
  if (!runs_in_ci())
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
               !hs_test::effects_tests::effects_full_suite()) {
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
 * @details Both directions are checked so a roster/expected-list drift in
 * either direction fails: every MODULES name must appear in argv, and every
 * argv name must be a roster module. Used by the CTest that pins
 * _hs_test_modules.
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

  hs_test::death_tests::self_exe() = (argc > 0) ? argv[0] : nullptr;

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
    hs_test::death_tests::run_child_case(death_case);
    return 0;
  }

  int failures = 0;
  if (argc > 1) {
    // Filtered run: execute only the named modules. An unknown name fails fast
    // (exit 2) so a typo never silently runs nothing.
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
  // Collapse to 0/1: a process exit status is only 8 bits on POSIX, so
  // returning a raw count would wrap (e.g. 256 failures -> 0 -> green CI).
  return (failures || hs_test::stats().failed > 0) ? 1 : 0;
}
