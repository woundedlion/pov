/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include <cstdio>
#include <cstring>

// Pull in the engine barrel first, exactly as a real target does, so geometry.h
// (which defines Orientation<CAP=4>) is ordered ahead of filter.h and the
// effects for the whole translation unit.
#include "core/engine/engine.h"

// One include per registered module, paired 1:1 with a row in
// HS_TEST_MODULE_LIST below; keep the two in sync when adding a module.
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
#include "tests/test_geometry.h"
#include "tests/test_spherical_field.h"
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
#include "tests/test_shader_workbench.h"
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

/**
 * @brief One entry in the test-module roster: a short name plus its entry
 * point.
 * @details An unfiltered run executes every module in array order; passing one
 * or more names on argv runs ONLY those modules, in the order given — the
 * iteration-speed counterpart to the HS_DEATH_CASE single-case dispatch below.
 * Modules are independent (each owns its own fixtures; the only cross-cutting
 * state, self_exe(), is set unconditionally in main), so any subset is safe to
 * run in isolation.
 */
struct TestModule {
  const char *name; /**< Short module name matched against argv. */
  int (*run)();     /**< Entry point; returns the module's failure count. */
};

// Single source of truth for the roster: expands into both MODULES[] and the
// derived HS_TEST_MODULE_COUNT below. Adding a module means an #include above
// AND one X(...) row here. Mirrors core/engine/effects.h's HS_EFFECT_LIST.
#define HS_TEST_MODULE_LIST(X)                                                 \
  X("3dmath", hs_test::math3d_tests::run_3dmath_tests)                         \
  X("concepts", hs_test::concepts_tests::run_concepts_tests)                   \
  X("memory", hs_test::memory_tests::run_memory_tests)                         \
  X("spatial", hs_test::spatial_tests::run_spatial_tests)                      \
  X("scb", hs_test::scb_tests::run_static_circular_buffer_tests)               \
  X("sdf", hs_test::sdf_tests::run_sdf_tests)                                  \
  X("conway", hs_test::conway_tests::run_conway_tests)                         \
  X("conway_morph", hs_test::conway_morph_tests::run_conway_morph_tests)       \
  X("conway_continuity",                                                       \
    hs_test::conway_continuity_tests::run_conway_continuity_tests)             \
  X("partition_seam", hs_test::partition_seam_tests::run_partition_seam_tests) \
  X("conway_soak", hs_test::conway_soak_tests::run_conway_soak_tests)          \
  X("opchain_probe", hs_test::opchain_probe_tests::run_opchain_probe_tests)    \
  X("opchain_arena_survey",                                                    \
    hs_test::opchain_arena_survey_tests::run_opchain_arena_survey_tests)       \
  X("hankin", hs_test::hankin_tests::run_hankin_tests)                         \
  X("hyper_lattice", hs_test::hyper_lattice_tests::run_hyper_lattice_tests)    \
  X("geometry", hs_test::geometry_tests::run_geometry_tests)                   \
  X("spherical_field",                                                         \
    hs_test::spherical_field_tests::run_spherical_field_tests)                 \
  X("mesh", hs_test::mesh_tests::run_mesh_tests)                               \
  X("solids", hs_test::solids_tests::run_solids_tests)                         \
  X("reaction_graph", hs_test::reaction_graph_tests::run_reaction_graph_tests) \
  X("color", hs_test::color_tests::run_color_tests)                            \
  X("palettes", hs_test::palettes_tests::run_palettes_tests)                   \
  X("easing_waves", hs_test::easing_waves_tests::run_easing_waves_tests)       \
  X("interpolate", hs_test::interpolate_tests::run_interpolate_tests)          \
  X("platform", hs_test::platform_tests::run_platform_tests)                   \
  X("profiling", hs_test::profiling_tests::run_profiling_tests)                \
  X("pullback", hs_test::pullback_tests::run_pullback_tests)                   \
  X("shader_chain", hs_test::shader_chain_tests::run_shader_chain_tests)       \
  X("filter", hs_test::filter_tests::run_filter_tests)                         \
  X("plot_scan", hs_test::plot_scan_tests::run_plot_scan_tests)                \
  X("canvas", hs_test::canvas_tests::run_canvas_tests)                         \
  X("scan", hs_test::scan_tests::run_scan_tests)                               \
  X("mesh_raster", hs_test::mesh_raster_tests::run_mesh_raster_tests)          \
  X("transformers", hs_test::transformers_tests::run_transformers_tests)       \
  X("noise", hs_test::noise_tests::run_noise_tests)                            \
  X("noise_field", hs_test::noise_field_tests::run_noise_field_tests)          \
  X("projections", hs_test::projections_tests::run_projections_tests)          \
  X("animation", hs_test::animation_tests::run_animation_tests)                \
  X("effects", hs_test::effects_tests::run_effects_tests)                      \
  X("mindsplatter", hs_test::mindsplatter_tests::run_mindsplatter_tests)       \
  X("effects_smoke", hs_test::effects_smoke_tests::run_effects_smoke_tests)    \
  X("effect_factory", hs_test::effect_factory_tests::run_effect_factory_tests) \
  X("shader_workbench",                                                        \
    hs_test::shader_workbench_tests::run_shader_workbench_tests)               \
  X("lattice_melt", hs_test::lattice_melt_tests::run_lattice_melt_tests)       \
  X("kaleidoscope_smooth",                                                     \
    hs_test::kaleidoscope_smooth_tests::run_kaleidoscope_smooth_tests)         \
  X("composed_effect",                                                         \
    hs_test::composed_effect_tests::run_composed_effect_tests)                 \
  X("shapeshifter_oracle",                                                     \
    hs_test::shapeshifter_oracle_tests::run_shapeshifter_oracle_tests)         \
  X("shapeshifter_tiles",                                                      \
    hs_test::shapeshifter_tiles_tests::run_shapeshifter_tiles_tests)           \
  X("dma_core", hs_test::dma_core_tests::run_dma_core_tests)                   \
  X("hd107s", hs_test::hd107s_tests::run_hd107s_tests)                         \
  X("dma_controller", hs_test::dma_controller_tests::run_dma_controller_tests) \
  X("pov_segmented", hs_test::pov_segmented_tests::run_pov_segmented_tests)    \
  X("pov_single", hs_test::pov_single_tests::run_pov_single_tests)             \
  X("pov_sync", hs_test::pov_sync_tests::run_pov_sync_tests)                   \
  X("param_marshal", hs_test::param_marshal_tests::run_param_marshal_tests)    \
  X("wasm_predicates",                                                         \
    hs_test::wasm_predicates_tests::run_wasm_predicates_tests)                 \
  X("led", hs_test::led_tests::run_led_tests)                                  \
  X("presets", hs_test::presets_tests::run_presets_tests)                      \
  X("styles", hs_test::styles_tests::run_styles_tests)                         \
  X("shading", hs_test::shading_tests::run_shading_tests)                      \
  X("death", hs_test::death_tests::run_death_tests)

#define HS_TEST_MODULE_ENTRY(name, fn) {name, fn},
static const TestModule MODULES[] = {HS_TEST_MODULE_LIST(HS_TEST_MODULE_ENTRY)};
#undef HS_TEST_MODULE_ENTRY

#define HS_TEST_MODULE_COUNT_ADD(name, fn) +1
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
 * effects_smoke or mindsplatter.
 */
static bool runs_effects(int argc, char **argv) {
  if (argc <= 1)
    return true;
  for (int i = 1; i < argc; ++i)
    if (std::strcmp(argv[i], "effects") == 0 ||
        std::strcmp(argv[i], "effects_smoke") == 0 ||
        std::strcmp(argv[i], "mindsplatter") == 0)
      return true;
  return false;
}

/**
 * @brief Smoke window every CI leg has to set.
 * @details Below this the frame-cyclic paths stay unreached and no preset
 * transition arms, so a shallower window reports green over untested code.
 */
constexpr int CI_MIN_SMOKE_FRAMES = 120;

/** @brief Reports whether CI-only runner depth levers must be enforced. */
static bool runs_in_ci() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *ci = std::getenv("CI");
#pragma clang diagnostic pop
  return ci && ci[0] != '\0';
}

/**
 * @brief Verifies the environment carries the depth levers a CI run must set.
 * @param argc Argument count as passed to main.
 * @param argv Argument vector as passed to main.
 * @return 0 when every lever this invocation needs is set, else 1.
 * @details Both levers default to the shallow local tier when unset. A workflow
 * that stops exporting one drops the deep smoke window and the
 * full-resolution roster passes while still reporting green. Under CI that is
 * a failure, the same stance the death harness takes on a suite it cannot run.
 * HS_EFFECTS_FULL is required only when an effects module actually runs, so
 * jobs that select other modules are untouched.
 */
static int check_ci_levers(int argc, char **argv) {
  if (!runs_in_ci())
    return 0;

  int missing = 0;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *frames = std::getenv("HS_SMOKE_FRAMES");
  const char *require_effects_full = std::getenv("HS_REQUIRE_EFFECTS_FULL");
#pragma clang diagnostic pop
  if (!frames || std::atoi(frames) < CI_MIN_SMOKE_FRAMES) {
    std::fprintf(stderr,
                 "run_tests: CI=on but HS_SMOKE_FRAMES is unset or below %d — "
                 "a shallower window skips frame-cyclic paths and arms no "
                 "preset transition. Set HS_SMOKE_FRAMES in the workflow "
                 "step's env.\n",
                 CI_MIN_SMOKE_FRAMES);
    ++missing;
  }
  if (require_effects_full && std::atoi(require_effects_full) > 0 &&
      runs_effects(argc, argv) &&
      !hs_test::effects_tests::effects_full_suite()) {
    std::fprintf(stderr,
                 "run_tests: CI=on but HS_EFFECTS_FULL is unset or 0 — the "
                 "effects modules would run the QUICK tier, dropping the "
                 "288x144 roster passes and the white-box block. Set "
                 "HS_EFFECTS_FULL=1 in the workflow step's env.\n");
    ++missing;
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
 * @return 0 on success, 1 if any test failed or a required CI depth lever is
 * missing, 2 on an unknown module name, 3 on a --check-modules divergence.
 * @details Dispatches the HS_DEATH_CASE child case if set, else runs the full
 * roster or only the modules named on argv.
 */
int main(int argc, char **argv) {
  // Unbuffered stdout so progress survives a trap/abort in a death-case child.
  std::setvbuf(stdout, nullptr, _IONBF, 0);

  hs_test::death_tests::self_exe() = (argc > 0) ? argv[0] : nullptr;

  // Child death-case dispatch: when HS_DEATH_CASE is set, run ONLY that single
  // trap-triggering case and exit — never the full suite (which would re-spawn
  // children recursively). The case is expected to __builtin_trap(); reaching
  // the return here means it did NOT trap, so we exit 0 and let the parent flag
  // the missing trap.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (const char *dc = std::getenv("HS_DEATH_CASE")) {
#pragma clang diagnostic pop
    if (dc[0] != '\0') {
      hs_test::death_tests::run_child_case(dc);
      return 0;
    }
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

  if (check_ci_levers(argc, argv))
    return 1;

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
  // Collapse to 0/1: a process exit status is only 8 bits on POSIX, so
  // returning a raw count would wrap (e.g. 256 failures -> 0 -> green CI).
  return failures ? 1 : 0;
}
