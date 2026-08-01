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
#include "tests/test_geometry.h"
#include "tests/test_spherical_field.h"
#include "tests/test_mesh.h"
#include "tests/test_solids.h"
#include "tests/test_reaction_graph.h"
#include "tests/test_color.h"
#include "tests/test_palettes.h"
#include "tests/test_easing_waves.h"
#include "tests/test_platform.h"
#include "tests/test_filter.h"
#include "tests/test_plot_scan.h"
#include "tests/test_canvas.h"
#include "tests/test_scan.h"
#include "tests/test_mesh_raster.h"
#include "tests/test_transformers.h"
#include "tests/test_noise.h"
#include "tests/test_generators.h"
#include "tests/test_animation.h"
#include "tests/test_effects.h"
#include "tests/test_shapeshifter_oracle.h"
#include "tests/test_dma_core.h"
#include "tests/test_hd107s_frame.h"
#include "tests/test_dma_controller.h"
#include "tests/test_pov_segmented.h"
#include "tests/test_pov_single.h"
#include "tests/test_pov_sync.h"
#include "tests/test_param_marshal.h"
#include "tests/test_wasm_predicates.h"
#include "tests/test_util.h"
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
  const char *name;   /**< Short module name matched against argv. */
  int (*run)();       /**< Entry point; returns the module's failure count. */
  int min_assertions; /**< Assertion floor; 0 means unfloored. */
};

// Single source of truth for the roster: expands into both MODULES[] and the
// derived HS_TEST_MODULE_COUNT below. Adding a module means an #include above
// AND one X(...) row here. Mirrors core/engine/effects.h's HS_EFFECT_LIST.
//
// The third column is the module's minimum assertion count, enforced by
// hs_test::end_module(). Cases are invoked by hand-written calls, so a case that
// is defined and never called compiles clean; the floor turns the resulting drop
// in assertions red.
//
// Every value is MEASURED, as the minimum over the configurations that can move
// a module's count: the smoke window (HS_SMOKE_FRAMES), Debug -O0 against
// RelWithDebInfo -O2, and the NDEBUG-gated cases in test_memory.h (7
// assertions, so memory's floor sits 7 under its Debug count). Adding cases
// raises a module's count and leaves its floor valid; lower a floor only after
// deliberately removing or cheapening a case. 0 leaves a module unfloored.
// death's floor covers its whole-suite skip exits because end_module exempts a
// module that reported a skip.

// Effects floors, one per depth tier (tests/test_effects.h
// effects_full_suite()). A single floor would have to be the QUICK count,
// leaving every FULL-tier case deletable without turning CI — the only runner of
// that tier — red. The roster row below selects the floor for the tier the
// environment picked.
constexpr int EFFECTS_QUICK_MIN_ASSERTIONS = 136043;
constexpr int EFFECTS_FULL_MIN_ASSERTIONS = 171001;

#define HS_TEST_MODULE_LIST(X)                                                 \
  X("3dmath", hs_test::math3d::run_3dmath_tests, 28826)                        \
  X("concepts", hs_test::concepts_tests::run_concepts_tests, 43273)            \
  X("memory", hs_test::mem::run_memory_tests, 271)                             \
  X("spatial", hs_test::spatial::run_spatial_tests, 179)                       \
  X("scb", hs_test::scb::run_static_circular_buffer_tests, 197)                \
  X("sdf", hs_test::sdf::run_sdf_tests, 271681)                                \
  X("conway", hs_test::conway_tests::run_conway_tests, 5982)                   \
  X("conway_morph", hs_test::conway_morph_tests::run_conway_morph_tests,       \
    1656015)                                                                   \
  X("conway_continuity",                                                       \
    hs_test::conway_continuity_tests::run_conway_continuity_tests, 49151)      \
  X("partition_seam", hs_test::partition_seam_tests::run_partition_seam_tests, \
    39)                                                                        \
  X("conway_soak", hs_test::conway_soak_tests::run_conway_soak_tests, 56)      \
  X("opchain_probe", hs_test::opchain_probe_tests::run_opchain_probe_tests,    \
    496265)                                                                    \
  X("opchain_arena_survey",                                                    \
    hs_test::opchain_arena_survey::run_opchain_arena_survey_tests, 21562)      \
  X("hankin", hs_test::hankin_tests::run_hankin_tests, 1948)                   \
  X("geometry", hs_test::geometry::run_geometry_tests, 4848)                   \
  X("spherical_field", hs_test::spherical_field::run_spherical_field_tests,    \
    231)                                                                       \
  X("mesh", hs_test::mesh_tests::run_mesh_tests, 42756)                        \
  X("solids", hs_test::solids_tests::run_solids_tests, 307007)                 \
  X("reaction_graph", hs_test::reaction_graph_tests::run_reaction_graph_tests, \
    983)                                                                       \
  X("color", hs_test::color_tests::run_color_tests, 630110)                    \
  X("palettes", hs_test::palettes_tests::run_palettes_tests, 95)               \
  X("easing_waves", hs_test::easing_waves_tests::run_easing_waves_tests, 6080) \
  X("platform", hs_test::platform_tests::run_platform_tests, 201119)           \
  X("filter", hs_test::filter_tests::run_filter_tests, 10297)                  \
  X("plot_scan", hs_test::plot_scan_tests::run_plot_scan_tests, 9109889)       \
  X("canvas", hs_test::canvas_tests::run_canvas_tests, 487)                    \
  X("scan", hs_test::scan_tests::run_scan_tests, 189869)                       \
  X("mesh_raster", hs_test::mesh_raster_tests::run_mesh_raster_tests, 2793)    \
  X("transformers", hs_test::transformers_tests::run_transformers_tests,       \
    216596)                                                                    \
  X("noise", hs_test::noise_tests::run_noise_tests, 201)                       \
  X("generators", hs_test::generators_tests::run_generators_tests, 59)         \
  X("animation", hs_test::animation_tests::run_animation_tests, 16488)         \
  X("effects", hs_test::effects_tests::run_effects_tests,                      \
    hs_test::effects_tests::effects_full_suite()                               \
        ? EFFECTS_FULL_MIN_ASSERTIONS                                          \
        : EFFECTS_QUICK_MIN_ASSERTIONS)                                        \
  X("shapeshifter_oracle",                                                   \
    hs_test::shapeshifter_oracle_tests::run_shapeshifter_oracle_tests, 109)   \
  X("dma_core", hs_test::dma_core::run_dma_core_tests, 12)                     \
  X("hd107s", hs_test::hd107s_tests::run_hd107s_tests, 280)                    \
  X("dma_controller", hs_test::dma_controller::run_dma_controller_tests, 67)   \
  X("pov_segmented", hs_test::pov_segmented_tests::run_pov_segmented_tests,    \
    263162)                                                                    \
  X("pov_single", hs_test::pov_single_tests::run_pov_single_tests, 8640)       \
  X("pov_sync", hs_test::pov_sync_tests::run_pov_sync_tests, 2015)             \
  X("param_marshal", hs_test::param_marshal_tests::run_param_marshal_tests,    \
    1158)                                                                      \
  X("wasm_predicates",                                                         \
    hs_test::wasm_predicates_tests::run_wasm_predicates_tests, 128)            \
  X("util", hs_test::util_tests::run_util_tests, 288)                          \
  X("led", hs_test::led_tests::run_led_tests, 12)                              \
  X("presets", hs_test::presets_tests::run_presets_tests, 23)                  \
  X("styles", hs_test::styles_tests::run_styles_tests, 17862)                  \
  X("shading", hs_test::shading_tests::run_shading_tests, 43)                  \
  X("death", hs_test::death::run_death_tests, 85)

#define HS_TEST_MODULE_ENTRY(name, fn, min_assertions)                         \
  {name, fn, min_assertions},
static const TestModule MODULES[] = {HS_TEST_MODULE_LIST(HS_TEST_MODULE_ENTRY)};
#undef HS_TEST_MODULE_ENTRY

#define HS_TEST_MODULE_COUNT_ADD(name, fn, min_assertions) +1
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
 * @brief Runs one roster module under its assertion floor.
 * @param m The roster entry to run.
 * @return The module's failure count, plus one if it fell below its floor.
 * @details Publishes the floor for the ModuleScope the module opens, so both the
 * full run and an argv-filtered subset enforce the same count.
 */
static int run_module(const TestModule &m) {
  hs_test::pending_min_assertions() = m.min_assertions;
  return m.run();
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
 * @return 0 on success, 1 if any test failed, 2 on an unknown module name, 3 on
 * a --check-modules divergence.
 * @details Dispatches the HS_DEATH_CASE child case if set, else runs the full
 * roster or only the modules named on argv.
 */
int main(int argc, char **argv) {
  // Unbuffered stdout so progress survives a trap/abort in a death-case child.
  std::setvbuf(stdout, nullptr, _IONBF, 0);

  hs_test::death::self_exe() = (argc > 0) ? argv[0] : nullptr;

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
      hs_test::death::run_child_case(dc);
      return 0;
    }
  }

  int failures = 0;
  if (argc > 1) {
    // Filtered run: execute only the named modules. An unknown name fails fast
    // (exit 2) so a typo never silently runs nothing.
    if (std::strcmp(argv[1], "--list") == 0 ||
        std::strcmp(argv[1], "-h") == 0 ||
        std::strcmp(argv[1], "--help") == 0) {
      std::printf("usage: run_tests [module...]\nmodules:\n");
      print_modules(stdout);
      return 0;
    }
    if (std::strcmp(argv[1], "--check-modules") == 0)
      return check_modules(argc - 2, argv + 2);
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
      failures += run_module(*match);
    }
  } else {
    for (const TestModule &m : MODULES)
      failures += run_module(m);
  }
  // Collapse to 0/1: a process exit status is only 8 bits on POSIX, so
  // returning a raw count would wrap (e.g. 256 failures -> 0 -> green CI).
  return failures ? 1 : 0;
}
