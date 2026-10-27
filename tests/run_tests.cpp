/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Complete native test-module roster.
 */
#include "tests/module_registry.h"

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
#define HS_TEST_MODULE_ENTRY(name, fn, effects_tier) int fn();
HS_TEST_MODULE_LIST(HS_TEST_MODULE_ENTRY)
#undef HS_TEST_MODULE_ENTRY
} // namespace hs_test::module_runners

#define HS_TEST_MODULE_ENTRY(name, fn, effects_tier)                           \
  {name, hs_test::module_runners::fn, effects_tier},

const hs_test::TestSuite &hs_test::test_suite() {
  static const TestModule MODULES[] = {
      HS_TEST_MODULE_LIST(HS_TEST_MODULE_ENTRY)};
  static const TestSuite SUITE{MODULES, module_runners::set_death_self_exe,
                               module_runners::run_death_child_case};
  return SUITE;
}
#undef HS_TEST_MODULE_ENTRY
