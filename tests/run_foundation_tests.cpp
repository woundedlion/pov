/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_3dmath.h"
#include "tests/test_concepts.h"
#include "tests/test_memory.h"
#include "tests/test_spatial.h"
#include "tests/test_static_circular_buffer.h"
#include "tests/test_palettes.h"
#include "tests/test_easing_waves.h"
#include "tests/test_interpolate.h"
#include "tests/test_platform.h"
#include "tests/test_profiling.h"
#include "tests/test_transformers.h"
#include "tests/test_noise.h"
#include "tests/test_noise_field.h"
#include "tests/test_projections.h"
#include "tests/test_animation.h"
#include "tests/test_wasm_predicates.h"
#include "tests/test_led.h"
#include "tests/test_presets.h"
#include "tests/test_styles.h"
#include "tests/test_shading.h"

namespace hs_test::module_runners {

HS_FORWARD_MODULE(math3d_tests, run_3dmath_tests)
HS_FORWARD_MODULE(concepts_tests, run_concepts_tests)
HS_FORWARD_MODULE(memory_tests, run_memory_tests)
HS_FORWARD_MODULE(spatial_tests, run_spatial_tests)
HS_FORWARD_MODULE(scb_tests, run_static_circular_buffer_tests)
HS_FORWARD_MODULE(palettes_tests, run_palettes_tests)
HS_FORWARD_MODULE(easing_waves_tests, run_easing_waves_tests)
HS_FORWARD_MODULE(interpolate_tests, run_interpolate_tests)
HS_FORWARD_MODULE(platform_tests, run_platform_tests)
HS_FORWARD_MODULE(profiling_tests, run_profiling_tests)
HS_FORWARD_MODULE(transformers_tests, run_transformers_tests)
HS_FORWARD_MODULE(noise_tests, run_noise_tests)
HS_FORWARD_MODULE(noise_field_tests, run_noise_field_tests)
HS_FORWARD_MODULE(projections_tests, run_projections_tests)
HS_FORWARD_MODULE(animation_tests, run_animation_tests)
HS_FORWARD_MODULE(wasm_predicates_tests, run_wasm_predicates_tests)
HS_FORWARD_MODULE(led_tests, run_led_tests)
HS_FORWARD_MODULE(presets_tests, run_presets_tests)
HS_FORWARD_MODULE(styles_tests, run_styles_tests)
HS_FORWARD_MODULE(shading_tests, run_shading_tests)

} // namespace hs_test::module_runners
