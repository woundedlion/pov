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
#include "tests/test_color.h"
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
#include "tests/test_param_marshal.h"
#include "tests/test_wasm_predicates.h"
#include "tests/test_led.h"
#include "tests/test_presets.h"
#include "tests/test_styles.h"
#include "tests/test_shading.h"

namespace hs_test::module_runners {

int run_3dmath_tests() { return hs_test::math3d_tests::run_3dmath_tests(); }
int run_concepts_tests() {
  return hs_test::concepts_tests::run_concepts_tests();
}
int run_memory_tests() { return hs_test::memory_tests::run_memory_tests(); }
int run_spatial_tests() { return hs_test::spatial_tests::run_spatial_tests(); }
int run_static_circular_buffer_tests() {
  return hs_test::scb_tests::run_static_circular_buffer_tests();
}
int run_color_tests() { return hs_test::color_tests::run_color_tests(); }
int run_palettes_tests() {
  return hs_test::palettes_tests::run_palettes_tests();
}
int run_easing_waves_tests() {
  return hs_test::easing_waves_tests::run_easing_waves_tests();
}
int run_interpolate_tests() {
  return hs_test::interpolate_tests::run_interpolate_tests();
}
int run_platform_tests() {
  return hs_test::platform_tests::run_platform_tests();
}
int run_profiling_tests() {
  return hs_test::profiling_tests::run_profiling_tests();
}
int run_transformers_tests() {
  return hs_test::transformers_tests::run_transformers_tests();
}
int run_noise_tests() { return hs_test::noise_tests::run_noise_tests(); }
int run_noise_field_tests() {
  return hs_test::noise_field_tests::run_noise_field_tests();
}
int run_projections_tests() {
  return hs_test::projections_tests::run_projections_tests();
}
int run_animation_tests() {
  return hs_test::animation_tests::run_animation_tests();
}
int run_param_marshal_tests() {
  return hs_test::param_marshal_tests::run_param_marshal_tests();
}
int run_wasm_predicates_tests() {
  return hs_test::wasm_predicates_tests::run_wasm_predicates_tests();
}
int run_led_tests() { return hs_test::led_tests::run_led_tests(); }
int run_presets_tests() { return hs_test::presets_tests::run_presets_tests(); }
int run_styles_tests() { return hs_test::styles_tests::run_styles_tests(); }
int run_shading_tests() { return hs_test::shading_tests::run_shading_tests(); }

} // namespace hs_test::module_runners
