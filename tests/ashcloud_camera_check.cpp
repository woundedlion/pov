/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file ashcloud_camera_check.cpp
 * @brief Runtime-geometry qualification of AshCloud camera columns. */
#include "effects/AshCloud.h"
#include "tests/composed_frame_fixture.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include <bit>

namespace hs_test::composed_effect_tests {
inline uint32_t bits(float value) { return std::bit_cast<uint32_t>(value); }
#include "tests/composed_effect/ashcloud_shared_color.h"
} // namespace hs_test::composed_effect_tests

int main() {
  hs_test::ModuleFixture fixture("AshCloud runtime camera");
  hs_test::composed_effect_tests::test_ash_cloud_camera_column_admission();
  HS_EXPECT_TRUE(
      math::set_display_geometry(.02f * math::PI_F, .98f * math::PI_F));
  math::init_geometry_luts<288, 144>();
  HS_EXPECT_FALSE((AshCloud<288, 144>::camera_column_admitted()));
  hs_test::composed_effect_tests::test_ash_cloud_shared_camera_qualification();
  hs_test::composed_effect_tests::test_ash_cloud_shared_color_identity();
  hs_test::composed_effect_tests::test_ash_cloud_shared_column_fallbacks();
  return fixture.result() ? 1 : 0;
}
