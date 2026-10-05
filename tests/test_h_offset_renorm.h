/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * South-pole Y-clip renormalization coverage for the legacy H_OFFSET mapping.
 * Compiled with HS_TEST_H_OFFSET=3; shipping targets use H_OFFSET == 0.
 */
#pragma once

#include "core/render/filter.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/render/plot.h"
#include "core/render/scan.h"
#include "tests/test_fixture.h"
#include "tests/pole_geometry_test_util.h"
#include "tests/test_harness.h"
#include "tests/test_pole_wrap.h"

namespace hs_test {
namespace h_offset_renorm {

// With HS_TEST_H_OFFSET=3, H_VIRT = H + 3 and the bottom physical row y=H-1
// lands SHORT of the south pole (sin(phi) > 0).
#ifdef HS_OFFSET_FULL_RESOLUTION
constexpr int W = 96;
constexpr int H = 20;
#else
constexpr int W = 32;
constexpr int H = 16;
#endif

inline void face_white(const math::Vector &, Fragment &fragment) {
  fragment.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
}

/** @brief Verifies Face row bounds and rasterization use the virtual height. */
inline void test_face_bounds_use_virtual_height() {
  constexpr float PHI_MIN = 1.2f;
  constexpr float PHI_MAX = math::PI_F * 0.5f;
  const auto point = [](float phi, float theta) {
    return math::Vector(sinf(phi) * cosf(theta), cosf(phi),
                        sinf(phi) * sinf(theta));
  };
  const math::Vector vertices[] = {point(PHI_MAX, -0.2f), point(PHI_MIN, 0.0f),
                                   point(PHI_MAX, 0.2f)};
  const uint16_t indices[] = {0, 1, 2};

  hs_test::StubEffect effect(W, H);
  math::PixelCoords centroid;
  {
    Canvas canvas(effect);
    Pipeline<W, H> pipeline;
    SDF::FaceScratchBuffer scratch;
    SDF::Face face(std::span<const math::Vector>(vertices, 3),
                   std::span<const uint16_t>(indices, 3), scratch,
                   H + hs::H_OFFSET, H, &canvas.clip());
    const SDF::Bounds actual = face.get_vertical_bounds<H>();
    const SDF::Bounds expected = SDF::phi_bounds_to_rows(
        PHI_MIN - SDF::BOUNDS_MARGIN, PHI_MAX + SDF::BOUNDS_MARGIN,
        H + hs::H_OFFSET, H);
    HS_EXPECT_EQ(actual.y_min, expected.y_min);
    HS_EXPECT_EQ(actual.y_max, expected.y_max);
    HS_EXPECT_LT(actual.y_max, H);
    centroid = math::vector_to_pixel<W, H>(face.center);
    Scan::rasterize_face<W, H>(pipeline, canvas, face, face_white);
  }
  effect.advance_display();
  const int x = static_cast<int>(std::lround(centroid.x)) % W;
  const int y = static_cast<int>(std::lround(centroid.y));
  HS_EXPECT_GE(y, 0);
  HS_EXPECT_LT(y, H);
  const Pixel &pixel = effect.get_pixel(x, y);
  HS_EXPECT_GT(static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b, 0U);
}

inline void test_face_rejects_only_the_virtual_south_rows() {
  const auto point = [](float row, float theta) {
    const float phi = row * math::PI_F / (H + hs::H_OFFSET - 1);
    return math::Vector(sinf(phi) * cosf(theta), cosf(phi),
                        sinf(phi) * sinf(theta));
  };
  const uint16_t indices[] = {0, 1, 2};
  for (bool virtual_only : {false, true}) {
    const float row = virtual_only ? H + 0.75f : H - 1.5f;
    const math::Vector vertices[] = {point(row + 0.5f, -0.2f), point(row, 0.0f),
                                     point(row + 0.5f, 0.2f)};
    SDF::FaceScratchBuffer scratch;
    SDF::Face face(std::span<const math::Vector>(vertices, 3),
                   std::span<const uint16_t>(indices, 3), scratch,
                   H + hs::H_OFFSET, H);
    const SDF::Bounds bounds = face.get_vertical_bounds<H>();
    if (virtual_only) {
      HS_EXPECT_EQ(face.count, 0);
      HS_EXPECT_GT(bounds.y_min, bounds.y_max);
      HS_EXPECT_TRUE(face.scratch_owner == nullptr);
    } else {
      HS_EXPECT_EQ(face.count, 3);
      HS_EXPECT_LE(bounds.y_min, H - 1);
      HS_EXPECT_EQ(bounds.y_max, H - 1);
      HS_EXPECT_TRUE(face.scratch_owner == &scratch);
    }
  }
}

/**
 * @brief Sums the tap alphas emitted by one AntiAlias::plot call.
 * @param aa The filter under test.
 * @param x Sub-pixel column coordinate.
 * @param y Sub-pixel row coordinate.
 * @param alpha Input blend alpha.
 * @param[out] tap_count Number of taps emitted (optional; nullptr to ignore).
 * @param taps_on_row Row to count taps landing on (optional; -1 to ignore).
 * @param[out] row_tap_count Count of taps that landed on taps_on_row.
 * @return Sum of the per-tap alphas (the deposited energy).
 */
inline float deposited_energy(Filter::Screen::AntiAlias<W, H> &aa, float x,
                              float y, float alpha, int *tap_count = nullptr,
                              int taps_on_row = -1,
                              int *row_tap_count = nullptr) {
  float sum = 0.0f;
  int count = 0;
  int on_row = 0;
  aa.plot(x, y, Pixel(1, 2, 3), 0.0f, alpha,
          [&](float, float ty, const Pixel &, float, float a) {
            sum += a;
            ++count;
            if (taps_on_row >= 0 && static_cast<int>(ty) == taps_on_row)
              ++on_row;
          });
  if (tap_count)
    *tap_count = count;
  if (row_tap_count)
    *row_tap_count = on_row;
  return sum;
}

/**
 * @brief Pins the offset and LUT used by the Scan, Plot, Feedback and Face cases.
 * @details H_VIRT is H + 3. The last physical row has sin(phi) > 0;
 *          the virtual bottom row samples the pole with float rounding.
 */
inline void test_offset_is_active_and_lut_nondegenerate() {
  using LUT = math::TrigLUT<W, H>;
  const int h_virt = LUT::H_VIRT;
  HS_EXPECT_EQ(hs::H_OFFSET, 3);
  HS_EXPECT_EQ(h_virt, H + 3);

  if (!LUT::initialized)
    LUT::init();

  // Last physical row sits short of the south pole: sin(phi) > 0.
  const float sin_last_phys = LUT::sin_phi[H - 1];
  HS_EXPECT_GT(sin_last_phys, 0.01f);

  // Final virtual row samples sin(PI_F), approximately zero.
  const float sin_virtual_pole = LUT::sin_phi[h_virt - 1];
  HS_EXPECT_NEAR(sin_virtual_pole, 0.0f, 1e-4f);
}

/**
 * @brief Energy conservation across a y-sweep straddling the south-pole clip.
 * @details For every sample whose center row is still on-image (y < H, so the
 *          y0 = floor(y) tap survives), the renorm must redistribute the clipped
 *          Y tap so the deposited alphas sum back to the input alpha. For a
 *          sample fully below the last row (y >= H) every tap is clipped and the
 *          deposited energy is zero — the LEDs stop short of the pole, so the
 *          image is clipped, not stretched. The boundary band [H-1, H) folds
 *          the clipped Y weight onto the surviving row, as in the ideal host profile.
 */
inline void test_energy_conserved_through_clip_boundary() {
  Filter::Screen::AntiAlias<W, H> aa;
  const float in_alpha = 0.8f;

  for (int step = 0; step <= 30; ++step) {
    const float y = static_cast<float>(H - 2) + 0.1f * static_cast<float>(step);
    int count = 0;
    float energy = deposited_energy(aa, 10.37f, y, in_alpha, &count);

    if (y < static_cast<float>(H)) {
      HS_EXPECT_NEAR(energy, in_alpha, 1e-4f);
      HS_EXPECT_GE(count, 1);
      HS_EXPECT_LE(count, 4);
    } else {
      HS_EXPECT_NEAR(energy, 0.0f, 1e-6f);
      HS_EXPECT_EQ(count, 0);
    }
  }
}

/**
 * @brief The boundary row splits across two columns and conserves alpha.
 * @details X weights depend on the framebuffer fraction at every latitude.
 *          The clipped y1 weight folds into y0 (wy0 == 1), as in the ideal host profile.
 */
inline void test_boundary_row_splits_two_columns_and_conserves() {
  Filter::Screen::AntiAlias<W, H> aa;
  const float in_alpha = 1.0f;

  // y in [H-1, H): y0 = H-1 survives, y1 = H is clipped -> the renorm fires.
  const float y = static_cast<float>(H - 1) + 0.4f;
  int count = 0, taps_on_boundary = 0;
  float energy = deposited_energy(aa, 7.5f, y, in_alpha, &count,
                                  /*taps_on_row=*/H - 1, &taps_on_boundary);

  HS_EXPECT_EQ(taps_on_boundary, 2);
  HS_EXPECT_EQ(count, 2);
  HS_EXPECT_NEAR(energy, in_alpha, 1e-4f);
}

/**
 * @brief Energy conservation at the boundary is independent of the X fraction.
 * @details Sweeps the X sub-pixel offset across a full cell at a fixed boundary
 *          y. The renorm sets wy0 = 1, so the deposited energy is wy0 * (sum of
 *          the X weights) = 1 * 1 for every X offset, regardless of how the
 *          quintic-eased split lands between the columns.
 */
inline void test_boundary_energy_independent_of_x_fraction() {
  Filter::Screen::AntiAlias<W, H> aa;
  const float in_alpha = 0.5f;
  const float y = static_cast<float>(H - 1) + 0.25f;

  for (int step = 1; step < 20; ++step) {
    const float xf = 0.05f * static_cast<float>(step);
    float energy = deposited_energy(aa, 12.0f + xf, y, in_alpha);
    HS_EXPECT_NEAR(energy, in_alpha, 1e-4f);
  }
}

/**
 * @brief The last physical row shades as a latitude ring under Scan, not a pole.
 * @details Scan::Shader reconstructs each pixel's direction through the
 *          H_VIRT-aware trig tables. At the legacy offset-3 mapping row H-1 sits at
 *          colatitude (H-1)*PI/(H_VIRT-1), so every column of that row shares
 *          one latitude well off the pole while its azimuth sweeps a full turn.
 *          An H_OFFSET == 0 build maps the row to the south pole up to float
 *          rounding.
 */
inline void test_scan_bottom_row_is_a_latitude_ring() {
  using LUT = math::TrigLUT<W, H>;
  if (!LUT::initialized)
    LUT::init();

  constexpr float FULL_SCALE = 60000.0f;
  hs_test::StubEffect fx(W, H);
  {
    Canvas c(fx);
    Scan::Shader::draw<W, H, 1>(c, [](const math::Vector &v) {
      float theta = std::atan2(v.z, v.x);
      if (theta < 0.0f)
        theta += 2.0f * math::PI_F;
      // Latitude in red, azimuth in green.
      return Color4(
          Pixel(static_cast<uint16_t>((v.y + 1.0f) * 0.5f * FULL_SCALE),
                static_cast<uint16_t>(theta * FULL_SCALE / (2.0f * math::PI_F)),
                0),
          1.0f);
    });
  }
  fx.advance_display();

  const float latitude = (LUT::cos_phi[H - 1] + 1.0f) * 0.5f * FULL_SCALE;
  HS_EXPECT_GT(latitude, 1000.0f);

  int azimuth_steps = 0;
  uint16_t previous = fx.get_pixel(W - 1, H - 1).g;
  for (int x = 0; x < W; ++x) {
    const Pixel &p = fx.get_pixel(x, H - 1);
    HS_EXPECT_NEAR(p.r, latitude, 2.0f);
    if (p.g != previous)
      ++azimuth_steps;
    previous = p.g;
  }
  HS_EXPECT_EQ(azimuth_steps, W);
}

/**
 * @brief Plot maps latitude through H_VIRT, so the sub-pole gap holds no data.
 * @details The LED ring stops H_OFFSET rows short of the south pole. A stroke
 *          past the last physical row therefore has no hardware to light and
 *          must be dropped, while one at the bottom row still draws. On an
 *          H_OFFSET == 0 build both colatitudes land on real rows and the
 *          distinction does not exist.
 */
inline void test_plot_below_last_row_is_clipped() {
  const float bottom_row_phi = math::y_to_phi<H>(static_cast<float>(H - 1));
  const float gap_phi = math::y_to_phi<H>(static_cast<float>(H + 1));

  HS_EXPECT_GT((hs_test::pole_geometry::plot_stroke<W, H>(bottom_row_phi)), 0);
  HS_EXPECT_EQ((hs_test::pole_geometry::plot_stroke<W, H>(gap_phi)), 0);
}

/**
 * @brief Verifies the feedback compositor drives the bottom row as a
 *        mid-latitude ring, not as a pole, at the legacy offset-3 mapping.
 * @details A rotation about +Y moves a mid-latitude ring by its angle and
 *          leaves a pole fixed. On the host (H_OFFSET == 0) row H-1 IS the
 *          south pole, so the compositor pins its warp origin and the row
 *          cannot move; test_filter.h covers that collapse. With the legacy offset-3
 *          mapping the LED ring stops short of the pole, so the same row has to
 *          carry the full longitude shift. Nothing else compiles the feedback
 *          pole path at this offset.
 */
inline void test_feedback_bottom_row_rotates_in_longitude() {
  using LUT = math::TrigLUT<W, H>;
  if (!LUT::initialized)
    LUT::init();
  HS_EXPECT_GT(LUT::sin_phi[H - 1], 0.1f);
  hs_test::pole_geometry::check_feedback_ring_centroid<W, H>(H - 1, 50000 / 4);
}

/**
 * @brief Runs the H_OFFSET renorm module.
 * @return The module's failure count.
 */
inline int run_h_offset_renorm_tests() {
  hs_test::ModuleFixture fixture("h_offset_renorm (HS_TEST_H_OFFSET=3)");
  test_offset_is_active_and_lut_nondegenerate();
  test_energy_conserved_through_clip_boundary();
  test_boundary_row_splits_two_columns_and_conserves();
  test_boundary_energy_independent_of_x_fraction();
  test_scan_bottom_row_is_a_latitude_ring();
  test_plot_below_last_row_is_clipped();
  test_feedback_bottom_row_rotates_in_longitude();
  test_face_bounds_use_virtual_height();
  test_face_rejects_only_the_virtual_south_rows();
  pole_wrap_tests::run_pole_wrap_cases();
  return fixture.result();
}

} // namespace h_offset_renorm
} // namespace hs_test
