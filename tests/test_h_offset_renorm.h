/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * South-pole Y-clip renormalization coverage for the device H_OFFSET path.
 *
 * Compiled with -DHS_TEST_H_OFFSET=3 so the entire pipeline (PhiLUT / TrigLUT /
 * the AntiAlias renorm) is built with the device geometry. On a host build
 * H_OFFSET == 0, the bottom row maps exactly to the south pole (sin(phi) == 0)
 * and the renorm path is never exercised.
 */
#pragma once

#include "core/render/filter.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace h_offset_renorm {

// With HS_TEST_H_OFFSET=3, H_VIRT = H + 3 and the bottom physical row y=H-1
// lands SHORT of the south pole (sin(phi) > 0).
constexpr int W = 32;
constexpr int H = 16;

/**
 * @brief Sums the tap alphas emitted by one AntiAlias::plot call.
 * @param aa The filter under test.
 * @param x Sub-pixel column coordinate.
 * @param y Sub-pixel row coordinate.
 * @param alpha Input blend alpha.
 * @param[out] tap_count Number of taps emitted (optional; nullptr to ignore).
 * @param[out] taps_on_row Row to count taps landing on (optional; -1 to ignore).
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
 * @brief Pins that this TU really compiled with the device offset, so the LUT
 *        the renorm reads is non-degenerate (the host-build precondition the
 *        whole file rests on).
 * @details H_VIRT must be H + 3, and the last PHYSICAL row (y = H-1) must carry
 *          sin(phi) > 0 — on an H_OFFSET == 0 build it would be sin(PI) == 0 and
 *          the two-column split below could never happen. The virtual bottom row
 *          (H_VIRT-1) still reaches the pole exactly (sin == 0).
 */
inline void test_offset_is_active_and_lut_nondegenerate() {
  using LUT = TrigLUT<W, H>;
  const int h_virt = LUT::H_VIRT;
  HS_EXPECT_EQ(hs::H_OFFSET, 3);
  HS_EXPECT_EQ(h_virt, H + 3);

  if (!LUT::initialized)
    LUT::init();

  // Last physical row sits short of the south pole: sin(phi) > 0.
  const float sin_last_phys = LUT::sin_phi[H - 1];
  HS_EXPECT_GT(sin_last_phys, 0.01f);

  // Final virtual row reaches the pole: sin(PI) == 0.
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
 *          image is clipped, not stretched. The boundary band [H-1, H) is where
 *          the device-only fold actually fires.
 */
inline void test_energy_conserved_through_clip_boundary() {
  Filter::Screen::AntiAlias<W, H> aa;
  const float in_alpha = 0.8f;

  for (float y = static_cast<float>(H - 2); y < static_cast<float>(H + 1) + 1e-3f;
       y += 0.1f) {
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
 * @brief The boundary row genuinely splits across two columns under the renorm.
 * @details This is the fact an H_OFFSET == 0 build cannot reproduce: there the
 *          last physical row has sin(phi) == 0, so the X fractional collapses to
 *          a single column (x_frac -> 0) and only one tap lands on the row. With
 *          the device offset the row has sin(phi) > 0, so a sample with a clear
 *          X fraction splits across columns x0 and x1 — and the two taps still
 *          sum to the input alpha because the clipped y1 row's weight was folded
 *          into the surviving y0 row (wy0 == 1).
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
 *          quintic-eased, sin(phi)-compensated split lands between the columns.
 */
inline void test_boundary_energy_independent_of_x_fraction() {
  Filter::Screen::AntiAlias<W, H> aa;
  const float in_alpha = 0.5f;
  const float y = static_cast<float>(H - 1) + 0.25f;

  for (float xf = 0.05f; xf < 1.0f; xf += 0.05f) {
    float energy = deposited_energy(aa, 12.0f + xf, y, in_alpha);
    HS_EXPECT_NEAR(energy, in_alpha, 1e-4f);
  }
}

/**
 * @brief Runs the H_OFFSET renorm module.
 * @return The module's failure count.
 */
inline int run_h_offset_renorm_tests() {
  hs_test::ModuleFixture fixture("h_offset_renorm");
  test_offset_is_active_and_lut_nondegenerate();
  test_energy_conserved_through_clip_boundary();
  test_boundary_row_splits_two_columns_and_conserves();
  test_boundary_energy_independent_of_x_fraction();
  return fixture.result();
}

} // namespace h_offset_renorm
} // namespace hs_test
