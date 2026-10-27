/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include "core/render/spacetime.h"
#include "tests/test_harness.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <optional>
#include <type_traits>
#include <vector>

namespace hs_test {
/** @brief Captures row-major linear RGB luminance from a displayed frame.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam E Effect type.
 * @param effect Effect whose displayed buffer is sampled.
 * @return Normalized mean RGB channel values. */
template <int W, int H, typename E>
std::vector<float> capture_luminance(const E &effect) {
  std::vector<float> luminance(static_cast<size_t>(W) * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const auto pixel = effect.get_pixel(x, y);
      luminance[static_cast<size_t>(y) * W + x] =
          (pixel.r + pixel.g + pixel.b) / (3.0f * 65535.0f);
    }
  return luminance;
}

/** @brief Displaces a reference vertex along its meridian without projection.
 * @tparam H Display height.
 * @param point Unit world-space vertex.
 * @param rows Signed latitude displacement in display rows.
 * @return Displaced unit vertex; zero displacement preserves every component. */
template <int H>
math::Vector displace_reference_rows(const math::Vector &point, float rows) {
  if (rows == 0.0f)
    return point;
  const double angle = double(rows) * math::RADIANS_PER_ROW<H>;
  const double sine = std::sin(angle), cosine = std::cos(angle);
  const double radius = std::hypot(double(point.x), double(point.z));
  const double displaced_radius = radius * cosine + double(point.y) * sine;
  const double y = double(point.y) * cosine - radius * sine;
  if (radius == 0.0)
    return {float(displaced_radius), float(y), 0.0f};
  const double scale = displaced_radius / radius;
  return {float(double(point.x) * scale), float(y),
          float(double(point.z) * scale)};
}

/** @brief Reference row shifts match analytic latitude and preserve identity. */
inline void test_reference_row_displacement_is_analytic() {
  constexpr int H = 144;
  for (float theta : {0.137f, 0.7f, 2.2f, 4.7f})
    for (float phi : {0.001f, 0.4f, 1.2f, 3.1405f}) {
      const math::Vector point = math::Vector::from_spherical(theta, phi);
      const auto identity = displace_reference_rows<H>(point, 0.0f);
      HS_EXPECT_EQ(identity.x, point.x);
      HS_EXPECT_EQ(identity.y, point.y);
      HS_EXPECT_EQ(identity.z, point.z);
      for (float rows : {-0.2f, -0.001f, 0.001f, 0.2f}) {
        const auto shifted = displace_reference_rows<H>(point, rows);
        const double latitude =
            double(phi) + double(rows) * math::RADIANS_PER_ROW<H>;
        HS_EXPECT_NEAR(shifted.x, std::sin(latitude) * std::cos(double(theta)),
                       2e-6);
        HS_EXPECT_NEAR(shifted.y, std::cos(latitude), 2e-6);
        HS_EXPECT_NEAR(shifted.z, std::sin(latitude) * std::sin(double(theta)),
                       2e-6);
        HS_EXPECT_NEAR(math::dot(shifted, shifted), 1.0f, 1e-6f);
      }
    }
}

/** @brief Tests the whole-frame temporal error against the endpoint budget.
 * @param error Summed luminance error over the frame.
 * @param endpoint_error Summed endpoint-render error over the frame.
 * @param share Fraction of the twentieth-error limit the frame may use.
 * @return Whether the temporal frame satisfies the scaled limit. */
inline bool temporal_frame_error_within_budget(float error,
                                               float endpoint_error,
                                               float share = 1.0f) {
  return error <= share * endpoint_error / 20.0f + 1e-5f;
}

/** @brief Compares columns and whole-frame error with frozen geometry.
 * @param render Timed or frozen render callback.
 * @param name Diagnostic label.
 * @param share Fraction of the twentieth-error frame limit allowed.
 * @param bidirectional Uses each column's larger error under positive and
 * negative uniform 0.2-pixel displacements. */
template <int W, int H, typename Render>
void seam_compare_render(Render render, const char *name, float share = 1.0f,
                         bool bidirectional = false) {
  const auto actual = render(true, 1.0f, 0.0f);
  const auto endpoint = render(false, 1.0f, 0.0f);
  float total = 0.0f, endpoint_total = 0.0f, displaced_total = 0.0f;
  float worst = 0.0f;
  int displacement_failures = 0;
  for (int column = 0; column < W / 2; ++column) {
    const float fraction = Spacetime::SamplingSurface<W>::time(column);
    const auto frozen = render(false, fraction, 0.0f);
    const auto shifted = render(false, fraction, 0.2f);
    std::optional<std::remove_cvref_t<decltype(shifted)>> opposite;
    if (bidirectional)
      opposite.emplace(render(false, fraction, -0.2f));
    for (int x : {column, column + W / 2}) {
      float error = 0.0f, endpoint_error = 0.0f, displaced_error = 0.0f;
      float opposite_error = 0.0f;
      for (int y = 0; y < H; ++y) {
        const size_t index = y * W + x;
        error += fabsf(actual[index] - frozen[index]);
        endpoint_error += fabsf(endpoint[index] - frozen[index]);
        displaced_error += fabsf(shifted[index] - frozen[index]);
        if (opposite)
          opposite_error += fabsf((*opposite)[index] - frozen[index]);
      }
      displaced_error = std::max(displaced_error, opposite_error);
      total += error;
      endpoint_total += endpoint_error;
      displaced_total += displaced_error;
      worst = std::max(worst, error);
      if (error > displaced_error + 1e-5f)
        std::printf("  frozen %s column %d: error=%.9f displaced=%.9f\n", name,
                    x, error, displaced_error);
      displacement_failures += error > displaced_error + 1e-5f;
    }
  }
  std::printf("  frozen %s: error=%.6f endpoint=%.6f displaced=%.6f "
              "worst_column=%.6f frame_limit=%.6f displacement_failures=%d\n",
              name, total, endpoint_total, displaced_total, worst,
              endpoint_total / 20.0f, displacement_failures);
  HS_EXPECT_GT(endpoint_total, 0.01f);
  HS_EXPECT(temporal_frame_error_within_budget(total, endpoint_total, share),
            "whole-frame error is within its share of one twentieth of "
            "endpoint error");
  HS_EXPECT_EQ(displacement_failures, 0);
}

/** @brief The endpoint budget is shared across columns and retains its scale. */
inline void test_temporal_oracle_uses_whole_frame_endpoint_budget() {
  HS_EXPECT(temporal_frame_error_within_budget(0.04f, 1.0f),
            "error within one twentieth");
  HS_EXPECT(!temporal_frame_error_within_budget(0.06f, 1.0f),
            "error above one twentieth");
  HS_EXPECT(!temporal_frame_error_within_budget(0.04f, 1.0f, 0.5f),
            "error above half of one twentieth");
  auto render = [](bool timed, float fraction, float displacement) {
    std::array<float, 4> pixels{};
    if (timed)
      pixels[0] = 0.02f;
    else if (displacement != 0.0f)
      pixels.fill(0.1f);
    else if (fraction == 1.0f)
      pixels[1] = 1.0f;
    return pixels;
  };
  seam_compare_render<4, 1>(render, "distributed endpoint budget");
}

} // namespace hs_test
