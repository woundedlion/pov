/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared framebuffer predicates for the suites that assert on a StubEffect's
 * pixels — canvas, filter, scan, and mesh raster.
 *
 * The two lit-pixel counters differ in the region they scan and are not
 * interchangeable: count_lit_canvas() takes the effect's own reported
 * dimensions, count_lit_region() a compile-time window that a test states
 * independently of the effect it built.
 */
#pragma once

#include "core/color/color.h"
#include "core/render/canvas.h"
#include "tests/test_fixture.h"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hs_test {

/**
 * @brief Tests whether a pixel is fully black (cleared-frame / unwritten).
 * @param p Pixel to inspect.
 * @return True when all of the pixel's RGB channels are zero.
 */
inline bool is_black(const Pixel &p) {
  return p.r == 0 && p.g == 0 && p.b == 0;
}

/**
 * @brief Tests whether a pixel has exactly the given channel values.
 * @param p Pixel to inspect.
 * @param r Expected red channel value.
 * @param g Expected green channel value.
 * @param b Expected blue channel value.
 * @return True if all three channels match exactly.
 */
inline bool pix_eq(const Pixel &p, uint16_t r, uint16_t g, uint16_t b) {
  return p.r == r && p.g == g && p.b == b;
}

#define HS_EXPECT_PIXEL(actual, red, green, blue)                              \
  HS_EXPECT_EQ((actual), Pixel((red), (green), (blue)))

/**
 * @brief Absolute gap between two 16-bit channel values.
 * @param a First channel value.
 * @param b Second channel value.
 * @return |a - b|.
 */
inline uint16_t channel_gap(uint16_t a, uint16_t b) {
  return a > b ? a - b : b - a;
}

/**
 * @brief Largest per-channel gap between two pixels.
 * @param a First pixel.
 * @param b Second pixel.
 * @return Maximum of the three RGB channel gaps, 16-bit scale.
 */
inline uint16_t max_channel_gap(const Pixel &a, const Pixel &b) {
  return std::max(
      {channel_gap(a.r, b.r), channel_gap(a.g, b.g), channel_gap(a.b, b.b)});
}

/**
 * @brief Asserts two shaded colors agree within a per-channel tolerance.
 * @param actual Color produced by the code under test.
 * @param expected Reference color.
 * @param max_gap Largest RGB channel gap allowed, 16-bit scale.
 * @param alpha_tolerance Absolute tolerance on alpha.
 */
inline void expect_color_within(const Color4 &actual, const Color4 &expected,
                                uint16_t max_gap, float alpha_tolerance) {
  HS_EXPECT_LE(max_channel_gap(actual.color, expected.color), max_gap);
  HS_EXPECT_NEAR(actual.alpha, expected.alpha, alpha_tolerance);
}

/**
 * @brief Running per-channel gap statistics over a sweep of pixel pairs.
 * @details Every RGB channel of every pair counts once; alpha is not folded in.
 */
struct ChannelError {
  uint16_t max = 0;
  uint64_t total = 0;
  uint64_t channels = 0;

  /**
   * @brief Folds one pixel pair's three channel gaps into the statistics.
   * @param a First pixel.
   * @param b Second pixel.
   */
  void add(const Pixel &a, const Pixel &b) {
    for (uint16_t gap : {channel_gap(a.r, b.r), channel_gap(a.g, b.g),
                         channel_gap(a.b, b.b)}) {
      max = std::max(max, gap);
      total += gap;
    }
    channels += 3;
  }

  /**
   * @brief Mean gap per channel, truncated toward zero.
   * @return total / channels.
   */
  uint64_t mean() const {
    HS_CHECK(channels > 0, "channel error mean requires samples");
    return total / channels;
  }
};

/**
 * @brief Counts the non-black pixels across the effect's reported canvas.
 * @param fx Effect whose framebuffer is scanned, fx.width() by fx.height().
 * @return Number of lit (non-black) pixels.
 */
inline size_t count_lit_canvas(const StubEffect &fx) {
  size_t n = 0;
  for (int y = 0; y < fx.height(); ++y)
    for (int x = 0; x < fx.width(); ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++n;
  return n;
}

/**
 * @brief Counts the non-black pixels across a fixed W by H canvas region.
 * @tparam W Region width in pixels.
 * @tparam H Region height in pixels.
 * @param fx Effect whose framebuffer is scanned.
 * @return Number of lit (non-black) pixels in the region.
 */
template <int W, int H> inline size_t count_lit_region(const StubEffect &fx) {
  size_t n = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++n;
  return n;
}

/**
 * @brief Copies an effect's displayed frame into a flat row-major buffer.
 * @tparam W Frame width in pixels.
 * @tparam H Frame height in pixels.
 * @tparam Fx Effect type exposing get_pixel(x, y) const.
 * @param fx Effect whose displayed frame is read.
 * @param out Destination, resized to W * H and indexed y * W + x.
 * @details Reads through get_pixel rather than the raw buffer, so an effect
 *          that overrides it with a per-pixel transform is captured as displayed.
 */
template <int W, int H, typename Fx>
inline void capture_frame(const Fx &fx, std::vector<Pixel> &out) {
  out.resize(static_cast<size_t>(W) * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      out[static_cast<size_t>(y) * W + x] = fx.get_pixel(x, y);
}

/**
 * @brief Sums every RGB channel across a fixed W by H framebuffer window.
 * @tparam W Window width in pixels.
 * @tparam H Window height in pixels.
 * @tparam Fx Effect type exposing get_pixel(x, y) const.
 * @param fx Effect whose displayed frame is read.
 * @return Channel sum; zero means the window is all black.
 */
template <int W, int H, typename Fx>
inline uint64_t frame_energy(const Fx &fx) {
  uint64_t energy = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &p = fx.get_pixel(x, y);
      energy += static_cast<uint64_t>(p.r) + p.g + p.b;
    }
  return energy;
}

/**
 * @brief Sums every RGB channel of a captured frame.
 * @param frame Pixels as capture_frame() lays them out.
 * @return Channel sum; zero means the frame is all black.
 */
inline uint64_t frame_energy(const std::vector<Pixel> &frame) {
  uint64_t energy = 0;
  for (const Pixel &p : frame)
    energy += static_cast<uint64_t>(p.r) + p.g + p.b;
  return energy;
}

} // namespace hs_test
