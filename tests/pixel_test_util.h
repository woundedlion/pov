/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared pixel comparison and framebuffer helpers for the native test suites.
 */
#pragma once

#include "core/color/color.h"
#include "core/render/canvas.h"
#include "core/render/spacetime.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <span>
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
 * @brief Counts the non-black pixels in a pixel span.
 * @param pixels Pixels to scan.
 * @return Number of lit (non-black) pixels.
 */
inline size_t count_lit(std::span<const Pixel> pixels) {
  size_t n = 0;
  for (const Pixel &p : pixels)
    n += !is_black(p);
  return n;
}

/**
 * @brief Counts the non-black pixels across the effect's reported canvas.
 * @param fx Effect whose framebuffer is scanned, fx.width() by fx.height().
 * @return Number of lit (non-black) pixels.
 */
template <typename Fx> inline size_t count_lit_canvas(const Fx &fx) {
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
template <int W, int H, typename Fx>
inline size_t count_lit_region(const Fx &fx) {
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
 * @details Reads through get_pixel, so a per-pixel display transform is
 *          captured.
 */
template <int W, int H, typename Fx>
inline void capture_frame(const Fx &fx, std::vector<Pixel> &out) {
  out.resize(static_cast<size_t>(W) * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      out[static_cast<size_t>(y) * W + x] = fx.get_pixel(x, y);
}

/** @brief Folds a displayed frame into a running channel hash.
 * @tparam W Frame width.
 * @tparam H Frame height.
 * @tparam Fx Effect type.
 * @param fx Effect to read.
 * @param hash Initial running hash.
 * @return Hash extended by every channel in row-major order. */
template <int W, int H, typename Fx>
inline uint64_t frame_hash(const Fx &fx, uint64_t hash = HASH_SEED64) {
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &pixel = fx.get_pixel(x, y);
      for (uint16_t channel : {pixel.r, pixel.g, pixel.b})
        hash = fnv1a64_channel(hash, channel);
    }
  return hash;
}

/**
 * @brief Whether committed pixel goldens apply to this build.
 * @details Goldens are captured by the native Windows Clang Debug build. Other
 * C libraries and floating-point modes round transcendental results
 * differently.
 */
#if defined(_WIN32) && !defined(__OPTIMIZE__) && !HS_TEST_FAST_MATH
inline constexpr bool PIXEL_GOLDENS = true;
#else
inline constexpr bool PIXEL_GOLDENS = false;
#endif

/**
 * @brief Compares a render hash with its committed golden where goldens apply.
 * @param name Scenario name printed beside the hash.
 * @param hash Hash of the rendered frames.
 * @param golden Committed hash of the same frames.
 */
inline void expect_pixel_golden(const char *name, uint64_t hash,
                                uint64_t golden) {
  std::printf("  %s: pixel hash 0x%016llx\n", name,
              static_cast<unsigned long long>(hash));
  if constexpr (PIXEL_GOLDENS)
    HS_EXPECT_EQ(hash, golden);
}

#if defined(__EMSCRIPTEN__) || HS_ENABLE_TEST_HOOKS
/** @brief Selects column-time sampling for a scope, then restores it. */
struct ScopedTemporalSampling {
  /** @brief Applies the mode.
   * @param enabled Whether columns use their emission times. */
  explicit ScopedTemporalSampling(bool enabled)
      : entering_flat(Spacetime::flat_surface()) {
    Spacetime::set_flat_surface(!enabled);
  }
  ~ScopedTemporalSampling() { Spacetime::set_flat_surface(entering_flat); }
  ScopedTemporalSampling(const ScopedTemporalSampling &) = delete;
  ScopedTemporalSampling &operator=(const ScopedTemporalSampling &) = delete;

private:
  bool entering_flat; /**< Flat-surface mode restored on scope exit. */
};
#endif

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

/** @brief Frame hash and channel energy of a pinned render sequence. */
struct RenderPin {
  uint64_t hash = HASH_SEED64; ///< frame_hash() chained over the frames.
  uint64_t energy = 0;         ///< frame_energy() summed over the frames.
};

/**
 * @brief Extends a render pin by an effect's displayed frame.
 * @tparam W Frame width.
 * @tparam H Frame height.
 * @tparam Fx Effect type.
 * @param pin Pin to extend.
 * @param fx Effect to read.
 */
template <int W, int H, typename Fx>
inline void fold_render_pin(RenderPin &pin, const Fx &fx) {
  pin.hash = frame_hash<W, H>(fx, pin.hash);
  pin.energy += frame_energy<W, H>(fx);
}

/**
 * @brief Checks a render against its pin, hash first.
 * @param rendered Pin folded from the render under test.
 * @param pin Recorded pin.
 * @param relative_band Energy band accepted when the hash moved.
 * @details A moved hash, as from a libm whose last bits differ from the
 * recording toolchain's, falls back to the channel energy within the band and
 * reports the fallback.
 */
inline void expect_render_pin(const RenderPin &rendered, const RenderPin &pin,
                              double relative_band = 1e-3) {
  HS_EXPECT_GT(pin.energy, uint64_t{0});
  if (rendered.hash == pin.hash) {
    HS_EXPECT_EQ(rendered.energy, pin.energy);
    return;
  }
  std::printf("  render pin hash moved (%llu, energy %llu); checking energy\n",
              static_cast<unsigned long long>(rendered.hash),
              static_cast<unsigned long long>(rendered.energy));
  HS_EXPECT_NEAR(static_cast<double>(rendered.energy),
                 static_cast<double>(pin.energy),
                 relative_band * static_cast<double>(pin.energy));
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

/** @brief Exact and aggregate channel-error statistics for two captures. */
struct FrameErrorStats {
  size_t compared_pixels = 0;
  size_t lit_reference_pixels = 0; /**< Reference pixels that are not black. */
  size_t coverage_mismatches = 0;  /**< Pixels black in exactly one capture. */
  size_t different_pixels = 0;
  size_t different_channels = 0;
  uint32_t max_absolute_error = 0;
  uint64_t total_absolute_error = 0;
  uint64_t total_squared_error = 0;
  int first_different_pixel = -1; /**< Index across all compared spans. */

  bool exact() const { return different_pixels == 0; }

  double mean_absolute_error() const {
    const size_t channels = compared_pixels * 3;
    return channels == 0 ? 0.0
                         : static_cast<double>(total_absolute_error) / channels;
  }

  double root_mean_squared_error() const {
    const size_t channels = compared_pixels * 3;
    return channels == 0
               ? 0.0
               : std::sqrt(static_cast<double>(total_squared_error) / channels);
  }
};

/**
 * @brief Compares two same-sized pixel spans channel by channel.
 * @param reference Reference pixels.
 * @param candidate Candidate pixels, same length as reference.
 * @param stats Statistics to accumulate into, for comparing several spans.
 * @return stats extended by this span pair.
 */
inline FrameErrorStats compare_buffers(std::span<const Pixel> reference,
                                       std::span<const Pixel> candidate,
                                       FrameErrorStats stats = {}) {
  HS_CHECK(reference.size() == candidate.size(), "captures differ in size");
  const size_t base = stats.compared_pixels;
  stats.compared_pixels += reference.size();
  for (size_t i = 0; i < reference.size(); ++i) {
    const Pixel &a = reference[i];
    const Pixel &b = candidate[i];
    stats.lit_reference_pixels += !is_black(a);
    stats.coverage_mismatches += is_black(a) != is_black(b);
    bool pixel_differs = false;
    for (const uint32_t error : {channel_gap(a.r, b.r), channel_gap(a.g, b.g),
                                 channel_gap(a.b, b.b)}) {
      if (error == 0)
        continue;
      pixel_differs = true;
      ++stats.different_channels;
      stats.max_absolute_error = std::max(stats.max_absolute_error, error);
      stats.total_absolute_error += error;
      stats.total_squared_error += static_cast<uint64_t>(error) * error;
    }
    if (!pixel_differs)
      continue;
    if (stats.first_different_pixel < 0)
      stats.first_different_pixel = static_cast<int>(base + i);
    ++stats.different_pixels;
  }
  return stats;
}

/**
 * @brief Checks two frames channel by channel within a bound.
 * @param label Draw named in failure context.
 * @param actual Rendered frame.
 * @param expected Oracle frame.
 * @param bound Largest channel difference allowed.
 * @return Number of lit oracle pixels, so a vacuous comparison is visible.
 */
inline int expect_frames_within(const char *label,
                                const std::vector<Pixel> &actual,
                                const std::vector<Pixel> &expected, int bound) {
  HS_CONTEXT(label);
  int maximum_error = 0;
  int lit = 0;
  for (size_t i = 0; i < expected.size(); ++i) {
    lit += !is_black(expected[i]);
    const Pixel &a = actual[i];
    const Pixel &b = expected[i];
    for (int error : {std::abs(int(a.r) - b.r), std::abs(int(a.g) - b.g),
                      std::abs(int(a.b) - b.b)})
      maximum_error = std::max(maximum_error, error);
  }
  HS_EXPECT_LE(maximum_error, bound);
  return lit;
}

/**
 * @brief Assembles a frame whose column pairs come from renders of the state
 *        at each pair's sweep time.
 * @tparam W Frame width.
 * @tparam H Frame height.
 * @param render Callable `std::vector<Pixel>(float u)` rendering the whole
 *        frame with every column at frame fraction `u`.
 * @return Frame whose columns x and x + W / 2 come from @p render at
 *         `u = (x % (W / 2)) / (W / 2)`.
 */
template <int W, int H, typename RenderFn>
inline std::vector<Pixel> column_time_oracle(RenderFn &&render) {
  std::vector<Pixel> out(static_cast<size_t>(W) * H);
  for (int column = 0; column < W / 2; ++column) {
    const std::vector<Pixel> frame =
        render(static_cast<float>(column) / (W / 2));
    for (int y = 0; y < H; ++y)
      for (int x : {column, column + W / 2})
        out[static_cast<size_t>(y) * W + x] =
            frame[static_cast<size_t>(y) * W + x];
  }
  return out;
}

/** @brief Largest share of a held frame's seam error a deformed-mesh frame
 *  keeps. */
inline constexpr double MOTION_SEAM_RATIO = 0.4;

/** @brief Smallest channel error a held frame must show at its seam columns
 *  for a seam comparison to bite. */
inline constexpr int MOTION_SEAM_FLOOR = 8192;

/**
 * @brief First-arm columns a seam check compares: the first and last columns
 *        of the arm, where held state is a whole frame off, and the middle.
 * @tparam W Frame width.
 */
template <int W>
inline constexpr std::array<int, 12> SEAM_COLUMNS = {
    0,     1,         2,         3,         W / 4 - 2, W / 4 - 1,
    W / 4, W / 4 + 1, W / 2 - 4, W / 2 - 3, W / 2 - 2, W / 2 - 1};

/**
 * @brief Checks that a frame whose state moves over the frame shows each
 *        seam column's state at the column's sweep time.
 * @tparam W Frame width.
 * @tparam H Frame height.
 * @tparam RenderFn Callable `std::vector<Pixel>(float u)` rendering the whole
 *         frame with the state frozen at frame fraction `u`.
 * @param label Draw named in failure context.
 * @param drawn The frame under test.
 * @param held The same frame drawn with its state held at the current
 *        endpoint.
 * @param render The independent time evaluator.
 * @param limit Largest drawn-to-held error ratio accepted.
 * @return The drawn frame's summed seam error over the held frame's.
 * @details Sums channel errors against the evaluator over SEAM_COLUMNS and
 * their opposing columns. The held frame's error grows to a whole frame of
 * motion at the first columns of each arm; the drawn frame must keep at most
 * @p limit of it.
 */
template <int W, int H, typename RenderFn>
inline double
expect_seam_resolved(const char *label, const std::vector<Pixel> &drawn,
                     const std::vector<Pixel> &held, RenderFn &&render,
                     double limit = MOTION_SEAM_RATIO) {
  HS_CONTEXT(label);
  double drawn_error = 0.0, held_error = 0.0;
  int held_peak = 0;
  for (const int column : SEAM_COLUMNS<W>) {
    const std::vector<Pixel> expected =
        render(static_cast<float>(column) / (W / 2));
    for (int y = 0; y < H; ++y)
      for (int x : {column, column + W / 2}) {
        const size_t i = static_cast<size_t>(y) * W + x;
        const auto error = [&](const Pixel &a) {
          return std::abs(int(a.r) - expected[i].r) +
                 std::abs(int(a.g) - expected[i].g) +
                 std::abs(int(a.b) - expected[i].b);
        };
        drawn_error += error(drawn[i]);
        held_error += error(held[i]);
        held_peak = std::max(held_peak, error(held[i]));
      }
  }
  const double ratio = held_error > 0.0 ? drawn_error / held_error : 1.0;
  std::printf("  [%s] seam error %.0f vs held %.0f (ratio %.3f)\n", label,
              drawn_error, held_error, ratio);
  HS_EXPECT_GT(held_peak, MOTION_SEAM_FLOOR);
  HS_EXPECT_LE(ratio, limit);
  return ratio;
}

} // namespace hs_test
