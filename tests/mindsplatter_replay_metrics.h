/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/canvas.h"

#include <algorithm>
#include <cstdint>

namespace mindsplatter_replay {

inline constexpr const char *SOURCE_REVISION = "msp-heavy-search-v3";

/** @brief Quadrant clips each searched frame is scored under. */
inline constexpr int SEARCH_CLIP_COUNT = 4;

/**
 * @brief Quadrant clip `index` of a W x H frame.
 * @param index Quadrant id in [0, SEARCH_CLIP_COUNT): bit 0 selects the right
 *        half, bit 1 the bottom half.
 * @return Clip region covering that quadrant.
 * @details Corpus::peak_clip indexes these, so the generator's search and the
 * replay derive the rectangle from this one definition.
 */
template <int W, int H> constexpr ClipRegion search_clip(int index) {
  ClipRegion clip;
  clip.y_start = (index & 2) ? H / 2 : 0;
  clip.y_end = (index & 2) ? H : H / 2;
  clip.x_start = (index & 1) ? W / 2 : 0;
  clip.x_end = (index & 1) ? W : W / 2;
  return clip;
}

/** @brief Per-channel difference metrics, written by every comparator. */
struct ReferenceStats {
  uint64_t framebuffer_hash = 1469598103934665603ull;
  uint64_t expected_hash = 1469598103934665603ull;
  uint64_t total_absolute_error = 0;
  uint32_t changed_pixels = 0;
  uint32_t changed_channels = 0;
  uint16_t max_channel_error = 0;
  ClipRegion clip{};
};

/**
 * @brief Difference metrics against a golden corpus frame.
 * @details Adds the luminance and lit-coverage terms, which need the corpus's
 * expanded sparse framebuffer and so are written only by compare_frame.
 */
struct FrameStats : ReferenceStats {
  uint64_t total_luminance_error = 0;
  uint64_t coverage_luminance_error = 0;
  int64_t luminance_bias = 0;
  uint32_t added_pixels = 0;
  uint32_t dropped_pixels = 0;
  uint16_t max_coverage_luminance = 0;
};

inline uint64_t hash_channel(uint64_t hash, uint16_t channel) {
  hash = (hash ^ static_cast<uint8_t>(channel)) * 1099511628211ull;
  return (hash ^ static_cast<uint8_t>(channel >> 8)) * 1099511628211ull;
}

inline uint16_t linear_luminance(uint16_t r, uint16_t g, uint16_t b) {
  return static_cast<uint16_t>((13933u * r + 46871u * g + 4732u * b + 32768u) >>
                               16);
}

template <int W, typename GoldenPixel>
FrameStats compare_frame(const Pixel *pixels, const ClipRegion &clip,
                         const GoldenPixel *expected, size_t expected_count) {
  FrameStats stats;
  stats.clip = clip;
  size_t expected_index = 0;
  for (int y = clip.y_start; y < clip.y_end; ++y) {
    for (int x = clip.x_start; x < clip.x_end; ++x) {
      const size_t pixel_index = static_cast<size_t>(y) * W + x;
      const Pixel &pixel = pixels[pixel_index];
      while (expected_index < expected_count &&
             expected[expected_index].index < pixel_index)
        ++expected_index;
      const bool expected_lit = expected_index < expected_count &&
                                expected[expected_index].index == pixel_index;
      const uint16_t actual[] = {pixel.r, pixel.g, pixel.b};
      const uint16_t reference[] = {
          expected_lit ? expected[expected_index].r : uint16_t{0},
          expected_lit ? expected[expected_index].g : uint16_t{0},
          expected_lit ? expected[expected_index].b : uint16_t{0},
      };
      const bool actual_lit = (pixel.r | pixel.g | pixel.b) != 0;
      if (actual_lit && !expected_lit)
        ++stats.added_pixels;
      if (!actual_lit && expected_lit)
        ++stats.dropped_pixels;
      const uint16_t actual_luminance =
          linear_luminance(actual[0], actual[1], actual[2]);
      const uint16_t reference_luminance =
          linear_luminance(reference[0], reference[1], reference[2]);
      stats.total_luminance_error +=
          actual_luminance > reference_luminance
              ? actual_luminance - reference_luminance
              : reference_luminance - actual_luminance;
      stats.luminance_bias += static_cast<int32_t>(actual_luminance) -
                              static_cast<int32_t>(reference_luminance);
      if (actual_lit != expected_lit) {
        const uint16_t coverage_error =
            actual_lit ? actual_luminance : reference_luminance;
        stats.coverage_luminance_error += coverage_error;
        stats.max_coverage_luminance =
            std::max(stats.max_coverage_luminance, coverage_error);
      }
      bool changed = false;
      for (int channel = 0; channel < 3; ++channel) {
        stats.framebuffer_hash =
            hash_channel(stats.framebuffer_hash, actual[channel]);
        stats.expected_hash =
            hash_channel(stats.expected_hash, reference[channel]);
        const uint16_t error = actual[channel] > reference[channel]
                                   ? actual[channel] - reference[channel]
                                   : reference[channel] - actual[channel];
        if (error == 0)
          continue;
        changed = true;
        ++stats.changed_channels;
        stats.total_absolute_error += error;
        stats.max_channel_error = std::max(stats.max_channel_error, error);
      }
      if (changed)
        ++stats.changed_pixels;
    }
  }
  return stats;
}

template <int W, typename GoldenPixel>
FrameStats compare_frame(Canvas &canvas, const ClipRegion &clip,
                         const GoldenPixel *expected, size_t expected_count) {
  return compare_frame<W>(&canvas(0, 0), clip, expected, expected_count);
}

inline size_t frame_pixel_count(const ClipRegion &clip) {
  return static_cast<size_t>(clip.x_end - clip.x_start) *
         static_cast<size_t>(clip.y_end - clip.y_start);
}

template <int W>
size_t capture_frame(const Pixel *pixels, const ClipRegion &clip, Pixel *output,
                     size_t capacity) {
  const size_t count = frame_pixel_count(clip);
  HS_CHECK(count <= capacity,
           "MindSplatter replay reference buffer is too small");
  size_t output_index = 0;
  for (int y = clip.y_start; y < clip.y_end; ++y)
    for (int x = clip.x_start; x < clip.x_end; ++x)
      output[output_index++] = pixels[static_cast<size_t>(y) * W + x];
  return count;
}

template <int W>
ReferenceStats
compare_frame_reference(const Pixel *pixels, const ClipRegion &clip,
                        const Pixel *reference, size_t reference_count) {
  HS_CHECK(reference_count == frame_pixel_count(clip),
           "MindSplatter replay reference dimensions differ");
  ReferenceStats stats;
  stats.clip = clip;
  size_t reference_index = 0;
  for (int y = clip.y_start; y < clip.y_end; ++y) {
    for (int x = clip.x_start; x < clip.x_end; ++x) {
      const Pixel &pixel = pixels[static_cast<size_t>(y) * W + x];
      const Pixel &expected = reference[reference_index++];
      const uint16_t actual[] = {pixel.r, pixel.g, pixel.b};
      const uint16_t baseline[] = {expected.r, expected.g, expected.b};
      bool changed = false;
      for (int channel = 0; channel < 3; ++channel) {
        stats.framebuffer_hash =
            hash_channel(stats.framebuffer_hash, actual[channel]);
        stats.expected_hash =
            hash_channel(stats.expected_hash, baseline[channel]);
        const uint16_t error = actual[channel] > baseline[channel]
                                   ? actual[channel] - baseline[channel]
                                   : baseline[channel] - actual[channel];
        if (error == 0)
          continue;
        changed = true;
        ++stats.changed_channels;
        stats.total_absolute_error += error;
        stats.max_channel_error = std::max(stats.max_channel_error, error);
      }
      if (changed)
        ++stats.changed_pixels;
    }
  }
  return stats;
}

template <int W>
ReferenceStats compare_frame_reference(Canvas &canvas, const ClipRegion &clip,
                                       const Pixel *reference,
                                       size_t reference_count) {
  return compare_frame_reference<W>(&canvas(0, 0), clip, reference,
                                    reference_count);
}

inline void clear_frame(Canvas &canvas, const ClipRegion &clip) {
  for (int y = clip.y_start; y < clip.y_end; ++y)
    for (int x = clip.x_start; x < clip.x_end; ++x)
      canvas(x, y) = Pixel(0, 0, 0);
}

} // namespace mindsplatter_replay
