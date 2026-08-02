/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/canvas.h"

#include <algorithm>
#include <cstdint>

namespace mindsplatter_replay {

struct FrameStats {
  uint64_t framebuffer_hash = 1469598103934665603ull;
  uint64_t expected_hash = 1469598103934665603ull;
  uint64_t total_absolute_error = 0;
  uint32_t changed_pixels = 0;
  uint32_t changed_channels = 0;
  uint16_t max_channel_error = 0;
  ClipRegion clip{};
};

inline uint64_t hash_channel(uint64_t hash, uint16_t channel) {
  hash = (hash ^ static_cast<uint8_t>(channel)) * 1099511628211ull;
  return (hash ^ static_cast<uint8_t>(channel >> 8)) * 1099511628211ull;
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
FrameStats compare_frame_reference(const Pixel *pixels, const ClipRegion &clip,
                                   const Pixel *reference,
                                   size_t reference_count) {
  HS_CHECK(reference_count == frame_pixel_count(clip),
           "MindSplatter replay reference dimensions differ");
  FrameStats stats;
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
FrameStats compare_frame_reference(Canvas &canvas, const ClipRegion &clip,
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
