/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke.h
 * @brief Four-stroke contour classification, presentation clock and linear-light resolve.
 */

#include "core/color/pixel.h"
#include "core/math/rotation_math.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>

namespace FourStroke {

/// Edge contrast as a function of normal alignment.
enum class Mode { REFERENCE_POLARITY, DIRECTIONAL_CONTRAST };

/// Authored material settings in output-pixel units.
struct Controls {
  float heading = 37; ///< Clockwise degrees from increasing column.
  float width = 1;    ///< Inner band width in output pixels.
  float strength = 1; ///< Straight linear RGB mix.
  Mode mode = Mode::DIRECTIONAL_CONTRAST; ///< Angular contrast law.
};

/// Unit direction in the output-pixel basis.
struct Direction {
  float x; ///< Column component.
  float y; ///< Row component.
};

/** @brief Prepares a periodic clockwise output heading.
 * @param degrees Finite heading in degrees.
 * @return Unit column/row vector. */
inline Direction heading(float degrees) {
  const float radians = std::remainder(degrees, 360.0f) * (math::PI_F / 180.0f);
  return {std::cos(radians), std::sin(radians)};
}

/** @brief Normalizes a material's output-direction override.
 * @param x Column component.
 * @param y Row component.
 * @param result Receives the normalized direction on success.
 * @return Whether the input is finite and nonzero. */
inline bool normalize_direction(float x, float y, Direction &result) {
  if (!std::isfinite(x) || !std::isfinite(y))
    return false;
  const float scale = std::max(std::abs(x), std::abs(y));
  if (scale == 0)
    return false;
  x /= scale;
  y /= scale;
  const float length = std::sqrt(x * x + y * y);
  result = {x / length, y / length};
  return true;
}

/** @brief Interpolates clockwise headings along the shortest arc.
 * @param a Initial degrees.
 * @param b Final degrees.
 * @param f Interpolation fraction.
 * @return Interpolated degrees modulo 360, with a positive half-turn tie. */
inline float interpolate_heading(float a, float b, float f) {
  a = std::remainder(a, 360.0f);
  float turn = std::remainder(std::remainder(b, 360.0f) - a, 360.0f);
  if (turn == -180)
    turn = 180;
  return std::fmod(a + turn * f + 360.0f, 360.0f);
}

/** @brief Projects an output normal and view-space arrow into the physical tangent metric.
 * @param p Unit view at the sample.
 * @param dx View derivative per output column.
 * @param dy View derivative per output row.
 * @param arrow Finite authored direction rotated into view space.
 * @param normal Outward pixel-space gradient.
 * @return Signed unit alignment, or zero at a degenerate sample. */
HS_HOT_FLASH_MEMBER inline float
tangent_alignment(math::Vector p, math::Vector dx, math::Vector dy,
                  math::Vector arrow, Direction normal) {
  const float scale =
      std::max({std::abs(arrow.x), std::abs(arrow.y), std::abs(arrow.z)});
  if (!(scale > 0) || !std::isfinite(scale))
    return 0;
  arrow = {arrow.x / scale, arrow.y / scale, arrow.z / scale};
  arrow = arrow / arrow.length();
  const auto tangent = arrow - p * math::dot(p, arrow);
  const float length = tangent.length();
  if (!(length > 1e-6f))
    return 0;

  const auto area = math::cross(dx, dy);
  const float a = math::dot(dx, dx), b = math::dot(dx, dy),
              c = math::dot(dy, dy);
  const float largest = 0.5f * (a + c + std::hypot(a - c, 2 * b));
  // |dx cross dy| / sigma_max^2 is sigma_min / sigma_max.
  if (!(area.length() > largest * 1e-6f))
    return 0;
  const auto physical =
      math::cross(dy, area) * normal.x + math::cross(area, dx) * normal.y;
  const float norm = physical.length();
  return norm > 0 ? std::clamp(math::dot(physical, tangent) / (length * norm),
                               -1.0f, 1.0f)
                  : 0;
}

/// Nearest exposed feature at a material sample.
struct Boundary {
  bool inside = false; ///< Membership at the query.
  float depth =
      std::numeric_limits<float>::infinity(); ///< Distance in output pixels.
  float alignment = 0;                        ///< Signed unit direction cosine.
};

/** @brief Classifies the opposing strips of a translated silhouette.
 * @param x,y Output sample coordinates.
 * @param shift Displacement of the silhouette in output pixels.
 * @param inside Coverage predicate in the same chart.
 * @return Coverage and strip polarity, with the primary strip painted last. */
inline Boundary translated_strips(float x, float y, Direction shift,
                                  const auto &inside) {
  if (!inside(x, y))
    return {};
  if (!inside(x + shift.x, y + shift.y))
    return {true, 0, 1};
  if (!inside(x - shift.x, y - shift.y))
    return {true, 0, -1};
  return {true};
}

/// Edge-neutral angular deadband.
inline constexpr float DEADBAND = 0.05f;

/** @brief Ranks edge polarity for equidistant features.
 * @param alignment Signed direction cosine.
 * @return Priority in ascending neutral, trailing, leading order. */
constexpr int priority(float alignment) {
  return alignment > DEADBAND ? 2 : alignment < -DEADBAND ? 1 : 0;
}

/** @brief Accumulates an exposed contour candidate, preserving membership.
 * @param nearest Current winning feature.
 * @param depth Candidate distance in output pixels.
 * @param alignment Candidate signed direction cosine. */
inline void consider(Boundary &nearest, float depth, float alignment) {
  if (depth < nearest.depth ||
      (depth == nearest.depth &&
       (priority(alignment) > priority(nearest.alignment) ||
        (priority(alignment) == priority(nearest.alignment) &&
         std::abs(alignment) > std::abs(nearest.alignment))))) {
    nearest.depth = depth;
    nearest.alignment = alignment;
  }
}

/** @brief Samples the shared two-arm presentation timeline.
 * @param t0 Frame start, elapsed microseconds.
 * @param t1 Frame end, elapsed microseconds, at or after the start.
 * @param x Wrapped logical output column.
 * @param width Positive even canvas width.
 * @param flat Whether all columns sample the frame end.
 * @return Elapsed microseconds at the column. */
constexpr uint64_t column_time(uint64_t t0, uint64_t t1, int x, int width,
                               bool flat) {
  if (flat)
    return t1;
  const uint64_t n = static_cast<unsigned>(width / 2);
  const uint64_t k = static_cast<unsigned>(x) % n;
  const uint64_t duration = t1 - t0;
  return t0 + duration / n * k + duration % n * k / n;
}

/// A beat's linear-light role colors.
struct Colors {
  Pixel fill;     ///< Interior and neutral edges.
  Pixel leading;  ///< Positive alignment.
  Pixel trailing; ///< Negative alignment.
};

/** @brief Decodes a packed sRGB color to linear RGB.
 * @param hex Packed 0xRRGGBB color.
 * @return Linear-light pixel. */
constexpr Pixel linear_rgb(uint32_t hex) {
  return {srgb_to_linear_lut[(hex >> 16) & 255],
          srgb_to_linear_lut[(hex >> 8) & 255], srgb_to_linear_lut[hex & 255]};
}

/// Opaque neutral reference backdrop.
inline constexpr Pixel GRAY = linear_rgb(0x808080);

/** @brief Looks up the linear-light roles for an exact presentation timestamp.
 * @param elapsed_us Microseconds since the shared effect epoch.
 * @return The active beat's palette. */
inline const Colors &colors(uint64_t elapsed_us) {
  static constexpr auto BEATS = [] {
    constexpr uint32_t HUES[] = {0x24b5b8, 0x5a78e6, 0xe06858, 0xd4a838};
    constexpr int CONTRAST_PERCENT = 45;
    constexpr int CHROMATIC_MIX = CONTRAST_PERCENT * 115;
    constexpr int NEUTRAL_MIX = CONTRAST_PERCENT * 72;
    const auto soften = [](uint32_t hex, int mix) {
      const auto channel = [=](int shift) -> uint32_t {
        const int value = (hex >> shift) & 255;
        return (128 * 10000 + (value - 128) * mix + 5000) / 10000;
      };
      return linear_rgb((channel(16) << 16) | (channel(8) << 8) | channel(0));
    };
    const Pixel dark = soften(0x222222, CHROMATIC_MIX);
    const Pixel light = soften(0xe6e6e6, CHROMATIC_MIX);
    const Pixel slate = soften(0x333b42, NEUTRAL_MIX);
    const Pixel mist = soften(0xd6d9dc, NEUTRAL_MIX);
    std::array<Colors, 2 * std::size(HUES)> beats{};
    for (size_t i = 0; i < std::size(HUES); i += 2) {
      const Pixel a = soften(HUES[i], CHROMATIC_MIX);
      const Pixel b = soften(HUES[i + 1], CHROMATIC_MIX);
      beats[2 * i] = {a, dark, light};
      beats[2 * i + 1] = {slate, a, b};
      beats[2 * i + 2] = {b, light, dark};
      beats[2 * i + 3] = {mist, b, a};
    }
    return beats;
  }();
  return BEATS[(elapsed_us / 62500) % BEATS.size()];
}

/** @brief Applies a four-stroke beat to straight material RGB.
 * @param input Incoming straight RGB.
 * @param boundary Membership and nearest outward feature.
 * @param beat Prepared role colors.
 * @param controls Authored strength, width and contrast law.
 * @return Shaded straight RGB. */
inline Pixel shade(Pixel input, Boundary boundary, const Colors &beat,
                   const Controls &controls) {
  if (!boundary.inside || controls.strength == 0)
    return input;
  Pixel role = beat.fill;
  if (boundary.depth < controls.width &&
      std::abs(boundary.alignment) > DEADBAND) {
    const Pixel edge = boundary.alignment > 0 ? beat.leading : beat.trailing;
    const float weight =
        (std::abs(boundary.alignment) - DEADBAND) / (1 - DEADBAND);
    role = controls.mode == Mode::REFERENCE_POLARITY
               ? edge
               : beat.fill.lerp16(edge, frac_to_q16(weight));
  }
  return input.lerp16(role, frac_to_q16(controls.strength));
}

/// Equal-area subpixel offsets on either output axis.
inline constexpr float OFFSETS[] = {-0.375f, -0.125f, 0.125f, 0.375f};
/// Subpixel taps resolve_pixel() averages.
inline constexpr int TAPS = int(std::size(OFFSETS) * std::size(OFFSETS));

/// Coverage-resolved linear color and alpha.
struct Resolved {
  Pixel premultiplied; ///< RGB after material alpha and coverage.
  float alpha;         ///< Mean material alpha over the footprint.
};

/** @brief Integrates material and coverage at TAPS subpixel centers.
 * @tparam Sample Straight-alpha material query.
 * @param sample Callable accepting center-relative column and row offsets.
 * @return Premultiplied RGB and alpha averaged over the footprint. */
template <typename Sample> inline Resolved resolve_pixel(Sample &&sample) {
  float r = 0, g = 0, b = 0, alpha = 0;
  for (float dy : OFFSETS) {
    for (float dx : OFFSETS) {
      const Color4 tap = sample(dx, dy);
      r += tap.color.r * tap.alpha;
      g += tap.color.g * tap.alpha;
      b += tap.color.b * tap.alpha;
      alpha += tap.alpha;
    }
  }
  return {{round_linear_channel(r / TAPS), round_linear_channel(g / TAPS),
           round_linear_channel(b / TAPS)},
          alpha / TAPS};
}

/** @brief Resolves premultiplied foreground over opaque linear RGB.
 * @param foreground Coverage-resolved material.
 * @param background Opaque backdrop.
 * @return Composited pixel. */
inline Pixel over(Resolved foreground, Pixel background) {
  const float remaining = 1 - foreground.alpha;
  return {round_linear_channel(foreground.premultiplied.r +
                               remaining * background.r),
          round_linear_channel(foreground.premultiplied.g +
                               remaining * background.g),
          round_linear_channel(foreground.premultiplied.b +
                               remaining * background.b)};
}

} // namespace FourStroke