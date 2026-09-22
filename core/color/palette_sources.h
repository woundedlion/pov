/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette_sources.h
 * @brief Gradient and procedural palette implementations.
 */

#include <array>
#include <initializer_list>
#include "color/palette.h"
#include "color/color_space.h"

/**
 * @brief A palette backed by a precomputed 256-entry linear-RGB lookup table,
 * filled at construction by interpolating between the color stops in OKLCH
 * perceptual space.
 */
class Gradient : public Palette {
public:
  /**
   * @brief Builds the 256-entry LUT by interpolating between color stops.
   * @param points Sorted-ascending (position in [0,1], color) stops.
   * @details Emptiness, stop bounds and ordering are trapped always-on
   * (construction is cold).
   */
  HS_FLASH_MEMBER
  Gradient(std::initializer_list<std::pair<float, CPixel>> points) : entries() {
    HS_CHECK(points.size() > 0, "Gradient requires at least one stop");

    float prev_check = -1.0f;
    for (const auto &stop : points) {
      HS_CHECK(stop.first >= 0.0f && stop.first <= 1.0f,
               "Gradient stop position out of [0,1]");
      HS_CHECK(stop.first >= prev_check,
               "Gradient stops must be sorted ascending");
      prev_check = stop.first;
    }

    auto it = points.begin();
    float prev_pos = it->first;
    CPixel prev_color = it->second;

    // Flat fills bake through the same OKLCH path as the segment endpoints, so
    // a flat region matches its adjacent segment at the shared stop.
    int first_stop = static_cast<int>(prev_pos * 255.0f + 0.5f);
    Pixel prev_solid =
        oklch_to_pixel(srgb_to_oklch(prev_color.r, prev_color.g, prev_color.b));
    for (int i = 0; i <= first_stop; i++)
      entries[i] = prev_solid;

    it++;
    while (it != points.end()) {
      float next_pos = it->first;
      CPixel next_color = it->second;

      int start = static_cast<int>(prev_pos * 255.0f + 0.5f);
      int end = static_cast<int>(next_pos * 255.0f + 0.5f);

      // end == start (two stops quantizing to the same index) is the intended
      // "hard stop" — an abrupt color boundary, not a dropped stop.
      if (end > start) {
        OKLCH a = srgb_to_oklch(prev_color.r, prev_color.g, prev_color.b);
        OKLCH b = srgb_to_oklch(next_color.r, next_color.g, next_color.b);
        for (int i = start; i < end; i++) {
          float t = static_cast<float>(i - start) / (end - start);
          entries[i] = oklch_to_pixel(lerp_oklch(a, b, t));
        }
      }
      prev_pos = next_pos;
      prev_color = next_color;
      it++;
    }

    int last_stop = static_cast<int>(prev_pos * 255.0f + 0.5f);
    Pixel last_solid =
        oklch_to_pixel(srgb_to_oklch(prev_color.r, prev_color.g, prev_color.b));
    for (int i = last_stop; i < 256; i++)
      entries[i] = last_solid;
  }

  /**
   * @brief LUT lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1].
   * @return The interpolated color (alpha 1.0).
   * @details Interpolated, not nearest-index, to avoid visible banding.
   */
  Color4 get(float t) const override {
    // Clamp before the int cast: t < 0 is float->int UB and NaN maps to the last
    // entry, both of which lut_sample_pixel requires the caller to have excluded.
    return Color4(
        lut_sample_pixel(entries, 256, hs::clamp(t, 0.0f, 1.0f) * 255.0f),
        1.0f);
  }

private:
  Pixel entries[256];
};

/**
 * @brief A palette defined by a mathematical cosine wave function.
 * C(t) = A + B * cos(2 * PI * (C * t + D))
 */
class ProceduralPalette : public Palette {
public:
  /**
   * @brief Default-constructs a palette with all-zero cosine coefficients.
   */
  constexpr ProceduralPalette()
      : a{0, 0, 0}, b{0, 0, 0}, c{0, 0, 0}, d{0, 0, 0} {}
  /**
   * @brief Constructs from the four cosine coefficient vectors.
   * @param a Bias term per RGB channel.
   * @param b Amplitude per RGB channel.
   * @param c Frequency per RGB channel.
   * @param d Phase per RGB channel.
   */
  constexpr ProceduralPalette(std::array<float, 3> a, std::array<float, 3> b,
                              std::array<float, 3> c, std::array<float, 3> d)
      : a(a), b(b), c(c), d(d) {}

  /**
   * @brief Evaluates the cosine palette at a coordinate.
   * @param t Lookup coordinate.
   * @return The color at t (alpha 1.0).
   * @details Computes color in float sRGB space, then converts to 16-bit linear
   * via the interpolated LUT, avoiding 8-bit quantization without a per-channel
   * powf.
   */
  Color4 get(float t) const override {
    float r_srgb =
        a[0] + b[0] * math::fast_cosf(2 * math::PI_F * (c[0] * t + d[0]));
    float g_srgb =
        a[1] + b[1] * math::fast_cosf(2 * math::PI_F * (c[1] * t + d[1]));
    float b_srgb =
        a[2] + b[2] * math::fast_cosf(2 * math::PI_F * (c[2] * t + d[2]));

    Pixel color(srgb_to_linear_interp(r_srgb), srgb_to_linear_interp(g_srgb),
                srgb_to_linear_interp(b_srgb));
    return Color4(color, 1.0f);
  }

  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~ProceduralPalette() {}

protected:
  std::array<float, 3> a, b, c, d;
};

/**
 * @brief A palette that allows continuous mutation between two procedural
 * palettes.
 */
class MutatingPalette : public ProceduralPalette {
public:
  /**
   * @brief Constructs from two endpoint cosine parameter sets.
   * @param a1 Start bias per channel.
   * @param b1 Start amplitude per channel.
   * @param c1 Start frequency per channel.
   * @param d1 Start phase per channel.
   * @param a2 End bias per channel.
   * @param b2 End amplitude per channel.
   * @param c2 End frequency per channel.
   * @param d2 End phase per channel.
   * @details Initializes the active parameters to the start set (mutate(0)).
   */
  MutatingPalette(std::array<float, 3> a1, std::array<float, 3> b1,
                  std::array<float, 3> c1, std::array<float, 3> d1,
                  std::array<float, 3> a2, std::array<float, 3> b2,
                  std::array<float, 3> c2, std::array<float, 3> d2)
      : ProceduralPalette(a1, b1, c1, d1), a1(a1), b1(b1), c1(c1), d1(d1),
        a2(a2), b2(b2), c2(c2), d2(d2) {
    mutate(0.0f);
  }

  /**
   * @brief Sets the active cosine parameters to the endpoint interpolation.
   * @param t Blend weight in [0, 1] between the start and end parameter sets.
   */
  void mutate(float t) {
    for (int i = 0; i < 3; ++i) {
      a[i] = hs::lerp(a1[i], a2[i], t);
      b[i] = hs::lerp(b1[i], b2[i], t);
      c[i] = hs::lerp(c1[i], c2[i], t);
      d[i] = hs::lerp(d1[i], d2[i], t);
    }
  }

private:
  std::array<float, 3> a1, b1, c1, d1;
  std::array<float, 3> a2, b2, c2, d2;

public:
  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~MutatingPalette() {}
};
