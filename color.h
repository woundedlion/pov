/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <variant>
#include <vector>
#include <FastLED.h>
#include "3dmath.h"
#include "util.h"

 // TODO: 3D Palettes

 ///////////////////////////////////////////////////////////////////////////////

 /**
  * @brief Represents a color with an alpha channel.
  */
struct Color4 {
  Pixel color;
  float alpha;

  Color4() : color(CRGB::Black), alpha(1.0f) {}
  Color4(const Pixel& p, float a = 1.0f) : color(p), alpha(a) {}
  Color4(uint8_t r, uint8_t g, uint8_t b, float a = 1.0f) : color(r, g, b), alpha(a) {}
};

/**
* @brief Scales a Pixel's RGB components by a floating-point scalar.
*/
Pixel operator*(const Pixel& p, float s) {
  return Pixel(
    static_cast<uint8_t>(std::clamp(static_cast<int>(p.r * s), 0, 255)),
    static_cast<uint8_t>(std::clamp(static_cast<int>(p.g * s), 0, 255)),
    static_cast<uint8_t>(std::clamp(static_cast<int>(p.b * s), 0, 255))
  );
}

/**
 * @brief Scales a Pixel's RGB components by a floating-point scalar in place.
 */
Pixel& operator*=(Pixel& p, float s) {
  p = p * s;
  return p;
}

/**
 * @brief Gamma brightness lookup table (LUT).
 * gamma = 2.20 steps = 256 range = 0-255
 */
const uint8_t gamma_lut[256] = {
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,
     1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   2,   2,
     3,   3,   3,   3,   3,   4,   4,   4,   4,   5,   5,   5,   5,   6,   6,   6,
     6,   7,   7,   7,   8,   8,   8,   9,   9,   9,  10,  10,  11,  11,  11,  12,
    12,  13,  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18,  19,  19,
    20,  20,  21,  22,  22,  23,  23,  24,  25,  25,  26,  26,  27,  28,  28,  29,
    30,  30,  31,  32,  33,  33,  34,  35,  35,  36,  37,  38,  39,  39,  40,  41,
    42,  43,  43,  44,  45,  46,  47,  48,  49,  49,  50,  51,  52,  53,  54,  55,
    56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,
    73,  74,  75,  76,  77,  78,  79,  81,  82,  83,  84,  85,  87,  88,  89,  90,
    91,  93,  94,  95,  97,  98,  99, 100, 102, 103, 105, 106, 107, 109, 110, 111,
   113, 114, 116, 117, 119, 120, 121, 123, 124, 126, 127, 129, 130, 132, 133, 135,
   137, 138, 140, 141, 143, 145, 146, 148, 149, 151, 153, 154, 156, 158, 159, 161,
   163, 165, 166, 168, 170, 172, 173, 175, 177, 179, 181, 182, 184, 186, 188, 190,
   192, 194, 196, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221,
   223, 225, 227, 229, 231, 234, 236, 238, 240, 242, 244, 246, 248, 251, 253, 255,
};

/**
 * @brief Applies gamma correction to a CRGB pixel color.
 */
Pixel gamma_correct(const Pixel& p) {
  return CRGB(
    gamma_lut[p.r],
    gamma_lut[p.g],
    gamma_lut[p.b]
  );
}

///////////////////////////////////////////////////////////////////////////////
// Blending Functions
///////////////////////////////////////////////////////////////////////////////

Pixel blend_max(const Pixel& c1, const Pixel& c2) {
  return Pixel(std::max(c1.r, c2.r), std::max(c1.g, c2.g), std::max(c1.b, c2.b));
}

Pixel blend_over(const Pixel& c1, const Pixel& c2) {
  return c2;
}

Pixel blend_under(const Pixel& c1, const Pixel& c2) {
  return c1;
}

Pixel blend_add(const Pixel& c1, const Pixel& c2) {
  return Pixel(qadd8(c1.r, c2.r), qadd8(c1.g, c2.g), qadd8(c1.b, c2.b));
}

auto blend_alpha(float a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      c1.r * (1.0f - a) + c2.r * a,
      c1.g * (1.0f - a) + c2.g * a,
      c1.b * (1.0f - a) + c2.b * a);
    };
}

auto blend_accumulate(float a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      qadd8(c1.r, c2.r * a),
      qadd8(c1.g, c2.g * a),
      qadd8(c1.b, c2.b * a));
    };
}

Pixel blend_over_max(const Pixel& c1, const Pixel& c2) {
  float m1 = sqrtf(c1.r * c1.r + c1.g * c1.g + c1.b * c1.b);
  float m2 = sqrtf(c2.r * c2.r + c2.g * c2.g + c2.b * c2.b);
  if (m2 == 0) return c1;
  float s = std::max(m1, m2) / m2;
  return c2 * s;
}

Pixel blend_over_min(const Pixel& c1, const Pixel& c2) {
  float m1 = sqrtf(c1.r * c1.r + c1.g * c1.g + c1.b * c1.b);
  float m2 = sqrtf(c2.r * c2.r + c2.g * c2.g + c2.b * c2.b);
  if (m2 == 0) return Pixel(0, 0, 0);
  float s = std::min(m1, m2) / m2;
  return c2 * s;
}

Pixel blend_mean(const Pixel& c1, const Pixel& c2) {
  return Pixel((c1.r + c2.r) / 2, (c1.g + c2.g) / 2, (c1.b + c2.b) / 2);
}

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Base class for all color palette systems.
 */
class Palette {
public:
  virtual Color4 get(float t) const = 0;
  virtual ~Palette() = default;
};

/**
 * @brief Defines types of color harmonies for generative palettes.
 */
enum class HarmonyType {
  TRIADIC,
  SPLIT_COMPLEMENTARY,
  COMPLEMENTARY,
  ANALOGOUS
};

/**
 * @brief Defines the visual shape or distribution of colors across the palette domain.
 */
enum class GradientShape {
  STRAIGHT,
  CIRCULAR,
  VIGNETTE,
  FALLOFF
};

/**
 * @brief Defines the overall brightness profile across the palette domain.
 */
enum class BrightnessProfile {
  ASCENDING,
  DESCENDING,
  FLAT,
  BELL,
  CUP
};

/**
 * @brief Defines the saturation profile.
 */
enum class SaturationProfile {
  PASTEL,
  MID,
  VIBRANT
};

/**
 * @brief Converts a float in the range [0.0, 1.0] to a 16-bit integer for FastLED lerping.
 */
uint16_t to_short(float zero_to_one) {
  return std::clamp(static_cast<int>(std::round(zero_to_one * 65535.0f)), 0, 65535);
}

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief A class representing a discrete color gradient/lookup table.
 */
class Gradient : public Palette {
public:
  /**
   * @brief Constructs a gradient from a list of {position, color} pairs.
   * @param points Initializer list of pairs, e.g. {{0.0f, CRGB::Red}, {1.0f, CRGB::Blue}}.
   */
  Gradient(std::initializer_list<std::pair<float, Pixel>> points) {
    if (points.size() == 0) return;

    auto it = points.begin();
    float last_pos = it->first;
    Pixel last_color = it->second;
    ++it;

    // Fill initial segment if first point > 0
    if (last_pos > 0) {
      fill_solid(entries, to_short(last_pos) / 256, last_color);
    }

    for (; it != points.end(); ++it) {
      float next_pos = it->first;
      Pixel next_color = it->second;

      int start_idx = static_cast<int>(last_pos * 255);
      int end_idx = static_cast<int>(next_pos * 255);

      if (end_idx > start_idx) {
        fill_gradient_RGB(entries, start_idx, last_color, end_idx, next_color);
      }

      last_pos = next_pos;
      last_color = next_color;
    }

    // Fill remaining
    if (last_pos < 1.0f) {
      int start_idx = static_cast<int>(last_pos * 255);
      fill_solid(entries + start_idx, 256 - start_idx, last_color);
    }
  }

  Color4 get(float t) const override {
    return Color4(ColorFromPalette(entries, static_cast<uint8_t>(t * 255), 255, NOBLEND), 1.0f);
  }

  CRGBPalette256 entries;
};

/**
 * @brief A palette that generates colors based on defined harmony and brightness profiles.
 */
class GenerativePalette : public Palette {
public:

  GenerativePalette(
    GradientShape gradient_shape,
    HarmonyType harmony_type,
    BrightnessProfile profile,
    SaturationProfile sat_profile = SaturationProfile::MID) :
    gradient_shape(gradient_shape),
    harmony_type(harmony_type),
    seed_hue(static_cast<uint8_t>(hs::rand_int(0, 256)))
  {
    uint8_t h1 = seed_hue;
    uint8_t h2, h3;

    // Advance seed for next generation
    seed_hue = static_cast<uint8_t>(
      (static_cast<uint32_t>(seed_hue) + static_cast<uint32_t>(G * 255.0f)) % 256);

    calc_hues(h1, h2, h3, harmony_type);

    uint8_t s1, s2, s3;
    switch (sat_profile) {
    case SaturationProfile::PASTEL:
      s1 = s2 = s3 = 100;
      break;
    case SaturationProfile::MID:
      s1 = hs::rand_int(153, 204);
      s2 = hs::rand_int(153, 204);
      s3 = hs::rand_int(153, 204);
      break;
    case SaturationProfile::VIBRANT:
      s1 = s2 = s3 = 255;
      break;
    }

    uint8_t v1, v2, v3;
    switch (profile) {
    case BrightnessProfile::ASCENDING:
      v1 = hs::rand_int(25, 76);
      v2 = hs::rand_int(127, 178);
      v3 = hs::rand_int(204, 255);
      break;
    case BrightnessProfile::DESCENDING:
      v1 = hs::rand_int(204, 255);
      v2 = hs::rand_int(127, 178);
      v3 = hs::rand_int(25, 76);
      break;
    case BrightnessProfile::FLAT:
      v1 = v2 = v3 = 255;
      break;
    case BrightnessProfile::BELL:
      v1 = hs::rand_int(51, 127);
      v2 = hs::rand_int(178, 255);
      v3 = v1;
      break;
    case BrightnessProfile::CUP:
      v1 = hs::rand_int(178, 255);
      v2 = hs::rand_int(51, 127);
      v3 = v1;
      break;
    }

    a = CHSV(h1, s1, v1);
    b = CHSV(h2, s2, v2);
    c = CHSV(h3, s3, v3);

    update_luts();
  }

  void update_luts() {
    const Pixel vignette_color(0, 0, 0);
    switch (gradient_shape) {
    case GradientShape::VIGNETTE:
      shape = { 0, 0.1f, 0.5f, 0.9f, 1.0f };
      colors = { vignette_color, a, b, c, vignette_color };
      size = 5;
      break;
    case GradientShape::STRAIGHT:
      shape = { 0, 0.5f, 1.0f };
      colors = { a, b, c };
      size = 3;
      break;
    case GradientShape::CIRCULAR:
      shape = { 0, 0.33f, 0.66f, 1.0f };
      colors = { a, b, c, a };
      size = 4;
      break;
    case GradientShape::FALLOFF:
      shape = { 0, 0.33f, 0.66f, 0.9f, 1.0f };
      colors = { a, b, c, vignette_color };
      size = 4;
      break;
    }
  }

  void lerp(const GenerativePalette& from, const GenerativePalette& to, float amount) {
    uint16_t fract = to_short(amount);
    a = from.a.lerp16(to.a, fract);
    b = from.b.lerp16(to.b, fract);
    c = from.c.lerp16(to.c, fract);
    update_luts();
  }

  Color4 get(float t) const override {
    int seg = -1;
    for (int i = 0; i < size - 1; ++i) {
      if (t >= shape[i] && t < shape[i + 1]) {
        seg = i;
        break;
      }
    }
    if (seg < 0) seg = size - 2;

    float start = shape[seg];
    float end = shape[seg + 1];
    Pixel c1 = colors[seg];
    Pixel c2 = colors[seg + 1];

    // Safe division check
    float dist = end - start;
    if (dist < 0.0001f) return Color4(c1, 1.0f);

    float p = (t - start) / dist;
    return Color4(c1.lerp16(c2, to_short(std::clamp(p, 0.0f, 1.0f))), 1.0f);
  }

private:
  uint8_t wrap_hue(int hue) const { return (hue % 256 + 256) % 256; }

  void calc_hues(uint8_t h1, uint8_t& h2, uint8_t& h3, HarmonyType harmony_type) const {
    const int h1_int = h1;
    switch (harmony_type) {
    case HarmonyType::TRIADIC:
      h2 = wrap_hue(h1_int + 85);
      h3 = wrap_hue(h1_int + 170);
      break;
    case HarmonyType::SPLIT_COMPLEMENTARY: {
      const int complement = wrap_hue(h1_int + 128);
      const int offset = 21;
      h2 = wrap_hue(complement - offset);
      h3 = wrap_hue(complement + offset);
      break;
    }
    case HarmonyType::COMPLEMENTARY: {
      h2 = wrap_hue(h1_int + 128);
      const int offset = hs::rand_int(-7, 7);
      h3 = wrap_hue(h1_int + offset);
      break;
    }
    case HarmonyType::ANALOGOUS:
    default: {
      const int dir = (hs::rand_int(0, 1) == 0) ? 1 : -1;
      const int offset1 = dir * hs::rand_int(11, 22);
      h2 = wrap_hue(h1_int + offset1);
      const int offset2 = dir * hs::rand_int(11, 22);
      h3 = wrap_hue(h2 + offset2);
      break;
    }
    }
  }

  GradientShape gradient_shape;
  HarmonyType harmony_type;
  uint8_t seed_hue;
  Pixel a, b, c;
  std::array<float, 5> shape;
  std::array<Pixel, 5> colors;
  int size = 0;
};

/**
 * @brief A palette defined by a mathematical cosine wave function.
 * C(t) = A + B * cos(2 * PI * (C * t + D))
 */
class ProceduralPalette : public Palette {
public:
  ProceduralPalette(
    std::array<float, 3> a,
    std::array<float, 3> b,
    std::array<float, 3> c,
    std::array<float, 3> d) : a(a), b(b), c(c), d(d) {
  }

  Color4 get(float t) const override {
    return Color4(Pixel(
      static_cast<uint8_t>(255 * std::clamp(a[0] + b[0] * cosf(2 * PI_F * (c[0] * t + d[0])), 0.0f, 1.0f)),
      static_cast<uint8_t>(255 * std::clamp(a[1] + b[1] * cosf(2 * PI_F * (c[1] * t + d[1])), 0.0f, 1.0f)),
      static_cast<uint8_t>(255 * std::clamp(a[2] + b[2] * cosf(2 * PI_F * (c[2] * t + d[2])), 0.0f, 1.0f))
    ), 1.0f);
  }

protected:
  std::array<float, 3> a, b, c, d;
};

/**
 * @brief A palette that allows continuous mutation between two procedural palettes.
 */
class MutatingPalette : public ProceduralPalette {
public:
  MutatingPalette(
    std::array<float, 3> a1, std::array<float, 3> b1, std::array<float, 3> c1, std::array<float, 3> d1,
    std::array<float, 3> a2, std::array<float, 3> b2, std::array<float, 3> c2, std::array<float, 3> d2
  ) : ProceduralPalette(a1, b1, c1, d1), // Initialize base with start values
    a1(a1), b1(b1), c1(c1), d1(d1),
    a2(a2), b2(b2), c2(c2), d2(d2)
  {
    mutate(0.0f);
  }

  void mutate(float t) {
    for (int i = 0; i < 3; ++i) {
      a[i] = lerp(a1[i], a2[i], t);
      b[i] = lerp(b1[i], b2[i], t);
      c[i] = lerp(c1[i], c2[i], t);
      d[i] = lerp(d1[i], d2[i], t);
    }
  }

private:
  float lerp(float x, float y, float t) { return x * (1.0f - t) + y * t; }
  std::array<float, 3> a1, b1, c1, d1;
  std::array<float, 3> a2, b2, c2, d2;
};

///////////////////////////////////////////////////////////////////////////////
// Palette Wrappers
///////////////////////////////////////////////////////////////////////////////

class ReversePalette : public Palette {
public:
  ReversePalette(const Palette& palette) : palette(palette) {}
  Color4 get(float t) const override {
    return palette.get().get(1.0f - t);
  }
private:
  std::reference_wrapper<const Palette> palette;
};

class VignettePalette : public Palette {
public:
  VignettePalette(const Palette& palette) : palette(palette) {}
  Color4 get(float t) const override {
    CRGB vignette_color(0, 0, 0);
    if (t < 0.2f) {
      return Color4(vignette_color.lerp16(palette.get().get(0.0f).color, to_short(t / 0.2f)), 1.0f);
    }
    else if (t >= 0.8f) {
      return Color4(palette.get().get(1.0f).color.lerp16(vignette_color, to_short((t - 0.8f) / 0.2f)), 1.0f);
    }
    else {
      return palette.get().get((t - 0.2f) / 0.6f);
    }
  }
private:
  std::reference_wrapper<const Palette> palette;
};

class TransparentVignette : public Palette {
public:
  TransparentVignette(const Palette& palette) : palette(palette) {}
  Color4 get(float t) const override {
    Color4 result;
    float factor = 1.0f;
    if (t < 0.2f) {
      result = palette.get().get(0.0f);
      factor = t / 0.2f;
    }
    else if (t >= 0.8f) {
      result = palette.get().get(1.0f);
      factor = (1.0f - (t - 0.8f) / 0.2f);
    }
    else {
      return palette.get().get((t - 0.2f) / 0.6f);
    }
    result.alpha *= factor;
    return result;
  }
private:
  std::reference_wrapper<const Palette> palette;
};

// Backwards compatibility for functional style
auto vignette(const Palette& palette) {
  return VignettePalette(palette);
}

/**
 * @brief Variant holding any supported palette type.
 */
using PaletteVariant = std::variant<
  std::monostate,
  GenerativePalette,
  ProceduralPalette,
  Gradient,
  MutatingPalette,
  ReversePalette,
  VignettePalette,
  TransparentVignette
>;

///////////////////////////////////////////////////////////////////////////////
// Predefined Palettes
///////////////////////////////////////////////////////////////////////////////

// Procedural Palettes
static const ProceduralPalette darkRainbow({ 0.367f, 0.367f, 0.367f }, { 0.500f, 0.500f, 0.500f }, { 1.000f, 1.000f, 1.000f }, { 0.000f, 0.330f, 0.670f });
static const ProceduralPalette bloodStream({ 0.169f, 0.169f, 0.169f }, { 0.313f, 0.313f, 0.313f }, { 0.231f, 0.231f, 0.231f }, { 0.036f, 0.366f, 0.706f });
static const ProceduralPalette vintageSunset({ 0.256f, 0.256f, 0.256f }, { 0.500f, 0.080f, 0.500f }, { 0.277f, 0.277f, 0.277f }, { 0.000f, 0.330f, 0.670f });
static const ProceduralPalette richSunset({ 0.309f, 0.500f, 0.500f }, { 1.000f, 1.000f, 0.500f }, { 0.149f, 0.148f, 0.149f }, { 0.132f, 0.222f, 0.521f });
static const ProceduralPalette undersea({ 0.000f, 0.000f, 0.000f }, { 0.500f, 0.276f, 0.423f }, { 0.296f, 0.296f, 0.296f }, { 0.374f, 0.941f, 0.000f });
static const ProceduralPalette lateSunset({ 0.337f, 0.500f, 0.096f }, { 0.500f, 1.000f, 0.176f }, { 0.261f, 0.261f, 0.261f }, { 0.153f, 0.483f, 0.773f });
static const ProceduralPalette mangoPeel({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.080f, 0.500f }, { 0.431f, 0.431f, 0.431f }, { 0.566f, 0.896f, 0.236f });
static const ProceduralPalette iceMelt({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.500f, 0.500f }, { 0.083f, 0.147f, 0.082f }, { 0.579f, 0.353f, 0.244f });
static const ProceduralPalette lemonLime({ 0.455f, 0.455f, 0.455f }, { 0.571f, 0.151f, 0.571f }, { 0.320f, 0.320f, 0.320f }, { 0.087f, 0.979f, 0.319f });
static const ProceduralPalette algae({ 0.210f, 0.210f, 0.210f }, { 0.500f, 1.000f, 0.021f }, { 0.086f, 0.086f, 0.075f }, { 0.419f, 0.213f, 0.436f });
static const ProceduralPalette embers({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.500f, 0.500f }, { 0.265f, 0.285f, 0.198f }, { 0.577f, 0.440f, 0.358f });
static const ProceduralPalette fireGlow({ 0.000f, 0.000f, 0.000f }, { 0.560f, 0.560f, 0.560f }, { 0.216f, 0.346f, 0.174f }, { 0.756f, 0.542f, 0.279f });
static const ProceduralPalette darkPrimary({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.610f, 0.500f }, { 0.746f, 0.347f, 0.000f }, { 0.187f, 0.417f, 0.670f });
static const ProceduralPalette mauveFade({ 0.583f, 0.000f, 0.583f }, { 1.000f, 0.000f, 1.000f }, { 0.191f, 0.348f, 0.191f }, { 0.175f, 0.045f, 0.150f });
static const ProceduralPalette lavenderLake({ 0.288f, 0.080f, 0.288f }, { 0.380f, 0.430f, 0.380f }, { 0.154f, 0.154f, 0.154f }, { 2.926f, 2.926f, 2.926f });

// Gradient Palettes
static const Gradient rainbow({
  {0.0f / 16, CRGB(0xFF0000)},
  {1.0f / 16, CRGB(0xD52A00)},
  {2.0f / 16, CRGB(0xAB5500)},
  {3.0f / 16, CRGB(0xAB7F00)},
  {4.0f / 16, CRGB(0xABAB00)},
  {5.0f / 16, CRGB(0x56D500)},
  {6.0f / 16, CRGB(0x00FF00)},
  {7.0f / 16, CRGB(0x00D52A)},
  {8.0f / 16, CRGB(0x00AB55)},
  {9.0f / 16, CRGB(0x0056AA)},
  {10.0f / 16, CRGB(0x0000FF)},
  {11.0f / 16, CRGB(0x2A00D5)},
  {12.0f / 16, CRGB(0x5500AB)},
  {13.0f / 16, CRGB(0x7F0081)},
  {14.0f / 16, CRGB(0xAB0055)},
  {15.0f / 16, CRGB(0xD5002B)},
  {16.0f / 16, CRGB(0xD5002B)}
  });

static const Gradient rainbowStripes({
  {0.0f / 16, CRGB(0xFF0000)},
  {1.0f / 16, CRGB(0x000000)},
  {2.0f / 16, CRGB(0xAB5500)},
  {3.0f / 16, CRGB(0x000000)},
  {4.0f / 16, CRGB(0xABAB00)},
  {5.0f / 16, CRGB(0x000000)},
  {6.0f / 16, CRGB(0x00FF00)},
  {7.0f / 16, CRGB(0x000000)},
  {8.0f / 16, CRGB(0x00AB55)},
  {9.0f / 16, CRGB(0x000000)},
  {10.0f / 16, CRGB(0x0000FF)},
  {11.0f / 16, CRGB(0x000000)},
  {12.0f / 16, CRGB(0x5500AB)},
  {13.0f / 16, CRGB(0x000000)},
  {14.0f / 16, CRGB(0xAB0055)},
  {15.0f / 16, CRGB(0x000000)},
  {16.0f / 16, CRGB(0xFF0000)}
  });

static const Gradient rainbowThinStripes({
  {0.0f, CRGB(0xFF0000)},
  {1.0f / 32, CRGB(0x000000)},
  {3.0f / 32, CRGB(0x000000)},
  {4.0f / 32, CRGB(0xAB5500)},
  {5.0f / 32, CRGB(0x000000)},
  {7.0f / 32, CRGB(0x000000)},
  {8.0f / 32, CRGB(0xABAB00)},
  {9.0f / 32, CRGB(0x000000)},
  {11.0f / 32, CRGB(0x000000)},
  {12.0f / 32, CRGB(0x00FF00)},
  {13.0f / 32, CRGB(0x000000)},
  {15.0f / 32, CRGB(0x000000)},
  {16.0f / 32, CRGB(0x00AB55)},
  {17.0f / 32, CRGB(0x000000)},
  {19.0f / 32, CRGB(0x000000)},
  {20.0f / 32, CRGB(0x0000FF)},
  {21.0f / 32, CRGB(0x000000)},
  {23.0f / 32, CRGB(0x000000)},
  {24.0f / 32, CRGB(0x5500AB)},
  {25.0f / 32, CRGB(0x000000)},
  {27.0f / 32, CRGB(0x000000)},
  {28.0f / 32, CRGB(0xAB0055)},
  {29.0f / 32, CRGB(0x000000)},
  {32.0f / 32, CRGB(0x000000)}
  });

static const Gradient grayToBlack({
  {0.0f, CRGB(0x888888)},
  {1.0f, CRGB(0x000000)}
  });

static const Gradient blueToBlack({
  {0.0f, CRGB(0xee00ee)},
  {1.0f, CRGB(0x000000)}
  });

static const Gradient emeraldForest({
  {0.0f, CRGB(0x004E64)},
  {0.2f, CRGB(0x0B6E4F)},
  {0.4f, CRGB(0x08A045)},
  {0.6f, CRGB(0x6BBF59)},
  {0.8f, CRGB(0x138086)},
  {1.0f, CRGB(0x000000)}
  });

static const Gradient g1({
  {0.0f, CRGB(0xffaa00)},
  {1.0f, CRGB(0xff0000)}
  });

static const Gradient g2({
  {0.0f, CRGB(0x0000ff)},
  {1.0f, CRGB(0x660099)}
  });

static const Gradient g3({
  {0.0f, CRGB(0xffff00)},
  {0.3f, CRGB(0xfc7200)},
  {0.8f, CRGB(0x06042f)},
  {1.0f, CRGB(0x000000)}
  });

static const Gradient g4({
  {0.0f, CRGB(0x0000ff)},
  {1.0f, CRGB(0x000000)}
  });