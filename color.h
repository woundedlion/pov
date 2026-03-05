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

#include "platform.h"
#include "3dmath.h"
#include "util.h"
#include "static_circular_buffer.h"

inline uint16_t srgb_to_linear(uint8_t srgb);

/**
 * @brief Represents a 16-bit per channel RGB color (Linear space).
 * @details Used for high-precision mixing and HDR rendering before
 * downsampling/tone-mapping to 8-bit output.
 */
struct Pixel16 {
  uint16_t r, g, b;

  constexpr Pixel16() : r(0), g(0), b(0) {}
  constexpr Pixel16(uint16_t _r, uint16_t _g, uint16_t _b)
      : r(_r), g(_g), b(_b) {}

  // Construct from HSV (converts to sRGB then Linear)
  Pixel16(const CHSV &hsv) {
    CRGB srgb(hsv);
    r = srgb_to_linear(srgb.r);
    g = srgb_to_linear(srgb.g);
    b = srgb_to_linear(srgb.b);
  }

  // Construct from CRGB (converts to Linear)
  Pixel16(const CRGB &c) {
    r = srgb_to_linear(c.r);
    g = srgb_to_linear(c.g);
    b = srgb_to_linear(c.b);
  }

  // Explicit forward declaration of conversion
  operator CRGB() const;

  // Basic arithmetic
  Pixel16 &operator+=(const Pixel16 &rhs) {
    r = (uint16_t)std::min((uint32_t)65535, (uint32_t)r + rhs.r);
    g = (uint16_t)std::min((uint32_t)65535, (uint32_t)g + rhs.g);
    b = (uint16_t)std::min((uint32_t)65535, (uint32_t)b + rhs.b);
    return *this;
  }

  Pixel16 operator*(float s) const {
    return Pixel16((uint16_t)std::clamp((int)(r * s), 0, 65535),
                   (uint16_t)std::clamp((int)(g * s), 0, 65535),
                   (uint16_t)std::clamp((int)(b * s), 0, 65535));
  }

  // Lerp (Linear interpolation 16-bit)
  Pixel16 lerp16(const Pixel16 &other, uint16_t frac) const {
    // frac 0..65535
    // result = this + (other - this) * frac
    // Exact: (a * (65535 - frac) + b * frac) / 65535
    uint32_t r32 =
        ((uint32_t)r * (65535 - frac) + (uint32_t)other.r * frac) / 65535;
    uint32_t g32 =
        ((uint32_t)g * (65535 - frac) + (uint32_t)other.g * frac) / 65535;
    uint32_t b32 =
        ((uint32_t)b * (65535 - frac) + (uint32_t)other.b * frac) / 65535;
    return Pixel16((uint16_t)r32, (uint16_t)g32, (uint16_t)b32);
  }

  bool operator==(const Pixel16 &rhs) const {
    return r == rhs.r && g == rhs.g && b == rhs.b;
  }

  bool operator!=(const Pixel16 &rhs) const { return !(*this == rhs); }

  bool operator==(const CHSV &rhs) const { return *this == Pixel16(rhs); }

  bool operator!=(const CHSV &rhs) const { return !(*this == rhs); }

  bool operator==(const CRGB &rhs) const { return *this == Pixel16(rhs); }

  bool operator!=(const CRGB &rhs) const { return !(*this == rhs); }
};

using Pixel = Pixel16;

/**
 * @brief Represents a color with an alpha channel.
 */
struct Color4 {
  Pixel color;
  float alpha;

  Color4() : color(Pixel(0, 0, 0)), alpha(1.0f) {}
  Color4(Pixel p, float a = 1.0f) : color(p), alpha(a) {}
  Color4(uint8_t r, uint8_t g, uint8_t b, float a = 1.0f)
      : color(Pixel(srgb_to_linear(r), srgb_to_linear(g), srgb_to_linear(b))),
        alpha(a) {}
  Color4(Color4 c, float a) : color(c.color), alpha(a) {}

  Color4 lerp(const Color4 &other, float t) const {
    uint16_t frac = static_cast<uint16_t>(hs::clamp(t, 0.0f, 1.0f) * 65535.0f);
    Pixel blended = color.lerp16(other.color, frac);
    float blended_a = alpha + (other.alpha - alpha) * t;
    return Color4(blended, blended_a);
  }

  operator CRGB() const { return color; }
};

#include "color_luts.h"

inline uint16_t srgb_to_linear(uint8_t srgb) {
  return srgb_to_linear_lut[srgb];
}

inline Pixel16::operator CRGB() const {
  return CRGB(linear_to_srgb_lut[r], linear_to_srgb_lut[g],
              linear_to_srgb_lut[b]);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Blending Functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * @brief Blends two pixels by taking the maximum of each component.
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return The max blend result.
 */
inline Pixel blend_max(const Pixel &c1, const Pixel &c2) {
  return Pixel(std::max(c1.r, c2.r), std::max(c1.g, c2.g),
               std::max(c1.b, c2.b));
}

/**
 * @brief Blends two pixels by overwriting (Painter's algorithm).
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return c2 (Source).
 */
inline Pixel blend_over(const Pixel &c1, const Pixel &c2) { return c2; }

/**
 * @brief Blends two pixels by keeping the destination (Under).
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return c1 (Destination).
 */
inline Pixel blend_under(const Pixel &c1, const Pixel &c2) { return c1; }

/**
 * @brief Blends two pixels by additive mixing (saturated).
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return c1 + c2 (clamped).
 */
inline Pixel blend_add(const Pixel &c1, const Pixel &c2) {
  // Saturated Add
  uint32_t r = (uint32_t)c1.r + c2.r;
  uint32_t g = (uint32_t)c1.g + c2.g;
  uint32_t b = (uint32_t)c1.b + c2.b;
  return Pixel((r > 65535) ? 65535 : (uint16_t)r,
               (g > 65535) ? 65535 : (uint16_t)g,
               (b > 65535) ? 65535 : (uint16_t)b);
}

inline auto blend_alpha(float a) {
  uint16_t ai = std::clamp((int)(a * 65535), 0, 65535);
  return [ai](const Pixel &c1, const Pixel &c2) { return c1.lerp16(c2, ai); };
}

inline auto blend_accumulate(float a) {
  return [a](const Pixel &c1, const Pixel &c2) {
    // c1 + c2 * a
    Pixel16 added = c2 * a;
    return blend_add(c1, added);
  };
}

inline Pixel blend_over_max(const Pixel &c1, const Pixel &c2) {
  float m1 = (float)c1.r * c1.r + (float)c1.g * c1.g + (float)c1.b * c1.b;
  float m2 = (float)c2.r * c2.r + (float)c2.g * c2.g + (float)c2.b * c2.b;
  if (m2 < 0.001f)
    return c1;

  if (m1 > m2) {
    float s = sqrtf(m1 / m2);
    return c2 * s;
  }
  return c2;
}

inline Pixel blend_over_min(const Pixel &c1, const Pixel &c2) {
  float m1 = (float)c1.r * c1.r + (float)c1.g * c1.g + (float)c1.b * c1.b;
  float m2 = (float)c2.r * c2.r + (float)c2.g * c2.g + (float)c2.b * c2.b;
  if (m2 < 0.001f)
    return Pixel(0, 0, 0);

  if (m1 < m2) {
    float s = sqrtf(m1 / m2);
    return c2 * s;
  }
  return c2;
}

/**
 * @brief Blends two pixels by averaging/mean.
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return (c1 + c2) / 2.
 */
inline Pixel blend_mean(const Pixel &c1, const Pixel &c2) {
  return Pixel((c1.r + c2.r) / 2, (c1.g + c2.g) / 2, (c1.b + c2.b) / 2);
}

/**
 * @brief Enumerates available blending modes for renderers.
 */
enum BlendMode : uint8_t { BLEND_OVER = 0, BLEND_ADD = 1, BLEND_MAX = 2 };

inline auto blend_add_alpha(float a) {
  return [a](const Pixel &dest, const Pixel &src) {
    return blend_add(dest, src * a);
  };
}

inline auto blend_max_alpha(float a) {
  return [a](const Pixel &dest, const Pixel &src) {
    Pixel16 s = src * a;
    return Pixel(std::max(dest.r, s.r), std::max(dest.g, s.g),
                 std::max(dest.b, s.b));
  };
}

///////////////////////////////////////////////////////////////////////////////

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
 * @brief Defines the visual shape or distribution of colors across the palette
 * domain.
 */
enum class GradientShape { STRAIGHT, CIRCULAR, VIGNETTE, FALLOFF };

/**
 * @brief Defines the overall brightness profile across the palette domain.
 */
enum class BrightnessProfile { ASCENDING, DESCENDING, FLAT, BELL, CUP };

/**
 * @brief Defines the saturation profile.
 */
enum class SaturationProfile { PASTEL, MID, VIBRANT };

/**
 * @brief Converts a float in the range [0.0, 1.0] to a 16-bit integer for
 * FastLED lerping.
 */
inline uint16_t to_short(float zero_to_one) {
  return std::clamp(static_cast<int>(std::round(zero_to_one * 65535.0f)), 0,
                    65535);
}

///////////////////////////////////////////////////////////////////////////////

constexpr uint8_t lerp8(uint8_t a, uint8_t b, float t) {
  return static_cast<uint8_t>(a + (b - a) * t);
}

/**
 * @brief A constexpr-compatible RGB pixel structure for Flash storage.
 * Layout compatible with CRGB but without non-constexpr constructors.
 */
struct CPixel {
  uint8_t r, g, b;
  constexpr CPixel() : r(0), g(0), b(0) {}
  constexpr CPixel(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
  constexpr CPixel(uint32_t hex)
      : r((hex >> 16) & 0xFF), g((hex >> 8) & 0xFF), b(hex & 0xFF) {}
  CPixel(const CRGB &c) : r(c.r), g(c.g), b(c.b) {}

  // Convert to FastLED CRGB (Pixel)
  operator Pixel() const { return CRGB(r, g, b); }
};

/**
 * @brief A class representing a discrete color gradient/lookup table.
 */
// Helper for high-precision conversion (sRGB float 0-1 -> Linear float 0-1)
constexpr float srgb_to_linear_float(float s) {
  return (s <= 0.04045f) ? s / 12.92f : powf((s + 0.055f) / 1.055f, 2.4f);
}

class Gradient {
public:
  Pixel entries[256];

  Gradient(std::initializer_list<std::pair<float, CPixel>> points) : entries() {
    Pixel black(0, 0, 0);
    for (int i = 0; i < 256; i++)
      entries[i] = black;

    if (points.size() == 0)
      return;

    auto it = points.begin();
    float prevPos = it->first;
    CPixel prevColor = it->second;

    // Fill start
    int firstStop = static_cast<int>(prevPos * 255);
    Pixel prevLinear(srgb_to_linear(prevColor.r), srgb_to_linear(prevColor.g),
                     srgb_to_linear(prevColor.b));
    for (int i = 0; i <= firstStop; i++)
      entries[i] = prevLinear;

    it++;
    while (it != points.end()) {
      float nextPos = it->first;
      CPixel nextColor = it->second;

      int start = static_cast<int>(prevPos * 255);
      int end = static_cast<int>(nextPos * 255);

      if (end > start) {
        for (int i = start; i <= end; i++) {
          float t = static_cast<float>(i - start) / (end - start);

          // Interpolate in sRGB float space
          float r = (float)prevColor.r * (1.0f - t) + (float)nextColor.r * t;
          float g = (float)prevColor.g * (1.0f - t) + (float)nextColor.g * t;
          float b = (float)prevColor.b * (1.0f - t) + (float)nextColor.b * t;

          entries[i] =
              Pixel(srgb_to_linear(r), srgb_to_linear(g), srgb_to_linear(b));
        }
      }
      prevPos = nextPos;
      prevColor = nextColor;
      it++;
    }

    // Fill end
    int lastStop = static_cast<int>(prevPos * 255);
    Pixel lastLinear(srgb_to_linear(prevColor.r), srgb_to_linear(prevColor.g),
                     srgb_to_linear(prevColor.b));
    for (int i = lastStop; i < 256; i++)
      entries[i] = lastLinear;
  }

  Color4 get(float t) const {
    uint8_t index = static_cast<uint8_t>(t * 255);
    return Color4(entries[index], 1.0f);
  }

  ~Gradient() {}
};

/**
 * @brief A palette that generates colors based on defined harmony and
 * brightness profiles.
 */
class GenerativePalette {
public:
  GenerativePalette()
      : GenerativePalette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS,
                          BrightnessProfile::FLAT) {}

  GenerativePalette(GradientShape gradient_shape, HarmonyType harmony_type,
                    BrightnessProfile profile,
                    SaturationProfile sat_profile = SaturationProfile::MID,
                    int manual_seed = -1)
      : gradient_shape(gradient_shape), harmony_type(harmony_type) {
    if (manual_seed != -1) {
      palette_hue = static_cast<uint8_t>(manual_seed);
    } else {
      palette_hue = g_hue_seed;
      // Advance global seed for next generation (Golden Ratio distribution)
      g_hue_seed = static_cast<uint8_t>((static_cast<uint32_t>(g_hue_seed) +
                                         static_cast<uint32_t>(G * 255.0f)) %
                                        256);
    }

    uint8_t h1 = palette_hue;
    uint8_t h2, h3;

    calc_hues(h1, h2, h3, harmony_type);

    uint8_t s1 = 0, s2 = 0, s3 = 0;
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

    uint8_t v1 = 0, v2 = 0, v3 = 0;
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

    a = CPixel(CRGB(CHSV(h1, s1, v1)));
    b = CPixel(CRGB(CHSV(h2, s2, v2)));
    c = CPixel(CRGB(CHSV(h3, s3, v3)));

    update_luts();
  }

  void update_luts() {
    const CPixel vignette_color(0, 0, 0);
    switch (gradient_shape) {
    case GradientShape::VIGNETTE:
      shape = {0, 0.1f, 0.5f, 0.9f, 1.0f};
      colors = {vignette_color, a, b, c, vignette_color};
      size = 5;
      break;
    case GradientShape::STRAIGHT:
      shape = {0, 0.5f, 1.0f};
      colors = {a, b, c};
      size = 3;
      break;
    case GradientShape::CIRCULAR:
      shape = {0, 0.33f, 0.66f, 1.0f};
      colors = {a, b, c, a};
      size = 4;
      break;
    case GradientShape::FALLOFF:
      shape = {0, 0.33f, 0.66f, 0.9f, 1.0f};
      colors = {a, b, c, vignette_color, vignette_color};
      size = 5;
      break;
    }
  }

  void lerp(const GenerativePalette &from, const GenerativePalette &to,
            float amount) {
    a = CPixel(lerp8(from.a.r, to.a.r, amount), lerp8(from.a.g, to.a.g, amount),
               lerp8(from.a.b, to.a.b, amount));
    b = CPixel(lerp8(from.b.r, to.b.r, amount), lerp8(from.b.g, to.b.g, amount),
               lerp8(from.b.b, to.b.b, amount));
    c = CPixel(lerp8(from.c.r, to.c.r, amount), lerp8(from.c.g, to.c.g, amount),
               lerp8(from.c.b, to.c.b, amount));
    update_luts();
  }

  Color4 get(float t) const {
    int seg = -1;
    for (int i = 0; i < size - 1; ++i) {
      if (t >= shape[i] && t < shape[i + 1]) {
        seg = i;
        break;
      }
    }
    if (seg < 0)
      seg = size - 2;

    float start = shape[seg];
    float end = shape[seg + 1];
    CPixel c1 = colors[seg];
    CPixel c2 = colors[seg + 1];

    float dist = end - start;
    if (dist < 0.0001f)
      return Color4(c1, 1.0f);

    float p = std::clamp((t - start) / dist, 0.0f, 1.0f);

    uint8_t r = lerp8(c1.r, c2.r, p);
    uint8_t g = lerp8(c1.g, c2.g, p);
    uint8_t b = lerp8(c1.b, c2.b, p);

    Pixel color(srgb_to_linear(r), srgb_to_linear(g), srgb_to_linear(b));
    return Color4(color, 1.0f);
  }

private:
  uint8_t wrap_hue(int hue) const { return (hue % 256 + 256) % 256; }

  void calc_hues(uint8_t h1, uint8_t &h2, uint8_t &h3,
                 HarmonyType harmony_type) const {
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
  uint8_t palette_hue;
  CPixel a, b, c;

  static inline uint8_t g_hue_seed = 0;
  std::array<float, 5> shape;
  std::array<CPixel, 5> colors;
  int size = 0;

public:
  constexpr ~GenerativePalette() {}
};

/**
 * @brief A palette defined by a mathematical cosine wave function.
 * C(t) = A + B * cos(2 * PI * (C * t + D))
 */
class ProceduralPalette {
public:
  constexpr ProceduralPalette(std::array<float, 3> a, std::array<float, 3> b,
                              std::array<float, 3> c, std::array<float, 3> d)
      : a(a), b(b), c(c), d(d) {}

  Color4 get(float t) const {
    // Determine color in high-precision float sRGB space first, then convert to
    // 16-bit Linear This avoids 8-bit quantization steps
    float r_srgb =
        hs::clamp(a[0] + b[0] * cosf(2 * PI_F * (c[0] * t + d[0])), 0.0f, 1.0f);
    float g_srgb =
        hs::clamp(a[1] + b[1] * cosf(2 * PI_F * (c[1] * t + d[1])), 0.0f, 1.0f);
    float b_srgb =
        hs::clamp(a[2] + b[2] * cosf(2 * PI_F * (c[2] * t + d[2])), 0.0f, 1.0f);

    Pixel color(srgb_to_linear(r_srgb * 255.0f),
                srgb_to_linear(g_srgb * 255.0f),
                srgb_to_linear(b_srgb * 255.0f));
    return Color4(color, 1.0f);
  }

protected:
  std::array<float, 3> a, b, c, d;

public:
  constexpr ~ProceduralPalette() {}
};

/**
 * @brief A palette that allows continuous mutation between two procedural
 * palettes.
 */
class MutatingPalette : public ProceduralPalette {
public:
  MutatingPalette(std::array<float, 3> a1, std::array<float, 3> b1,
                  std::array<float, 3> c1, std::array<float, 3> d1,
                  std::array<float, 3> a2, std::array<float, 3> b2,
                  std::array<float, 3> c2, std::array<float, 3> d2)
      : ProceduralPalette(a1, b1, c1, d1), // Initialize base with start values
        a1(a1), b1(b1), c1(c1), d1(d1), a2(a2), b2(b2), c2(c2), d2(d2) {
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

public:
  constexpr ~MutatingPalette() {}
};

///////////////////////////////////////////////////////////////////////////////
// Palette Modifiers
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Linearly cycles the palette coordinate.
 */
struct CycleModifier {
  const float *offset;

  CycleModifier(const float *driver_offset = nullptr) : offset(driver_offset) {}

  float modify(float t) const { return offset ? t + *offset : t; }
};

/**
 * @brief Oscillates the palette coordinate (Breathing).
 */
struct BreatheModifier {
  const float *phase;
  float amplitude;

  BreatheModifier(const float *driver_phase, float amp = 0.1f)
      : phase(driver_phase), amplitude(amp) {}

  float modify(float t) const {
    if (!phase)
      return t;
    return t + sinf(*phase) * amplitude;
  }
};

/**
 * @brief Distorts the palette spatially with a sine wave, creating a liquid
 * ripple effect. When applied to a spatial coordinate, colors will compress and
 * expand like waves.
 */
struct RippleModifier {
  const float *phase;
  float frequency;
  float amplitude;

  RippleModifier(const float *phase, float freq = 3.0f, float amp = 0.1f)
      : phase(phase), frequency(freq), amplitude(amp) {}

  float modify(float t) const {
    if (!phase)
      return t;
    // Calculate a local distortion based on the current t position
    return t + sinf(t * frequency * PI_F * 2.0f + *phase) * amplitude;
  }
};

/**
 * @brief Folds the palette back and forth like a kaleidoscope.
 * A folds value of 2.0 maps [0...1] to [0 -> 1 -> 0 -> 1 -> 0].
 */
struct FoldModifier {
  const float *phase;
  float folds;

  FoldModifier(float folds = 2.0f, const float *phase = nullptr)
      : phase(phase), folds(folds) {}

  float modify(float t) const {
    float shift = phase ? *phase : 0.0f;
    float scaled = (t * folds) + shift;

    // Triangle wave formula: creates symmetrical, continuous bouncing between 0
    // and 1
    return fabsf(fmodf(scaled, 2.0f) - 1.0f);
  }
};

/**
 * @brief Pinches or expands the center of the palette.
 * positive tension pulls colors toward the center, negative pushes them to the
 * edges. Creates an elastic, tension-release visual dynamic.
 */
struct PinchModifier {
  const float *tension; // Expects roughly -0.9 to 0.9

  PinchModifier(const float *t) : tension(t) {}

  float modify(float t) const {
    if (!tension)
      return t;

    // Shift t to -1.0 to 1.0 range based on a 0.0 to 1.0 domain block
    float wrapped_t = wrap_t(t);
    float centered = wrapped_t * 2.0f - 1.0f;
    float sign = centered < 0.0f ? -1.0f : 1.0f;

    // Apply power curve
    float amount = std::clamp(*tension, -0.99f, 0.99f);
    float power = (amount < 0.0f) ? (1.0f / (1.0f + std::abs(amount)))
                                  : (1.0f + amount * 3.0f);

    centered = sign * powf(std::abs(centered), power);

    // Map back to 0.0 to 1.0 and add to integer spatial frame
    return floorf(t) + ((centered + 1.0f) * 0.5f);
  }
};

/**
 * @brief Snaps smooth gradients into harsh, distinct bands (Posterization).
 * By animating the step count or offsetting before quantizing, you get a
 * retro/glitch aesthetic.
 */
struct QuantizeModifier {
  const float *dynamic_steps;
  float base_steps;

  QuantizeModifier(float steps, const float *d_steps = nullptr)
      : dynamic_steps(d_steps), base_steps(steps) {}

  float modify(float t) const {
    float s = dynamic_steps ? *dynamic_steps : base_steps;
    if (s < 1.0f)
      s = 1.0f;

    // Round to nearest step in the infinite domain
    return roundf(t * s) / s;
  }
};

/**
 * @brief Multiplies the palette coordinate, increasing the frequency
 * so the palette repeats multiple times across the domain.
 */
struct ScaleModifier {
  const float *dynamic_scale;
  float base_scale;

  ScaleModifier(float s = 1.0f, const float *d_scale = nullptr)
      : dynamic_scale(d_scale), base_scale(s) {}

  float modify(float t) const {
    return t * (dynamic_scale ? *dynamic_scale : base_scale);
  }
};

/**
 * @brief Variant holding any supported modifier type.
 */
using ModifierVariant =
    std::variant<CycleModifier, BreatheModifier, RippleModifier, FoldModifier,
                 PinchModifier, QuantizeModifier, ScaleModifier>;

class AnimatedPalette;
class CircularPalette;
class ReversePalette;
class VignettePalette;
class TransparentVignette;
class AlphaFalloffPalette;
class SolidColorPalette;

///////////////////////////////////////////////////////////////////////////////
// Palette Variants
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Variant holding any supported palette type.
 */
using PaletteVariant =
    std::variant<std::monostate, GenerativePalette, ProceduralPalette, Gradient,
                 MutatingPalette, ReversePalette, VignettePalette,
                 TransparentVignette, AlphaFalloffPalette, SolidColorPalette,
                 CircularPalette, AnimatedPalette>;

inline Color4 get_color(const PaletteVariant &pv, float t);

/**
 * @brief A wrapper that applies a chain of modifiers to the palette parameter
 * 't'.
 */
class AnimatedPalette {
public:
  AnimatedPalette(const PaletteVariant *source) : source(source) {}

  /**
   * Explicitly connect a modifier to this palette.
   * Stores a value.
   */
  AnimatedPalette &add(const ModifierVariant &modifier) {
    if (!modifiers.is_full()) {
      modifiers.push_back(modifier);
    }
    return *this;
  }

  Color4 get(float t) const;

private:
  const PaletteVariant *source;
  StaticCircularBuffer<ModifierVariant, 4> modifiers;

public:
  ~AnimatedPalette() = default;
};

/**
 * @brief Palette wrapper that makes the source palette symmetrical/circular.
 * Maps range [0, 1] -> [0, 1, 0].
 */
class CircularPalette {
public:
  CircularPalette(const PaletteVariant *palette) : palette(palette) {}
  Color4 get(float t) const;

private:
  const PaletteVariant *palette;
};

/**
 * @brief Palette wrapper that reverses the source palette.
 * Maps inputs t -> 1-t.
 */
class ReversePalette {
public:
  ReversePalette(const PaletteVariant *palette) : palette(palette) {}
  Color4 get(float t) const;

private:
  const PaletteVariant *palette;
};

/**
 * @brief Palette wrapper that applies black vignetting at the edges.
 * fades to black at t < 0.2 and t > 0.8.
 */
class VignettePalette {
public:
  VignettePalette(const PaletteVariant *palette) : palette(palette) {}

  Color4 get(float t) const;

private:
  const PaletteVariant *palette;
};

/**
 * @brief Palette wrapper that applies alpha transparency vignetting at the
 * edges.
 */
class TransparentVignette {
public:
  TransparentVignette(const PaletteVariant *palette) : palette(palette) {}

  Color4 get(float t) const;

private:
  const PaletteVariant *palette;
};

class AlphaFalloffPalette {
public:
  using FalloffFunction = float (*)(float);
  AlphaFalloffPalette(FalloffFunction fn, const PaletteVariant *source)
      : fn(fn), source(source) {}
  Color4 get(float t) const;

private:
  FalloffFunction fn;
  const PaletteVariant *source;
};

class SolidColorPalette {
public:
  SolidColorPalette(const Color4 &color) : color(color) {}
  Color4 get(float t) const { return color; }
  Color4 color;
};

inline Color4 get_color(const PaletteVariant &pv, float t) {
  return std::visit(
      [t](auto &&arg) -> Color4 {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, std::monostate>) {
          return Color4(CRGB(0, 0, 0), 0.0f);
        } else {
          return arg.get(t);
        }
      },
      pv);
}

// Implementations that require 'get_color' to be fully aware of the variants
inline Color4 AnimatedPalette::get(float t) const {
  float final_t = t;
  for (size_t i = 0; i < modifiers.size(); ++i) {
    final_t = std::visit([final_t](auto &&arg) { return arg.modify(final_t); },
                         modifiers[i]);
  }

  // Fast wrap the final coordinate right before querying the source palette
  final_t = wrap_t(final_t);

  return source ? get_color(*source, final_t) : Color4(CRGB(0, 0, 0), 0.0f);
}

inline Color4 CircularPalette::get(float t) const {
  return palette ? get_color(*palette, 1.0f - std::abs(2.0f * t - 1.0f))
                 : Color4(CRGB(0, 0, 0), 0.0f);
}

inline Color4 ReversePalette::get(float t) const {
  return palette ? get_color(*palette, 1.0f - t) : Color4(CRGB(0, 0, 0), 0.0f);
}

inline Color4 VignettePalette::get(float t) const {
  CRGB vignette_color(0, 0, 0);
  if (!palette)
    return Color4(vignette_color, 0.0f);

  if (t < 0.2f) {
    return Color4(vignette_color.lerp16(get_color(*palette, 0.0f).color,
                                        to_short(quintic_kernel(t / 0.2f))),
                  1.0f);
  } else if (t >= 0.8f) {
    return Color4(get_color(*palette, 1.0f)
                      .color.lerp16(vignette_color, to_short(quintic_kernel(
                                                        (t - 0.8f) / 0.2f))),
                  1.0f);
  } else {
    return get_color(*palette, (t - 0.2f) / 0.6f);
  }
}

inline Color4 TransparentVignette::get(float t) const {
  Color4 result;
  float factor = 1.0f;
  if (!palette)
    return result;

  if (t < 0.2f) {
    result = get_color(*palette, 0.0f);
    factor = quintic_kernel(t / 0.2f);
  } else if (t >= 0.8f) {
    result = get_color(*palette, 1.0f);
    factor = quintic_kernel(1.0f - (t - 0.8f) / 0.2f);
  } else {
    return get_color(*palette, (t - 0.2f) / 0.6f);
  }
  result.alpha *= factor;
  return result;
}

inline Color4 AlphaFalloffPalette::get(float t) const {
  if (!source)
    return Color4(CRGB(0, 0, 0), 0.0f);
  Color4 c = get_color(*source, t);
  c.alpha *= fn(t);
  return c;
}
