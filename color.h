/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

#include "platform.h"
#include "3dmath.h"
#include "util.h"
#include "static_circular_buffer.h"

#include "memory.h"

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

  Color4 &operator+=(const Color4 &rhs) {
    color += rhs.color;
    alpha += rhs.alpha;
    return *this;
  }

  Color4 &operator*=(float s) {
    color = color * s;
    alpha *= s;
    return *this;
  }

  operator CRGB() const { return color; }
};

/**
 * @brief Rotates the hue of a Color4 by a fractional amount (0..1 = full turn).
 * Uses Rodrigues rotation around the (1,1,1) gray axis in linear RGB space.
 */
inline Color4 hue_rotate(const Color4 &c, float amount) {
  float angle = amount * (2.0f * PI_F);
  float cos_a = cosf(angle);
  float sin_a = sinf(angle);
  float one_third = 1.0f / 3.0f;
  float sqrt3_inv = 0.57735026919f; // 1/sqrt(3)

  float r = c.color.r, g = c.color.g, b = c.color.b;
  float avg = (r + g + b) * one_third;
  float s = sin_a * sqrt3_inv;

  float nr = avg + (r - avg) * cos_a + (g - b) * s;
  float ng = avg + (g - avg) * cos_a + (b - r) * s;
  float nb = avg + (b - avg) * cos_a + (r - g) * s;

  Color4 result = c;
  result.color.r = static_cast<uint16_t>(std::clamp(nr, 0.0f, 65535.0f));
  result.color.g = static_cast<uint16_t>(std::clamp(ng, 0.0f, 65535.0f));
  result.color.b = static_cast<uint16_t>(std::clamp(nb, 0.0f, 65535.0f));
  return result;
}

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

// Inverse: Linear float 0-1 -> sRGB float 0-1
constexpr float linear_to_srgb_float(float l) {
  return (l <= 0.0031308f) ? l * 12.92f : 1.055f * powf(l, 1.0f / 2.4f) - 0.055f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// OKLab / OKLCH Color Space (Björn Ottosson, 2020)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct OKLab { float L, a, b; };
struct OKLCH { float L, C, h; };

/// Linear RGB [0,1] -> OKLab
inline OKLab linear_rgb_to_oklab(float r, float g, float b) {
  float l = 0.4122214708f * r + 0.5363325363f * g + 0.0514459929f * b;
  float m = 0.2119034982f * r + 0.6806995451f * g + 0.1073969566f * b;
  float s = 0.0883024619f * r + 0.2817188376f * g + 0.6299787005f * b;

  float l_ = cbrtf(l), m_ = cbrtf(m), s_ = cbrtf(s);

  return {0.2104542553f * l_ + 0.7936177850f * m_ - 0.0040720468f * s_,
          1.9779984951f * l_ - 2.4285922050f * m_ + 0.4505937099f * s_,
          0.0259040371f * l_ + 0.7827717662f * m_ - 0.8086757660f * s_};
}

/// OKLab -> Linear RGB [0,1]
inline void oklab_to_linear_rgb(OKLab lab, float &r, float &g, float &b) {
  float l_ = lab.L + 0.3963377774f * lab.a + 0.2158037573f * lab.b;
  float m_ = lab.L - 0.1055613458f * lab.a - 0.0638541728f * lab.b;
  float s_ = lab.L - 0.0894841775f * lab.a - 1.2914855480f * lab.b;

  float l = l_ * l_ * l_, m = m_ * m_ * m_, s = s_ * s_ * s_;

  r = +4.0767416621f * l - 3.3077115913f * m + 0.2309699292f * s;
  g = -1.2684380046f * l + 2.6097574011f * m - 0.3413193965f * s;
  b = -0.0041960863f * l - 0.7034186147f * m + 1.7076147010f * s;
}

inline OKLCH oklab_to_oklch(OKLab lab) {
  float C = sqrtf(lab.a * lab.a + lab.b * lab.b);
  float h = atan2f(lab.b, lab.a);
  return {lab.L, C, h};
}

inline OKLab oklch_to_oklab(OKLCH lch) {
  return {lch.L, lch.C * cosf(lch.h), lch.C * sinf(lch.h)};
}

/// sRGB [0-255] -> OKLCH (convenience)
inline OKLCH srgb_to_oklch(uint8_t r, uint8_t g, uint8_t b) {
  float rf = srgb_to_linear_float(r / 255.0f);
  float gf = srgb_to_linear_float(g / 255.0f);
  float bf = srgb_to_linear_float(b / 255.0f);
  return oklab_to_oklch(linear_rgb_to_oklab(rf, gf, bf));
}

/// OKLCH -> 16-bit linear Pixel (gamut-clamped)
inline Pixel oklch_to_pixel(OKLCH lch) {
  float r, g, b;
  oklab_to_linear_rgb(oklch_to_oklab(lch), r, g, b);
  return Pixel(static_cast<uint16_t>(hs::clamp(r, 0.0f, 1.0f) * 65535.0f),
               static_cast<uint16_t>(hs::clamp(g, 0.0f, 1.0f) * 65535.0f),
               static_cast<uint16_t>(hs::clamp(b, 0.0f, 1.0f) * 65535.0f));
}

/// Interpolate two OKLCH colors (shortest-arc hue)
inline OKLCH lerp_oklch(OKLCH a, OKLCH b, float t) {
  // Handle achromatic cases (near-zero chroma)
  float h;
  if (a.C < 1e-4f && b.C < 1e-4f) {
    h = 0.0f;
  } else if (a.C < 1e-4f) {
    h = b.h;
  } else if (b.C < 1e-4f) {
    h = a.h;
  } else {
    // Shortest arc interpolation
    float dh = b.h - a.h;
    if (dh > PI_F) dh -= 2.0f * PI_F;
    if (dh < -PI_F) dh += 2.0f * PI_F;
    h = a.h + dh * t;
  }
  return {a.L + (b.L - a.L) * t, a.C + (b.C - a.C) * t, h};
}

/// Interpolate two sRGB CPixels in OKLCH space -> 16-bit linear Pixel
inline Pixel lerp_oklch(const CPixel &c1, const CPixel &c2, float t) {
  OKLCH a = srgb_to_oklch(c1.r, c1.g, c1.b);
  OKLCH b = srgb_to_oklch(c2.r, c2.g, c2.b);
  return oklch_to_pixel(lerp_oklch(a, b, t));
}

/// Interpolate two sRGB CPixels in OKLCH space -> sRGB CPixel
inline CPixel lerp_oklch_srgb(const CPixel &c1, const CPixel &c2, float t) {
  OKLCH a = srgb_to_oklch(c1.r, c1.g, c1.b);
  OKLCH b = srgb_to_oklch(c2.r, c2.g, c2.b);
  OKLCH result = lerp_oklch(a, b, t);
  float r, g, b_val;
  oklab_to_linear_rgb(oklch_to_oklab(result), r, g, b_val);
  return CPixel(
    static_cast<uint8_t>(hs::clamp(linear_to_srgb_float(hs::clamp(r, 0.0f, 1.0f)) * 255.0f + 0.5f, 0.0f, 255.0f)),
    static_cast<uint8_t>(hs::clamp(linear_to_srgb_float(hs::clamp(g, 0.0f, 1.0f)) * 255.0f + 0.5f, 0.0f, 255.0f)),
    static_cast<uint8_t>(hs::clamp(linear_to_srgb_float(hs::clamp(b_val, 0.0f, 1.0f)) * 255.0f + 0.5f, 0.0f, 255.0f)));
}

/// Abstract base for all palettes. Provides a uniform interface for color
/// lookup, replacing std::variant + std::visit with a single vtable pointer.
class Palette {
public:
  virtual Color4 get(float t) const = 0;
  virtual ~Palette() = default;
};

/// Abstract base for all palette modifiers.
class Modifier {
public:
  virtual float modify(float t) const = 0;
  virtual ~Modifier() = default;
};

class Gradient : public Palette {
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
        // Pre-convert stop colors to linear float
        float pr = srgb_to_linear_float(prevColor.r / 255.0f);
        float pg = srgb_to_linear_float(prevColor.g / 255.0f);
        float pb = srgb_to_linear_float(prevColor.b / 255.0f);
        float nr = srgb_to_linear_float(nextColor.r / 255.0f);
        float ng = srgb_to_linear_float(nextColor.g / 255.0f);
        float nb = srgb_to_linear_float(nextColor.b / 255.0f);
        for (int i = start; i <= end; i++) {
          float t = static_cast<float>(i - start) / (end - start);

          // Interpolate in linear space
          float r_lin = pr * (1.0f - t) + nr * t;
          float g_lin = pg * (1.0f - t) + ng * t;
          float b_lin = pb * (1.0f - t) + nb * t;

          entries[i] = Pixel(
            static_cast<uint16_t>(r_lin * 65535.0f),
            static_cast<uint16_t>(g_lin * 65535.0f),
            static_cast<uint16_t>(b_lin * 65535.0f));
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

  Color4 get(float t) const override {
    uint8_t index = static_cast<uint8_t>(t * 255);
    return Color4(entries[index], 1.0f);
  }

  ~Gradient() {}
};

/**
 * @brief A palette that generates colors based on defined harmony and
 * brightness profiles.
 */
class GenerativePalette : public Palette {
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

  /// Lightweight snapshot of just the 3 color keys (9 bytes).
  struct Snapshot {
    CPixel a, b, c;
  };

  Snapshot snapshot() const { return {a, b, c}; }

  void lerp(const GenerativePalette &from, const GenerativePalette &to,
            float amount) {
    lerp(from.snapshot(), to.snapshot(), amount);
  }

  void lerp(const Snapshot &from, const Snapshot &to, float amount) {
    a = lerp_oklch_srgb(from.a, to.a, amount);
    b = lerp_oklch_srgb(from.b, to.b, amount);
    c = lerp_oklch_srgb(from.c, to.c, amount);
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

    // Interpolate in OKLCH for perceptually uniform gradients
    return Color4(lerp_oklch(c1, c2, p), 1.0f);
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
class ProceduralPalette : public Palette {
public:
  constexpr ProceduralPalette()
      : a{0, 0, 0}, b{0, 0, 0}, c{0, 0, 0}, d{0, 0, 0} {}
  constexpr ProceduralPalette(std::array<float, 3> a, std::array<float, 3> b,
                              std::array<float, 3> c, std::array<float, 3> d)
      : a(a), b(b), c(c), d(d) {}

  Color4 get(float t) const override {
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
struct CycleModifier : public Modifier {
  const float *offset;

  CycleModifier(const float *driver_offset = nullptr) : offset(driver_offset) {}

  float modify(float t) const override { return offset ? t + *offset : t; }
};

/**
 * @brief Oscillates the palette coordinate (Breathing).
 */
struct BreatheModifier : public Modifier {
  const float *phase;
  float amplitude;

  BreatheModifier(const float *driver_phase, float amp = 0.1f)
      : phase(driver_phase), amplitude(amp) {}

  float modify(float t) const override {
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
struct RippleModifier : public Modifier {
  const float *phase;
  float frequency;
  float amplitude;

  RippleModifier(const float *phase, float freq = 3.0f, float amp = 0.1f)
      : phase(phase), frequency(freq), amplitude(amp) {}

  float modify(float t) const override {
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
struct FoldModifier : public Modifier {
  const float *phase;
  float folds;

  FoldModifier(float folds = 2.0f, const float *phase = nullptr)
      : phase(phase), folds(folds) {}

  float modify(float t) const override {
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
struct PinchModifier : public Modifier {
  const float *tension; // Expects roughly -0.9 to 0.9

  PinchModifier(const float *t) : tension(t) {}

  float modify(float t) const override {
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
struct QuantizeModifier : public Modifier {
  const float *dynamic_steps;
  float base_steps;

  QuantizeModifier(float steps, const float *d_steps = nullptr)
      : dynamic_steps(d_steps), base_steps(steps) {}

  float modify(float t) const override {
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
struct ScaleModifier : public Modifier {
  const float *dynamic_scale;
  float base_scale;

  ScaleModifier(float s = 1.0f, const float *d_scale = nullptr)
      : dynamic_scale(d_scale), base_scale(s) {}

  float modify(float t) const override {
    return t * (dynamic_scale ? *dynamic_scale : base_scale);
  }
};

///////////////////////////////////////////////////////////////////////////////
// Compile-Time Palette Composition
///////////////////////////////////////////////////////////////////////////////

/// Structural concept: any type with a const modify(float)->float method.
/// All existing Modifier subclasses satisfy this automatically.
template <typename T>
concept PaletteModifier = requires(const T m, float t) {
  { m.modify(t) } -> std::convertible_to<float>;
};

/**
 * @brief A compile-time composition of a Source palette and N modifiers.
 * @details Default constructible to allow safe wiring in init().
 * Uses a fold expression to inline the modifier chain with zero runtime
 * overhead. Follows the ArenaVector idiom: default construct, then bind().
 */
template <typename Source, PaletteModifier... Mods>
class StaticPalette {
public:
  StaticPalette() = default;

  void bind(const Source *src, const Mods *...ms) {
    source_ = src;
    mods_ = std::make_tuple(ms...);
  }

  Color4 get(float t) const {
    assert(source_ != nullptr && "StaticPalette used before bind()!");

    float ft = t;
    std::apply([&](const auto *...m) { ((ft = m->modify(ft)), ...); }, mods_);

    return source_->get(wrap_t(ft));
  }

private:
  const Source *source_ = nullptr;
  std::tuple<const Mods *...> mods_{};
};

class AnimatedPalette;
class CircularPalette;
class ReversePalette;
class VignettePalette;
class TransparentVignette;
class AlphaFalloffPalette;
class SolidColorPalette;

inline Color4 get_color(const Palette &p, float t) { return p.get(t); }
inline Color4 get_color(const Palette *p, float t) {
  return p ? p->get(t) : Color4(CRGB(0, 0, 0), 0.0f);
}

/**
 * @brief A wrapper that applies a chain of modifiers to the palette parameter
 * 't'.
 */
class AnimatedPalette : public Palette {
public:
  AnimatedPalette(const Palette *source) : source(source) {}

  /**
   * Explicitly connect a modifier to this palette.
   * Stores a pointer (modifier must outlive this palette).
   */
  AnimatedPalette &add(const Modifier &modifier) {
    if (num_modifiers < MAX_MODIFIERS) {
      modifiers[num_modifiers++] = &modifier;
    }
    return *this;
  }

  Color4 get(float t) const override;

private:
  static constexpr int MAX_MODIFIERS = 4;
  const Palette *source;
  const Modifier *modifiers[MAX_MODIFIERS] = {};
  int num_modifiers = 0;
};

/**
 * @brief Palette wrapper that makes the source palette symmetrical/circular.
 * Maps range [0, 1] -> [0, 1, 0].
 */
class CircularPalette : public Palette {
public:
  CircularPalette(const Palette *palette) : palette(palette) {}
  Color4 get(float t) const override;

private:
  const Palette *palette;
};

/**
 * @brief Palette wrapper that reverses the source palette.
 * Maps inputs t -> 1-t.
 */
class ReversePalette : public Palette {
public:
  ReversePalette(const Palette *palette) : palette(palette) {}
  Color4 get(float t) const override;

private:
  const Palette *palette;
};

/**
 * @brief Palette wrapper that applies black vignetting at the edges.
 * fades to black at t < 0.2 and t > 0.8.
 */
class VignettePalette : public Palette {
public:
  VignettePalette(const Palette *palette) : palette(palette) {}

  Color4 get(float t) const override;

private:
  const Palette *palette;
};

/**
 * @brief Palette wrapper that applies alpha transparency vignetting at the
 * edges.
 */
class TransparentVignette : public Palette {
public:
  TransparentVignette(const Palette *palette) : palette(palette) {}

  Color4 get(float t) const override;

private:
  const Palette *palette;
};

class AlphaFalloffPalette : public Palette {
public:
  using FalloffFunction = float (*)(float);
  AlphaFalloffPalette(FalloffFunction fn, const Palette *source)
      : fn(fn), source(source) {}
  Color4 get(float t) const override;

private:
  FalloffFunction fn;
  const Palette *source;
};

class SolidColorPalette : public Palette {
public:
  SolidColorPalette(const Color4 &color) : color(color) {}
  Color4 get(float t) const override { return color; }
  Color4 color;
};

/// Pre-baked 256-entry Pixel16 LUT allocated in an arena.
/// Converts any Palette into a fast table lookup with lerp interpolation.
/// Not a Palette subclass — call get(t) directly for zero-overhead lookups.
class BakedPalette {
public:
  static constexpr int LUT_SIZE = 256;

  BakedPalette() = default;

  /// Bake a palette into a 256-entry LUT in the given arena.
  void bake(Arena &arena, const Palette &source) {
    lut_ = static_cast<Pixel16 *>(
        arena.allocate(LUT_SIZE * sizeof(Pixel16), alignof(Pixel16)));
    for (int i = 0; i < LUT_SIZE; ++i) {
      float t = static_cast<float>(i) / (LUT_SIZE - 1);
      Color4 c = source.get(t);
      lut_[i] = c.color;
    }
  }

  /// Fast lookup with linear interpolation between adjacent entries.
  Color4 get(float t) const {
    float idx = t * (LUT_SIZE - 1);
    int lo = static_cast<int>(idx);
    if (lo >= LUT_SIZE - 1) return Color4(lut_[LUT_SIZE - 1], 1.0f);
    if (lo < 0) return Color4(lut_[0], 1.0f);
    uint16_t frac = static_cast<uint16_t>((idx - lo) * 65535.0f);
    return Color4(lut_[lo].lerp16(lut_[lo + 1], frac), 1.0f);
  }

private:
  Pixel16 *lut_ = nullptr;
};

// Implementations that require all palette types to be complete
inline Color4 AnimatedPalette::get(float t) const {
  float final_t = t;
  for (int i = 0; i < num_modifiers; ++i) {
    final_t = modifiers[i]->modify(final_t);
  }

  // Fast wrap the final coordinate right before querying the source palette
  final_t = wrap_t(final_t);

  return get_color(source, final_t);
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
